#
# The Chronus Quantum (ChronusQ) software package is high-performace 
# computational chemistry software with a strong emphasis on explicitly 
# time-dependent and post-SCF quantum mechanical methods.
# 
# Copyright (C) 2014-2015 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
# 
#
import os,sys
import libpythonapi as chronusQ
from libpythonapi import CErrMsg
from parseBasis import parseBasis
from meta.knownKeywords import requiredKeywords
from meta.knownJobs import *
from meta.enumMaps import sdrMethodMap
from meta.enumMaps import aointAlg
from meta.enumMaps import guessMap
from meta.enumMaps import exchMap 
from meta.enumMaps import corrMap 
from meta.enumMaps import kernelMap 
from meta.enumMaps import gridMap 
from meta.enumMaps import dftWeightScheme 
from meta.enumMaps import envMap

#
#  Parse the QM section of the input file and populate
#  the CQ::SingleSlater<T> and CQ::BasisSet objects
#  as well as any subsequent job objects:
#
#    - RealTime<T>
#    - SDResponse<T>
#
#
#  ** Note that this does DOES NOT allocate the SingleSlater **
#  ** object, but does allocate the BasisSet object          **
#
#  Input:
#    workers      -     the workers array of CQ objects
#    secDict      -     total parsed section dictionary
#
#  Output:
#    JobStr       -     string containing a map to standard
#                       job functions
#
def parseQM(workers,secDict): 
#  print 'Parsing QM Information'
  ssSettings = secDict["QM"]

#
# Check for unknown keywords in the QM section
#
# knownKeywords = [ 'reference', 'basis', 'job' ]
# for i in ssSettings:
#   if i not in knownKeywords:
#     print "Keyword QM."+ str(i) +" not recognized"

#
#  Check that all of the required keywords for Molecule
#  object are found
#
  for i in requiredKeywords['QM']:
    if i not in ssSettings:
      msg = 'Required keyword QM.' + str(i) + ' not found'
      CErrMsg(workers['CQFileIO'],msg)

#
# Try to set the reference for CQ::SingleSlater
#
  handleReference(workers,ssSettings)

#
# Try to set the basis for the QM Job
#
  parseBasis(workers,ssSettings)


#
# Set integral algorithm
#

  try:
    workers['CQAOIntegrals'].setAlgorithm(aointAlg[str(ssSettings['INTS'])])
  except KeyError:
    pass

  try:
    parseSCF(workers,secDict['SCF'])
  except KeyError:
    pass

  if str(ssSettings['JOB']) in knownJobs:
    if ssSettings['JOB'] in ('RT'):
      parseRT(workers,secDict['RT']) 
    elif ssSettings['JOB'] in ('RPA','CIS','STAB','PPRPA','PPATDA','PPCTDA'):
      if chronusQ.getSize() > 1:
        msg = "Response cannot run with >1 MPI Processes"
        CErrMsg(workers['CQFileIO'],msg)
      else:
        parseSDR(workers,secDict)
  else:
    msg = 'QM.Job ' + str(ssSettings['JOB']) + ' not recognized'
    CErrMsg(workers['CQFileIO'],str(msg))

#
# Set global print level
#
  try:
    workers['CQMolecule'].setPrintLevel(ssSettings['PRINT'])
    workers['CQBasisSet'].setPrintLevel(ssSettings['PRINT'])
    workers['CQSingleSlater'].setPrintLevel(ssSettings['PRINT'])
    workers['CQRealTime'].setPrintLevel(ssSettings['PRINT'])
  except KeyError:
    pass

  return str(ssSettings['JOB'])

#  # Space filler to pasify error 
#  workers["CQSingleSlaterDouble"].communicate(
#    workers["CQMolecule"], workers["CQBasisSet"], workers["CQAOIntegrals"],
#    workers["CQFileIO"], workers["CQControls"]
#  )
#  workers["CQSingleSlaterDouble"].initMeta()
#  workers["CQSingleSlaterDouble"].genMethString()
#  workers["CQSingleSlaterDouble"].alloc()
#  runSCF(workers)  

#
# Set the reference for CQ::SingleSlater
#
def handleReference(workers,settings):
  mult = workers["CQMolecule"].multip()
  ref = settings['REFERENCE']
  ref = ref.split()
  # Decide if reference is complex or not (defaults to real if not specified
  if 'COMPLEX' in ref:
    workers["CQSingleSlater"] = workers["CQSingleSlaterComplex"]
  else:
    workers["CQSingleSlater"] = workers["CQSingleSlaterDouble"]

  # Set SS Reference
  if 'HF' in ref:
    # Smartly figure out of reference is R/U
    if mult == 1:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
      workers["CQSingleSlater"].isClosedShell = True
    else:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.UHF)
  elif 'RHF' in ref:
    # Force RHF
    if mult != 1:
      msg = 'Non-singlet multiplicity is not suitable for RHF'
      CErrMsg(workers['CQFileIO'],str(msg))
 
    workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
    workers["CQSingleSlater"].isClosedShell = True
  elif 'UHF' in ref:
    # Forch UHF
    workers["CQSingleSlater"].setRef(chronusQ.Reference.UHF)
    workers["CQSingleSlater"].isClosedShell = False
  elif 'CUHF' in ref:
    # Use Constrained UHF (not complex UHF) 
    workers["CQSingleSlater"].setRef(chronusQ.Reference.CUHF)
  elif 'GHF' in ref:
    # Do GHF
    workers["CQSingleSlater"].setRef(chronusQ.Reference.TCS)
  elif 'KS' in ref:
    # Smartly figure out of reference is R/U
    if mult == 1:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
      workers["CQSingleSlater"].isClosedShell = True
    else:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.UHF)

    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['USERDEFINED'])
  
    corrKernel = 0
    exchKernel = 0
    try:
      corrKernel = settings['CORR']
      exchKernel = settings['EXCHANGE']
    except KeyError:
      msg = "Must specify both Correlation and Exchange Kernel\n"
      msg = msg + " for user defined QM.KS reference"
      CErrMsg(workers['CQFileIO'],str(msg))


    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])

    try:
      workers["CQSingleSlater"].setCorrKernel(corrMap[corrKernel])
    except KeyError:
      msg = "Specified Correlation Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    try:
      workers["CQSingleSlater"].setExchKernel(exchMap[exchKernel])
    except KeyError:
      msg = "Specified Exchange Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    workers["CQSingleSlater"].checkDFTType()
  elif 'RKS' in ref:
    # Force RKS
    if mult != 1:
      msg = 'Non-singlet multiplicity is not suitable for RKS'
      CErrMsg(workers['CQFileIO'],str(msg))
 
    workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
    workers["CQSingleSlater"].isClosedShell = True

    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['USERDEFINED'])
  
    corrKernel = 0
    exchKernel = 0
    try:
      corrKernel = settings['CORR']
      exchKernel = settings['EXCHANGE']
    except KeyError:
      msg = "Must specify both Correlation and Exchange Kernel\n"
      msg = msg + " for user defined QM.KS reference"
      CErrMsg(workers['CQFileIO'],str(msg))


    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])


    try:
      workers["CQSingleSlater"].setCorrKernel(corrMap[corrKernel])
    except KeyError:
      msg = "Specified Correlation Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    try:
      workers["CQSingleSlater"].setExchKernel(exchMap[exchKernel])
    except KeyError:
      msg = "Specified Exchange Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    workers["CQSingleSlater"].checkDFTType()
  elif 'UKS' in ref:
    # forch uhf
    workers["CQSingleSlater"].setRef(chronusQ.Reference.UHF)
    workers["CQSingleSlater"].isClosedShell = False


    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['USERDEFINED'])
  
    corrKernel = 0
    exchKernel = 0
    try:
      corrKernel = settings['CORR']
      exchKernel = settings['EXCHANGE']
    except KeyError:
      msg = "Must specify both Correlation and Exchange Kernel\n"
      msg = msg + " for user defined QM.KS reference"
      CErrMsg(workers['CQFileIO'],str(msg))

    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])

    try:
      workers["CQSingleSlater"].setCorrKernel(corrMap[corrKernel])
    except KeyError:
      msg = "Specified Correlation Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    try:
      workers["CQSingleSlater"].setExchKernel(exchMap[exchKernel])
    except KeyError:
      msg = "Specified Exchange Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    workers["CQSingleSlater"].checkDFTType()
  elif 'CUKS' in ref:
    # Use Constrained UHF (not complex UHF) 
    workers["CQSingleSlater"].setRef(chronusQ.Reference.CUHF)
    workers["CQSingleSlater"].isClosedShell = False

    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['USERDEFINED'])
  
    corrKernel = 0
    exchKernel = 0
    try:
      corrKernel = settings['CORR']
      exchKernel = settings['EXCHANGE']
    except KeyError:
      msg = "Must specify both Correlation and Exchange Kernel\n"
      msg = msg + " for user defined QM.KS reference"
      CErrMsg(workers['CQFileIO'],str(msg))


    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])

    try:
      workers["CQSingleSlater"].setCorrKernel(corrMap[corrKernel])
    except KeyError:
      msg = "Specified Correlation Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    try:
      workers["CQSingleSlater"].setExchKernel(exchMap[exchKernel])
    except KeyError:
      msg = "Specified Exchange Kernel is not Defined"
      CErrMsg(workers['CQFileIO'],str(msg))

    workers["CQSingleSlater"].checkDFTType()
  elif 'LSDA' in ref:
    # Smartly figure out of reference is R/U
    if mult == 1:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
      workers["CQSingleSlater"].isClosedShell = True
    else:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.UHF)

    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['LSDA'])
    workers["CQSingleSlater"].setCorrKernel(corrMap['VWN3'])
    workers["CQSingleSlater"].setExchKernel(exchMap['SLATER'])

    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])
  elif 'RLSDA' in ref:
    # Force RKS
    if mult != 1:
      msg = 'Non-singlet multiplicity is not suitable for RKS'
      CErrMsg(workers['CQFileIO'],str(msg))
 
    workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
    workers["CQSingleSlater"].isClosedShell = True

    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['LSDA'])
    workers["CQSingleSlater"].setCorrKernel(corrMap['VWN3'])
    workers["CQSingleSlater"].setExchKernel(exchMap['SLATER'])

    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])
  elif 'ULSDA' in ref:
    # Force RKS
    if mult != 1:
      msg = 'Non-singlet multiplicity is not suitable for RKS'
      CErrMsg(workers['CQFileIO'],str(msg))
 
    workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
    workers["CQSingleSlater"].isClosedShell = True

    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['LSDA'])
    workers["CQSingleSlater"].setCorrKernel(corrMap['VWN3'])
    workers["CQSingleSlater"].setExchKernel(exchMap['SLATER'])

    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])
  elif 'CULSDA' in ref:
    # Use Constrained UHF (not complex UHF) 
    workers["CQSingleSlater"].setRef(chronusQ.Reference.CUHF)

    workers["CQSingleSlater"].isDFT = True
    workers["CQSingleSlater"].isHF  = False
    workers["CQSingleSlater"].setDFTKernel(kernelMap['LSDA'])
    workers["CQSingleSlater"].setCorrKernel(corrMap['VWN3'])
    workers["CQSingleSlater"].setExchKernel(exchMap['SLATER'])

    if 'DFT_GRID' in settings:
      workers["CQSingleSlater"].setDFTGrid(gridMap[settings['DFT_GRID']])
    if 'DFT_WEIGHTS' in settings:
      workers["CQSingleSlater"].setDFTWeightScheme(
        dftWeightScheme[settings['DFT_WEIGHTS']])

    if 'DFT_SCREEN' in settings:
      if settings['DFT_SCREEN']:
        pass
      else:
        workers["CQSingleSlater"].turnOffDFTScreening()

    if 'DFT_SCRTOL' in settings:
      workers["CQSingleSlater"].setDFTScreenTol(settings['DFT_SCRTOL'])

    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])

    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])


  else:
    msg = 'Reference ' + str(sum(ref)) + ' not able to be parsed'
    CErrMsg(workers['CQFileIO'],str(msg))

  # Check if reference is 2-Component
  TCMethods = [chronusQ.Reference.TCS, chronusQ.Reference.GKS]
  if workers["CQSingleSlater"].Ref() in TCMethods:
    workers["CQSingleSlater"].setNTCS(2)
    workers["CQAOIntegrals"].setNTCS(2)

def parseRT(workers,settings):

  # Set RT object based on SS reference
  if workers['CQSingleSlater'] == workers['CQSingleSlaterDouble']:
    workers['CQRealTime'] = workers['CQRealTimeDouble']
  else:
    workers['CQRealTime'] = workers['CQRealTimeComplex']

  # Define a map from keyword to function
  optMap = {
    'MAXSTEP'  :workers['CQRealTime'].setMaxSteps ,
    'TIMESTEP' :workers['CQRealTime'].setStepSize ,
    'EDFIELD'  :workers['CQRealTime'].setFieldAmp ,
    'TIME_ON'  :workers['CQRealTime'].setTOn      ,
    'TIME_OFF' :workers['CQRealTime'].setTOff     ,
    'FREQUENCY':workers['CQRealTime'].setFreq     ,
    'PHASE'    :workers['CQRealTime'].setPhase    ,
    'SIGMA'    :workers['CQRealTime'].setSigma    ,
    'ENVELOPE' :workers['CQRealTime'].setEnvelope ,
    'ORTHO'    :workers['CQRealTime'].setOrthoTyp ,
    'INIDEN'   :workers['CQRealTime'].setInitDen  ,
    'UPROP'    :workers['CQRealTime'].setFormU    ,
    'ELL_POL'  :workers['CQRealTime'].setEllPol     
  }


  # Loop over optional keywords, set options accordingly
  # note that because these are optional, if the keyword
  # is not found in setings, no error is thrown and the
  # next keyword is processed
# for i in optMap:
#   try:
#     if i not in ('EDFIELD'):
#       optMap[i](settings[i])
#     else:
#       optMap[i](settings[i][0],settings[i][1],settings[i][2])
#   except KeyError:
#     continue
  for i in optMap:
    if (i in settings) and (i not in ('EDFIELD')):
      optMap[i](settings[i])
    elif (i in settings) and (i in('EDFIELD')):
      optMap[i](settings[i][0],settings[i][1],settings[i][2])


  #if 'EDFIELD' not in settings:
  #  optMap['EDFIELD'](0.0,0.0,0.0)
    
  if 'TARCSVS' in settings:
    if not settings['TARCSVS']:
      workers['CQRealTime'].doNotTarCSV()

  # Idiot Checks

  if 'ENVELOPE' in settings:
    env = settings['ENVELOPE']

    if (env == envMap['PW']) and ('EDFIELD' not in settings):
      msg = "Must specify EDFIELD with PW envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))
 
    if (env == envMap['LINRAMP']) and ('EDFIELD' not in settings):
      msg = "Must specify EDFIELD with LINRAMP envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))
 
    if (env == envMap['LINRAMP']) and ('FREQUENCY' not in settings):
      msg = "Must specify FREQUENCY with LINRAMP envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env == envMap['GAUSSIAN']) and ('EDFIELD' not in settings):
      msg = "Must specify EDFIELD with GAUSSIAN envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env == envMap['GAUSSIAN']) and ('SIGMA' not in settings):
      msg = "Must specify SIGMA with GAUSSIAN envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env == envMap['STEP']) and ('TIME_ON' not in settings):
      msg = "Must specify TIME_ON with STEP envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env == envMap['STEP']) and ('TIME_OFF' not in settings):
      msg = "Must specify TIME_OFF with STEP envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env == envMap['SINSQ']):
      msg = "SINSQ envelope NYI"
      CErrMsg(workers['CQFileIO'],str(msg))



def parseSDR(workers,secDict):
  jobSettings = {}
  JOB = secDict['QM']['JOB']
  try:
    jobSettings = secDict[JOB]
  except KeyError:
    if JOB in ('STAB'):
      pass
    else: 
      msg = "Must specify an options sections for " + JOB
      CErrMsg(workers['CQFileIO'],str(msg))

   
  # Set SDR object based on SS reference
  if workers['CQSingleSlater'] == workers['CQSingleSlaterDouble']:
    workers['CQSDResponse']  = workers['CQSDResponseDouble']
    workers['CQMOIntegrals'] = workers['CQMOIntegralsDouble']
  elif workers['CQSingleSlater'] == workers['CQSingleSlaterComplex']:
    msg = "Wave Function Response using a Complex Reference is NYI"
    CErrMsg(workers['CQFileIO'],str(msg))
    workers['CQSDResponse']  = workers['CQSDResponseComplex']
    workers['CQMOIntegrals'] = workers['CQMOIntegralsComplex']
  else:
    msg = "Error in Reference Recieved by SDResponse"
    CErrMsg(workers['CQFileIO'],str(msg))

  if workers['CQSingleSlater'].nTCS() == 2 and workers['CQAOIntegrals'].integralAlgorithm != aointAlg['INCORE']:
    msg = "Wave Function Response using a Two-Component Reference is\n"
    msg = msg + "Only Implemented using INCORE integrals (QM.ints = INCORE)"
    CErrMsg(workers['CQFileIO'],str(msg))
    

  try:
    workers['CQSDResponse'].setNSek(jobSettings['NSTATES'])
  except KeyError:
    if JOB in ('STAB'):
      workers['CQSDResponse'].setNSek(3)
    else: 
      msg = "Must specify number of desired roots for " + JOB
      CErrMsg(workers['CQFileIO'],str(msg))
     
  workers['CQSDResponse'].setMeth(sdrMethodMap[str(JOB)])

def parseSCF(workers,scfSettings):
  optMap = {
    'SCFDENTOL' :workers['CQSingleSlater'].setSCFDenTol,
    'SCFENETOL' :workers['CQSingleSlater'].setSCFEneTol,
    'SCFMAXITER':workers['CQSingleSlater'].setSCFMaxIter,
    'FIELD'     :workers['CQSingleSlater'].setField,
    'GUESS'     :workers['CQSingleSlater'].setGuess ,
    'PRINT'     :workers['CQSingleSlater'].setPrintLevel
  }
  # Loop over optional keywords, set options accordingly
  # note that because these are optional, if the keyword
  # is not found in setings, no error is thrown and the
  # next keyword is processed
  for i in optMap:
    try:
      if i not in ('FIELD'):
        optMap[i](scfSettings[i])
      else:
        optMap[i](scfSettings[i][0],scfSettings[i][1],scfSettings[i][2])
    except KeyError:
      continue

  try:
    if scfSettings['GUESS'] == guessMap['READ']:
      workers['CQFileIO'].doRestart = True
  except KeyError:
    pass

  try:
    if scfSettings['DIIS']:
      workers['CQSingleSlater'].doDIIS = True
    else:
      workers['CQSingleSlater'].doDIIS = False
  except KeyError:
    pass

