#
# This file is part of the Chronus Quantum (ChronusQ) software package
# 
# Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
#from meta.enumMaps import sdrMethodMap
#from meta.enumMaps import aointAlg
#from meta.enumMaps import guessMap
#from meta.enumMaps import exchMap 
#from meta.enumMaps import corrMap 
#from meta.enumMaps import kernelMap 
#from meta.enumMaps import gridMap 
#from meta.enumMaps import dftWeightScheme 
#from meta.enumMaps import envMap

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

#
#  Check that all of the required keywords for QM
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
# Default to CORE guess for 1 atom
#
  if workers['CQMolecule'].nAtoms() == 1:
    workers['CQSingleSlater'].setGuess(str('CORE'))

#  if str(ssSettings['JOB']) in knownJobs:
#    if ssSettings['JOB'] in ('RT'):
#      pass
##      parseRT(workers,secDict['RT']) 
#    elif ssSettings['JOB'] in ('RPA','CIS','STAB','PPRPA','PPATDA','PPCTDA'):
#      if chronusQ.getSize() > 1:
#        msg = "Response cannot run with >1 MPI Processes"
#        CErrMsg(workers['CQFileIO'],msg)
##     else:
##       parseSDR(workers,secDict)
#  else:
#    msg = 'QM.Job ' + str(ssSettings['JOB']) + ' not recognized'
#    CErrMsg(workers['CQFileIO'],str(msg))

  if str(ssSettings['JOB']) not in knownJobs:
    msg = 'QM.Job ' + str(ssSettings['JOB']) + ' not recognized'
    CErrMsg(workers['CQFileIO'],str(msg))

  return str(ssSettings['JOB'])

#
# Set the reference for CQ::SingleSlater
#
def handleReference(workers,settings):
  ref = settings['REFERENCE']
  ref = ref.split()
  # Decide if reference is complex or not (defaults to real if not specified
  if 'COMPLEX' in ref or 'X2C' in ref:
    workers["CQSingleSlater"] = chronusQ.SingleSlater_complex()
  else:
    workers["CQSingleSlater"] = chronusQ.SingleSlater_double()

  isDFT = False

  if len(ref) > 1:
    workers["CQSingleSlater"].setRef(str(ref[1]))
    isDFT = 'HF' not in ref[1] and 'X2C' not in ref[1]
  else:
    workers["CQSingleSlater"].setRef(str(ref[0]))
    isDFT = 'HF' not in ref[0] and 'X2C' not in ref[0]


  if isDFT:
    if 'DFT_NRAD' in settings:
      workers["CQSingleSlater"].setDFTNRad(settings['DFT_NRAD'])
    if 'DFT_NANG' in settings:
      workers["CQSingleSlater"].setDFTNAng(settings['DFT_NANG'])

    # No 2C DFT Yet
    if isDFT and workers["CQSingleSlater"].nTCS() == 2:
      msg = 'Two-Component Kohn-Sham NYI'
      CErrMsg(workers['CQFileIO'],str(msg))


def parseRT(workers,secDict):


  if workers["CQSingleSlater"].nTCS() == 2:
    msg = 'Two-Component Real-Time Simulations NYI'
    CErrMsg(workers['CQFileIO'],str(msg))

  ref = secDict['QM']['REFERENCE']
  ref = ref.split()
  
  rtSettings = secDict['RT']

  # Set RT object based on SS reference
  if len(ref) == 1 or ref[0] == 'REAL':
    workers["CQRealTime"] = chronusQ.RealTime_double();
  else:
    workers["CQRealTime"] = chronusQ.RealTime_complex();

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
    'ENVELOPE' :workers['CQRealTime'].setEnvlp    ,
#    'ORTHO'    :workers['CQRealTime'].setOrthoTyp ,
#    'INIDEN'   :workers['CQRealTime'].setInitDen  ,
#    'UPROP'    :workers['CQRealTime'].setFormU    ,
    'IRSTRT'   :workers['CQRealTime'].setIRstrt   ,
#    'ELL_POL'  :workers['CQRealTime'].setEllPol     
    'INTSCHEME':workers['CQRealTime'].setIntScheme,
    'MMUTSCHEME':workers['CQRealTime'].setMMUTRstScheme,
    'MATEXP'   :workers['CQRealTime'].setPropMeth,
    'NPOLYMAX' :workers['CQRealTime'].setNPolyExpMax,
    'POLYEPS'  :workers['CQRealTime'].setPolyEps
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
    if (i in rtSettings) and (i not in ('EDFIELD')):
      optMap[i](rtSettings[i])
    elif (i in rtSettings) and (i in('EDFIELD')):
      optMap[i](rtSettings[i][0],rtSettings[i][1],rtSettings[i][2])


  if 'EDFIELD' not in rtSettings:
    optMap['EDFIELD'](0.0,0.0,0.0)
    
  if 'TARCSVS' in rtSettings:
    if not rtSettings['TARCSVS']:
      workers['CQRealTime'].doNotTarCSV()

  # Idiot Checks

  if 'ENVELOPE' in rtSettings:

    env = rtSettings['ENVELOPE']

    if env != "DELTA":
      msg = 'Non-Delta Field Envelopes NYI'
      CErrMsg(workers['CQFileIO'],str(msg))

    needsEDField = ['PLANEWAVE','CONSTANT','GAUSSIAN','DELTA',"STEP"]
    needsFreq    = ['LINRAMP','GAUSSIAN']
    needsSigma   = ['GAUSSIAN']
    needsTConst  = ['STEP']
    NYI          = ['SINSQ']


    if (env in needsEDField) and ('EDFIELD' not in rtSettings):
      msg = "Must specify EDFIELD with " + env + " envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env in needsFreq ) and ('FREQUENCY' not in rtSettings):
      msg = "Must specify FREQUENCY with " + env + "envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env in needsSigma ) and ('SIGMA' not in rtSettings):
      msg = "Must specify SIGMA with " + env + "envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env in needsTConst ) and ('TIME_ON' not in rtSettings or 'TIME_OFF' not in rtSettings):
      msg = "Must specify both TIME_ON and TIME_OFF with " + env + "envelope for RT"
      CErrMsg(workers['CQFileIO'],str(msg))

    if (env in NYI):
      msg = env + "NYI"
      CErrMsg(workers['CQFileIO'],str(msg))

#    if (env == envMap['PW']) and ('EDFIELD' not in rtSettings):
#      msg = "Must specify EDFIELD with PW envelope for RT"
#      CErrMsg(workers['CQFileIO'],str(msg))
# 
#    if (env == envMap['LINRAMP']) and ('EDFIELD' not in rtSettings):
#      msg = "Must specify EDFIELD with LINRAMP envelope for RT"
#      CErrMsg(workers['CQFileIO'],str(msg))
 
#    if (env == envMap['LINRAMP']) and ('FREQUENCY' not in rtSettings):
#      msg = "Must specify FREQUENCY with LINRAMP envelope for RT"
#      CErrMsg(workers['CQFileIO'],str(msg))

#    if (env == envMap['GAUSSIAN']) and ('EDFIELD' not in rtSettings):
#      msg = "Must specify EDFIELD with GAUSSIAN envelope for RT"
#      CErrMsg(workers['CQFileIO'],str(msg))

#    if (env == envMap['GAUSSIAN']) and ('SIGMA' not in rtSettings):
#      msg = "Must specify SIGMA with GAUSSIAN envelope for RT"
#      CErrMsg(workers['CQFileIO'],str(msg))

#    if (env == envMap['STEP']) and ('TIME_ON' not in rtSettings):
#      msg = "Must specify TIME_ON with STEP envelope for RT"
#      CErrMsg(workers['CQFileIO'],str(msg))
#
#    if (env == envMap['STEP']) and ('TIME_OFF' not in rtSettings):
#      msg = "Must specify TIME_OFF with STEP envelope for RT"
#      CErrMsg(workers['CQFileIO'],str(msg))
#
#    if (env == envMap['SINSQ']):
#      msg = "SINSQ envelope NYI"
#      CErrMsg(workers['CQFileIO'],str(msg))



#def parseSDR(workers,secDict):
#  jobSettings = {}
#  JOB = secDict['QM']['JOB']
#  try:
#    jobSettings = secDict[JOB]
#  except KeyError:
#    if JOB in ('STAB'):
#      pass
#    else: 
#      msg = "Must specify an options sections for " + JOB
#      CErrMsg(workers['CQFileIO'],str(msg))
#
#   
#  # Set SDR object based on SS reference
#  if workers['CQSingleSlater'] == workers['CQSingleSlaterDouble']:
#    workers['CQSDResponse']  = workers['CQSDResponseDouble']
#    workers['CQMOIntegrals'] = workers['CQMOIntegralsDouble']
#  elif workers['CQSingleSlater'] == workers['CQSingleSlaterComplex']:
#    msg = "Wave Function Response using a Complex Reference is NYI"
#    CErrMsg(workers['CQFileIO'],str(msg))
#    workers['CQSDResponse']  = workers['CQSDResponseComplex']
#    workers['CQMOIntegrals'] = workers['CQMOIntegralsComplex']
#  else:
#    msg = "Error in Reference Recieved by SDResponse"
#    CErrMsg(workers['CQFileIO'],str(msg))
#
#  if workers['CQSingleSlater'].nTCS() == 2 and workers['CQAOIntegrals'].integralAlgorithm != aointAlg['INCORE']:
#    msg = "Wave Function Response using a Two-Component Reference is\n"
#    msg = msg + "Only Implemented using INCORE integrals (QM.ints = INCORE)"
#    CErrMsg(workers['CQFileIO'],str(msg))
#    
#
#  try:
#    workers['CQSDResponse'].setNSek(jobSettings['NSTATES'])
#  except KeyError:
#    if JOB in ('STAB'):
#      workers['CQSDResponse'].setNSek(3)
#    else: 
#      msg = "Must specify number of desired roots for " + JOB
#      CErrMsg(workers['CQFileIO'],str(msg))
#     
#  workers['CQSDResponse'].setMeth(sdrMethodMap[str(JOB)])

def parseSCF(workers,scfSettings):
  optMap = {
    'DENTOL' :workers['CQSingleSlater'].setSCFDenTol,
    'ENETOL' :workers['CQSingleSlater'].setSCFEneTol,
    'MAXITER':workers['CQSingleSlater'].setSCFMaxIter,
    'DT'        :workers['CQSingleSlater'].setITPdt,
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
        if i in ('GUESS'):
          optMap[i](str(scfSettings[i]))
        else:
          optMap[i](scfSettings[i])
      else:
        optMap[i](scfSettings[i][0],scfSettings[i][1],scfSettings[i][2])
    except KeyError:
      continue

  if workers['CQMolecule'].nAtoms() == 1:
    optMap['GUESS'](str('CORE'))

  if 'GUESS' in scfSettings:
    if scfSettings['GUESS'] == 'READ':
      workers['CQFileIO'].doRestart = True

  try:
    if scfSettings['DIIS']:
      workers['CQSingleSlater'].doDIIS = True
    else:
      workers['CQSingleSlater'].doDIIS = False
  except KeyError:
    pass

  try:
    if scfSettings['ITP']:
      workers['CQSingleSlater'].doITP = True
    else:
      workers['CQSingleSlater'].doITP = False
  except KeyError:
    pass

