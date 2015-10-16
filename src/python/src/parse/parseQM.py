import os,sys
import libpythonapi as chronusQ
from parseBasis import parseBasis
#from standardJobs import *

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
def parseQM(workers,secDict): 
  print 'Parsing QM Information'
  ssSettings = secDict["qm"]

#
# Check for unknown keywords in the QM section
#
  knownKeywords = [ 'reference', 'basis', 'job' ]
  for i in ssSettings:
    if i not in knownKeywords:
      print "Keyword QM."+ str(i) +" not recognized"

#
# Try to set the reference for CQ::SingleSlater
#
  try:
    handleReference(workers,ssSettings['reference'])
  except KeyError:
    print 'No Reference keyword found'

#
# Try to set the basis for the QM Job
#
  try:
    parseBasis(workers,ssSettings['basis'])
  except KeyError:
    print 'No Basis Set keyword found'

  if ssSettings['job'] in ('RT'):
    parseRT(workers,secDict['rt']) 

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
def handleReference(workers,ref):
  mult = workers["CQMolecule"].multip()
  if ref in ('HF','hf'):
    if mult == 1:
      workers["CQSingleSlaterDouble"].setRef(chronusQ.Reference.RHF)
      workers["CQSingleSlaterDouble"].isClosedShell = True
    else:
      workers["CQSingleSlaterDouble"].setRef(chronusQ.Reference.UHF)
  elif ref in ('RHF','rhf'):
    workers["CQSingleSlaterDouble"].setRef(chronusQ.Reference.RHF)
    workers["CQSingleSlaterDouble"].isClosedShell = True
  elif ref in ('UHF','uhf'):
    workers["CQSingleSlaterDouble"].setRef(chronusQ.Reference.UHF)
  elif ref in ('CUHF','cuhf'):
    workers["CQSingleSlaterDouble"].setRef(chronusQ.Reference.CUHF)
  elif ref in ('GHF','ghf'):
    workers["CQSingleSlaterDouble"].setRef(chronusQ.Reference.TCS)

  TCMethods = [chronusQ.Reference.TCS, chronusQ.Reference.GKS]
  if workers["CQSingleSlaterDouble"].Ref() in TCMethods:
    workers["CQSingleSlaterDouble"].setNTCS(2)

def parseRT(workers,settings):
  requiredKeywords = []
  optionalKeywords = [ 'maxstep' , 'timestep' , 'edfield' , 'time_on',
                       'time_off', 'frequency', 'phase'   , 'sigma'  ,
                       'envelope', 'ortho'    , 'iniden'  , 'uprop'  ]
  knownKeywords = requiredKeywords + optionalKeywords

  reqMap = {}
  optMap = {   'maxstep':workers['CQRealTime'].setMaxSteps ,
              'timestep':workers['CQRealTime'].setStepSize ,
               'edfield':workers['CQRealTime'].setFieldAmp ,
               'time_on':workers['CQRealTime'].setTOn      ,
              'time_off':workers['CQRealTime'].setTOff     ,
             'frequency':workers['CQRealTime'].setFreq     ,
                 'phase':workers['CQRealTime'].setPhase    ,
                 'sigma':workers['CQRealTime'].setSigma    ,
              'envelope':workers['CQRealTime'].setEnvelope ,
                 'ortho':workers['CQRealTime'].setOrthoTyp ,
                'iniden':workers['CQRealTime'].setInitDen  ,
                 'uprop':workers['CQRealTime'].setFormU    }

  orthoMap = {    'lowdin':chronusQ.RealTime_ORTHO.Lowdin    ,
                'cholesky':chronusQ.RealTime_ORTHO.Cholesky  ,
               'canonical':chronusQ.RealTime_ORTHO.Canonical }

  formUMap = { 'eigendecomp':chronusQ.RealTime_FORM_U.EigenDecomp ,
                    'taylor':chronusQ.RealTime_FORM_U.Taylor      }

  envMap   = {       'pw':chronusQ.RealTime_ENVELOPE.Constant ,
                'linramp':chronusQ.RealTime_ENVELOPE.LinRamp  ,
               'gaussian':chronusQ.RealTime_ENVELOPE.Gaussian ,
                   'step':chronusQ.RealTime_ENVELOPE.Step     ,
                  'sinsq':chronusQ.RealTime_ENVELOPE.SinSq    }

  for i in settings:
    if i not in knownKeywords:
      print "Keyword RealTime."+ str(i) +" not recognized"
    elif i in ('maxstep','MaxStep','MAXSTEP','iniden','IniDen','INIDEN'):
      settings[i] = int(settings[i])
    elif i in ('ortho','Ortho','ORTHO'):
      settings[i] = orthoMap[i]
    elif i in ('envelope','Envelope','ENVELOPE'):
      settings[i] = envMap[i]
    elif i in ('uprop','Uprop','UProp','UPROP'):
      settings[i] = formUMap[i]
    elif i in ('edfield','EDfield','EDField'):
      settings[i] = settings[i].split()
      for j in range(len(settings[i])): settings[i][j] = float(settings[i][j])
    else:
      settings[i] = float(settings[i])


  for i in optMap:
    try:
      if i not in ('edfield','EDfield','EDField'):
        optMap[i](settings[i])
      else:
        optMap[i](settings[i][0],settings[i][1],settings[i][2])
        
    except KeyError:
      continue


  

