import os,sys
import libpythonapi as chronusQ
from parseBasis import parseBasis
from meta.knownKeywords import requiredKeywords
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
      print 'Required keyword QM.' + str(i) + ' not found'
      exit(1)

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

  if ssSettings['job'] in ('rt'):
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
  ref = ref.split()
  if 'complex' in ref:
    workers["CQSingleSlater"] = workers["CQSingleSlaterComplex"]
  else:
    workers["CQSingleSlater"] = workers["CQSingleSlaterDouble"]

  print workers["CQSingleSlater"]
  print workers["CQSingleSlaterDouble"]
  if 'hf' in ref:
    if mult == 1:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
      workers["CQSingleSlater"].isClosedShell = True
    else:
      workers["CQSingleSlater"].setRef(chronusQ.Reference.UHF)
  elif 'rhf' in ref:
    workers["CQSingleSlater"].setRef(chronusQ.Reference.RHF)
    workers["CQSingleSlater"].isClosedShell = True
  elif 'uhf' in ref:
    workers["CQSingleSlater"].setRef(chronusQ.Reference.UHF)
  elif 'cuhf' in ref:
    workers["CQSingleSlater"].setRef(chronusQ.Reference.CUHF)
  elif 'ghf' in ref:
    workers["CQSingleSlater"].setRef(chronusQ.Reference.TCS)

  TCMethods = [chronusQ.Reference.TCS, chronusQ.Reference.GKS]
  if workers["CQSingleSlater"].Ref() in TCMethods:
    workers["CQSingleSlater"].setNTCS(2)

def parseRT(workers,settings):
# requiredKeywords = []
# optionalKeywords = [ 'maxstep' , 'timestep' , 'edfield' , 'time_on',
#                      'time_off', 'frequency', 'phase'   , 'sigma'  ,
#                      'envelope', 'ortho'    , 'iniden'  , 'uprop'  ]
# knownKeywords = requiredKeywords + optionalKeywords

# reqMap = {}
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



  for i in optMap:
    try:
      if i not in ('edfield'):
        optMap[i](settings[i])
      else:
        optMap[i](settings[i][0],settings[i][1],settings[i][2])
        
    except KeyError:
      continue


  

