import os,sys
#sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
#sys.path.append('/home/dbwy/git_repo/chronusq/src/python')
import libpythonapi as chronusQ
from parseBasis import parseBasis
#from standardJobs import *

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

#stdJobs = [ 'SCF', 'CIS', 'TDA', 'RPA', 'RT', 'STAB' ]
#stdJobFunctions = {
#  'SCF':runSCF
#}
def parseQM(workers,settings): 
  print 'Parsing QM Information'
  try:
    handleReference(workers,settings['reference'])
  except KeyError:
    print 'No Reference keyword found'

  try:
    parseBasis(workers,settings['basis'])
  except KeyError:
    print 'No Basis Set keyword found'

  # Space filler to pasify error 
  workers["CQSingleSlaterDouble"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQAOIntegrals"],
    workers["CQFileIO"], workers["CQControls"]
  )
  workers["CQSingleSlaterDouble"].initMeta()
  workers["CQSingleSlaterDouble"].genMethString()
  workers["CQSingleSlaterDouble"].alloc()
  
