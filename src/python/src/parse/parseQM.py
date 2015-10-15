import os,sys
import libpythonapi as chronusQ
from parseBasis import parseBasis
#from standardJobs import *

#
#  Parse the QM section of the input file and populate
#  the CQ::SingleSlater<T> and CQ::BasisSet objects
#
#  ** Note that this does DOES NOT allocate the SingleSlater **
#  ** object, but does allocate the BasisSet object          **
#
#  Input:
#    workers      -     the workers array of CQ objects
#    settings     -     the settings parsed by configparser for 
#                       the QM section
#
def parseQM(workers,settings): 
  print 'Parsing QM Information'

#
# Check for unknown keywords in the QM section
#
  knownKeywords = [ 'reference', 'basis' ]
  for i in settings:
    if i not in knownKeywords:
      print "Keyword Molecule."+ str(i) +" not recognized"

#
# Try to set the reference for CQ::SingleSlater
#
  try:
    handleReference(workers,settings['reference'])
  except KeyError:
    print 'No Reference keyword found'

#
# Try to set the basis for the QM Job
#
  try:
    parseBasis(workers,settings['basis'])
  except KeyError:
    print 'No Basis Set keyword found'

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

