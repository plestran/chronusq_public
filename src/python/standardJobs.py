import os,sys
import libpythonapi as chronusQ

def communicate(workers):
  workers["CQAOIntegrals"].communicate(
    workers["CQMolecule"] workers["CQBasisSet"], workers["CQControls"], 
    workers["CQDFBasisSet"]
  )
  workers["CQSingleSlaterDouble"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQAOIntegrals"],
    workers["CQFileIO"], workers["CQControls"]
  )

def initialize(workers):
  # Set up Wavefunction Information
  workers["CQSingleSlaterDouble"].initMeta()
  workers["CQSingleSlaterDouble"].genMethString()
  
def alloc(workers):
  # Allocate Space for AO Integrals
  workers["CQAOIntegrals"].alloc()

  # Allocate Space for Wavefunction Information
  workers["CQSingleSlaterDouble"].alloc()

def runSCF(workers):
  communicate(workers)
  initialize(workers)
  alloc(workers)

  workers["CQMolecule"].printInfo(workers["CQFileIO"])
