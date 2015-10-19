import os,sys
import libpythonapi as chronusQ

def communicate(workers):
  workers["CQAOIntegrals"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQFileIO"], 
    workers["CQControls"]
  )
  workers["CQSingleSlater"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQAOIntegrals"],
    workers["CQFileIO"], workers["CQControls"]
  )
  workers["CQRealTime"].communicate(workers["CQFileIO"],workers["CQControls"],
    workers["CQAOIntegrals"],workers["CQSingleSlater"])

def initialize(workers):
  # Set Up AOIntegrals Metadata
  workers["CQAOIntegrals"].initMeta()

  # Set up Wavefunction Information
  workers["CQSingleSlater"].initMeta()
  workers["CQSingleSlater"].genMethString()

  # RT
  
def alloc(workers):
  # Allocate Space for AO Integrals
  workers["CQAOIntegrals"].alloc()

  # Allocate Space for Wavefunction Information
  workers["CQSingleSlater"].alloc()

  # RT

def runSCF(workers):
  communicate(workers)
  initialize(workers)
  alloc(workers)

  workers["CQMolecule"].printInfo(workers["CQFileIO"])
  workers["CQBasisSet"].printInfo();

  workers["CQSingleSlater"].formGuess()
  workers["CQSingleSlater"].formFock()
  workers["CQSingleSlater"].computeEnergy()
  workers["CQSingleSlater"].SCF()
  workers["CQSingleSlater"].computeMultipole()
  workers["CQSingleSlater"].printMultipole()
  workers["CQRealTime"].initMeta()
  workers["CQRealTime"].alloc()
  workers["CQRealTime"].iniDensity()
  workers["CQRealTime"].doPropagation()
  
  
