import os,sys
import libpythonapi as chronusQ

def communicate(workers):
  workers["CQAOIntegrals"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQFileIO"], 
    workers["CQControls"]
  )
  workers["CQSingleSlaterDouble"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQAOIntegrals"],
    workers["CQFileIO"], workers["CQControls"]
  )
  workers["CQRealTime"].communicate(workers["CQFileIO"],workers["CQControls"],
    workers["CQAOIntegrals"],workers["CQSingleSlaterDouble"])

def initialize(workers):
  # Set Up AOIntegrals Metadata
  workers["CQAOIntegrals"].initMeta()

  # Set up Wavefunction Information
  workers["CQSingleSlaterDouble"].initMeta()
  workers["CQSingleSlaterDouble"].genMethString()

  # RT
  workers["CQSingleSlaterDouble"].initMeta()
  
def alloc(workers):
  # Allocate Space for AO Integrals
  workers["CQAOIntegrals"].alloc()

  # Allocate Space for Wavefunction Information
  workers["CQSingleSlaterDouble"].alloc()

  # RT
  workers["CQSingleSlaterDouble"].alloc()

def runSCF(workers):
  communicate(workers)
  initialize(workers)
  alloc(workers)

  workers["CQMolecule"].printInfo(workers["CQFileIO"])
  workers["CQBasisSet"].printInfo();

  workers["CQSingleSlaterDouble"].formGuess()
  workers["CQSingleSlaterDouble"].formFock()
  workers["CQSingleSlaterDouble"].computeEnergy()
  workers["CQSingleSlaterDouble"].SCF()
  workers["CQSingleSlaterDouble"].computeMultipole()
  workers["CQSingleSlaterDouble"].printMultipole()
  workers["CQRealTime"].iniDensity()
  workers["CQRealTime"].doPropagation()
  
  
