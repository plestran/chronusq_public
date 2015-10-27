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
  try:
    workers["CQRealTime"].communicate(workers["CQFileIO"],workers["CQControls"],
      workers["CQAOIntegrals"],workers["CQSingleSlater"])
  except KeyError:
    pass

  try:
    workers["CQSDResponse"].communicate(workers["CQMolecule"],
      workers["CQBasisSet"],workers["CQSingleSlater"],workers["CQMOIntegrals"],
      workers["CQFileIO"],workers["CQControls"])
  except KeyError:
    pass

  try:
    workers["CQMOIntegrals"].communicate(workers["CQMolecule"],
      workers["CQBasisSet"],workers["CQFileIO"],workers["CQControls"],
      workers["CQAOIntegrals"],workers["CQSingleSlater"])
  except KeyError:
    pass



def runSCF(workers,meta):
  # Make the classes know about eachother
  communicate(workers)

  # Print some information pertaining to the job
  # FIXME: These two are general and should always be printed
  #        regardless of job
  workers["CQMolecule"].printInfo(workers["CQFileIO"])
  workers["CQBasisSet"].printInfo();

  
  # Set Up AOIntegrals Metadata
  workers["CQAOIntegrals"].initMeta()

  # Set up Wavefunction Metadata
  workers["CQSingleSlater"].initMeta()
  workers["CQSingleSlater"].genMethString()

  # Allocate Space for AO Integrals
  workers["CQAOIntegrals"].alloc()

  # Allocate Space for Wavefunction Information
  workers["CQSingleSlater"].alloc()

  workers["CQSingleSlater"].formGuess()
  workers["CQSingleSlater"].formFock()
  workers["CQSingleSlater"].computeEnergy()
  workers["CQSingleSlater"].SCF()
  workers["CQSingleSlater"].computeMultipole()
  workers["CQSingleSlater"].printMultipole()

def runRT(workers,meta):
  runSCF(workers)
  workers["CQRealTime"].initMeta()
  workers["CQRealTime"].alloc()
  workers["CQRealTime"].iniDensity()
  workers["CQRealTime"].doPropagation()

def runSDR(workers,meta):
  runSCF(workers)
  workers["CQMOIntegrals"].initMeta()
  workers["CQSDResponse"].initMeta()
  workers["CQSDResponse"].initMeth()
  workers["CQSDResponse"].alloc()
  workers["CQSDResponse"].IterativeRPA()
