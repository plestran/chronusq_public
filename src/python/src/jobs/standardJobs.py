#
# The Chronus Quantum (ChronusQ) software package is high-performace 
# computational chemistry software with a strong emphasis on explicitly 
# time-dependent and post-SCF quantum mechanical methods.
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

def communicate(workers):
  workers["CQAOIntegrals"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQFileIO"],
    workers["CQMemManager"]
  )
  workers["CQSingleSlater"].communicate(
    workers["CQMolecule"], workers["CQBasisSet"], workers["CQAOIntegrals"],
    workers["CQFileIO"], workers["CQMemManager"]
  )
  try:
    workers["CQRealTime"].communicate(workers["CQSingleSlater"])
  except KeyError:
    pass

#  try:
#    workers["CQSDResponse"].communicate(workers["CQMolecule"],
#      workers["CQBasisSet"],workers["CQSingleSlater"],workers["CQMOIntegrals"],
#      workers["CQFileIO"],workers["CQControls"])
#  except KeyError:
#    pass
#
#  try:
#    workers["CQMOIntegrals"].communicate(workers["CQMolecule"],
#      workers["CQBasisSet"],workers["CQFileIO"],
#      workers["CQAOIntegrals"],workers["CQSingleSlater"])
#  except KeyError:
#    pass

def runCQJob(job,workers,meta):
  workers["CQFileIO"].iniH5Files()

  # Make the classes know about eachother
  communicate(workers)

  # Print some information pertaining to the job
  workers["CQMolecule"].printInfo(workers["CQFileIO"])
  workers["CQBasisSet"].printInfo();

  # Set Up AOIntegrals Metadata
  workers["CQAOIntegrals"].initMeta()

  # Set up Wavefunction Metadata
  workers["CQSingleSlater"].initMeta()

  # Allocate Space for AO Integrals
  workers["CQAOIntegrals"].alloc()
 
  # Allocate Space for Wavefunction Information
  workers["CQSingleSlater"].alloc()

  # Always Form Guess
  workers["CQSingleSlater"].formGuess()

#  if job == 'SCF':
  workers["CQSingleSlater"].SCF() 

  workers["CQSingleSlater"].computeProperties()
  
  workers["CQSingleSlater"].printProperties()

  if job == 'RT':
    workers["CQRealTime"].alloc()
    workers["CQRealTime"].doPropagation()


  if job == 'SCF':
    meta.E          = workers["CQSingleSlater"].totalEnergy()
    meta.dipole     = workers["CQSingleSlater"].dipole()
    meta.quadrupole = workers["CQSingleSlater"].quadrupole()
    meta.octupole   = workers["CQSingleSlater"].octupole()
  elif job == 'RT':
    pass
    

def runSCF(workers):
  # Make the classes know about eachother
  communicate(workers)

  # Print some information pertaining to the job
  # FIXME: These two are general and should always be printed
  #        regardless of job
  if chronusQ.getRank() == 0:
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

  #workers["CQSingleSlater"].formFock()
  #workers["CQSingleSlater"].computeEnergy()

  workers["CQSingleSlater"].SCF()

  workers["CQSingleSlater"].computeProperties()
  if chronusQ.getRank() == 0:
    workers["CQSingleSlater"].printProperties()


  #meta.E          = workers["CQSingleSlater"].totalEnergy
  #meta.scfIters   = workers["CQSingleSlater"].nSCFIter
  #meta.dipole     = workers["CQSingleSlater"].dipole()
  #meta.quadrupole = workers["CQSingleSlater"].quadrupole()
  #meta.octupole   = workers["CQSingleSlater"].octupole()

#def runRT(workers,meta):
#  runSCF(workers,meta)
#  workers["CQRealTime"].initMeta()
#
#  workers["CQRealTime"].alloc()
#  if chronusQ.getRank() == 0:
#    workers["CQRealTime"].iniDensity()
#  workers["CQRealTime"].doPropagation()
#
#  meta.lastDipole = workers['CQRealTime'].lastDipole()
#  meta.lastEnergy = workers['CQRealTime'].lastEnergy()
#
#def runSDR(workers,meta):
#
#  runSCF(workers,meta)
#  workers["CQMOIntegrals"].initMeta()
#  workers["CQSDResponse"].initMeta()
#  workers["CQSDResponse"].setPPRPA(1)
#  if chronusQ.getRank() == 0:
#    workers["CQSDResponse"].initMeth()
#    workers["CQSDResponse"].alloc()
#    workers["CQSDResponse"].IterativeRESP()
#
#  meta.davIters = workers["CQSDResponse"].nIter
#  meta.excEne   = workers["CQSDResponse"].excitationEnergies()
#  meta.oscStr   = workers["CQSDResponse"].oscStrengths()
