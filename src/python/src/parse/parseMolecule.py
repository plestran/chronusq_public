#
# The Chronus Quantum (ChronusQ) software package is high-performace 
# computational chemistry software with a strong emphasis on explicitly 
# time-dependent and post-SCF quantum mechanical methods.
# 
# Copyright (C) 2014-2015 Li Research Group (University of Washington)
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
from meta.knownKeywords import requiredKeywords

#
#  Parse the Molecule section of the input file and populate
#  the CQ::Molecule object
#
#  ** Note that this includes allocation of the Molecules object **
#
#  Input:
#    workers      -     the workers array of CQ objects
#    settings     -     the settings parsed by configparser for 
#                       the Molecules section
#
def parseMolecule(workers,settings): 
#  print 'Parsing Molecular Information'

#
#  Check that all of the required keywords for Molecule
#  object are found
#
  for i in requiredKeywords['MOLECULE']:
    if i not in settings:
      msg = 'Required keyword Molecule.' + str(i) + ' not found'
      CErrMsg(workers['CQFileIO'],str(msg))
#
# Grab charge and multiplicity
#
  charge = settings['CHARGE']
  mult   = settings['MULT'  ]
#
# Populate charge and multiplicity of the
# CQ Molecules object
#
  workers["CQMolecule"].setCharge(charge)
  workers["CQMolecule"].setMultip(mult  )
#
# Read the geometry from input file
#
  readGeom(workers,settings)
#
# 1) Convert coordinates from Angstrom to Bohr
# 2) Compute Nuclear Repulsion Energy
# 3) Compute Intranuclear distance tensor
# 4) Compute Moment of Inertia
#
  workers["CQMolecule"].convBohr()      # 1
  workers["CQMolecule"].computeNucRep() # 2 
  workers["CQMolecule"].computeRij()    # 3
  workers["CQMolecule"].computeI()      # 4

  try:
    workers["CQMolecule"].setPrintLevel(settings['PRINT'])
  except KeyError:
    pass


#
# Read the molecular geometry from input file
#
# ** Note that the entire workers array is passed **
# ** so that the Molecule object can be allocated **
# ** in reference to the CQ::FileIO object        **
#
def readGeom(workers,settings):
#
# Grab the string that contains the molecular geometry
# and split up by line ends
#
  geomStr = settings['GEOM']
  geomStr = geomStr.split('\n')
#
# Determine how many atoms we have and try to remove
# blank lines as parsed by ConfigParser
#
# FIXME: Need a more robust way to remove blank lines
#
  nAtoms = 0
  for i in geomStr:
    if i not in ('', ' '):
      nAtoms += 1
  oldLen = len(geomStr)
  for i in range(oldLen - nAtoms):
    geomStr.remove('')
#
#  Populate the CQ::Molecule object with nAtoms and
#  allocate the nessacary memory
#
  workers["CQMolecule"].setNAtoms(nAtoms)
  workers["CQMolecule"].alloc(workers["CQFileIO"])
  
#
# Parse the geometry string and populate the Cartesian
# position tensor of the CQ::Molecule object
#
# Two methods of geometry specification are valid:
#   
#  1) <ATM> X Y Z
#  2) <ATM> <ISO> X Y Z
#
#  where <ATM> is the atomic symbol and <ISO> is some
#  isotope identifier (** not yet fully implemented **)
#
# This routine also populates:
#   CQ::Molecule::index    -  a map from atom number to index in CQ::atoms array
#   CQ::Molecule::nTotalE  -  total Number of electrons 
#                               (modulated by CQ::Molecule::charge)
#
#  nTotalE = 0
  nTotalE = - settings['CHARGE']
  for i in range(len(geomStr)):
    line = geomStr[i]
    lineSplit = line.split()
    indx = -1
    if len(lineSplit) != 5 and len(lineSplit) != 4:
      msg = 'Input Error: Invalid Geometry Specification'
      CErrMsg(workers['CQFileIO'],str(msg))
    if any(char.isdigit() for char in lineSplit[0]):
      if len(lineSplit) == 5:
        indx = chronusQ.HashZ(int(lineSplit[0]),int(lineSplit[1]))
      elif len(lineSplit) == 4:
        indx = chronusQ.HashZ(int(lineSplit[0]),0)
    else:
      if len(lineSplit) == 5:
        indx = chronusQ.HashAtom(str(lineSplit[0]),int(lineSplit[1]))
      elif len(lineSplit) == 4:
        indx = chronusQ.HashAtom(str(lineSplit[0]),0)

    if indx != -1:
      workers["CQMolecule"].setIndex(i,indx)
    else:
      msg = 'Input Error: Invalid Atomic Symbol or Mass Number'
      CErrMsg(workers['CQFileIO'],str(msg))

    nTotalE += chronusQ.getAtomicNumber(indx)

    if len(lineSplit) == 5:
      x = float(lineSplit[2])
      y = float(lineSplit[3])
      z = float(lineSplit[4])
    elif len(lineSplit) == 4:
      x = float(lineSplit[1])
      y = float(lineSplit[2])
      z = float(lineSplit[3])

    workers["CQMolecule"].setCart(i,x,y,z)
  workers["CQMolecule"].setNTotalE(nTotalE)
