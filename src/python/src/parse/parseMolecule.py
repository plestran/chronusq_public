import os,sys
#sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
#sys.path.append('/home/dbwy/git_repo/chronusq/src/python')
import libpythonapi as chronusQ

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
  print 'Parsing Molecular Information'

#
# Check for unknown keywords in the Molecules section
#
# knownKeywords = [ 'charge', 'mult', 'geom' ]
# for i in settings:
#   if i not in knownKeywords:
#     print "Keyword Molecule."+ str(i) +" not recognized"
#
# Grab charge and multiplicity
#
# FIXME: Need a check if these keywords aren't there
#        or an ugly python error occurs
#
  charge = settings['charge']
  mult   = settings['mult'  ]
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
  geomStr = settings['geom']
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
#   CQ::Molecule::index    -   a map from atom number to index in CQ::atoms array
#   CQ::Molecule::nTotalE  -   total Number of electrons 
#                                (modulated by CQ::Molecule::charge)
#
  nTotalE = 0
  for i in range(len(geomStr)):
    line = geomStr[i]
    lineSplit = line.split()
    indx = -1
    if len(lineSplit) == 5:
      indx = chronusQ.HashAtom(str(lineSplit[0]),int(lineSplit[1]))
    elif len(lineSplit) == 4:
      indx = chronusQ.HashAtom(str(lineSplit[0]),0)
    else:
      print 'Input Error: Invalid Geometry Specification'

    if indx != -1:
      workers["CQMolecule"].setIndex(i,indx)
    else:
      print 'Input Error: Invalid Atomic Symbol or Mass Number'

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
