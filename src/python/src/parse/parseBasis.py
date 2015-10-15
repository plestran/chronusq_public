import os,sys
import libpythonapi as chronusQ

#
# Parse Basis Set file
#
# ** Note that this allocates the CQ::BasisSet Object **
#
def parseBasis(workers,basis):
  nTCS = workers["CQSingleSlaterDouble"].nTCS()

  # Make CQ::BasisSet aware of CQ::FileIO
  workers["CQBasisSet"].communicate(workers["CQFileIO"])
  
#
# 1) Try to find the basis set file in standard locations
# 2) Allocate and populate a global basis set definition
# 3) Construct a local subset (viz. CQ::Molecule) of the basis set
# 4) Construct the varios maps incolving the basis set
# 5) Renormalize the LibInt2::Shell's (this is outdated in newer Libint)
#
  workers["CQBasisSet"].findBasisFile(str(basis))             # 1
  workers["CQBasisSet"].parseGlobal()                         # 2
  workers["CQBasisSet"].constructLocal(workers["CQMolecule"]) # 3
  workers["CQBasisSet"].makeMaps(nTCS,workers["CQMolecule"])  # 4
  workers["CQBasisSet"].renormShells()                        # 5
  
