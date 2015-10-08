import os,sys
#sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
#sys.path.append('/home/dbwy/git_repo/chronusq/src/python')
import libpythonapi as chronusQ

def parseBasis(workers,basis):
  nTCS = workers["CQSingleSlaterDouble"].nTCS()
  workers["CQBasisSet"].communicate(workers["CQFileIO"])
  
  workers["CQBasisSet"].findBasisFile(str(basis))
  workers["CQBasisSet"].parseGlobal()
  workers["CQBasisSet"].constructLocal(workers["CQMolecule"])
  workers["CQBasisSet"].makeMaps(nTCS,workers["CQMolecule"])
  workers["CQBasisSet"].renormShells()
  
