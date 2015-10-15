import os,sys
sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
#sys.path.append('/home/dbwy/git_repo/chronusq/src/python')
import parse.parseInput as PI
import libpythonapi as chronusQ
from standardJobs import *

fname = sys.argv[1]

mol        = chronusQ.Molecule()
basisSet   = chronusQ.BasisSet()
DFbasisSet = chronusQ.BasisSet()
controls   = chronusQ.Controls()
aoints     = chronusQ.AOIntegrals()
hf         = chronusQ.SingleSlater_double()
out        = chronusQ.FileIO(fname)

workers = {"CQMolecule":mol,
           "CQBasisSet":basisSet,
           "CQDFBasisSet":DFbasisSet,
           "CQControls":controls,
           "CQAOIntegrals":aoints,
           "CQSingleSlaterDouble":hf,
           "CQFileIO":out}


controls.iniControls()

PI.parseInput(workers,fname+".inp")
runSCF(workers)

#controls.printSettings()
#mol.printInfo(out,controls)
#basisSet.printInfo()
