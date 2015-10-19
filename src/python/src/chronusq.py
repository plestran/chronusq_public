import os,sys
sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
import libpythonapi as chronusQ
from parse.meta.standardJobs import *
from parse.parseInput import parseInput
from parse.parseMolecule import parseMolecule
from parse.parseQM import parseQM

fname = sys.argv[1]

mol        = chronusQ.Molecule()
basisSet   = chronusQ.BasisSet()
DFbasisSet = chronusQ.BasisSet()
controls   = chronusQ.Controls()
aoints     = chronusQ.AOIntegrals()

hf_double  = chronusQ.SingleSlater_double()
rt_double  = chronusQ.RealTime_double()

out        = chronusQ.FileIO(fname)

workers = {"CQMolecule":mol,
           "CQBasisSet":basisSet,
           "CQDFBasisSet":DFbasisSet,
           "CQControls":controls,
           "CQAOIntegrals":aoints,
           "CQSingleSlaterDouble":hf_double,
           "CQRealTimeDouble":rt_double,
           "CQFileIO":out}


controls.iniControls()

secDict = parseInput(workers,fname+".inp")

parseMolecule(workers,secDict["MOLECULE"])
parseQM(workers,secDict)

runSCF(workers)

#controls.printSettings()
#mol.printInfo(out,controls)
#basisSet.printInfo()
