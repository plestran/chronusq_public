import os,sys
sys.path.append('/home/jjgoings/chronusq/build_gcc_libint_openmp/src/python')
import libpythonapi as chronusQ
from jobs.jobMap import *
from parse.parseInput import parseInput
from parse.parseMolecule import parseMolecule
from parse.parseQM import parseQM
from parse.parseMisc import parseMisc

try:
    fname = sys.argv[1]
except IndexError:
    sys.exit("Please specify the input filename")

mol        = chronusQ.Molecule()
basisSet   = chronusQ.BasisSet()
DFbasisSet = chronusQ.BasisSet()
controls   = chronusQ.Controls()
aoints     = chronusQ.AOIntegrals()

hf_double     = chronusQ.SingleSlater_double()
rt_double     = chronusQ.RealTime_double()
sdr_double    = chronusQ.SDResponse_double()
moints_double = chronusQ.MOIntegrals_double()

out        = chronusQ.FileIO(fname)

workers = {
    "CQMolecule"          :mol,
    "CQBasisSet"          :basisSet,
    "CQDFBasisSet"        :DFbasisSet,
    "CQControls"          :controls,
    "CQAOIntegrals"       :aoints,
    "CQSingleSlaterDouble":hf_double,
    "CQRealTimeDouble"    :rt_double,
    "CQSDResponseDouble"  :sdr_double,
    "CQMOIntegrals"       :moints_double,
    "CQFileIO"            :out
}

header = """
   ______ __                                      ____ 
  / ____// /_   _____ ____   ____   __  __ _____ / __ \ 
 / /    / __ \ / ___// __ \ / __ \ / / / // ___// / / /
/ /___ / / / // /   / /_/ // / / // /_/ /(__  )/ /_/ / 
\____//_/ /_//_/    \____//_/ /_/ \__,_//____/ \___\_\ 
                                                       
"""

out.write(header)


chronusQ.initCQ()
controls.iniControls()

secDict = parseInput(workers,fname+".inp")
print secDict

parseMisc(workers,secDict["MISC"])
parseMolecule(workers,secDict["MOLECULE"])
jobStr = parseQM(workers,secDict)

jobMap[jobStr](workers)

chronusQ.finalizeCQ()

#controls.printSettings()
#mol.printInfo(out,controls)
#basisSet.printInfo()
