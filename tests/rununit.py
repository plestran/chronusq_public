import sys,os
import getopt
from refval import *
from chronusq import *

##############################
#      Class Definitions     #
##############################

class UnitTest:
    def __init__(self,infile,testClass,desc):
        self.infile    = infile
        self.testClass = testClass
        self.desc      = desc

##############################
#     Routine Definitions    #
##############################
#--------------------------------------------------------------------
def appendSummary(fname,errors,jobtype):
# find where to insert the new results
  outf = open(summaryf,'r')
  contents = outf.readlines()
  find = False
  index  = 0
  for line in contents:
    index += 1
    if jobtype in line: find = True
    if line in ['\n', '\r\n'] and find :
      insert = index
      find = False
  outf.close()
#
# insert the new results
  if jobtype == 'SCF':
    string = '%-28s %.4E %6s %.4E %10s %.4E %8s %.4E' % (fname,errors[0],"",errors[1],"",errors[2],"",errors[3])
    if errors[0] < 5E-10 and errors[1] < 5E-5 and errors[2] < 5E-5 and errors[3] < 5E-5: string = string+"  YES\n"
    else: string = string+"  ** NO **\n"
  elif jobtype == 'RESP':
    nstates = len(ref[fname[:8]].w)
    string = '%-28s %.4E %5s %.4E %6s %3d' % (fname,errors[0],"",errors[1],"",nstates)
    if errors[0] < 1E-7 and errors[1] < 1E-7: string = string+"  YES\n"
    else: string = string+"  ** NO **\n"
  elif jobtype == 'RT':
    string = '%-32s %.4E %10s %.4E' % (fname,errors[0],"",errors[1])
    if errors[0] < 1E-10 and errors[1] < 1E-6: string = string+"  YES\n"
    else: string = string+"  ** NO **\n"
  else:
    print "Job type not recognized: "+jobtype
    sys.exit()
  contents.insert(insert-1,string)
  outf = open(summaryf, "w")
  for line in contents: outf.write(line)
  outf.close()
#
# kill the run if a job fails and doKill is set to True
  if "** NO **" in string: 
    if doPrint: print "\n** WARNING: Job "+fname+" did not pass**\n"
    if doKill: sys.exit(1)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def findFile(name,path):
# Checks whether file "name" exists in "path"
  found = False
  for root, dirs, files in os.walk(path):
    if name in files:
      found = True
      break
  return found
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def genTable():
# Reads test.index to see which tests to run
  start = False
  f = open("test.index",'r')
  table = []
  for line in f:
    strx = line.split()
    if ".inp" in line and '#' not in line:

#     check if the type falls within the set we want to run
      runTest = False
      for sub in subTests:
        if sub in line.lower():
          runTest = True
          continue

#     run all the tests not commented out
      if testType == 'all': runTest = True

#     checks for exclusive options
      if testSize == 'large' and "small" in line.lower(): runTest = False
      if testSize == 'small' and "large" in line.lower(): runTest = False
      if testPar  == 'off' and "openmp" in line.lower(): runTest = False
      if testPar  == 'on'  and "serial" in line.lower(): runTest = False
      if testInts == 'incore' and "direct" in line.lower(): runTest = False
      if testInts == 'direct' and "incore" in line.lower(): runTest = False
      if "dfield" in subTests and "dfield" not in line.lower(): runTest = False
      if testComp == 'no' and "complex" in line.lower(): runTest = False
      if testComp == 'yes'and "real"    in line.lower(): runTest = False
      if testBasis != 'all':
        for basis in subBasis:
          if basis not in line.lower(): 
            runTest = False
            continue

#     add the tests to the list
      if runTest:
        desc = ' '.join(strx[2:])
        entry = UnitTest(strx[0],strx[1],desc)
        table.append(entry)

  f.close()
  return table
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def initSummary():
# Opens the output file and adds headers for the different tests
  global summaryf
  summaryf = "summary.txt"
  exists = os.path.exists(summaryf)
  if exists: os.remove(summaryf)
  outf = open(summaryf,'w')

  SCFheader = "SCF Test Job                  |dEnergy|    max(|dDipole|)    max(|dQuadrupole|)    max(|dOctupole|)  Passed\n--------------------------  -----------  ----------------  --------------------  ------------------  --------\n\n"
  RESPheader = "RESP Test Job                 max(|df|)    max(|domega|)    NStates  Passed\n--------------------------  -----------  ---------------  ---------  --------\n\n"
  RTheader = "RT Test Job                   |dLastEnergy|    max(|dLastDipole|)  Passed\n--------------------------  ---------------  --------------------  --------\n\n"
  outf.write(SCFheader+RESPheader+RTheader)
  outf.close()
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def runUnit(doPrint):
# Runs the unit tests
  global errors
  refvalues()
  tests = [None for x in xrange(1000)] 
  testtable = genTable()

  k = 0
  for i in testtable:
    errors  = []
    if findFile(i.infile,"."):
#
#     run chronus
      if doPrint: print "running job: "+i.infile
      tests[k] = runCQ(i.infile,'')
#
#     test SCF values
      if 'SCF' in ref[i.infile[:8]].typ:
        testSCF(ref[i.infile[:8]],tests[k])
        appendSummary(i.infile,errors,'SCF')
#
#     test RESP values
      elif 'RESP' in ref[i.infile[:8]].typ:
        testRESP(ref[i.infile[:8]],tests[k])
        appendSummary(i.infile,errors,'RESP')
#
#     test RT values 
      elif 'RT' in ref[i.infile[:8]].typ:
        testRT(ref[i.infile[:8]],tests[k])
        appendSummary(i.infile,errors,'RT')

      else:
        print "Job type not recognized: "+ref[i.infile[:8]].typ
        sys.exit()

    k += 1
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def testRESP(ref,tests):
  auToeV = 27.2113961

# test excitation energies
  maxerr = 0.0
  for i in range(0,len(ref.w)):
    abserr = abs(ref.w[i]   - tests.excEne[i]*auToeV)
    if abserr > maxerr:
      maxerr = abserr
  errors.append(abserr)

#  test oscillator strengths
  maxerr = 0.0
  for i in range(0,len(ref.w)):
    abserr = abs(ref.osc[i] - tests.oscStr[i])
    if abserr > maxerr:
      maxerr = abserr
  errors.append(abserr)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def testRT(ref,tests):
  auToD   = 0.3934303070
  auToAng = 0.5291772083

# test RT energy
  abserr = abs(ref.eng - tests.lastEnergy)
  errors.append(abserr)

# test time evolving molecular dipoles
  maxerr = 0.0
  for i in range(4):
    abserr = abs(abs(ref.dip[i]) - abs(tests.lastDipole[i]))
    if abserr > maxerr:
      maxerr = abserr
  errors.append(maxerr)
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def testSCF(ref,tests):
  auToD   = 0.3934303070
  auToAng = 0.5291772083

# test SCF energy
  abserr = abs(ref.eng - tests.E)
  errors.append(abserr)

# test molecular dipoles
  maxerr = 0.0
  for i in range(3):
    abserr = abs(ref.dip[i] - tests.dipole[i]/auToD)
    if abserr > maxerr:
      maxerr = abserr
  errors.append(maxerr)

# test molecular quadrupole
  k = 0
  maxerr = 0.0
  for i in range(3):
    for j in range(i,3):
      abserr = abs(ref.quad[k] - tests.quadrupole[i][j]*auToAng/auToD)
      if abserr > maxerr:
        maxerr = abserr
      k += 1
  errors.append(maxerr)

# test molecular octupole
  l = 0
  maxerr = 0.0
  for i in range(3):
    for j in range(i,3):
      for k in range(j,3):
        abserr = abs(ref.octu[l] - tests.octupole[i][j][k]*auToAng*auToAng/auToD)
        if abserr > maxerr:
          maxerr = abserr
        l += 1
  errors.append(maxerr)
#--------------------------------------------------------------------

##############################
#        Main Program        #
##############################

if __name__ in "__main__":
  chronusQ.initCQ(len(sys.argv),sys.argv)

  if chronusQ.getSize() > 1:
    if chronusQ.getRank() == 0:
      print "Cannot run tests jobs with more than one MPI process"
    chronusQ.finalizeCQ()
    sys.exit(1)
# parse user options
  msg = """
  python runtests.py [-o --option=]

  Options:
    -h, --help        Print usage instructions
    -s, --silent      Disable Print
    -k, --kill        Stop testing if a job fails
    --type=           Determines types of tests to run. Multiple options
                      can be specified by separating with a comma.
                      3 classes of tests  = [SCF,RESP,RT] 
                      Specify References  = [RHF,UHF,CUHF,GHF]
                                            [RKS,UKS,SLATER,LSDA,SVWN5]
                      Reference and Type  = [(R|U|CU)HF-SCF,HF-CIS,HF-RPA]
                                            [(R|U)KS-SCF,SCF-LSDA] 
                      Dipole Field        = [DField]
    --integrals=      Integral evaluation = [incore] or [direct]
    --parallel=       Whether to run parallel jobs = [on] or [off]
    --size=           Size of jobs to run = [small] or [large] or [both]
                      [small] is the default
    --complex=        Complex Jobs = [yes] or [no] or [both]
                      [both] is the default
    --basis=          Only run tests for this basis set
                      [STO-3G,6-31G,cc-pVDZ,def2-SVPD]
"""
  doPrint  = True
  doKill   = False
  testType = 'all'
  testInts = ''
  testPar  = ''
  testSize = 'small'
  testComp = ''
  testBasis = 'all'
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hks",["help","silent","kill","type=","integrals=","parallel=","size=","complex=","basis="])
  except getopt.GetoptError:
    print msg
    sys.exit(2)
  for opt, arg in opts:
    if opt in ('-h',"--help"):
      print msg
      sys.exit()
    elif opt in ('-s',"--silent"):
      doPrint = False
    elif opt in ('-k',"--kill"):
      doKill = True
    elif opt in ("--type"):
      testType = arg.lower()
    elif opt in ("--integrals"):
      testInts = arg.lower()
    elif opt in ("--parallel"):
      testPar  = arg.lower()
    elif opt in ("--size"):
      testSize = arg.lower()
    elif opt in ("--complex"):
      testComp = arg.lower()
    elif opt in ("--basis"):
      testBasis = arg.lower()

  subTests = testType.split(',')
  subBasis = testBasis.split(',')

# run unit tests
  initSummary()
  runUnit(doPrint)    


  chronusQ.finalizeCQ()
#--------------------------------------------------------------------

