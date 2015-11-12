import sys,os
import getopt
from refval import *
from chronusq import *
from tabulate import tabulate

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
def findFile(name,path):
#  Checks whether file "name" exists in "path"
  found = False
  for root, dirs, files in os.walk(path):
    if name in files:
      found = True
      break
  return found
#--------------------------------------------------------------------

#--------------------------------------------------------------------
def genSummary(testtable,summary):
#  Prints the results in summary.txt
  outf = open('summary.txt','w')

# SCF output    
  headers = ["SCF Test Job","|dEnergy|","max(|dDipole|)","max(|dQuadrupole|)","max(|dOctupole|)","Passed"]
  sumrytable = []
  j = 0
  for i in testtable:
    if 'SCF' in ref[i.infile[:8]].typ:
      entry = []
      entry.append(testtable[j].infile)
      for k in range(4):
        entry.append(summary[j][k])
      if summary[j][0] < 1E-10 and summary[j][1] < 1E-8 and summary[j][2] < 1E-8 and summary[j][3] < 1E-8:
        entry.append('YES')
      else:
        entry.append('** NO **')
      sumrytable.append(entry)
    j += 1
  outf.write(tabulate(sumrytable,headers,tablefmt="simple",floatfmt=".4E"))

# RESP output
  headers = ["RESP Test Job","max(|f|)","max(|omega|)","NStates","Passed"]
  sumrytable = []
  j = 0
  for i in testtable:
    if 'RESP' in ref[i.infile[:8]].typ:
      entry = []
      entry.append(testtable[j].infile)
      entry.append(summary[j][0])
      entry.append(summary[j][1])
      entry.append(len(ref[i.infile[:8]].w))
      if summary[j][0] < 1E-7 and summary[j][1] < 1E-7:
        entry.append('YES')
      else:
        entry.append('** NO **')
      sumrytable.append(entry)
    j += 1
  outf.write("\n\n")
  outf.write(tabulate(sumrytable,headers,tablefmt="simple",floatfmt=".4E"))

# RT output    
  headers = ["RT Test Job","|dLastEnergy|","max(|dLastDipole|)","Passed"]
  sumrytable = []
  j = 0
  for i in testtable:
    if 'RT' in ref[i.infile[:8]].typ:
      entry = []
      entry.append(testtable[j].infile)
      entry.append(summary[j][0])
      entry.append(summary[j][1])
      if summary[j][0] < 1E-10 and summary[j][1] < 1E-6:
        entry.append('YES')
      else:
        entry.append('** NO **')
      sumrytable.append(entry)
    j += 1
  outf.write("\n\n")
  outf.write(tabulate(sumrytable,headers,tablefmt="simple",floatfmt=".4E"))
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
def runUnit(doPrint):
# Runs the unit tests
  global errors, summary
  refvalues()
  tests = [None for x in xrange(1000)] 
  testtable = genTable()
  summary = []

  k = 0
  for i in testtable:
    errors  = []
    if findFile(i.infile,"."):
#
#     run chronus
      if doPrint:
        print "running job: "+i.infile
      tests[k] = runCQ(i.infile,'')
#
#     test SCF values
      if 'SCF' in ref[i.infile[:8]].typ:
        testSCF(ref[i.infile[:8]],tests[k])
        summary.append(errors)

#     test RESP values
      elif 'RESP' in ref[i.infile[:8]].typ:
        testRESP(ref[i.infile[:8]],tests[k])
        summary.append(errors)

      elif 'RT' in ref[i.infile[:8]].typ:
        testRT(ref[i.infile[:8]],tests[k])
        summary.append(errors)

      else:
        print "Not recognized job type for ", ref[i.infile[:8]].typ
        sys.exit()

    k += 1
  genSummary(testtable,summary)
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

# parse user options
  msg = """
  python runtests.py [-o --option=]

  Options:
    -h, --help        Print usage instructions
    -s, --silent      Disable Print
    --type=XXX        Determines types of tests to run. Multiple options
                      can be specified by separating with a comma.
                      3 classes of tests  = [SCF,RESP,RT] 
                      Specify References  = [RHF,UHF,CUHF,GHF]
                                            [RKS,UKS,SLATER,LSDA,SVWN5]
                      Reference and Type  = [(R|U|CU)HF-SCF,HF-CIS,HF-RPA]
                                            [(R|U)KS-SCF,SCF-LSDA] 
                      Dipole Field        = [DField]
    --integrals=XXX   Integral evaluation = [incore] or [direct]
    --parallel=XXX    Whether to run parallel jobs = [on] or [off]
    --size=XXX        Size of jobs to run = [small] or [large] or [both]
                      [small] is the default
    --complex=XXX     Complex Jobs = [yes] or [no] or [both]
                      [both] is the default
    --basis=XXX       Only run tests for this basis set
                      [STO-3G,6-31G,cc-pVDZ,def2-SVPD]
"""
  doPrint  = True
  testType = 'all'
  testInts = ''
  testPar  = ''
  testSize = 'small'
  testComp = ''
  testBasis = 'all'
  try:
    opts, args = getopt.getopt(sys.argv[1:],"hs",["help","silent","type=","integrals=","parallel=","size=","complex=","basis="])
  except getopt.GetoptError:
    print msg
    sys.exit(2)
  for opt, arg in opts:
    if opt in ('-h',"--help"):
      print msg
      sys.exit()
    elif opt in ('-s',"--silent"):
      doPrint = False
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
  runUnit(doPrint)    
#--------------------------------------------------------------------

