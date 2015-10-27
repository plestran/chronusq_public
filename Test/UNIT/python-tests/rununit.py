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

def findFile(name,path):
#  Checks whether file "name" exists in "path"
	found = False
	for root, dirs, files in os.walk(path):
		if name in files:
			found = True
			break
	return found

def genSummary(testtable,summary):
#  Prints the results in summary.txt
	outf = open('summary.txt','w')

# SCF output    
	headers = ["Test Job","|dEnergy|","max(|dDipole|)","max(|dQuadrupole|)","max(|dOctupole|)","Passed"]
	sumrytable = []
	j = 0
	for i in testtable:
		if 'SCF' in ref[j].typ:
			entry = []
			entry.append(testtable[j].infile.replace(".inp",''))
			for k in range(4):
				entry.append(summary[j][k])
			if summary[j][0] < 1E-7 and summary[j][1] < 1E-4 and summary[j][2] < 1E-7 and summary[j][3] < 1E-7:
				entry.append('YES')
			else:
				entry.append('** NO **')
			sumrytable.append(entry)
		j += 1
	outf.write(tabulate(sumrytable,headers,tablefmt="simple",floatfmt=".4E"))

def genTable():
# Reads test.index to see which tests to run
	start = False
	f = open("test.index",'r')
	table = []
	for line in f:
		strx = line.split()
		if '-----' in line:
			start = True
			continue
		if start:    
			if ".inp" in line and '#' not in line:
				desc = ' '.join(strx[2:])
				entry = UnitTest(strx[0],strx[1],desc)
				table.append(entry)
 	f.close()
 	return table

def runUnit(doKill,doPrint):
# Runs the unit tests
	global errors, summary
	refvalues()
	tests = [[None for x in xrange(100000)] for y in xrange(10)]
	testtable = genTable()
	summary = []
	errors  = []

	k = 0
	for i in testtable:
		if findFile(i.infile,"."):
#
#			run chronus
			if doPrint:
				print "running file: "+i.infile.replace(".inp",'')
			tests[k][0] = runCQ(i.infile.replace(".inp",''))
			print tests[k][0].E
#
#			test SCF values
			if 'SCF' in ref[k].typ:
				testSCF(ref[k],tests[k][0])
				summary.append(errors)

		k += 1
	genSummary(testtable,summary)

def testSCF(ref,tests):

	auToD = 0.393456
	auToAng = 0.529177
# test SCF energy
	abserr = abs(ref.scf - tests.E)
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


##############################
#        Main Program        #
##############################

if __name__ in "__main__":
	msg = """python runtests.py [-o --option]

  Options:
    -h, --help        Print usage instructions
    --enable-kill    Enable script termination on Unit Test failure
    --enable-travisci    Enable options for Travis-CI run
    -s, --silent        Disable Print"""
	doKill = False
	doPrint = True
	try:
		opts, args = getopt.getopt(sys.argv[1:],"hs",["enable-travisci","enable-kill","help","silent"])
	except getopt.GetoptError:
		print msg
		sys.exit(2)
	for opt, arg in opts:
		if opt in ('-h',"--help"):
			print msg
			sys.exit()
		elif opt in ("--enable-travisci","--enable-kill"):
			doKill = True
		elif opt in ('-s',"--silent"):
			doPrint = False

	runUnit(doKill,doPrint)    

