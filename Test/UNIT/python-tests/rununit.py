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

def findFile(name,path):
#  Checks whether file "name" exists in "path"
	found = False
	for root, dirs, files in os.walk(path):
		if name in files:
			found = True
			break
	return found

def runUnit(doKill,doPrint):
# Runs the unit tests
	tests = [[None for x in xrange(100000)] for y in xrange(10)]
	testtable = genTable()
	refvalues()

	k = 0
	for i in testtable:
		if findFile(i.infile,"."):
			print "file = "+i.infile.replace(".inp",'')
#			print ref[k].scf
			tests[k][0] = runCQ(i.infile.replace(".inp",''))
			print tests[k][0].E
		k += 1

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


# Load the reference values into ref
#refvalues()
#print ref[0].scf
#
#fname = "test0001_serial_incore"
#tests = [[None for x in xrange(100000)] for y in xrange(10)]
#tests[0][0] = runCQ(fname)
#print tests[0][0].E



