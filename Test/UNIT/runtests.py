import os,sys
import GAUREF.genref as GR
from tabulate import tabulate

class UnitTest:
	def __init__(self,infile,testClass,desc):
		self.infile    = infile
		self.testClass = testClass
		self.desc      = desc

class RespData:
	def __init__(self,energy,f):
		self.energy = energy
		self.f      = f

class MaxErrorExcedeed(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

#
#  Generate a Table of Unit Test input files
#  and their descriptions from "test.index"
#
#  This routine ignores lines with '#':
#  This can act as a convinient way to toggle
#  which tests actually run
#
def genTable():
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
#				table.append([strx[0],desc])
				table.append(entry)
	f.close()
	return table
#
#  Checks whether file "name" exists in "path"
#
def findFile(name,path):
	found = False
	for root, dirs, files in os.walk(path):
		if name in files:
			found = True
			break
	return found

#
#  Generate a dictionary that maps test number string
#  of form "testXXXX" to array that holds reference
#  values.
#
#  No assumption is made of the class of test
#
def genRefDict():
	ref = open("GAUREF/ref.val",'r')
	refdict = {}
	for line in ref:
		strx=line.split('/')
		refdict[strx[0]] = strx[1:]
	
	for i in refdict:
		for j in range(len(refdict[i])):
			try:
				refdict[i][j] = float(refdict[i][j])
			except ValueError:
				if ',' in refdict[i][j]:
			 		sp = refdict[i][j].split(',')
					refdict[i][j] = RespData(float(sp[0]),float(sp[1]))
				else:
					raise
			
	ref.close()
	return refdict

def genSummary(testtable,summary):
	outf = open('summary.txt','w')
	
	headers = ["Test Job","|dEnergy|","max(|dDipole|)","max(|dQuadrupole|)","max(|dOctupole|)","Passed"]
	sumrytable = []
	for i in range(len(testtable)):
		entry = []
		entry.append(testtable[i].infile.replace(".inp",''))
		entry.append(summary[i][0])
		entry.append(max(summary[i][1:4]))
		entry.append(max(summary[i][5:11]))
		entry.append(max(summary[i][12:22]))
		if summary[i][0] < 1E-08 and max(summary[i][1:4]) < 1E-8 and max(summary[i][5:11]) < 1E-8 and max(summary[i][12:22]) < 1E-8:
			entry.append('YES')
		else:
			entry.append('** NO **')
		sumrytable.append(entry)
	
	outf.write(tabulate(sumrytable,headers,tablefmt="simple",floatfmt=".4E"))



#
#  Runs Unit Tests
#
def runUnit():
	testtable = genTable()
	refdict = genRefDict()
		
	summary = []
	doKill = False
	for i in testtable:
		if findFile(i.infile,"."):
			print "../../chronusQ "+i.infile.replace(".inp",'')
			os.system("../../chronusQ "+i.infile.replace(".inp",'')+" > tmp")
			tmp = open("tmp",'r')
			line = tmp.readline()
			vals = line.split('/')
			line = tmp.readline()
			strx = line.split()
			time = float(strx[6])
			err = []
			for j in range(len(refdict[i.infile[:8]])):
				try:
					if refdict[i.infile[:8]][j] == 0:
						abserr = abs(float(vals[j]) - refdict[i.infile[:8]][j])
					else:
						abserr = abs(float(vals[j]) - refdict[i.infile[:8]][j])/refdict[i.infile[:8]][j]
					err.append(abs(abserr))
					if abs(abserr) > 1E-8:
						raise MaxErrorExcedeed(abs(abserr))
				except MaxErrorExcedeed:
					if doKill:
						raise
					else:
						continue
			summary.append(err)
	genSummary(testtable,summary)


if __name__ in "__main__":
	runUnit()	
