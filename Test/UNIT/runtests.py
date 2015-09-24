import os,sys
import GAUREF.genref as GR
from tabulate import tabulate

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
				desc = ' '.join(strx[1:])
				table.append([strx[0],desc])
	return table

def findFile(name,path):
	found = False
	for root, dirs, files in os.walk(path):
		if name in files:
			found = True
			break
	return found

testtable = genTable()
#for i in testtable:
#	if findFile(i[0],"."):
#		print "../../chronusQ "+i[0].replace(".inp",'')
	
ref = open("GAUREF/ref.val",'r')
refdict = {}
for line in ref:
	strx=line.split('/')
	refdict[strx[0]] = strx[1:]

for i in refdict:
	for j in range(len(refdict[i])):
		refdict[i][j] = float(refdict[i][j])

summary = []
for i in testtable:
	if findFile(i[0],"."):
 		print "../../chronusQ "+i[0].replace(".inp",'')
		os.system("../../chronusQ "+i[0].replace(".inp",'')+" > tmp")
		tmp = open("tmp",'r')
		line = tmp.readline()
		vals = line.split('/')
		line = tmp.readline()
		strx = line.split()
		time = float(strx[6])	
		err = []
		for j in range(len(refdict[i[0][:8]])):
			if refdict[i[0][:8]][j] == 0:
				abserr = abs(float(vals[j]) - refdict[i[0][:8]][j])
			else:
				abserr = abs(float(vals[j]) - refdict[i[0][:8]][j])/refdict[i[0][:8]][j]
			err.append(abs(abserr))
		summary.append(err)

outf = open('summary.txt','w')

headers = ["Test Job","|dEnergy|","max(|dDipole|)","max(|dQuadrupole|)","max(|dOctupole|)","Passed"]
sumrytable = []
for i in range(len(testtable)):
	entry = []
	entry.append(testtable[i][0].replace(".inp",''))
	entry.append(summary[i][0])
	entry.append(max(summary[i][1:4]))
	entry.append(max(summary[i][5:11]))
	entry.append(max(summary[i][12:22]))
	if summary[i][0] < 1E-08 and max(summary[i][1:4]) < 1E-8 and max(summary[i][5:11]) < 1E-8 and max(summary[i][12:22]) < 1E-8:
		entry.append('XXXX')
	else:
		entry.append('')
	sumrytable.append(entry)

outf.write(tabulate(sumrytable,headers,tablefmt="simple",floatfmt=".4E"))



		
			
