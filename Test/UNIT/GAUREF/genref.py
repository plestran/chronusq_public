import sys,os
from grabinfo_SCF import *

def genTable():
	start = False
	f = open("ref.index",'r')
	table = []
	for line in f:
		strx = line.split()
		if '-----' in line:
			start = True
			continue
		if start:
			table.append([strx[0],strx[1].replace('.com','.log')])
	return table

def genRef():
	table = genTable()
	
	of = open("ref.val",'w')
	for i in range(len(table)):
		of.write(table[i][0]+'/'+grabinfoSCF(table[i][1])+'\n')
