import sys,os
from grabinfo_SCF import *
from grabinfo_RESP import *

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
			table.append([strx[0],strx[1],strx[2].replace('.com','.log')])
	return table

class UnitTestNotFound(Exception):
	def __init__(self,value):
		self.value = value
	def __str__(self):
		return repr(self.value)

def genRef():
	table = genTable()
	
	of = open("ref.val",'w')
	for i in range(len(table)):
		of.write(table[i][0]+'/')
		if table[i][1] == 'SCF':
			of.write(grabinfoSCF(table[i][2]))
		elif table[i][1] == 'RESP':
			of.write(grabinfoRESP(table[i][2]))
		else:
			raise UnitTestNotFound(table[i][1])	
		of.write('\n')

if __name__ in "__main__":
	genRef()
