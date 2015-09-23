import os,sys

fname = sys.argv[1]
f = open(fname,'r')

es = []
osc = []

for line in f:
	strx = line.split()
	if 'Excited State  ' in line:
		es.append(float(strx[4]))
		osc.append(float(strx[8][2:]))

stry = ''
for i in range(len(es)):
	stry += str(es[i])+','+str(osc[i])+'/'
print stry
