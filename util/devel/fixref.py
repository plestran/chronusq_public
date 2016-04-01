import os,sys

refFileName = sys.argv[1]

refFile = open(refFileName)

refVals = []
for line in refFile:
  splt = line.split('/')
  refVals.append(splt)

uniqueVals = []
for i in range(len(refVals)):
  if i == 0:
    uniqueVals.append(refVals[i])
  elif refVals[i][0] != refVals[i-1][0]:
    uniqueVals.append(refVals[i])
  else:
    for j in range(len(refVals[i])):
      if j != 0 and j!= (len(refVals[i])-1):
        if abs(float(refVals[i][j])-float(refVals[i-1][j])) > 1e-6:
	  print 'Warning! ',refVals[i][0],i,j


print uniqueVals
fixedFName = refFileName+".fixed"
fixedF = open(fixedFName,'w')

for i in range(len(uniqueVals)):
  for j in range(len(uniqueVals[i])):
    if j != (len(uniqueVals[i]) - 1):
      print 'here'
      fixedF.write(uniqueVals[i][j]+'/')
    else:
      fixedF.write(uniqueVals[i][j]+'\n')
