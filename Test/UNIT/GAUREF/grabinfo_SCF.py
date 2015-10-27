import os,sys
def grabinfoSCF(fname):
	f = open(fname,'r')
	
	E = 0.0
	dX = 0.0
	dY = 0.0
	dZ = 0.0
	doTQ = False
	dXX = 0.0
	dXY = 0.0
	dXZ = 0.0
	dYY = 0.0
	dYZ = 0.0
	dZZ = 0.0
	dTQXX = 0.0
	dTQXY = 0.0
	dTQXZ = 0.0
	dTQYY = 0.0
	dTQYZ = 0.0
	dTQZZ = 0.0
	dXXX = 0.0
	dYYY = 0.0
	dZZZ = 0.0
	dXYY = 0.0
	dXXY = 0.0
	dXXZ = 0.0
	dXZZ = 0.0
	dYZZ = 0.0
	dYYZ = 0.0
	dXYZ = 0.0
	for line in f:
		strx = line.split()
		if 'SCF Done' in line:
			E = float(strx[4])
		if ' X=' in line:
			dX = float(strx[1])
			dY = float(strx[3])
			dZ = float(strx[5])
		if ' XX=' in line and not doTQ:
			dXX = float(strx[1])
			dYY = float(strx[3])
			dZZ = float(strx[5])
		if ' XY=' in line and not doTQ:
			dXY = float(strx[1])
			dXZ = float(strx[3])
			dYZ = float(strx[5])
		if 'Traceless Quadrupole' in line:
			doTQ = True
		if ' XX=' in line and doTQ:
			dTQXX = float(strx[1])
			dTQYY = float(strx[3])
			dTQZZ = float(strx[5])
		if ' XY=' in line and doTQ:
			dTQXY = float(strx[1])
			dTQXZ = float(strx[3])
			dTQYZ = float(strx[5])
		if ' XXX=' in line:
			dXXX = float(strx[1])
			dYYY = float(strx[3])
			dZZZ = float(strx[5])
			dXYY = float(strx[7])
		if ' XXY=' in line:
			dXXY = float(strx[1])
			dXXZ = float(strx[3])
			dXZZ = float(strx[5])
			dYZZ = float(strx[7])
		if ' YYZ=' in line:
			dYYZ = float(strx[1])
			dXYZ = float(strx[3])
	return '{0:.10f}/{1:.4f}/{2:.4f}/{3:.4f}/{4:.4f}/{5:.4f}/{6:.4f}/{7:.4f}/{8:.4f}/{9:.4f}/{10:.4f}/{11:.4f}/{12:.4f}/{13:.4f}/{14:.4f}/{15:.4f}/{16:.4f}/{17:.4f}/{18:.4f}/{19:.4f}'.format(E,dX,dY,dZ,dXX,dXY,dXZ,dYY,dYZ,dZZ,dXXX,dXXY,dXXZ,dXYY,dXYZ,dXZZ,dYYY,dYYZ,dYZZ,dZZZ)
