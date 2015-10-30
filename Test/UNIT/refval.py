import sys,os
#
# Collects all the data from the reference calculations
#	in the data class
#
class data:
	def __init__(self,eng,dip,quad,octu,osc,w,typ):
		self.eng   = eng
		self.dip   = dip
		self.quad  = quad
		self.octu  = octu
		self.osc   = osc
		self.w     = w
		self.typ   = typ

#
#	Grab the reference values and store them
#
ref = {}
with open("chronus-ref.val") as f:
	for line in f:
		dipole     = []
		quadrupole = []
		octupole   = []
		oscstr     = []
		omega      = []
		val        = line.split('/')

		if val[-1].rstrip() == 'SCF':
			for i in range(2,5):
				dipole.append(float(val[i]))
			for i in range(5,11):
				quadrupole.append(float(val[i]))
			for i in range(11,20):
				octupole.append(float(val[i]))
			octupole.append(float(val[20].rstrip()))
			ref[val[0]] = data(float(val[1]),dipole,quadrupole,octupole,oscstr,omega,'SCF'); 

		elif val[-1].rstrip() == 'RESP':
			for i in range(1,len(val)-1):
				if i % 2 == 0:
					oscstr.append(float(val[i]))
				else:
					omega.append(float(val[i]))
			ref[val[0]] = data(0.,dipole,quadrupole,octupole,oscstr,omega,'RESP'); 

		elif val[-1].rstrip() == 'RT':
			for i in range(2,6):
				dipole.append(float(val[i]))
			ref[val[0]] = data(float(val[1]),dipole,quadrupole,octupole,oscstr,omega,'RT'); 

		else:
			print "Unrecognized job type"
			print line
			sys.exit()

def refvalues():
	return ref

