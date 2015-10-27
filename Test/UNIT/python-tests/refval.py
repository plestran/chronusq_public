import sys,os
#
# Collects all the data from the reference calculations
#	in the data class
#
class data:
#	def __init__(self,scf,dip,quad,octu,osc,w,typ):
	def __init__(self,scf,dip,quad,octu,typ):
		self.scf   = scf
		self.dip   = dip
		self.quad  = quad
		self.octu  = octu
#		self.osc   = osc
#		self.w     = w
		self.typ   = typ

#
#	Grab the reference values and store them
#
ref = {}
with open("ref.val") as f:
	for line in f:
		dipole     = []
		quadrupole = []
		octupole   = []
		val        = line.split('/')
		for i in range(2,5):
			dipole.append(float(val[i]))
		for i in range(5,11):
			quadrupole.append(float(val[i]))
		for i in range(11,20):
			octupole.append(float(val[i]))
		octupole.append(float(val[20].rstrip()))
		ref[val[0]] = data(float(val[1]),dipole,quadrupole,octupole,'SCF'); 

def refvalues():
	return ref


