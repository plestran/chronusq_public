import sys,os
#
# Collects all the data from the reference calculations
#	in the data class
#
class data:
#	def __init__(self,scf,dip,quad,octu,osc,w):
	def __init__(self,scf,dip,quad,octu):
		self.scf   = scf
		self.dip   = dip
		self.quad  = quad
		self.octu  = octu
#		self.osc   = osc
#		self.w     = w

ref = [None]*100000
dipole = [0.0000,1.5340,0.0000]
quadrupole = [-4.3168,-5.2054,-6.2076,0.0000,0.0000,0.0000]
octupole = [0.0000,-0.9572,0.0000,0.0000,0.3916,0.0000,0.0000,-0.4476,0.0000,0.0000]
ref[0] = data(-74.9420799245,dipole,quadrupole,octupole);

def refvalues():
	return ref


