#
# The Chronus Quantum (ChronusQ) software package is high-performace 
# computational chemistry software with a strong emphasis on explicitly 
# time-dependent and post-SCF quantum mechanical methods.
# 
# Copyright (C) 2014-2015 Li Research Group (University of Washington)
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 
# Contact the Developers:
#   E-Mail: xsli@uw.edu
# 
#
import os,sys
import libpythonapi as chronusQ

#
# Parse Basis Set file
#
# ** Note that this allocates the CQ::BasisSet Object **
#
def parseBasis(workers,basis):
  nTCS = workers["CQSingleSlater"].nTCS()

  # Make CQ::BasisSet aware of CQ::FileIO
  workers["CQBasisSet"].communicate(workers["CQFileIO"])
  
#
# 1) Try to find the basis set file in standard locations
# 2) Allocate and populate a global basis set definition
# 3) Construct a local subset (viz. CQ::Molecule) of the basis set
# 4) Construct the varios maps incolving the basis set
# 5) Renormalize the LibInt2::Shell's (this is outdated in newer Libint)
#
  # FIXME: BasisSet Keywords?? str(basis).lower() is a hack
  workers["CQBasisSet"].findBasisFile(str(basis).lower())     # 1
  workers["CQBasisSet"].parseGlobal()                         # 2
  workers["CQBasisSet"].constructLocal(workers["CQMolecule"]) # 3
  workers["CQBasisSet"].makeMaps(nTCS,workers["CQMolecule"])  # 4
  workers["CQBasisSet"].renormShells()                        # 5
  
