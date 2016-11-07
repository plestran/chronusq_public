#
# The Chronus Quantum (ChronusQ) software package is high-performace 
# computational chemistry software with a strong emphasis on explicitly 
# time-dependent and post-SCF quantum mechanical methods.
# 
# Copyright (C) 2014-2016 Li Research Group (University of Washington)
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
from meta.knownKeywords import requiredKeywords
from libpythonapi import CErrMsg

#
# Parse Basis Set file
#
# ** Note that this allocates the CQ::BasisSet Object **
#
def parseBasis(workers,basisSettings):
#
#  Check that all of the required keywords for Basis
#  object are found
#
  for i in requiredKeywords['BASIS']:
    if i not in basisSettings:
      msg = 'Required keyword BasisSet.' + str(i) + ' not found'
      CErrMsg(workers['CQFileIO'],msg)
  
#  workers["CQBasisSet"] = chronusQ.BasisSet()
  basis =  basisSettings['BASIS']

  # Make CQ::BasisSet aware of CQ::FileIO
  workers["CQBasisSet"].communicate(workers["CQFileIO"])

  # Check to see if we're forcing cartesian functions
  forceCart = False
  if 'FORCECART' in basisSettings:
    if basisSettings['FORCECART']:
      workers["CQBasisSet"].forceCart()
  
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
  workers["CQBasisSet"].makeMaps(workers["CQMolecule"])  # 4
  workers["CQBasisSet"].renormShells()                        # 5

