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

def parseMisc(workers,settings):
#  print 'Parsing Misc Information'

#
#  Check that all of the required keywords for Molecule
#  object are found
#
  for i in requiredKeywords['MISC']:
    if i not in settings:
      msg = 'Required keyword Misc.' + str(i) + ' not found'
      CErrMsg(workers['CQFileIO'],str(msg))

  optMap = {
    'NSMP':chronusQ.CQSetNumThreads
  }

  # Loop over optional keywords, set options accordingly
  # note that because these are optional, if the keyword
  # is not found in setings, no error is thrown and the
  # next keyword is processed
  for i in optMap:
    try:
      optMap[i](settings[i])
    except KeyError:
      continue

  if 'PRINT' in settings:
    workers['CQMolecule'].setPrintLevel(settings['PRINT'])
    workers['CQBasisSet'].setPrintLevel(settings['PRINT'])
    workers['CQSingleSlater'].setPrintLevel(settings['PRINT'])
    workers['CQAOIntegrals'].setPrintLevel(settings['PRINT'])
    if 'CQRealTime' in workers:
      workers['CQRealTime'].setPrintLevel(settings['PRINT'])
    
