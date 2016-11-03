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

#
# Class to store metadata about a keyword
# (self explanitory)
#
class CQKeyword:
  def __init__(self,name,typ,req):
    self.name = name
    self.typ  = typ
    self.req  = req

# Dictionary of known sections
knownSections  = [ 
  "MOLECULE", 
  "QM", 
  "MISC",
  "BASIS",
  "INTS",
  "SCF",
  "RT"
#  "CIS",
#  "RPA",
#  "STAB",
#  "PPRPA",
#  "PPATDA",
#  "PPCTDA"
]

#
# Initialize an empty dictionary to store known
# keywords for each known section
#
knownKeywords = {}

# Dictionary for known keywords in the MOLECULE input section
knownKeywords['MOLECULE'] = {
  'CHARGE':CQKeyword('CHARGE','I',True),
  'MULT'  :CQKeyword('MULT'  ,'I',True),
  'GEOM'  :CQKeyword('GEOM'  ,'S',True),
  'PRINT' :CQKeyword('PRINT' ,'I',False)
} 

# Dictionary for known keywords in the QM input section
knownKeywords['QM'] = {
  'REFERENCE'   :CQKeyword('REFERENCE'  ,'S',True ),
  'JOB'         :CQKeyword('JOB'        ,'S',True ),
#  'INTS'        :CQKeyword('INTS'       ,'S',False),
  'EXCHANGE'    :CQKeyword('EXCHANGE'   ,'S',False),
  'CORR'        :CQKeyword('CORR'       ,'S',False),
#  'DFT_GRID'    :CQKeyword('DFT_GRID'   ,'S',False),
#  'DFT_WEIGHTS' :CQKeyword('DFT_WEIGHTS','S',False),
  'DFT_NRAD'    :CQKeyword('DFT_NRAD'   ,'I',False),
  'DFT_NANG'    :CQKeyword('DFT_NANG'   ,'I',False),
  'DFT_SCREEN'  :CQKeyword('DFT_SCREEN' ,'B',False),
  'DFT_SCRTOL'  :CQKeyword('DFT_SCRTOL' ,'D',False),
  'PRINT'       :CQKeyword('PRINT'      ,'I',False)
}

knownKeywords['INTS'] = {
  'ALG' : CQKeyword("ALG",'S',False)
}

# Dictionary for known keywords in the BASIS section
knownKeywords['BASIS'] = {
  'BASIS'       :CQKeyword('BASIS'      ,'S',True ),
  'FORCECART'   :CQKeyword('FORCECART'  ,'B',False)
}

# Dictionary for known keywords in the RT input section
knownKeywords['RT'] = {
  'MAXSTEP'  :CQKeyword('MAXSTEP'  ,'I' ,False),    
  'TIMESTEP' :CQKeyword('TIMESTEP' ,'D' ,False), 
  "IRSTRT"   :CQKeyword('IRSRT'    ,'I' ,False),
  'EDFIELD'  :CQKeyword('EDFIELD'  ,'D3',False),
  'TIME_ON'  :CQKeyword('TIME_ON'  ,'D' ,False),
  'TIME_OFF' :CQKeyword('TIME_OFF' ,'D' ,False),
  'FREQUENCY':CQKeyword('FREQUENCY','D' ,False), 
  'PHASE'    :CQKeyword('PHASE'    ,'D' ,False),
  'SIGMA'    :CQKeyword('SIGMA'    ,'D' ,False),
  'ENVELOPE' :CQKeyword('ENVELOPE' ,'S' ,False),
  'ORTHO'    :CQKeyword('ORTHO'    ,'S' ,False), 
  'ELL_POL'  :CQKeyword('ELL_POL'  ,'S' ,False), 
  'INIDEN'   :CQKeyword('INIDEN'   ,'I' ,False),
  'UPROP'    :CQKeyword('UPROP'    ,'S' ,False),
  'PRINT'    :CQKeyword('PRINT'    ,'I' ,False),
  'TARCSVS'  :CQKeyword('TARCSVS'  ,'B' ,False)
}

knownKeywords['MISC'] = {
  'NSMP'     :CQKeyword('NSMP'    ,'I',False),
  'UNITTEST' :CQKeyword('UNITTEST','S',False)
}

knownKeywords['SCF'] = {
  'DENTOL' :CQKeyword('DENTOL' ,'D'    ,False),
  'ENETOL' :CQKeyword('ENETOL' ,'D'    ,False),
  'MAXITER':CQKeyword('MAXITER','I'    ,False),
  'FIELD'     :CQKeyword('FIELD'     ,'D3'   ,False),
  'GUESS'     :CQKeyword('GUESS'     ,'S'    ,False),
  'DIIS'      :CQKeyword('DIIS'      ,'B'    ,False),
  'ITP'       :CQKeyword('ITP'       ,'B'    ,False),
  'DT'        :CQKeyword('DT'        ,'D'    ,False),
  'PRINT'     :CQKeyword('PRINT'     ,'I'    ,False)
}

#knownKeywords['CIS'] = {
#  'NSTATES':CQKeyword("NSTATES",'I',True)
#}
#
#knownKeywords['RPA'] = {
#  'NSTATES':CQKeyword("NSTATES",'I',True)
#}
#
#knownKeywords['PPRPA'] = {
#  'NSTATES':CQKeyword("NSTATES",'I',True)
#}
#knownKeywords['PPATDA'] = {
#  'NSTATES':CQKeyword("NSTATES",'I',True)
#}
#knownKeywords['PPCTDA'] = {
#  'NSTATES':CQKeyword("NSTATES",'I',True)
#}
#
#knownKeywords['STAB'] = {
#  'NSTATES':CQKeyword("NSTATES",'I',False)
#}



# Create a dictionary of required keywords
requiredKeywords = {}
#requiredKeywords['MOLECULE'] = []
#requiredKeywords['QM'] = []
#requiredKeywords['RT'] = []
#requiredKeywords['SCF'] = []
#requiredKeywords['MISC'] = []
#requiredKeywords['CIS'] = []
#requiredKeywords['RPA'] = []
#requiredKeywords['STAB'] = []

for sec in knownKeywords:
  requiredKeywords[sec] = []
  for keyWord in knownKeywords[sec]:
    if knownKeywords[sec][keyWord].req:
      requiredKeywords[sec].append(keyWord)
