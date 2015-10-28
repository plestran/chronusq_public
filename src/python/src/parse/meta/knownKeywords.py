
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
  "SCF",
  "RT",
  "CIS",
  "RPA",
  "STAB" 
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
  'GEOM'  :CQKeyword('GEOM'  ,'S',True)
} 

# Dictionary for known keywords in the QM input section
knownKeywords['QM'] = {
  'REFERENCE':CQKeyword('REFERENCE','S',True ),
  'BASIS'    :CQKeyword('BASIS'    ,'S',True ),
  'JOB'      :CQKeyword('JOB'      ,'S',True ),
  'INTS'     :CQKeyword('INTS'     ,'S',False),
  'EXCHANGE' :CQKeyword('EXCHANGE' ,'S',False),
  'CORR'     :CQKeyword('CORR'     ,'S',False)
}

# Dictionary for known keywords in the RT input section
knownKeywords['RT'] = {
  'MAXSTEP'  :CQKeyword('MAXSTEP'  ,'I'      ,False),    
  'TIMESTEP' :CQKeyword('TIMESTEP' ,'D'      ,False), 
  'EDFIELD'  :CQKeyword('EDFIELD'  ,'D3'     ,False),
  'TIME_ON'  :CQKeyword('TIME_ON'  ,'D'      ,False),
  'TIME_OFF' :CQKeyword('TIME_OFF' ,'D'      ,False),
  'FREQUENCY':CQKeyword('FREQUENCY','D'      ,False), 
  'PHASE'    :CQKeyword('PHASE'    ,'D'      ,False),
  'SIGMA'    :CQKeyword('SIGMA'    ,'D'      ,False),
  'ENVELOPE' :CQKeyword('ENVELOPE' ,'O-ENV'  ,False),
  'ORTHO'    :CQKeyword('ORTHO'    ,'O-ORTH' ,False), 
  'INIDEN'   :CQKeyword('INIDEN'   ,'I'      ,False),
  'UPROP'    :CQKeyword('UPROP'    ,'O-FORMU',False)
}

knownKeywords['MISC'] = {
  'NSMP'     :CQKeyword('NSMP'    ,'I',False),
  'UNITTEST' :CQKeyword('UNITTEST','S',False)
}

knownKeywords['SCF'] = {
  'SCFDENTOL' :CQKeyword('SCFDENTOL' ,'D'    ,False),
  'SCFENETOL' :CQKeyword('SCFENETOL' ,'D'    ,False),
  'SCFMAXITER':CQKeyword('SCFMAXITER','I'    ,False),
  'FIELD'     :CQKeyword('FIELD'     ,'D3'   ,False),
  'GUESS'     :CQKeyword('GUESS'     ,'O-GS' ,False),
  'DIIS'      :CQKeyword('DIIS'      ,'B'    ,False),
  'PRINT'     :CQKeyword('PRINT'     ,'I'    ,False)
}

knownKeywords['CIS'] = {
  'NSTATES':CQKeyword("NSTATES",'I',True)
}

knownKeywords['RPA'] = {
  'NSTATES':CQKeyword("NSTATES",'I',True)
}

knownKeywords['STAB'] = {
  'NSTATES':CQKeyword("NSTATES",'I',False)
}



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
