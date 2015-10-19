
class CQKeyword:
  def __init__(self,name,typ,req):
    self.name = name
    self.typ  = typ
    self.req  = req

knownKeywords = {}
knownKeywords['MOLECULE'] = {
  'CHARGE':CQKeyword('CHARGE','I',True),
  'MULT'  :CQKeyword('MULT'  ,'I',True),
  'GEOM'  :CQKeyword('GEOM'  ,'S',True)
} 

knownKeywords['QM'] = {
  'REFERENCE':CQKeyword('REFERENCE','S',True),
  'BASIS'    :CQKeyword('BASIS'    ,'S',True),
  'JOB'      :CQKeyword('JOB'      ,'S',True)
}

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


requiredKeywords = {}
requiredKeywords['MOLECULE'] = []
requiredKeywords['QM'] = []
requiredKeywords['RT'] = []

for sec in knownKeywords:
  for keyWord in knownKeywords[sec]:
    if knownKeywords[sec][keyWord].req:
      requiredKeywords[sec].append(keyWord)
