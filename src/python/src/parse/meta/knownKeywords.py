
class CQKeyword:
  def __init__(self,name,typ,req):
    self.name = name
    self.typ  = typ
    self.req  = req

knownKeywords = {}
knownKeywords['Molecule'] = {
  'charge':CQKeyword('charge','i',True),
  'mult'  :CQKeyword('mult'  ,'i',True),
  'geom'  :CQKeyword('geom'  ,'s',True)
} 

knownKeywords['QM'] = {
  'reference':CQKeyword('reference','s',True),
  'basis'    :CQKeyword('basis'    ,'s',True),
  'job'      :CQKeyword('job'      ,'s',True)
}

knownKeywords['RT'] = {
  'maxstep'  :CQKeyword('maxstep'  ,'i'      ,False),    
  'timestep' :CQKeyword('timestep' ,'d'      ,False), 
  'edfield'  :CQKeyword('edfield'  ,'d3'     ,False),
  'time_on'  :CQKeyword('time_on'  ,'d'      ,False),
  'time_off' :CQKeyword('time_off' ,'d'      ,False),
  'frequency':CQKeyword('frequency','d'      ,False), 
  'phase'    :CQKeyword('phase'    ,'d'      ,False),
  'sigma'    :CQKeyword('sigma'    ,'d'      ,False),
  'envelope' :CQKeyword('envelope' ,'o-env'  ,False),
  'ortho'    :CQKeyword('ortho'    ,'o-orth' ,False), 
  'iniden'   :CQKeyword('iniden'   ,'i'      ,False),
  'uprop'    :CQKeyword('uprop'    ,'o-formu',False)
}


requiredKeywords = {}
requiredKeywords['Molecule'] = []
requiredKeywords['QM'] = []
requiredKeywords['RT'] = []

for sec in knownKeywords:
  for keyWord in knownKeywords[sec]:
    if knownKeywords[sec][keyWord].req:
      requiredKeywords[sec].append(keyWord)
