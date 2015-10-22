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
