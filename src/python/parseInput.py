import os,sys
#sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
#sys.path.append('/home/dbwy/git_repo/chronusq/src/python')
import configparser
import libpythonapi as chronusQ
from parseMolecule import parseMolecule
from parseQM import parseQM

class CQError(Exception):
  def __init__(self,msg):
    self.msg = msg
  def __str__(self):
    return msg


def genSecDict(parser,section):
  dict1 = {}
  opts  = parser.options(section)
  
  for opt in opts:
    dict1[opt] = parser.get(section,opt)

  return dict1


def parseMisc(workers,settings): 
  print 'hello'

def parseInput(workers,iFileName):
  inputParser = configparser.ConfigParser()
  inputParser.read(iFileName)

  secDict = {}

  knownKeywords  = [ "molecule", "qm", "misc"  ]
  parseFunctions = { "molecule":parseMolecule ,
                     "qm":parseQM             ,
                     "misc":parseMisc          }

  for section in inputParser.sections():
    secStr = str(section).lower()
    if secStr not in knownKeywords:
      print "Keyword " + secStr + " not recognized"
      continue
    secDict[secStr] = genSecDict(inputParser,section) 



# print secDict['molecule']['charge']
# print secDict['molecule']['mult']
# print secDict['molecule']['geom']
# for i in secDict: 
#   parseFunctions[i](workers,secDict[i])

  parseFunctions["molecule"](workers,secDict["molecule"])
  parseFunctions["qm"](workers,secDict["qm"])
