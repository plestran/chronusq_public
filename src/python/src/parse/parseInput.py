import os,sys
#sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
#sys.path.append('/home/dbwy/git_repo/chronusq/src/python')
import configparser
import libpythonapi as chronusQ
from parseMolecule import parseMolecule
from parseQM import parseQM
from meta.knownKeywords import knownKeywords
import meta.enumMaps as enumMaps

def parserOrth(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).upper()
  return enumMaps.orthoMap[strX]

def parserEnv(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).upper()
  return enumMaps.envMap[strX]

def parserFormU(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).upper()
  return enumMaps.formUMap[strX]

def parserDoubleArray(parser,section,opt):
  strX = parser.get(section,opt).split()
  for i in range(len(strX)):
    strX[i] = float(strX[i])
  return strX

def genSecDict(parser,section):
  dict1 = {}
  opts  = parser.options(section)

  parseMap = { 
    'I'      :parser.getint,
    'D'      :parser.getfloat,
    'S'      :parser.get,
    'O-ENV'  :parserEnv,
    'O-ORTH' :parserOrth,
    'O-FORMU':parserFormU,
    'D3'     :parserDoubleArray
  }

  keydict = ''
  try: 
    keydict = knownKeywords[section.upper()]
  except KeyError:
    print 'No dictionary of known keywords is known for ' + section
    return

  for opt in opts:
    keyword = ''
    readtyp = ''

    try: 
      keyword = knownKeywords[section.upper()][opt.upper()]
    except KeyError:
      print 'Keyword: '+str(opt)+' not recognized for section '+section
      continue
     
    readTyp = keyword.typ

    if readTyp in ('I','D'):
      dict1[opt.upper()] = parseMap[readTyp](section,opt)
    elif readTyp in ('S'):
      dict1[opt.upper()] = parseMap[readTyp](section,opt).upper()
    elif readTyp in ('O-ORTH','O-ENV','O-FORMU','D3'):
      dict1[opt.upper()] = parseMap[readTyp](parser,section,opt)
    else:
      print 'Option data type not recognized'

  return dict1


def parseMisc(workers,settings): 
  print 'hello'

def parseInput(workers,iFileName):
  inputParser = configparser.ConfigParser()
  inputParser.read(iFileName)

  secDict = {}

  knownSections  = [ "MOLECULE", "QM", "MISC","RT"  ]
# parseFunctions = { "molecule":parseMolecule ,
#                    "qm":parseQM             ,
#                    "misc":parseMisc          }

  for section in inputParser.sections():
    secStr = str(section).upper()
    if secStr not in knownSections:
      print "Keyword section " + secStr + " not recognized"
      continue
    secDict[secStr] = genSecDict(inputParser,section) 

  return secDict

# parseFunctions["molecule"](workers,secDict["molecule"])
# parseFunctions["qm"](workers,secDict)
