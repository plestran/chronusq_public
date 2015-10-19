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
  strX = str(strX).lower()
  return enumMaps.orthoMap[strX]

def parserEnv(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).lower()
  return enumMaps.envMap[strX]

def parserFormU(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).lower()
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
    'i'      :parser.getint,
    'd'      :parser.getfloat,
    's'      :parser.get,
    'o-env'  :parserEnv,
    'o-orth' :parserOrth,
    'o-formu':parserFormU,
    'd3'     :parserDoubleArray
  }
  
# for opt in opts:
#    dict1[opt] = parser.get(section,opt)
#    try:
#      readTyp = knownKeywords[section][opt.lower()].typ
##      print opt
#      if readTyp in ('i','d','s'):
#        print 'reading I, D, or S: '+ repr(parseMap[readTyp])
#      elif readTyp in ('o-orth','o-env','o-formu'):
#        print 'reading something from RT: '+ repr(parseMap[readTyp])
#        dict1[opt.lower()] = parseMap[readTyp](parser,section,opt)
#      elif 'd' in readTyp and len(readTyp) > 1:
#        print 'reading a double array of dim ' + readTyp[1:] + ": " + repr(parseMap[readTyp])
#      else:
#        print 'reading something else: '+str(opt)
#     if readTyp in ('i','d','s'):
#       dict1[opt.lower()] = parseMap[readTyp](section,opt)
#     elif readTyp in ('o-orth','o-env','o-formu','d3'):
#       dict1[opt.lower()] = parseMap[readTyp](parser,section,opt)
#     else:
#       print 'reading something else: '+str(opt)
      
#    except KeyError:
#     print 'unrecognized keyword: '+str(opt)
#    else:
#      continue

  try: 
    keydict = knownKeywords[section]
  except KeyError:
    print 'No dictionary of known keywords is known for ' + section
    return

  for opt in opts:
    keydict = ''
    keyword = ''
    readtyp = ''

    try: 
      keyword = knownKeywords[section][opt.lower()]
    except KeyError:
      print 'Keyword: '+str(opt)+' not recognized for section '+section
      continue
     
    readTyp = keyword.typ

    if readTyp in ('i','s','d'):
      dict1[opt.lower()] = parseMap[readTyp](section,opt)
    elif readTyp in ('o-orth','o-env','o-formu','d3'):
      dict1[opt.lower()] = parseMap[readTyp](parser,section,opt)
    else:
      print 'Option data type not recognized'
  return dict1


def parseMisc(workers,settings): 
  print 'hello'

def parseInput(workers,iFileName):
  inputParser = configparser.ConfigParser()
  inputParser.read(iFileName)

  secDict = {}

  knownSections  = [ "molecule", "qm", "misc","rt"  ]
# parseFunctions = { "molecule":parseMolecule ,
#                    "qm":parseQM             ,
#                    "misc":parseMisc          }

  for section in inputParser.sections():
    secStr = str(section).lower()
    if secStr not in knownSections:
      print "Keyword " + secStr + " not recognized"
      continue
    secDict[secStr] = genSecDict(inputParser,section) 

  return secDict

# parseFunctions["molecule"](workers,secDict["molecule"])
# parseFunctions["qm"](workers,secDict)
