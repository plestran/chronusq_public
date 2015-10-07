import os,sys
sys.path.append('/home/dbwy/git_repo/chronusq/build_gcc_libint_openmp/src/python')
sys.path.append('/home/dbwy/git_repo/chronusq/src/python')
import configparser
import libpythonapi as chronusQ

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

def readGeom(workers,settings):
  geomStr = settings['geom']
  geomStr = geomStr.split('\n')
  nAtoms = 0
  for i in geomStr:
    if i not in ('', ' '):
      nAtoms += 1
  oldLen = len(geomStr)
  for i in range(oldLen - nAtoms):
    geomStr.remove('')
  print oldLen,len(geomStr)
  workers["CQMolecule"].setNAtoms(nAtoms)
  workers["CQMolecule"].alloc(workers["CQFileIO"])
  
  nTotalE = 0
  for i in range(len(geomStr)):
    line = geomStr[i]
    lineSplit = line.split()
    indx = -1
    if len(lineSplit) == 5:
      indx = chronusQ.HashAtom(str(lineSplit[0]),int(lineSplit[1]))
    elif len(lineSplit) == 4:
      indx = chronusQ.HashAtom(str(lineSplit[0]),0)
    else:
      print 'Input Error: Invalid Geometry Specification'

    if indx != -1:
      workers["CQMolecule"].setIndex(i,indx)
    else:
      print 'Input Error: Invalid Atomic Symbol or Mass Number'

    nTotalE += chronusQ.getAtomicNumber(indx)

    if len(lineSplit) == 5:
      x = float(lineSplit[2])
      y = float(lineSplit[3])
      z = float(lineSplit[4])
    elif len(lineSplit) == 4:
      x = float(lineSplit[1])
      y = float(lineSplit[2])
      z = float(lineSplit[3])

    workers["CQMolecule"].setCart(i,x,y,z)
  workers["CQMolecule"].setNTotalE(nTotalE)
  workers["CQMolecule"].computeNucRep()
  workers["CQMolecule"].convBohr()
  workers["CQMolecule"].computeRij()
  workers["CQMolecule"].computeI()
    
    

def parseMolecule(workers,settings): 
  workers["CQMolecule"].setCharge(int(settings['charge']))
  workers["CQMolecule"].setMultip(int(settings['mult'  ]))
  readGeom(workers,settings)

def parseQM(workers,settings): 
  print 'hello'
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
  for i in secDict: 
    parseFunctions[i](workers,secDict[i])

