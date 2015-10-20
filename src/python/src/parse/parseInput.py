import os,sys
import configparser
import libpythonapi as chronusQ
from parseMolecule import parseMolecule
from parseQM import parseQM
from meta.knownKeywords import knownKeywords
import meta.enumMaps as enumMaps

#
# Convert string obtained from ConfigParser
# for RT.ortho into RT.ORTHO enum
#
def parserOrth(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).upper()
  return enumMaps.orthoMap[strX]

#
# Convert string obtained from ConfigParser
# for RT.envelope into RT.ENVELOPE enum
#
def parserEnv(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).upper()
  return enumMaps.envMap[strX]

#
# Convert string obtained from ConfigParser
# for RT.uprop into RT.FORM_U enum
#
def parserFormU(parser,section,opt):
  strX = parser.get(section,opt)
  strX = str(strX).upper()
  return enumMaps.formUMap[strX]

#
# Convert string obtained from ConfigParser
# into an N dimenstional float array
#
def parserDoubleArray(parser,section,opt):
  strX = parser.get(section,opt).split()
  for i in range(len(strX)):
    strX[i] = float(strX[i])
  return strX

#
# Generate a dictionary for a given section of
# the input file
#
def genSecDict(parser,section):

  # Initialize an empy dictionary
  dict1 = {}

  # Grab options for given section
  opts  = parser.options(section)

  # Define a map from data type to parser function
  # (Includes user defined wrappers for special types)
  parseMap = { 
    'I'      :parser.getint,
    'D'      :parser.getfloat,
    'S'      :parser.get,
    'O-ENV'  :parserEnv,
    'O-ORTH' :parserOrth,
    'O-FORMU':parserFormU,
    'D3'     :parserDoubleArray
  }

  # Check if a dictionary exists for this section
  # FIXME: Need to kill the job if a dict can't be found
  keydict = ''
  try: 
    keydict = knownKeywords[section.upper()]
  except KeyError:
    print 'No dictionary of known keywords is known for ' + section
    return

  # Loop for options for section
  for opt in opts:
    keyword = ''
    readtyp = ''

    # Check if opt is a known keyword for the section
    # FIXME: Need to kill the job if keyword is not recognized
    try: 
      keyword = knownKeywords[section.upper()][opt.upper()]
    except KeyError:
      print 'Keyword: '+str(opt)+' not recognized for section '+section
      continue
     
    readTyp = keyword.typ

    # Call the proper function for the type being read
    # FIXME: Need to kill the job if the type is not recognized
    if readTyp in ('I','D'):
      dict1[opt.upper()] = parseMap[readTyp](section,opt)
    elif readTyp in ('S'):
      dict1[opt.upper()] = parseMap[readTyp](section,opt).upper()
    elif readTyp in ('O-ORTH','O-ENV','O-FORMU','D3'):
      dict1[opt.upper()] = parseMap[readTyp](parser,section,opt)
    else:
      print 'Option data type not recognized'
  # END LOOP

  # Return the dictionary
  return dict1


def parseMisc(workers,settings): 
  print 'hello'

#
# Parse the input file
#
def parseInput(workers,iFileName):
  # Read the entire input file using ConfigParser
  inputParser = configparser.ConfigParser()
  inputParser.read(iFileName)

  secDict = {}

  # Dictionary of known sections
  # FIXME: This should be in meta somehow
  knownSections  = [ "MOLECULE", "QM", "MISC","RT"  ]

  # Loop over all of the sections parser by ConfigParse
  # and generate a dictionary for each of them to populate
  # a master dictionary for the input file
  for section in inputParser.sections():
    secStr = str(section).upper()
    # Check if section is a known section
    # FIXME: Need to kill the job if not recognized
    if secStr not in knownSections:
      print "Keyword section " + secStr + " not recognized"
      continue
    secDict[secStr] = genSecDict(inputParser,section) 
  # END LOOP

  # Return the dictionary
  return secDict

