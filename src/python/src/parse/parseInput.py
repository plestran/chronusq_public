#
# This file is part of the Chronus Quantum (ChronusQ) software package
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
import os,sys
import configparser
import libpythonapi as chronusQ
from libpythonapi import CErrMsg
from libpythonapi import CErr
from meta.knownKeywords import knownKeywords
from meta.knownKeywords import knownSections
#import meta.enumMaps as enumMaps

#### THIS IS OLD RT CODE, ADD BACK IF NEEDED IN FUTURE
##
## Convert string obtained from ConfigParser
## for RT.ortho into RT.ORTHO enum
##
#def parserOrth(parser,section,opt):
#  strX = parser.get(section,opt)
#  strX = str(strX).upper()
#  return enumMaps.orthoMap[strX]
#
##
## Convert string obtained from ConfigParser
## for RT.ell_pol into RT.ELL_POL enum
##
#def parserEllPol(parser,section,opt):
#  strX = parser.get(section,opt)
#  strX = str(strX).upper()
#  return enumMaps.ellPolMap[strX]
#
##
## Convert string obtained from ConfigParser
## for RT.envelope into RT.ENVELOPE enum
##
#def parserEnv(parser,section,opt):
#  strX = parser.get(section,opt)
#  strX = str(strX).upper()
#  return enumMaps.envMap[strX]
#
##
## Convert string obtained from ConfigParser
## for RT.uprop into RT.FORM_U enum
##
#def parserFormU(parser,section,opt):
#  strX = parser.get(section,opt)
#  strX = str(strX).upper()
#  return enumMaps.formUMap[strX]

#
# Convert string obtained from ConfigParser
# for SCF.guess into Guess enum
#
#def parserGuess(parser,section,opt):
#  strX = parser.get(section,opt)
#  strX = str(strX).upper()
#  return enumMaps.guessMap[strX]

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
def genSecDict(workers,parser,section):

  # Initialize an empy dictionary
  dict1 = {}

  # Grab options for given section
  opts  = parser.options(section)

  # Define a map from data type to parser function
  # (Includes user defined wrappers for special types)
  parseMap = { 
    'I'      :parser.getint,
    'D'      :parser.getfloat,
    'B'      :parser.getboolean,
    'S'      :parser.get,
#    'O-ENV'  :parserEnv,
#    'O-ORTH' :parserOrth,
#    'O-FORMU':parserFormU,
#    'O-GS'   :parserGuess,
    'D3'     :parserDoubleArray
#    'O-ELL'  :parserEllPol
  }

  # Check if a dictionary exists for this section
  keydict = ''
  try: 
    keydict = knownKeywords[section.upper()]
  except KeyError:
    msg = 'No dictionary of known keywords is known for ' + section
    CErrMsg(workers['CQFileIO'],str(msg))
    return

  # Loop for options for section
  for opt in opts:
    keyword = ''
    readtyp = ''

    # Check if opt is a known keyword for the section
    try: 
      keyword = knownKeywords[section.upper()][opt.upper()]
    except KeyError:
      msg = 'Keyword: ' + str(opt) + ' not recognized for section ' + section
      CErrMsg(workers['CQFileIO'],str(msg))
      continue
     
    readTyp = keyword.typ

    # Call the proper function for the type being read
    if readTyp in ('I','D','B'):
      dict1[opt.upper()] = parseMap[readTyp](section,opt)
    elif readTyp in ('S'):
      dict1[opt.upper()] = str(parseMap[readTyp](section,opt).upper())
    elif readTyp in ('O-ORTH','O-ENV','O-FORMU','O-GS','D3','O-ELL'):
      dict1[opt.upper()] = parseMap[readTyp](parser,section,opt)
    else:
      msg = 'Option data ' + str(readTyp) + ' type not recognized'
      CErrMsg(workers['CQFileIO'],str(msg))
  # END LOOP

  # Return the dictionary
  return dict1


#
# Parse the input file
#
def parseInput(workers,iFileName):
  # Read the entire input file using ConfigParser
  inputParser = configparser.ConfigParser()
  inputParser.read(iFileName)

  secDict = {}

  # Loop over all of the sections parser by ConfigParse
  # and generate a dictionary for each of them to populate
  # a master dictionary for the input file
  for section in inputParser.sections():
    secStr = str(section).upper()
    # Check if section is a known section
    if secStr not in knownSections:
      msg = "Keyword section " + secStr + " not recognized"
      CErrMsg(workers['CQFileIO'],str(msg))
      continue
    secDict[secStr] = genSecDict(workers,inputParser,section) 
  # END LOOP

  # Return the dictionary
  return secDict

