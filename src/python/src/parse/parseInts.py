import os,sys
from meta.knownKeywords import requiredKeywords
import libpythonapi as chronusQ
from libpythonapi import CErrMsg

def parseInts(workers,intSettings):
  workers['CQAOIntegrals'] = chronusQ.AOIntegrals()
  if 'ALG' in intSettings:
    workers['CQAOIntegrals'].setAlgorithm(str(intSettings['ALG']))
