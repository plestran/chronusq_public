import sys,os
import getopt
from refval import *
from chronusq import *
from tabulate import tabulate

# Load the reference values into ref
refvalues()
print ref[0].scf

fname = "test0001_serial_incore"
tests = [[None for x in xrange(100000)] for y in xrange(10)]
tests[0][0] = runCQ(fname)
print tests[0][0].E



