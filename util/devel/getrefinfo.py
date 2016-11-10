import numpy as np
import h5py

import sys,os

fname = sys.argv[1]
dataset = sys.argv[2]
typ = sys.argv[3]

with h5py.File(fname,'r') as hf:
  if typ == 'd':
    if "Output" in dataset:
      tmp = np.array(hf.get(dataset))
      for x in tmp:
        try:
          sys.stdout.write(unichr(x))
        except ValueError:
          sys.stdout.write(' ')
    else:
      print np.array(hf.get(dataset))
  elif typ == 'g':
    gp = hf.get(dataset)
    for itm in gp.items():
      print itm[0] + ":"
      if "Output" in itm[0]:
        tmp = np.array(gp.get(itm[0]))
        for x in tmp:
          try:
            sys.stdout.write(unichr(x))
          except ValueError:
            sys.stdout.write(' ')
      else:
        print np.array(gp.get(itm[0]))
      print '\n\n'
    
