#!/bin/sh
make -j8
./src/numdiff/numdiff     `# Test` \
  0.000 0.002 0.001       `# Scan X parameters` \
  1.617436 1.619436 0.001 `# Scan Y parameters` \
  cc-pVDZ                 `# Basis Set` \
  CIS                     `# Response Type` \
  0                       `# Charge` \
  4                       `# NStates` \
  3 6                     `# Precision X / Y` \
  LiH                     `# File Prefix` \
  > test6                 `# Output File`
