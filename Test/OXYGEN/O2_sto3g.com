#p uhf/sto-3g opt=tight output=RawMatrixElement

O2 SCF

0 3
O 0.000000    0.000000    0.608586
O 0.000000    0.000000   -0.608586

O2_sto3g.matel
