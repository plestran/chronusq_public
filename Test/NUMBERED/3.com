%kjob l303
%mem=1000mb
%subst l302 .
#p hf/sto-3g nosymm 6d iop(3/1=123456)

d

0 1
 C                  0.00000000    0.00000000    0.00000000
 C                  1.50000000    0.00000000    0.00000000
 C                  3.00000000    0.00000000    0.00000000
 C                  4.50000000    0.00000000    0.00000000
 C                  6.00000000    0.00000000    0.00000000
 C                  7.50000000    0.00000000    0.00000000
 C                  9.00000000    0.00000000    0.00000000
 C                 10.50000000    0.00000000    0.00000000
 C                 12.00000000    0.00000000    0.00000000
 C                 13.50000000    0.00000000    0.00000000
 C                  0.00000000    1.50000000    0.00000000
 C                  0.00000000    3.00000000    0.00000000
 C                  0.00000000    4.50000000    0.00000000
 C                  0.00000000    6.00000000    0.00000000
 C                  0.00000000    7.50000000    0.00000000
 C                  0.00000000    9.00000000    0.00000000
 C                  0.00000000   10.50000000    0.00000000
 C                  0.00000000   12.00000000    0.00000000
 C                  0.00000000   13.50000000    0.00000000
 C                  1.50000000    1.50000000    0.00000000
 C                  1.50000000    3.00000000    0.00000000
 C                  1.50000000    4.50000000    0.00000000
 C                  1.50000000    6.00000000    0.00000000
 C                  1.50000000    7.50000000    0.00000000
 C                  1.50000000    9.00000000    0.00000000
 C                  1.50000000   10.50000000    0.00000000
 C                  1.50000000   12.00000000    0.00000000
 C                  1.50000000   13.50000000    0.00000000
 C                  3.00000000    1.50000000    0.00000000
 C                  3.00000000    3.00000000    0.00000000
 C                  3.00000000    4.50000000    0.00000000
 C                  3.00000000    6.00000000    0.00000000
 C                  3.00000000    7.50000000    0.00000000
 C                  3.00000000    9.00000000    0.00000000
 C                  3.00000000   10.50000000    0.00000000
 C                  3.00000000   12.00000000    0.00000000
 C                  3.00000000   13.50000000    0.00000000
 C                  4.50000000    1.50000000    0.00000000
 C                  4.50000000    3.00000000    0.00000000
 C                  4.50000000    4.50000000    0.00000000
 C                  4.50000000    6.00000000    0.00000000
 C                  4.50000000    7.50000000    0.00000000
 C                  4.50000000    9.00000000    0.00000000
 C                  4.50000000   10.50000000    0.00000000
 C                  4.50000000   12.00000000    0.00000000
 C                  4.50000000   13.50000000    0.00000000
 C                  6.00000000    1.50000000    0.00000000
 C                  6.00000000    3.00000000    0.00000000
 C                  6.00000000    4.50000000    0.00000000
 C                  6.00000000    6.00000000    0.00000000
 C                  6.00000000    7.50000000    0.00000000
 C                  6.00000000    9.00000000    0.00000000
 C                  6.00000000   10.50000000    0.00000000
 C                  6.00000000   12.00000000    0.00000000
 C                  6.00000000   13.50000000    0.00000000
 C                  7.50000000    1.50000000    0.00000000
 C                  7.50000000    3.00000000    0.00000000
 C                  7.50000000    4.50000000    0.00000000
 C                  7.50000000    6.00000000    0.00000000
 C                  7.50000000    7.50000000    0.00000000
 C                  7.50000000    9.00000000    0.00000000
 C                  7.50000000   10.50000000    0.00000000
 C                  7.50000000   12.00000000    0.00000000
 C                  7.50000000   13.50000000    0.00000000
 C                  9.00000000    1.50000000    0.00000000
 C                  9.00000000    3.00000000    0.00000000
 C                  9.00000000    4.50000000    0.00000000
 C                  9.00000000    6.00000000    0.00000000
 C                  9.00000000    7.50000000    0.00000000
 C                  9.00000000    9.00000000    0.00000000
 C                  9.00000000   10.50000000    0.00000000
 C                  9.00000000   12.00000000    0.00000000
 C                  9.00000000   13.50000000    0.00000000
 C                 10.50000000    1.50000000    0.00000000
 C                 10.50000000    3.00000000    0.00000000
 C                 10.50000000    4.50000000    0.00000000
 C                 10.50000000    6.00000000    0.00000000
 C                 10.50000000    7.50000000    0.00000000
 C                 10.50000000    9.00000000    0.00000000
 C                 10.50000000   10.50000000    0.00000000
 C                 10.50000000   12.00000000    0.00000000
 C                 10.50000000   13.50000000    0.00000000
 C                 12.00000000    1.50000000    0.00000000
 C                 12.00000000    3.00000000    0.00000000
 C                 12.00000000    4.50000000    0.00000000
 C                 12.00000000    6.00000000    0.00000000
 C                 12.00000000    7.50000000    0.00000000
 C                 12.00000000    9.00000000    0.00000000
 C                 12.00000000   10.50000000    0.00000000
 C                 12.00000000   12.00000000    0.00000000
 C                 12.00000000   13.50000000    0.00000000
 C                 13.50000000    1.50000000    0.00000000
 C                 13.50000000    3.00000000    0.00000000
 C                 13.50000000    4.50000000    0.00000000
 C                 13.50000000    6.00000000    0.00000000
 C                 13.50000000    7.50000000    0.00000000
 C                 13.50000000    9.00000000    0.00000000
 C                 13.50000000   10.50000000    0.00000000
 C                 13.50000000   12.00000000    0.00000000
 C                 13.50000000   13.50000000    0.00000000
 C                  0.00000000    0.00000000    1.50000000
 C                  0.00000000    0.00000000    3.00000000
 C                  0.00000000    0.00000000    4.50000000
 C                  0.00000000    0.00000000    6.00000000
 C                  0.00000000    0.00000000    7.50000000
 C                  0.00000000    0.00000000    9.00000000
 C                  0.00000000    0.00000000   10.50000000
 C                  0.00000000    0.00000000   12.00000000
 C                  0.00000000    0.00000000   13.50000000
 C                  1.50000000    0.00000000    1.50000000
 C                  1.50000000    0.00000000    3.00000000
 C                  1.50000000    0.00000000    4.50000000
 C                  1.50000000    0.00000000    6.00000000
 C                  1.50000000    0.00000000    7.50000000
 C                  1.50000000    0.00000000    9.00000000
 C                  1.50000000    0.00000000   10.50000000
 C                  1.50000000    0.00000000   12.00000000
 C                  1.50000000    0.00000000   13.50000000
 C                  3.00000000    0.00000000    1.50000000
 C                  3.00000000    0.00000000    3.00000000
 C                  3.00000000    0.00000000    4.50000000
 C                  3.00000000    0.00000000    6.00000000
 C                  3.00000000    0.00000000    7.50000000
 C                  3.00000000    0.00000000    9.00000000
 C                  3.00000000    0.00000000   10.50000000
 C                  3.00000000    0.00000000   12.00000000
 C                  3.00000000    0.00000000   13.50000000
 C                  4.50000000    0.00000000    1.50000000
 C                  4.50000000    0.00000000    3.00000000
 C                  4.50000000    0.00000000    4.50000000
 C                  4.50000000    0.00000000    6.00000000
 C                  4.50000000    0.00000000    7.50000000
 C                  4.50000000    0.00000000    9.00000000
 C                  4.50000000    0.00000000   10.50000000
 C                  4.50000000    0.00000000   12.00000000
 C                  4.50000000    0.00000000   13.50000000
 C                  6.00000000    0.00000000    1.50000000
 C                  6.00000000    0.00000000    3.00000000
 C                  6.00000000    0.00000000    4.50000000
 C                  6.00000000    0.00000000    6.00000000
 C                  6.00000000    0.00000000    7.50000000
 C                  6.00000000    0.00000000    9.00000000
 C                  6.00000000    0.00000000   10.50000000
 C                  6.00000000    0.00000000   12.00000000
 C                  6.00000000    0.00000000   13.50000000
 C                  7.50000000    0.00000000    1.50000000
 C                  7.50000000    0.00000000    3.00000000
 C                  7.50000000    0.00000000    4.50000000
 C                  7.50000000    0.00000000    6.00000000
 C                  7.50000000    0.00000000    7.50000000
 C                  7.50000000    0.00000000    9.00000000
 C                  7.50000000    0.00000000   10.50000000
 C                  7.50000000    0.00000000   12.00000000
 C                  7.50000000    0.00000000   13.50000000
 C                  9.00000000    0.00000000    1.50000000
 C                  9.00000000    0.00000000    3.00000000
 C                  9.00000000    0.00000000    4.50000000
 C                  9.00000000    0.00000000    6.00000000
 C                  9.00000000    0.00000000    7.50000000
 C                  9.00000000    0.00000000    9.00000000
 C                  9.00000000    0.00000000   10.50000000
 C                  9.00000000    0.00000000   12.00000000
 C                  9.00000000    0.00000000   13.50000000
 C                 10.50000000    0.00000000    1.50000000
 C                 10.50000000    0.00000000    3.00000000
 C                 10.50000000    0.00000000    4.50000000
 C                 10.50000000    0.00000000    6.00000000
 C                 10.50000000    0.00000000    7.50000000
 C                 10.50000000    0.00000000    9.00000000
 C                 10.50000000    0.00000000   10.50000000
 C                 10.50000000    0.00000000   12.00000000
 C                 10.50000000    0.00000000   13.50000000
 C                 12.00000000    0.00000000    1.50000000
 C                 12.00000000    0.00000000    3.00000000
 C                 12.00000000    0.00000000    4.50000000
 C                 12.00000000    0.00000000    6.00000000
 C                 12.00000000    0.00000000    7.50000000
 C                 12.00000000    0.00000000    9.00000000
 C                 12.00000000    0.00000000   10.50000000
 C                 12.00000000    0.00000000   12.00000000
 C                 12.00000000    0.00000000   13.50000000
 C                 13.50000000    0.00000000    1.50000000
 C                 13.50000000    0.00000000    3.00000000
 C                 13.50000000    0.00000000    4.50000000
 C                 13.50000000    0.00000000    6.00000000
 C                 13.50000000    0.00000000    7.50000000
 C                 13.50000000    0.00000000    9.00000000
 C                 13.50000000    0.00000000   10.50000000
 C                 13.50000000    0.00000000   12.00000000
 C                 13.50000000    0.00000000   13.50000000
 C                  0.00000000    1.50000000    1.50000000
 C                  0.00000000    1.50000000    3.00000000
 C                  0.00000000    1.50000000    4.50000000
 C                  0.00000000    1.50000000    6.00000000
 C                  0.00000000    1.50000000    7.50000000
 C                  0.00000000    1.50000000    9.00000000
 C                  0.00000000    1.50000000   10.50000000
 C                  0.00000000    1.50000000   12.00000000
 C                  0.00000000    1.50000000   13.50000000
 C                  0.00000000    3.00000000    1.50000000
 C                  0.00000000    3.00000000    3.00000000
 C                  0.00000000    3.00000000    4.50000000
 C                  0.00000000    3.00000000    6.00000000
 C                  0.00000000    3.00000000    7.50000000
 C                  0.00000000    3.00000000    9.00000000
 C                  0.00000000    3.00000000   10.50000000
 C                  0.00000000    3.00000000   12.00000000
 C                  0.00000000    3.00000000   13.50000000
 C                  0.00000000    4.50000000    1.50000000
 C                  0.00000000    4.50000000    3.00000000
 C                  0.00000000    4.50000000    4.50000000
 C                  0.00000000    4.50000000    6.00000000
 C                  0.00000000    4.50000000    7.50000000
 C                  0.00000000    4.50000000    9.00000000
 C                  0.00000000    4.50000000   10.50000000
 C                  0.00000000    4.50000000   12.00000000
 C                  0.00000000    4.50000000   13.50000000
 C                  0.00000000    6.00000000    1.50000000
 C                  0.00000000    6.00000000    3.00000000
 C                  0.00000000    6.00000000    4.50000000
 C                  0.00000000    6.00000000    6.00000000
 C                  0.00000000    6.00000000    7.50000000
 C                  0.00000000    6.00000000    9.00000000
 C                  0.00000000    6.00000000   10.50000000
 C                  0.00000000    6.00000000   12.00000000
 C                  0.00000000    6.00000000   13.50000000
 C                  0.00000000    7.50000000    1.50000000
 C                  0.00000000    7.50000000    3.00000000
 C                  0.00000000    7.50000000    4.50000000
 C                  0.00000000    7.50000000    6.00000000
 C                  0.00000000    7.50000000    7.50000000
 C                  0.00000000    7.50000000    9.00000000
 C                  0.00000000    7.50000000   10.50000000
 C                  0.00000000    7.50000000   12.00000000
 C                  0.00000000    7.50000000   13.50000000
 C                  0.00000000    9.00000000    1.50000000
 C                  0.00000000    9.00000000    3.00000000
 C                  0.00000000    9.00000000    4.50000000
 C                  0.00000000    9.00000000    6.00000000
 C                  0.00000000    9.00000000    7.50000000
 C                  0.00000000    9.00000000    9.00000000
 C                  0.00000000    9.00000000   10.50000000
 C                  0.00000000    9.00000000   12.00000000
 C                  0.00000000    9.00000000   13.50000000
 C                  0.00000000   10.50000000    1.50000000
 C                  0.00000000   10.50000000    3.00000000
 C                  0.00000000   10.50000000    4.50000000
 C                  0.00000000   10.50000000    6.00000000
 C                  0.00000000   10.50000000    7.50000000
 C                  0.00000000   10.50000000    9.00000000
 C                  0.00000000   10.50000000   10.50000000
 C                  0.00000000   10.50000000   12.00000000
 C                  0.00000000   10.50000000   13.50000000
 C                  0.00000000   12.00000000    1.50000000
 C                  0.00000000   12.00000000    3.00000000
 C                  0.00000000   12.00000000    4.50000000
 C                  0.00000000   12.00000000    6.00000000
 C                  0.00000000   12.00000000    7.50000000
 C                  0.00000000   12.00000000    9.00000000
 C                  0.00000000   12.00000000   10.50000000
 C                  0.00000000   12.00000000   12.00000000
 C                  0.00000000   12.00000000   13.50000000
 C                  0.00000000   13.50000000    1.50000000
 C                  0.00000000   13.50000000    3.00000000
 C                  0.00000000   13.50000000    4.50000000
 C                  0.00000000   13.50000000    6.00000000
 C                  0.00000000   13.50000000    7.50000000
 C                  0.00000000   13.50000000    9.00000000
 C                  0.00000000   13.50000000   10.50000000
 C                  0.00000000   13.50000000   12.00000000
 C                  0.00000000   13.50000000   13.50000000
 C                  1.50000000    1.50000000    1.50000000
 C                  1.50000000    1.50000000    3.00000000
 C                  1.50000000    1.50000000    4.50000000
 C                  1.50000000    1.50000000    6.00000000
 C                  1.50000000    1.50000000    7.50000000
 C                  1.50000000    1.50000000    9.00000000
 C                  1.50000000    1.50000000   10.50000000
 C                  1.50000000    1.50000000   12.00000000
 C                  1.50000000    1.50000000   13.50000000
 C                  1.50000000    3.00000000    1.50000000
 C                  1.50000000    3.00000000    3.00000000
 C                  1.50000000    3.00000000    4.50000000
 C                  1.50000000    3.00000000    6.00000000
 C                  1.50000000    3.00000000    7.50000000
 C                  1.50000000    3.00000000    9.00000000
 C                  1.50000000    3.00000000   10.50000000
 C                  1.50000000    3.00000000   12.00000000
 C                  1.50000000    3.00000000   13.50000000
 C                  1.50000000    4.50000000    1.50000000
 C                  1.50000000    4.50000000    3.00000000
 C                  1.50000000    4.50000000    4.50000000
 C                  1.50000000    4.50000000    6.00000000
 C                  1.50000000    4.50000000    7.50000000
 C                  1.50000000    4.50000000    9.00000000
 C                  1.50000000    4.50000000   10.50000000
 C                  1.50000000    4.50000000   12.00000000
 C                  1.50000000    4.50000000   13.50000000
 C                  1.50000000    6.00000000    1.50000000
 C                  1.50000000    6.00000000    3.00000000
 C                  1.50000000    6.00000000    4.50000000
 C                  1.50000000    6.00000000    6.00000000
 C                  1.50000000    6.00000000    7.50000000
 C                  1.50000000    6.00000000    9.00000000
 C                  1.50000000    6.00000000   10.50000000
 C                  1.50000000    6.00000000   12.00000000
 C                  1.50000000    6.00000000   13.50000000
 C                  1.50000000    7.50000000    1.50000000
 C                  1.50000000    7.50000000    3.00000000
 C                  1.50000000    7.50000000    4.50000000
 C                  1.50000000    7.50000000    6.00000000
 C                  1.50000000    7.50000000    7.50000000
 C                  1.50000000    7.50000000    9.00000000
 C                  1.50000000    7.50000000   10.50000000
 C                  1.50000000    7.50000000   12.00000000
 C                  1.50000000    7.50000000   13.50000000
 C                  1.50000000    9.00000000    1.50000000
 C                  1.50000000    9.00000000    3.00000000
 C                  1.50000000    9.00000000    4.50000000
 C                  1.50000000    9.00000000    6.00000000
 C                  1.50000000    9.00000000    7.50000000
 C                  1.50000000    9.00000000    9.00000000
 C                  1.50000000    9.00000000   10.50000000
 C                  1.50000000    9.00000000   12.00000000
 C                  1.50000000    9.00000000   13.50000000
 C                  1.50000000   10.50000000    1.50000000
 C                  1.50000000   10.50000000    3.00000000
 C                  1.50000000   10.50000000    4.50000000
 C                  1.50000000   10.50000000    6.00000000
 C                  1.50000000   10.50000000    7.50000000
 C                  1.50000000   10.50000000    9.00000000
 C                  1.50000000   10.50000000   10.50000000
 C                  1.50000000   10.50000000   12.00000000
 C                  1.50000000   10.50000000   13.50000000
 C                  1.50000000   12.00000000    1.50000000
 C                  1.50000000   12.00000000    3.00000000
 C                  1.50000000   12.00000000    4.50000000
 C                  1.50000000   12.00000000    6.00000000
 C                  1.50000000   12.00000000    7.50000000
 C                  1.50000000   12.00000000    9.00000000
 C                  1.50000000   12.00000000   10.50000000
 C                  1.50000000   12.00000000   12.00000000
 C                  1.50000000   12.00000000   13.50000000
 C                  1.50000000   13.50000000    1.50000000
 C                  1.50000000   13.50000000    3.00000000
 C                  1.50000000   13.50000000    4.50000000
 C                  1.50000000   13.50000000    6.00000000
 C                  1.50000000   13.50000000    7.50000000
 C                  1.50000000   13.50000000    9.00000000
 C                  1.50000000   13.50000000   10.50000000
 C                  1.50000000   13.50000000   12.00000000
 C                  1.50000000   13.50000000   13.50000000
 C                  3.00000000    1.50000000    1.50000000
 C                  3.00000000    1.50000000    3.00000000
 C                  3.00000000    1.50000000    4.50000000
 C                  3.00000000    1.50000000    6.00000000
 C                  3.00000000    1.50000000    7.50000000
 C                  3.00000000    1.50000000    9.00000000
 C                  3.00000000    1.50000000   10.50000000
 C                  3.00000000    1.50000000   12.00000000
 C                  3.00000000    1.50000000   13.50000000
 C                  3.00000000    3.00000000    1.50000000
 C                  3.00000000    3.00000000    3.00000000
 C                  3.00000000    3.00000000    4.50000000
 C                  3.00000000    3.00000000    6.00000000
 C                  3.00000000    3.00000000    7.50000000
 C                  3.00000000    3.00000000    9.00000000
 C                  3.00000000    3.00000000   10.50000000
 C                  3.00000000    3.00000000   12.00000000
 C                  3.00000000    3.00000000   13.50000000
 C                  3.00000000    4.50000000    1.50000000
 C                  3.00000000    4.50000000    3.00000000
 C                  3.00000000    4.50000000    4.50000000
 C                  3.00000000    4.50000000    6.00000000
 C                  3.00000000    4.50000000    7.50000000
 C                  3.00000000    4.50000000    9.00000000
 C                  3.00000000    4.50000000   10.50000000
 C                  3.00000000    4.50000000   12.00000000
 C                  3.00000000    4.50000000   13.50000000
 C                  3.00000000    6.00000000    1.50000000
 C                  3.00000000    6.00000000    3.00000000
 C                  3.00000000    6.00000000    4.50000000
 C                  3.00000000    6.00000000    6.00000000
 C                  3.00000000    6.00000000    7.50000000
 C                  3.00000000    6.00000000    9.00000000
 C                  3.00000000    6.00000000   10.50000000
 C                  3.00000000    6.00000000   12.00000000
 C                  3.00000000    6.00000000   13.50000000
 C                  3.00000000    7.50000000    1.50000000
 C                  3.00000000    7.50000000    3.00000000
 C                  3.00000000    7.50000000    4.50000000
 C                  3.00000000    7.50000000    6.00000000
 C                  3.00000000    7.50000000    7.50000000
 C                  3.00000000    7.50000000    9.00000000
 C                  3.00000000    7.50000000   10.50000000
 C                  3.00000000    7.50000000   12.00000000
 C                  3.00000000    7.50000000   13.50000000
 C                  3.00000000    9.00000000    1.50000000
 C                  3.00000000    9.00000000    3.00000000
 C                  3.00000000    9.00000000    4.50000000
 C                  3.00000000    9.00000000    6.00000000
 C                  3.00000000    9.00000000    7.50000000
 C                  3.00000000    9.00000000    9.00000000
 C                  3.00000000    9.00000000   10.50000000
 C                  3.00000000    9.00000000   12.00000000
 C                  3.00000000    9.00000000   13.50000000
 C                  3.00000000   10.50000000    1.50000000
 C                  3.00000000   10.50000000    3.00000000
 C                  3.00000000   10.50000000    4.50000000
 C                  3.00000000   10.50000000    6.00000000
 C                  3.00000000   10.50000000    7.50000000
 C                  3.00000000   10.50000000    9.00000000
 C                  3.00000000   10.50000000   10.50000000
 C                  3.00000000   10.50000000   12.00000000
 C                  3.00000000   10.50000000   13.50000000
 C                  3.00000000   12.00000000    1.50000000
 C                  3.00000000   12.00000000    3.00000000
 C                  3.00000000   12.00000000    4.50000000
 C                  3.00000000   12.00000000    6.00000000
 C                  3.00000000   12.00000000    7.50000000
 C                  3.00000000   12.00000000    9.00000000
 C                  3.00000000   12.00000000   10.50000000
 C                  3.00000000   12.00000000   12.00000000
 C                  3.00000000   12.00000000   13.50000000
 C                  3.00000000   13.50000000    1.50000000
 C                  3.00000000   13.50000000    3.00000000
 C                  3.00000000   13.50000000    4.50000000
 C                  3.00000000   13.50000000    6.00000000
 C                  3.00000000   13.50000000    7.50000000
 C                  3.00000000   13.50000000    9.00000000
 C                  3.00000000   13.50000000   10.50000000
 C                  3.00000000   13.50000000   12.00000000
 C                  3.00000000   13.50000000   13.50000000
 C                  4.50000000    1.50000000    1.50000000
 C                  4.50000000    1.50000000    3.00000000
 C                  4.50000000    1.50000000    4.50000000
 C                  4.50000000    1.50000000    6.00000000
 C                  4.50000000    1.50000000    7.50000000
 C                  4.50000000    1.50000000    9.00000000
 C                  4.50000000    1.50000000   10.50000000
 C                  4.50000000    1.50000000   12.00000000
 C                  4.50000000    1.50000000   13.50000000
 C                  4.50000000    3.00000000    1.50000000
 C                  4.50000000    3.00000000    3.00000000
 C                  4.50000000    3.00000000    4.50000000
 C                  4.50000000    3.00000000    6.00000000
 C                  4.50000000    3.00000000    7.50000000
 C                  4.50000000    3.00000000    9.00000000
 C                  4.50000000    3.00000000   10.50000000
 C                  4.50000000    3.00000000   12.00000000
 C                  4.50000000    3.00000000   13.50000000
 C                  4.50000000    4.50000000    1.50000000
 C                  4.50000000    4.50000000    3.00000000
 C                  4.50000000    4.50000000    4.50000000
 C                  4.50000000    4.50000000    6.00000000
 C                  4.50000000    4.50000000    7.50000000
 C                  4.50000000    4.50000000    9.00000000
 C                  4.50000000    4.50000000   10.50000000
 C                  4.50000000    4.50000000   12.00000000
 C                  4.50000000    4.50000000   13.50000000
 C                  4.50000000    6.00000000    1.50000000
 C                  4.50000000    6.00000000    3.00000000
 C                  4.50000000    6.00000000    4.50000000
 C                  4.50000000    6.00000000    6.00000000
 C                  4.50000000    6.00000000    7.50000000
 C                  4.50000000    6.00000000    9.00000000
 C                  4.50000000    6.00000000   10.50000000
 C                  4.50000000    6.00000000   12.00000000
 C                  4.50000000    6.00000000   13.50000000
 C                  4.50000000    7.50000000    1.50000000
 C                  4.50000000    7.50000000    3.00000000
 C                  4.50000000    7.50000000    4.50000000
 C                  4.50000000    7.50000000    6.00000000
 C                  4.50000000    7.50000000    7.50000000
 C                  4.50000000    7.50000000    9.00000000
 C                  4.50000000    7.50000000   10.50000000
 C                  4.50000000    7.50000000   12.00000000
 C                  4.50000000    7.50000000   13.50000000
 C                  4.50000000    9.00000000    1.50000000
 C                  4.50000000    9.00000000    3.00000000
 C                  4.50000000    9.00000000    4.50000000
 C                  4.50000000    9.00000000    6.00000000
 C                  4.50000000    9.00000000    7.50000000
 C                  4.50000000    9.00000000    9.00000000
 C                  4.50000000    9.00000000   10.50000000
 C                  4.50000000    9.00000000   12.00000000
 C                  4.50000000    9.00000000   13.50000000
 C                  4.50000000   10.50000000    1.50000000
 C                  4.50000000   10.50000000    3.00000000
 C                  4.50000000   10.50000000    4.50000000
 C                  4.50000000   10.50000000    6.00000000
 C                  4.50000000   10.50000000    7.50000000
 C                  4.50000000   10.50000000    9.00000000
 C                  4.50000000   10.50000000   10.50000000
 C                  4.50000000   10.50000000   12.00000000
 C                  4.50000000   10.50000000   13.50000000
 C                  4.50000000   12.00000000    1.50000000
 C                  4.50000000   12.00000000    3.00000000
 C                  4.50000000   12.00000000    4.50000000
 C                  4.50000000   12.00000000    6.00000000

C 0
S    3 1.00
 0.7161683735e+02  0.1543289673e+00
 0.1304509632e+02  0.5353281423e+00
 0.3530512160e+01  0.4446345422e+00
SP   3 1.00
 0.2941249355e+01 -0.9996722919e-01  0.1559162750e+00
 0.6834830964e+00  0.3995128261e+00  0.6076837186e+00
 0.2222899159e+00  0.7001154689e+00  0.3919573931e+00
D   3  1.00
       30.8534100         0.919990500E-01
       8.26498500         0.398502100
       2.49533200         0.691789700
****
