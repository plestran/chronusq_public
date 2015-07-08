%chk=5.chk
%subst l302 .
%subst l502 .
#p hf/sto-3g iop(3/1=123456,3/33=5,4/33=5,5/33=10) GFInput nosymm scf=incore 6d scf=tight guess=read

d

0 1
c  0.0  0.1 0.0

