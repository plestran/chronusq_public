C
C The Chronus Quantum (ChronusQ) software package is high-performace 
C computational chemistry software with a strong emphasis on explictly 
C time-dependent and post-SCF quantum mechanical methods.
C 
C Copyright (C) 2014-2015 Li Research Group (University of Washington)
C 
C This program is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation; either version 2 of the License, or
C (at your option) any later version.
C 
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C 
C You should have received a copy of the GNU General Public License along
C with this program; if not, write to the Free Software Foundation, Inc.,
C 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
C 
C Contact the Developers:
C   E-Mail: xsli@uw.edu
C 
C
*Deck EigSrt
      Subroutine EigSrt(IOp,WR,WI,VR,VL,N,M)
      Implicit Real*8 (A-H,O-Z)
C
      Dimension WR(N), WI(M), VR(N,N), VL(M,M)
C
C     If(IOp.eq.1.and.N.ne.M) then
C       Write(*,*) 'Illegal Options in EigSrt'
C       Stop
C       endIf
C
      Do 10 I = 1,N-1
        Do 20 J = I+1,N
          If(WR(I).gt.WR(J)) then
            Temp = WR(I)
            WR(I) = WR(J)
            WR(J) = Temp
            Do 30 K = 1,N
              Scr = VR(K,I)
              VR(K,I) = VR(K,J)
              VR(K,J) = Scr
   30         Continue
            endIf
   20     Continue
   10   Continue
      Return
      End 

