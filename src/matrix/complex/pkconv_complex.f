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
*Deck PkConv
      Subroutine PkConv(IOp,UPLO,N,A,LDA,AP)
C
C     DBWY Packed Storage Conversion (2014)
C     Yay!
C
      Implicit Real*8(A-H,O-Z)
C
      Logical LSAME
      External LSAME
      Character*1 UPLO
      Dimension A(LDA,*), AP(*)
C
      Do 10 I = 1,N
        If(LSAME(UPLO,'U')) then
          Do 20 J = I,N
            If(IOp.eq.1) then
              AP(I+J*(J-1)/2) = A(I,J)
            else If(IOp.eq.2) then
              A(I,J) = AP(I+J*(J-1)/2)
              If(I.ne.J) A(J,I) = AP(I+J*(J-1)/2)
              endIf
 20         Continue
        else If(LSAME(UPLO,'L')) then
          Do 30 J = 1,I
            If(IOp.eq.1) then
              AP(I+(2*N-J)*(J-1)/2) = A(I,J)
            else If(IOp.eq.2) then
              A(I,J) = AP(I+(2*N-J)*(J-1)/2)
              If(I.ne.J) A(J,I) = AP(I+(2*N-J)*(J-1)/2)
              endIf
 30         Continue
          endIf
 10     Continue
C     Call MPrint(1,N*(N+1)/2,AP)
      Return
      End
       
