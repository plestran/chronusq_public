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
       
