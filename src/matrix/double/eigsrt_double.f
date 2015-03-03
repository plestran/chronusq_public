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

