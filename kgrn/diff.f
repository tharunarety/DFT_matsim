      SUBROUTINE diff(dx,n,it,ita,is)
C   ******************************************************************
C   *                                                                *
C   *   Differentiate the potential for use in CDIRAC.               *
C   *                                                                *
C   ******************************************************************
      USE potential
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: a1 = 6.6666666667d0 
      REAL(KIND=8), PARAMETER :: a2 = 0.16666666667d0
      REAL(KIND=8) :: dx
      INTEGER      :: n, it, ita, is, nm2, i
      NM2=N-2
      VP(1,ita,IT,IS)=((6.D0*V(2,ita,IT,IS)+a1*V(4,ita,IT,IS)
     .           +1.2D0*V(6,ita,IT,IS))-(2.45D0*V(1,ita,IT,IS)
     .           +7.5D0*V(3,ita,IT,IS)+3.75D0*V(5,ita,IT,IS)
     .           +a2*V(7,ita,IT,IS)))/DX
      VP(2,ita,IT,IS)=((6.D0*V(3,ita,IT,IS)+a1*V(5,ita,IT,IS)
     .           +1.2D0*V(7,ita,IT,IS))-(2.45D0*V(2,ita,IT,IS)
     .           +7.5D0*V(4,ita,IT,IS)+3.75D0*V(6,ita,IT,IS)
     .           +a2*V(8,ita,IT,IS)))/DX
      DO 20 I=3,NM2
   20 VP(I,ita,IT,IS)=((V(I-2,ita,IT,IS)+8.D0*V(I+1,ita,IT,IS))
     .           -(8.D0*V(I-1,ita,IT,IS)+V(I+2,ita,IT,IS)))/12.D0/DX
      VP(N-1,ita,IT,IS)=(V(N,ita,IT,IS)-V(N-2,ita,IT,IS))/2.D0/DX
      VP(N,ita,IT,IS)=(V(N-2,ita,IT,IS)/2.D0-2.D0*V(N-1,ita,IT,IS)
     .          +3.D0/2.D0*V(N,ita,IT,IS))/DX
      RETURN
      END
