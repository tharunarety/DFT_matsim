      SUBROUTINE SIMPN(F,D,N,FI)
C     ***************************************************************
C     *                                                             *
C     *    Integrates F from X(1) to X(n) where X is Louck's mesh.  *
C     *                                                             *
C     ***************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION F(N)
  100 FORMAT(' SIMPN:**  N =',I4,'. Must be even')
      N2=N-1
      IF(MOD(N,2).EQ.0) THEN
         WRITE(6,100) N
         STOP
      ENDIF
      FF=0.D0
      DO 20 J=2,N2,2
      FF=FF+(F(J-1)+4.D0*F(J)+F(J+1))
   20 CONTINUE
      FI=D*FF/3.D0
      RETURN
      END
