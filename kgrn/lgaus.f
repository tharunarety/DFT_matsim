      SUBROUTINE lgaus(a,b,n,eps,ags,wgs)
C   ******************************************************************
C   *                                                                *
C   *    Generate logarithmic mesh for use in the Gaussian integra-  *
C   *    tion close to (EF,0.) for A negative or close to (EF,A)     *
C   *    for A positive but less than EPS.                           *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      INTEGER       :: n, i
      REAL(KIND=8), DIMENSION(n) :: ags, wgs
      REAL(KIND=8)  :: a, b, eps, xa, xb, x, expx
C
      xa=LOG(1.d0-a/eps)
      xb=LOG(1.d0-b/eps)
      CALL wagaus(xa,xb,ags,wgs,n)
      DO 20 i=1,n
      x=ags(i)
      expx=EXP(x)
      wgs(i)=-eps*expx*wgs(i)
   20 ags(i)=eps*(1.d0-expx)
      RETURN
      END
