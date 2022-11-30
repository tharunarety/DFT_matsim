      SUBROUTINE setcst
C   ******************************************************************
C   *                                                                *
C   *     Set constants.                                             *
C   *                                                                *
C   ******************************************************************
      USE csts
      IMPLICIT NONE
C
      pi=ACOS(-one)
      sqrtpi=SQRT(pi)
      twopi=2.d0*pi
      fourpi=4.d0*pi
      sqrt2=DSQRT(2.d0)
      r2h=SQRT2/2.d0
      sqfpi=DSQRT(4.d0*pi)
C
      RETURN
      END

C

