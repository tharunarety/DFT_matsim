      SUBROUTINE SETCST
C   ******************************************************************
C   *                                                                *
C   *     Set constants.                                             *
C   *                                                                *
C   ******************************************************************
      USE csts
      IMPLICIT NONE
C
      PI=ACOS(-ONE)
      TWOPI=2.D0*PI
      FOURPI=4.D0*PI
      SQRTPI=SQRT(PI)
      SQFPI=DSQRT(4.D0*PI)
      SQRT2=DSQRT(2.D0)
      R2H=SQRT2/2.D0
C
      RETURN
      END

