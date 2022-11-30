      SUBROUTINE setpsn(dx)
C   ******************************************************************
C   *                                                                *
C   *    Set constants for POISSON.                                  *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE poissonparam
      IMPLICIT NONE
      INTEGER :: lp
      REAL(KIND=8) :: dx, dx2, aa, b
C
      ALLOCATE(a(nl),c2(nl),f1(nl))
      dx2=dx*dx
      edl=DEXP(dx/2.d0)
      c=dx2/6.d0
      DO 20 lp=1,nl
      aa=0.25d0+(lp-1)*lp
      a(lp)=1.d0-dx2*aa/12.d0
      b=-2.d0-5.d0*dx2*aa/6.d0
      c2(lp)=-b/a(lp)
      f1(lp)=DEXP(dx*DSQRT(aa))
   20 CONTINUE
      RETURN
      END
