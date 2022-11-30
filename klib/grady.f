      SUBROUTINE grady(ylm,nlmhp,lmaxh,uz,gt,gp,nlmh)
C   ******************************************************************
C   *                                                                *
C   * Calculate the grad_theta Ylm and grad_phi Ylm .                *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      INTEGER :: lmaxh, nlmh, nlmhp, l, m, lm, lp, lmp, lmm
      REAL(KIND=8), DIMENSION(nlmhp) :: ylm
      REAL(KIND=8), DIMENSION(nlmh)  :: gt, gp
      REAL(KIND=8)  :: fl, facl, uz
      REAL(KIND=8), INTRINSIC :: SQRT
C
      DO 20 l =0,lmaxh
      fl      =(2.d0*l+1.d0)/(2.d0*l+3.d0)
      DO 20 m =-l,l
      lm      =l*l+l+m+1
      lp      =l+1
      facl    =SQRT(fl*(lp*lp-m*m))
      lmp     =lp*lp+lp+m+1
      gt(lm)  =facl*ylm(lmp)-lp*uz*ylm(lm)
      lmm     =l*l+l-m+1
      gp(lm)  =-m*ylm(lmm)
   20 CONTINUE
C
      RETURN
      END
