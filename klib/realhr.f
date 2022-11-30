      SUBROUTINE realhr(ylm,lmax,x,y,z)
C   ******************************************************************
C   *                                                                *
C   *    Calculate real harmonics. The realhinit has to be called    *
C   *    first to calculate normalization factors ylnorm.            *
C   *                                                                *
C   ******************************************************************
      USE realhr_norm
      IMPLICIT NONE
      INTEGER :: lmax, l, m, lm
      REAL(KIND=8) :: x, y, z, rxy, r, rinv, sig, yr, yi
      REAL(KIND=8) :: sinth, costh, sinfi, cosfi, sinmfi, cosmfi 
      REAL(KIND=8) :: plm, plm1m, plp1m, plp1mp1
      REAL(KIND=8), DIMENSION(*) :: ylm
C
      if(lmax.gt.lmaxi) write(*,*) 'REALHR:!!! ',lmax,lmaxi
C
      IF(lmax.LT.0) RETURN
      rxy=DSQRT(x*x+y*y)
      r=DSQRT(x*x+y*y+z*z)
      IF(r.LT.small) THEN
         costh=1.d0
         sinth=0.d0
      ELSE
         rinv=1.d0/r
         costh=z*rinv
         sinth=rxy*rinv
      ENDIF
      IF(rxy.LT.small) THEN
         cosfi=1.d0
         sinfi=0.d0
      ELSE
         rinv=1.d0/rxy
         cosfi=x*rinv
         sinfi=y*rinv
      ENDIF
C
      cosmfi=1.d0
      sinmfi=0.d0
      plm=1.d0
      sig=1.d0
C
C     m = 0 case
C
      plp1m=plm*costh
      plp1mp1=-plm*sinth
      yr=plm*ylnorm(0,0)*cosmfi
      ylm(1)= yr
C
      IF(lmax.eq.0) RETURN
C
      yr=plp1m*ylnorm(0,1)*cosmfi
      ylm(3)= yr
      DO 11 l=1,lmax-1
      plm1m=plm
      plm=plp1m
      plp1m=(tlp1(l)*plm*costh-al(l)*plm1m)*alinv(l+1)
      yr=plp1m*ylnorm(0,l+1)*cosmfi
      lm=(l+1)*(l+1)+l+1+1
      ylm(lm)= yr
   11 CONTINUE
      plm=plp1mp1
      sig=-sig
      r=cosmfi*cosfi-sinmfi*sinfi
      sinmfi=sinmfi*cosfi+cosmfi*sinfi
      cosmfi=r
C
C     m > 0 case
C
      DO 20 m=1,lmax-1
      plp1m=tlp1(m)*plm*costh
      plp1mp1=-tlp1(m)*plm*sinth
      yr=plm*ylnorm(m,m)*cosmfi
      yi=plm*ylnorm(m,m)*sinmfi
      lm=m*m+m+1
      ylm(lm-m)=yi*sqr2*sig
      ylm(lm+m)=yr*sqr2*sig
      yr=plp1m*ylnorm(m,m+1)*cosmfi
      yi=plp1m*ylnorm(m,m+1)*sinmfi
      lm=(m+1)*(m+1)+m+1+1
      ylm(lm-m)=yi*sqr2*sig
      ylm(lm+m)=yr*sqr2*sig
      DO 21 l=m+1,lmax-1
      plm1m=plm
      plm=plp1m
      plp1m=(tlp1(l)*plm*costh-al(l+m)*plm1m)*alinv(l-m+1)
      yr=plp1m*ylnorm(m,l+1)*cosmfi
      yi=plp1m*ylnorm(m,l+1)*sinmfi
      lm=(l+1)*(l+1)+l+1+1
      ylm(lm-m)=yi*sqr2*sig
      ylm(lm+m)=yr*sqr2*sig
   21 CONTINUE
      plm=plp1mp1
      sig=-sig
      r=cosmfi*cosfi-sinmfi*sinfi
      sinmfi=sinmfi*cosfi+cosmfi*sinfi
      cosmfi=r
   20 CONTINUE
C
      yr=plm*ylnorm(lmax,lmax)*cosmfi
      yi=plm*ylnorm(lmax,lmax)*sinmfi
C
      lm=lmax*lmax+lmax+1
      ylm(lm-lmax)=yi*sqr2*sig
      ylm(lm+lmax)=yr*sqr2*sig
C
      RETURN
      END
