      SUBROUTINE madl3m
C   ******************************************************************
C   *                                                                *
C   *    Transform the Madelung matrix calculated in SMTRX by the    *
C   *    Ewald technique fron the spherical harmonic to the real     *
C   *    harmonic representation. (For L =1,..,NL-1 )                *
C   *                                                                *
C   *    NPRN  : Print vectors if = 1 or 3.                          *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE control_data
      USE madelung_matrix
      USE message
      USE s_matrix
      IMPLICIT NONE
      COMPLEX(KIND=8) :: sum, spherm, fty, cilpl, cim 
      INTEGER :: lmax, iq, jq
      INTEGER :: iqlm, lp, nu, lpnu, lqnu, mp, lpmp, lqmp,
     .           jqlm,  l, mu,  lmu, lqmu,  m,   lm,  lqm
C
      cim=CMPLX(0.d0,1.d0,8)
      lmax=nl-1
C
      DO 10 iq=1,nq
      iqlm=(iq-1)*nlm
      DO 10 jq=1,nq
      jqlm=(jq-1)*nlm
C
C     RMad(l'nu;lmu) = sum(m';m){a(nu;m')*Mad(l'm';lm)*a(m;mu)
C
      DO 20 lp=0,lmax
      DO 20 nu=-lp,lp
      lpnu=lp*lp+lp+nu+1
      lqnu=iqlm+lpnu
      DO 20 l=0,lmax
      cilpl=((-1.d0)**l)*cim**(lp+l)
      DO 20 mu=-l,l
      lmu=l*l+l+mu+1
      lqmu=jqlm+lmu
C
      sum=CMPLX(0.d0,0.d0,8)
C
      DO 21 mp=-lp,lp
      lpmp=lp*lp+lp+mp+1
      lqmp=iqlm+lpmp
      DO 21 m=-l,l
      lm=l*l+l+m+1
      lqm=jqlm+lm
      spherm=CMPLX(srlmq(lqmp,lqm),silmq(lqmp,lqm),8)
   21 sum=sum+CONJG(fty(nu,mp))*spherm*fty(mu,m)
      sum=cilpl*sum
C
      IF(ABS(AIMAG(sum)).GT.1.d-8) THEN
         WRITE(m6,'(''SUM must bee real '',2I6,2E15.6)')
     .             lpnu,lmu,sum
c         STOP
      ENDIF
   20 vmad(lqnu,lqmu)=-REAL(sum,8)
   10 CONTINUE
C
      RETURN
      END
