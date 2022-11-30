      SUBROUTINE xchstr
C   ******************************************************************
C   *                                                                *
C   *    Calculate the exchange-correlation force components         *
C   *    from the stress field sigma_xc = n(eps_xc-mu_xc).           *
C   *                                                                *
C   *    The calculation is done for r = ri(jsi), which is the       *
C   *    closest mesh point to the inscribed sphere Si >~ ri(jsi).   *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    Fy,Fz,Fx: the m=-1,0,1 components of the force.             *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text ; USE csts
      USE density ; USE force ; USE gaussi ; USE message
      USE radialmesh   ; USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: denl
      REAL(KIND=8) :: tol = 1.d-16
      REAL(KIND=8) :: teta, cost, sint, phi, cosf, sinf, ux, uy, uz
      REAL(KIND=8) :: si, rho1, rho2, rho, exc, v1, v2
      REAL(KIND=8) :: weight, giz
      INTEGER :: iq, it, ita, jsi, l, m, lm, ifi, ith, is
      INTEGER :: lmaxhp, nlmhp
C
C     Loop for sites
C
      IF(ixc.GE.8) THEN
         WRITE(m6,100) ixc
         STOP
      ENDIF
      WRITE(m6,110)
      IF(msgl.NE.0) WRITE(msgio,110)
C
      lmaxhp=lmaxh+1
      nlmhp=(lmaxhp+1)*(lmaxhp+1)
C
      ALLOCATE(denl(nlmh,2))
      ALLOCATE(ylm(nlmhp))
C
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
C
C     Set up the average Si
C
      jsi=jwsi(ita,it)-shf
      si=ri(jsi,ita,it)
C
      DO 21 is=1,ns
      DO 21 l=0,lmaxh
      DO 21 m=-l,l
      lm=l*l+l+m+1
   21 denl(lm,is)=chdl(ita,iq,is,lm,jsi)/si/si
      IF(ns.EQ.1) THEN
         denl(1:nlmh,1)=denl(1:nlmh,1)/2.d0
         denl(1:nlmh,2)=denl(1:nlmh,1)
      ENDIF
C
C     Integrate on the inscribed sphere (si >~ ri(jsi))
C
      DO 22 ith=1,nth
      teta    =thvec(ith)
      cost    =DCOS(teta)
      sint    =DSIN(teta)
      DO 22 ifi=1,nfi
      phi     =fivec(ifi)
      cosf    =DCOS(phi)
      sinf    =DSIN(phi)
      ux      =cosf*sint
      uy      =sinf*sint
      uz      =cost
      CALL realhr(ylm,lmaxhp,ux,uy,uz)
      weight=wth(ith)*sint*wfi(ifi)
C
C     Average over the sites
C
      weight=conc(ita,it)*weight
C
C     Set up the exchange-correlation stress tensor sigma(r) for r(ux,uy,uz)
C     (only (r,r), (r,theta) and (r,phi) components)
C
      rho1=SUM(denl(1:nlmh,1)*ylm(1:nlmh))
      rho2=SUM(denl(1:nlmh,2)*ylm(1:nlmh))
      IF(rho1.LT.tol) rho1=0.d0
      IF(rho2.LT.tol) rho2=0.d0
      rho=rho1+rho2
      IF(rho.LT.tol) GO TO 22
C
      rhod(1:2)=0.d0
      rhodd(1:2)=0.d0
      CALL xcpot(ixc,rho1,rho2,rho,rhod,rhodd,si,v1,v2,exc)
      giz=rho*exc-rho1*v1-rho2*v2
      fxx(iq)=fxx(iq)+weight*giz*ux
      fxy(iq)=fxy(iq)+weight*giz*uy
      fxz(iq)=fxz(iq)+weight*giz*uz
C
   22 CONTINUE
C
      fxy(iq)=fxy(iq)*si*si
      fxz(iq)=fxz(iq)*si*si
      fxx(iq)=fxx(iq)*si*si
c
   20 CONTINUE
      DEALLOCATE(ylm,denl)
C
  100 FORMAT(/,' XCHSTR: IXC =',i3,' is not implemented')
  110 FORMAT(/,' XCHSTR: Set up the exchange-correlation stress tensor')
      RETURN
      END
