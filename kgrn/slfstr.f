      SUBROUTINE slfstr
C   ******************************************************************
C   *                                                                *
C   *    Calculate the self-interaction force components             *
C   *    from the stress field.                                      *
C   *                                                                *
C   *    The calculation is done for r = ri(jsi), which is the       *
C   *    closest mesh point to the inscribed sphere Si >~ ri(jsi).   *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    Fy,Fz,Fx: the m=-1,0,1 components of the force.             *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE csts ; USE force
      USE gaussi ; USE lattice ; USE message ; USE radialmesh
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: zpot, zpotd
      REAL(KIND=8), DIMENSION(nv,4) :: rsl
      REAL(KIND=8), DIMENSION(nv)   :: zsl
      REAL(KIND=8) :: teta, cost, sint, phi, cosf, sinf, ux, uy, uz
      REAL(KIND=8) :: si, sisl, rl, fl, ziv
      REAL(KIND=8) :: weight, srr, srt, srp, erp, etp, epp
      INTEGER, DIMENSION(nv) :: iqsl
      INTEGER :: iq, it, ita, jsi, l, m, lm, ifi, ith, iv
      INTEGER :: lmaxhp, nlmhp, nvsl
C
C     Loop for sites
C
      WRITE(m6,100)
      IF(msgl.NE.0) WRITE(msgio,100)
C
      lmaxhp=lmaxh+1
      nlmhp=(lmaxhp+1)*(lmaxhp+1)
C
      ALLOCATE(ylm(nlmhp),gt(nlmh),gp(nlmh))
C
      DO 20 iq=1,nq
      it=itq(iq)
C
C     Set up the average Si
C
      si=0.d0
      DO ita=1,nta(it)
      jsi=jwsi(ita,it)-shf
      si=si+conc(ita,it)*ri(jsi,ita,it)
      ENDDO
C
      CALL slfcls(iq,nvsl,zsl,rsl)
      ALLOCATE(zpot(nlmh,nvsl),zpotd(nlmh,nvsl))
C
      DO 21 iv=1,nvsl
      ux=rsl(iv,1)
      uy=rsl(iv,2)
      uz=rsl(iv,3)
      rl=rsl(iv,4)
      CALL realhr(ylm,lmaxh,ux,uy,uz)
C
      DO 22 l=0,lmaxh
      fl=4.d0*pi/(2.d0*l+1.d0)
      sisl=(si/rl)**l
      DO 22 m=-l,l
      lm=l*l+l+m+1
C
      zpot(lm,iv)=fl*sisl/rl*ylm(lm)
   22 zpotd(lm,iv)=l*zpot(lm,iv)/si
   21 CONTINUE
C
C     Integrate on the inscribed sphere (si >~ ri(jsi))
C
      DO 23 ith=1,nth
      teta    =thvec(ith)
      cost    =DCOS(teta)
      sint    =DSIN(teta)
      DO 23 ifi=1,nfi
      phi     =fivec(ifi)
      cosf    =DCOS(phi)
      sinf    =DSIN(phi)
      ux      =cosf*sint
      uy      =sinf*sint
      uz      =cost
      CALL realhr(ylm,lmaxhp,ux,uy,uz)
      CALL grady(ylm,nlmhp,lmaxh,uz,gt,gp,nlmh)
      gt=gt/sint
      gp=gp/sint
      weight=wth(ith)*sint*wfi(ifi)
C
C     Set up the electrostatic stress tensor sigma(r) for r(ux,uy,uz)
C     (only (r,r), (r,theta) and (r,phi) components)
C
      DO 24 iv=1,nvsl
      ziv=zsl(iv)*zsl(iv)
      erp=-SUM(zpotd(1:nlmh,iv)*ylm(1:nlmh))
      etp=-SUM(zpot(1:nlmh,iv)*gt(1:nlmh))/si
      epp=-SUM(zpot(1:nlmh,iv)*gp(1:nlmh))/si
C
      srr=ziv*(erp*erp-etp*etp-epp*epp)/2.d0/fourpi
      srt=ziv*erp*etp/fourpi
      srp=ziv*erp*epp/fourpi
C
      fiy(iq)=fiy(iq)+weight*(srr*uy+srt*sinf*cost+srp*cosf)
      fiz(iq)=fiz(iq)+weight*(srr*uz-srt*sint)
      fix(iq)=fix(iq)+weight*(srr*ux+srt*cosf*cost-srp*sinf)
C
   24 CONTINUE
   23 CONTINUE
C
      fiy(iq)=2.d0*fiy(iq)*si*si
      fiz(iq)=2.d0*fiz(iq)*si*si
      fix(iq)=2.d0*fix(iq)*si*si
C
      DEALLOCATE(zpot,zpotd)
   20 CONTINUE
      DEALLOCATE(ylm,gt,gp)
C
  100 FORMAT(/,' SLFSTR: Set up the self-interaction stress tensor')
      RETURN
      END
