      SUBROUTINE elestr
C   ******************************************************************
C   *                                                                *
C   *    Calculate the electrostatic force components from the       *
C   *    electrostatic stress field.                                 *
C   *                                                                *
C   *    The calculation is done for r = ri(jsi), which is the       *
C   *    closest mesh point to the inscribed sphere Si >~ ri(jsi).   *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    Fy,Fz,Fx: the m=-1,0,1 components of the force.             *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE csts    ; USE density
      USE force      ; USE gaussi       ; USE message ; USE moments
      USE radialmesh ; USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: denl
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: vm, vmd, vi, vid
      REAL(KIND=8) :: teta, cost, sint, phi, cosf, sinf, ux, uy, uz
      REAL(KIND=8) :: si, sisl, r, qint, pint, smm, fl
      REAL(KIND=8) :: weight, srr, srt, srp, erp, etp, epp
      INTEGER :: iq, it, ita, jri, jws, jsi, l, m, lm, ifi, ith, ir, jr
      INTEGER :: lmaxhp, nlmhp, nlmd
C
C     Loop for sites
C
      WRITE(m6,100)
      IF(msgl.NE.0) WRITE(msgio,100)
C
      lmaxhp=lmaxh+1
      nlmhp=(lmaxhp+1)*(lmaxhp+1)
C
      ALLOCATE(vm(nlmh),vmd(nlmh),denl(nlmh,dimr))
      ALLOCATE(vi(nlmh),vid(nlmh))
      ALLOCATE(ylm(nlmhp),gt(nlmh),gp(nlmh))
C
      nlmd=MIN0(nlmh,nlmmad)
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      jws=jwss(ita,it)
      jri=jws+2
      jsi=jwsi(ita,it)-shf
      si=ri(jsi,ita,it)
C
C     Madelung contribution
C
      DO 21 l=0,lmaxh
      sisl=(si/sws)**l
      sisl=sqfpi*sisl/(2.d0*l+1.d0)
      DO 21 m=-l,l
      lm=l*l+l+m+1
      IF(lm.LE.nlmd) THEN
         smm=SUM(vmad(iq,1:nq,lm,1:nlmd)*qlmn(1:nq,1:nlmd))
      ELSE
         smm=0.d0
      ENDIF
      vm(lm)=-sisl*smm/(2.d0*sws)
   21 vmd(lm)=l*vm(lm)/si
C
C     Intracell contribution
C
      DO 22 ir=1,jri
      r=ri(ir,ita,it)
   22 denl(1,ir)=SUM(chde(ita,it,1:ns,ir))*r/sqfpi
      DO 23 l=1,lmaxh
      DO 23 m=-l,l
      lm=l*l+l+m+1
      DO 23 ir=1,jri
      r=ri(ir,ita,it)
   23 denl(lm,ir)=SUM(chdl(ita,iq,1:ns,lm,ir))*r
C
      DO 24 l=0,lmaxh
      fl=4.d0*pi/(2.d0*l+1.d0)
      sisl=si**l
      DO 24 m=-l,l
      lm=l*l+l+m+1
      DO 25 ir=1,jsi
      r=ri(ir,ita,it)
   25 fi(ir)=denl(lm,ir)*(r**l)
      CALL gensim(fi,dx,jsi,qint)
      jr=0
      DO 26 ir=jsi,jri
      jr=jr+1
      r=ri(ir,ita,it)
   26 fi(jr)=denl(lm,ir)/(r**(l+1.))
      jr=jws-jsi+1
      CALL gensim(fi,dx,jr,pint)
      vi(lm)=-fl*(qint/sisl/si+pint*sisl)
      vid(lm)=-fl*((-l-1.d0)*qint/sisl/si+l*pint*sisl)/si
      IF(l.EQ.0) THEN
         vi(lm)=vi(lm)+sqfpi*nz(ita,it)/si
         vid(lm)=vid(lm)-sqfpi*nz(ita,it)/si/si
      ENDIF
   24 CONTINUE
      vi=vi+vm
      vid=vid+vmd
C
C     Integrate on the inscribed sphere (si >~ ri(jsi))
C
      DO 30 ith=1,nth
      teta    =thvec(ith)
      cost    =DCOS(teta)
      sint    =DSIN(teta)
      DO 30 ifi=1,nfi
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
      erp=-SUM(vid(1:nlmh)*ylm(1:nlmh))
      etp=-SUM(vi(1:nlmh)*gt(1:nlmh))/si
      epp=-SUM(vi(1:nlmh)*gp(1:nlmh))/si
C
      srr=(erp*erp-etp*etp-epp*epp)/2.d0/fourpi
      srt=erp*etp/fourpi
      srp=erp*epp/fourpi
C
      fey(iq)=fey(iq)+weight*(srr*uy+srt*sinf*cost+srp*cosf)
      fez(iq)=fez(iq)+weight*(srr*uz-srt*sint)
      fex(iq)=fex(iq)+weight*(srr*ux+srt*cosf*cost-srp*sinf)
C
   30 CONTINUE
C
      fey(iq)=2.d0*fey(iq)*si*si
      fez(iq)=2.d0*fez(iq)*si*si
      fex(iq)=2.d0*fex(iq)*si*si
C
   20 CONTINUE
      DEALLOCATE(ylm,gt,gp)
      DEALLOCATE(vm,vmd,vi,vid,denl)
C
  100 FORMAT(/,' ELESTR: Set up the electrostatic stress tensor')
      RETURN
      END
