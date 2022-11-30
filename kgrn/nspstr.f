      SUBROUTINE nspstr
C   ******************************************************************
C   *                                                                *
C   *    Calculate the non-spherical force components from the       *
C   *    non-spherical potential.                                    *
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
      USE pota       ; USE radialmesh   ; USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: denl, denld
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: vnsxc
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: vnsh, denl0
      REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dvr, dvt
      REAL(KIND=8) :: teta, cost, sint, phi, cosf, sinf, ux, uy, uz
      REAL(KIND=8) :: si, sisl, r, rp, qint, pint, smm, fl, fint
      REAL(KIND=8) :: weight, snn, srr, srt, srp
      REAL(KIND=8) :: rho0, exc, excd, excdd, excddd, muxcd, facns
      INTEGER :: iq, it, ita, jri, jws, jsi, ifi, ith, ir, jr, kr
      INTEGER :: lmaxhp, nlmhp, is, lmp, lm, l, m
C
      IF(ixc.GE.8.OR.ns.EQ.2) THEN
         WRITE(m6,*) 'Warning !!! in NSPSTR'
         RETURN
      ENDIF
      facns=spinfc/2.d0
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
      ALLOCATE(denl(ns,nlmh,dimr),denld(ns,nlmh,dimr))
      ALLOCATE(vnsh(2:nlmh,dimr),vnsxc(ns,2:nlmh,dimr))
      ALLOCATE(dvr(nlmh,nlmh),dvt(nlmh,nlmh),denl0(ns,dimr))
C
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      jws=jwss(ita,it)
      jri=jws+2
      jsi=jwsi(ita,it)-shf
      si=ri(jsi,ita,it)
C
C     Madelung contribution to the non-spherical potential
C
      DO 21 l=1,lmaxh
      DO 21 m=-l,l
      lm=l*l+l+m+1
      IF(lm.LE.diml) THEN
         smm=SUM(vmad(iq,1:nq,lm,1:diml)*qlmn(1:nq,1:diml))
      ELSE
         smm=0.d0
      ENDIF
      DO 21 ir=1,jws
      r=ri(ir,ita,it)
      sisl=(r/sws)**l
      sisl=sqfpi*sisl/(2.d0*l+1.d0)
   21 vnsh(lm,ir)=-sisl*smm/(2.d0*sws)
C
C     Intracell contribution to the non-spherical potential
C
      DO 22 l=1,lmaxh
      fl=4.d0*pi/(2.d0*l+1.d0)
      DO 22 m=-l,l
      lm=l*l+l+m+1
C
      DO 23 ir=1,jws
      r=ri(ir,ita,it)
      sisl=r**l
C
      DO 24 jr=1,ir+1
      rp=ri(jr,ita,it)**(l+1.)
   24 fi(jr)=SUM(chdl(ita,iq,1:ns,lm,jr))*rp
      CALL gensim(fi,dx,ir,qint)
C
      kr=0
      DO 25 jr=ir,jri
      kr=kr+1
      rp=ri(jr,ita,it)**l
   25 fi(kr)=SUM(chdl(ita,iq,1:ns,lm,jr))/rp
      kr=jws-ir+1
      CALL gensim(fi,dx,kr,pint)
C
   23 vnsh(lm,ir)=vnsh(lm,ir)-fl*(qint/sisl/r+pint*sisl)
   22 CONTINUE
      vnsh=-vnsh
C
C     Exchange-correlation contribution to the non-spherical potential
C
      DO 26 is=1,ns
      denl0(is,1:jri)=chde(ita,it,is,1:jri)/ri2(1:jri,ita,it)/sqfpi
      fi(1:jri)=denl0(is,1:jri)
      fi(1:jws)=fi(1:jws)
     .          -facns*cor(ita,it,1:jws)/ri2(1:jws,ita,it)/sqfpi
      CALL diffn(fi,rhop,rhopp,jri,dx)
      denl(is,1,1:jri)=fi(1:jri)
      denld(is,1,1:jri)=rhop(1:jri)/ri(1:jri,ita,it)
      DO 26 l=1,lmaxh
      DO 26 m=-l,l
      lm=l*l+l+m+1
      fi(1:jri)=chdl(ita,iq,is,lm,1:jri)/ri2(1:jri,ita,it)
      CALL diffn(fi,rhop,rhopp,jri,dx)
      denl(is,lm,1:jri)=fi(1:jri)
   26 denld(is,lm,1:jri)=rhop(1:jri)/ri(1:jri,ita,it)
C
C     Only the nonpolarized LDA is implemented !!
C
      DO 27 ir=1,jri
      DO 27 is=1,ns
      rho0=denl0(is,ir)/sqfpi
      CALL pwlda(rho0,exc,excd,excdd,excddd)
      muxcd=(2.d0*excd+rho0*excdd)/2.d0
      vnsxc(is,2:nlmh,ir)=muxcd*denl(is,2:nlmh,ir)
   27 CONTINUE
C
C     Integrate the radial part
C
      DO 28 lm=1,nlmh
      DO 28 lmp=2,nlmh
      DO 29 ir=1,jri
      r=ri(ir,ita,it)
   29 fi(ir)=(SUM(denld(1:ns,lm,ir))*vnsh(lmp,ir)+
     .        SUM(denld(1:ns,lm,ir)*vnsxc(1:ns,lmp,ir)))*r*r*r
      CALL gensim(fi,dx,jsi,fint)
      dvr(lm,lmp)=fint
      DO 30 ir=1,jri
      r=ri(ir,ita,it)
   30 fi(ir)=(SUM(denl(1:ns,lm,ir))*vnsh(lmp,ir)+
     .        SUM(denl(1:ns,lm,ir)*vnsxc(1:ns,lmp,ir)))*r*r
      CALL gensim(fi,dx,jsi,fint)
      dvt(lm,lmp)=fint
   28 CONTINUE
C
C     Performe the surface integral
C
      DO 31 ith=1,nth
      teta    =thvec(ith)
      cost    =DCOS(teta)
      sint    =DSIN(teta)
      DO 31 ifi=1,nfi
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
      snn=0.d0
      DO 32 lm=1,nlmh
      DO 32 lmp=2,nlmh
   32 snn=snn+(SUM(denl(1:ns,lm,jsi))*vnsh(lmp,jsi)+
     .   SUM(denl(1:ns,lm,jsi)*vnsxc(1:ns,lmp,jsi)))*ylm(lm)*ylm(lmp)
C
      srr=0.d0
      srt=0.d0
      srp=0.d0
      DO 33 lm=1,nlmh
      DO 33 lmp=2,nlmh
      srr=srr+dvr(lm,lmp)*ylm(lm)*ylm(lmp)
      srt=srt+dvt(lm,lmp)*gt(lm)*ylm(lmp)
   33 srp=srp+dvt(lm,lmp)*gp(lm)*ylm(lmp)
C
      fqy(iq)=fqy(iq)+weight*snn*uy
      fqz(iq)=fqz(iq)+weight*snn*uz
      fqx(iq)=fqx(iq)+weight*snn*ux
C
      fcy(iq)=fcy(iq)+weight*(srr*uy+srt*sinf*cost+srp*cosf)
      fcz(iq)=fcz(iq)+weight*(srr*uz-srt*sint)
      fcx(iq)=fcx(iq)+weight*(srr*ux+srt*cosf*cost-srp*sinf)
C
   31 CONTINUE
C
      fqy(iq)=2.d0*fqy(iq)*si*si
      fqz(iq)=2.d0*fqz(iq)*si*si
      fqx(iq)=2.d0*fqx(iq)*si*si
C
      fcy(iq)=-2.d0*fcy(iq)
      fcz(iq)=-2.d0*fcz(iq)
      fcx(iq)=-2.d0*fcx(iq)
C
   20 CONTINUE
      DEALLOCATE(ylm,gt,gp)
      DEALLOCATE(denl,denld)
      DEALLOCATE(vnsh,vnsxc)
      DEALLOCATE(dvr,dvt,denl0)
C
  100 FORMAT(/,' NSPSTR: Set up the non-spherical force')
      RETURN
      END
