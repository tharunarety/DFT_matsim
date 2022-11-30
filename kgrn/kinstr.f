      SUBROUTINE kinstr(nzm1,nzm2,gah,zm,wgm,nzm)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the kinetic force components from the stress      *
C   *    field sigma_ki = -2*sum_i grad(phi_i) o grad(phi_i).        *
C   *                                                                *
C   *    The calculation is done for r = ri(jsi), which is the       *
C   *    closest mesh point to the inscribed sphere Si >~ ri(jsi).   *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    Fy,Fz,Fx: the m=-1,0,1 components of the force.             *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE control_text ; USE csts ; USE density
      USE force        ; USE gaussi       ; USE logderivative
      USE message      ; USE partialwaves ; USE radialmesh
      USE slope        ; USE symmetry     ; USE temporary
      IMPLICIT NONE
      INTEGER         :: nzm, nzm1, nzm2
      COMPLEX(KIND=8), DIMENSION(nzm,nq,ns,nlmh,nlmh) :: gah
      COMPLEX(KIND=8), DIMENSION(nzm)       :: zm, wgm
      COMPLEX(KIND=8), DIMENSION(nzm1:nzm2) :: path
      COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: z, zd
      COMPLEX(KIND=8) :: path0
      REAL(KIND=8),    DIMENSION(:,:),   ALLOCATABLE :: zdzd, zdz
      REAL(KIND=8),    DIMENSION(:),     ALLOCATABLE :: chdb
      REAL(KIND=8),    DIMENSION(nq,2:4) :: ksf, ksf0, knf, knf0
      REAL(KIND=8) :: teta, cost, sint, phi, cosf, sinf, ux, uy, uz
      REAL(KIND=8) :: si, weight, snn, srr, srt, srp, smm, wq
      INTEGER :: iq, it, ita, jsi, jrn, ith, ifi, lmaxhp, nlmhp, is, lz
      INTEGER :: l, m, lm, lp, mp, lmp, irot, iqrot, i, ir
C
C     Loop for sites
C
      WRITE(m6,100)
      IF(msgl.NE.0) WRITE(msgio,100)
C
      lmaxhp=lmaxh+1
      nlmhp=(lmaxhp+1)*(lmaxhp+1)
C
      ALLOCATE(zdzd(nlmh,nlmh),zdz(nlmh,nlmh))
      ALLOCATE(chdb(nlmh))
C
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      jrn=jrsm(ita,it)
      IF(fcd.EQ.'Y') jrn=MAX0(jrn,jwsc(ita,it))
      jsi=jwsi(ita,it)-shf
      si=ri(jsi,ita,it)
C
      ALLOCATE(z(0:lmaxh,ns,nzm1:nzm2),zd(0:lmaxh,ns,nzm1:nzm2))
C
C     Set up the gradient of the wave function
C
      DO 21 l=0,lmaxh
      DO 21 is=1,ns
      DO 21 lz=nzm1,nzm2
      z(l,is,lz)=cf(jsi,l,ita,it,is,lz)/si
      fi(1:jrn) =REAL(cf(1:jrn,l,ita,it,is,lz),8)/ri(1:jrn,ita,it)
      fip(1:jrn)=AIMAG(cf(1:jrn,l,ita,it,is,lz))/ri(1:jrn,ita,it)
      CALL diffn(fi,rhxcd1,rhxcdd1,jrn,dx)
      CALL diffn(fip,rhxcd2,rhxcdd2,jrn,dx)
   21 zd(l,is,lz)=CMPLX(rhxcd1(jsi),rhxcd2(jsi),8)
C
      zdzd=0.d0
      zdz=0.d0
      DO 22 is=1,ns
      DO 22 lp=0,lmaxh
      DO 22 mp=-lp,lp
      lmp     =lp*lp+lp+mp+1
      DO 22 l =0,lmaxh
      DO 22 m =-l,l
      lm      =l*l+l+m+1
      path(nzm1:nzm2)=gah(nzm1:nzm2,iq,is,lmp,lm)
      IF(lm.EQ.lmp.AND.l.LE.lmax) THEN
         path(nzm1:nzm2)=path(nzm1:nzm2)+
     .   zone/dfi(nzm1:nzm2,ita,it,l,is,1)
      ENDIF
      path(nzm1:nzm2)=wgm(nzm1:nzm2)*path(nzm1:nzm2)
      path(nzm1:nzm2)=zd(lp,is,nzm1:nzm2)*path(nzm1:nzm2)
      zdzd(lmp,lm)=zdzd(lmp,lm)+SUM(AIMAG(path*zd(l,is,nzm1:nzm2)))
      zdz (lmp,lm)=zdz (lmp,lm)+SUM(AIMAG(path*z (l,is,nzm1:nzm2)))
   22 CONTINUE
C
      DO 23 is=1,ns
      DO 23 l=0,lmax
      DO 23 i=1,necr(l,ita,it,is)
      DO 23 lz=nzm1,nzm2
      path0=wgm(lz)*nocr(l,ita,it,is,i)/(zm(lz)-ecr(l,ita,it,is,i))
      DO 24 m=-l,l
      lm=l*l+l+m+1
      zdzd(lm,lm)=zdzd(lm,lm)+cfrr(l,ita,it,is,i)*AIMAG(path0)
   24 zdz (lm,lm)=zdz (lm,lm)+cfrt(l,ita,it,is,i)*AIMAG(path0)
   23 CONTINUE
C
      zdzd=zdzd/pi
      zdz=zdz/pi
C
      DEALLOCATE(z,zd)
C
C     Set up the gradien^2 n(r) term
C
      DO 25 l=0,lmaxh
      DO 25 m=-l,l
      lm=l*l+l+m+1
      DO 26 ir=1,jrn
   26 fi(ir)=SUM(chdl(ita,iq,1:ns,lm,ir))/ri2(ir,ita,it)
      CALL diffn(fi,rhop,rhopp,jrn,dx)
   25 chdb(lm)=rhopp(jsi)+rhop(jsi)-l*(l+1.d0)*fi(jsi)
C
      ALLOCATE(ylm(nlmhp),gt(nlmh),gp(nlmh))
C
C     Integrate on the inscribed sphere (si >~ ri(jsi))
C
      DO 27 ith=1,nth
      teta    =thvec(ith)
      cost    =DCOS(teta)
      sint    =DSIN(teta)
      DO 27 ifi=1,nfi
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
C     Set up the kinetic stress tensor sigma(r) for r(ux,uy,uz)
C     (only (r,r), (r,theta) and (r,phi) components)
C
      snn=0.d0
      srr=0.d0
      srt=0.d0
      srp=0.d0
      DO 28 lmp=1,nlmh
      snn=snn+chdb(lmp)*ylm(lmp)
      DO 28 lm=1,nlmh
      srr=srr+zdzd(lmp,lm)*ylm(lmp)*ylm(lm)
      srt=srt+zdz (lmp,lm)*ylm(lmp)*gt (lm)
      srp=srp+zdz (lmp,lm)*ylm(lmp)*gp (lm)
   28 CONTINUE
C
      fny(iq)=fny(iq)+weight*snn*uy
      fnz(iq)=fnz(iq)+weight*snn*uz
      fnx(iq)=fnx(iq)+weight*snn*ux
C
      fsy(iq)=fsy(iq)+weight*(srr*uy+srt*sinf*cost+srp*cosf)
      fsz(iq)=fsz(iq)+weight*(srr*uz-srt*sint)
      fsx(iq)=fsx(iq)+weight*(srr*ux+srt*cosf*cost-srp*sinf)
C
   27 CONTINUE
      fny(iq) =spinfc*fny(iq)/2.d0/pi
      fnz(iq) =spinfc*fnz(iq)/2.d0/pi
      fnx(iq) =spinfc*fnx(iq)/2.d0/pi
C
      fsy(iq) =-2.d0*fsy(iq)
      fsz(iq) =-2.d0*fsz(iq)
      fsx(iq) =-2.d0*fsx(iq)
c
      DEALLOCATE(ylm,gt,gp)
   20 CONTINUE
      DEALLOCATE(zdzd,zdz)
      DEALLOCATE(chdb)
C
C     Integrate in the full-BZ
C
      IF(fllbz.NE.'Y') THEN
         knf(1:nq,2)=fny(1:nq)
         knf(1:nq,3)=fnz(1:nq)
         knf(1:nq,4)=fnx(1:nq)
         knf0=0.d0
         ksf(1:nq,2)=fsy(1:nq)
         ksf(1:nq,3)=fsz(1:nq)
         ksf(1:nq,4)=fsx(1:nq)
         ksf0=0.d0
         DO 30 iq=1,nq
         wq=wqst(iq)
         l=1
         DO 31 m=-l,l
         lm=l*l+l+m+1
         DO 32 irot=1,nrot
         IF(ibzr(irot,iq).NE.0) THEN
            smm=0.d0
            snn=0.d0
            DO 33 mp=-l,l
            lmp=l*l+l+mp+1
            snn=snn+ugam(irot,l,m,mp)*knf(iq,lmp)
   33       smm=smm+ugam(irot,l,m,mp)*ksf(iq,lmp)
            iqrot=iprmt(irot,iq)
            knf0(iqrot,lm)=knf0(iqrot,lm)+snn*wq
            ksf0(iqrot,lm)=ksf0(iqrot,lm)+smm*wq
         ENDIF
   32    CONTINUE
   31    CONTINUE
   30    CONTINUE
         fny(1:nq)=knf0(1:nq,2)
         fnz(1:nq)=knf0(1:nq,3)
         fnx(1:nq)=knf0(1:nq,4)
         fsy(1:nq)=ksf0(1:nq,2)
         fsz(1:nq)=ksf0(1:nq,3)
         fsx(1:nq)=ksf0(1:nq,4)
      ENDIF
C
  100 FORMAT(/,' KINSTR: Set up the kinetic stress tensor')
      RETURN
      END
