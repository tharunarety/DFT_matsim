      SUBROUTINE arrays(all)
C   ******************************************************************
C   *                                                                *
C   * Allocate different arrays.                                     *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE density    ; USE dosmom       ; USE energymesh
      USE greenfunc  ; USE kinkmatrix   ; USE logderivative
      USE moments    ; USE partialwaves ; USE poissonparam
      USE potential  ; USE potparam     ; USE radialmesh
      USE realgaunt  ; USE slope        ; USE softcore
      USE temporary  ; USE totalenergy
      IMPLICIT NONE
      INTEGER :: all, i
C
      IF(all.EQ.1) THEN
         ALLOCATE(rp(dimr),rq(dimr))
         rp=0.d0
         rq=0.d0
         ALLOCATE(chdn(mnta,nt,ns,dimr),chdo(mnta,nt,ns,dimr))
         ALLOCATE(chde(mnta,nt,ns,dimr))
         chdn=0.d0
         chdo=0.d0
         chde=0.d0
         ALLOCATE(chdh(mnta,nt,ns,dimr),chdr(mnta,nt,dimr))
         chdh=0.d0
         chdr=0.d0
         ALLOCATE(fi(dimr),fip(dimr))
         fi=0.d0
         fip=0.d0
         ALLOCATE(tnos(mnta,nt,0:lmax,ns),emom(mnta,nt,0:lmax,ns))
         ALLOCATE(tdos(mnta,nt,0:lmax,ns),entr(mnta,nt,0:lmax,ns))
         tnos=0.d0
         tdos=0.d0
         emom=0.d0
         entr=0.d0
         ALLOCATE(amag(mnta,nt))
         IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
            ALLOCATE(tnos0(mnta,nt,0:lmax,ns),emom0(mnta,nt,0:lmax,ns))
            tnos0=0.d0
            emom0=0.d0
         ENDIF
         amag=0.d0
         ALLOCATE(gak(nlmq,nlmq),work(nlmq),worl(nlm))
         gak=zero
         work=0.d0
         worl=0.d0
         ALLOCATE(kinkm(nlmq,nlmq),unit(nlmq,nlmq),unil(nlm,nlm))
         kinkm=zero
         unit=zero
         unil=zero
         DO 10 i=1,nlmq
   10    unit(i,i)=CMPLX(1.d0,0.d0,8)
         DO 11 i=1,nlm
   11    unil(i,i)=CMPLX(1.d0,0.d0,8)
C
         ALLOCATE(sgm(0:lmax,mnta,nt))
         sgm=0.d0
         ALLOCATE(rhoa(dimr),rhov(dimr),wc(dimr),chdt(dimr))
         rhoa=0.d0
         rhov=0.d0
         wc=0.d0
         chdt=0.d0
         ALLOCATE(rhop(dimr),rhopp(dimr),rhosp1(dimr),rhosp2(dimr))
         rhop=0.d0
         rhopp=0.d0
         rhosp1=0.d0
         rhosp2=0.d0
         ALLOCATE(rhxcd1(dimr),rhxcdd1(dimr),rhxcd2(dimr),rhxcdd2(dimr))
         rhxcd1=0.d0
         rhxcdd1=0.d0
         rhxcd2=0.d0
         rhxcdd2=0.d0
         ALLOCATE(ecr(0:lmax,mnta,nt,ns,dimecr),necr(0:lmax,mnta,nt,ns))
         ALLOCATE(nocr(0:lmax,mnta,nt,ns,dimecr))
         ecr=0.d0
         nocr=0.d0
         necr=0
         ALLOCATE(p(dimr),q(dimr))
         p=zero
         q=zero
         ALLOCATE(cfrr(0:lmax,mnta,nt,ns,dimecr))
         ALLOCATE(cfrt(0:lmax,mnta,nt,ns,dimecr))
         cfrr=0.d0
         cfrt=0.d0
         ALLOCATE(core(mnta,nt,dimr),atmc(mnta,nt,dimr))
         core=0.d0
         atmc=0.d0
         ALLOCATE(v(dimr,mnta,nt,ns),vp(dimr,mnta,nt,ns))
         ALLOCATE(fullp(dimr,mnta,nt,ns),potm(mnta,nt))
         v=0.d0
         vp=0.d0
         potm=0.d0
         fullp=0.d0
         ALLOCATE(qs(mnta,nt),potmc(mnta,nt))
         qs=0.d0
         potmc=0.d0
         ALLOCATE(qsca(nt),qcpa(mnta,nt))
         qsca=0.d0
         qcpa=0.d0
         ALLOCATE(eone(mnta,nt),emadl(mnta,nt),vint(mnta,nt))
         ALLOCATE(ents(mnta,nt))
         ents=0.d0
         eone=0.d0
         emadl=0.d0
         vint=0.d0
         ALLOCATE(vintc(mnta,nt),enuc(mnta,nt),ecor(mnta,nt))
         vintc=0.d0
         enuc=0.d0
         ecor=0.d0
         ALLOCATE(eval(mnta,nt),exct(mnta,nt),excc(mnta,nt))
         eval=0.d0
         exct=0.d0
         excc=0.d0
         ALLOCATE(ekin(mnta,nt),etot(mnta,nt),okae(mnta,nt))
         ekin=0.d0
         etot=0.d0
         okae=0.d0
         ALLOCATE(pots(mnta,nt,ns),potw(mnta,nt,ns))
         pots=0.d0
         potw=0.d0
         ALLOCATE(e(dimr),f(dimr))
         e=0.d0
         f=0.d0
         ALLOCATE(qlmn(nq,diml),qlmo(nq,diml),qlm0(nq,diml))
         qlmn=0.d0
         qlmo=0.d0
         qlm0=0.d0
         IF(func.NE.'ASA') THEN
            ALLOCATE(grnfm(nq,ns,nlm,nlm,0:lmax2))
            grnfm=0.d0
         ENDIF
      ELSE
         IF(func.NE.'ASA') THEN
            DEALLOCATE(gnt,lmpg,lmg,lmppg,lppg)
            DEALLOCATE(grnfm)
         ENDIF
         DEALLOCATE(rp,rq)
         DEALLOCATE(qs,potmc,qsca)
         DEALLOCATE(chdn,chdr)
         DEALLOCATE(rhoa,rhov,wc,chdt,rhosp1,rhosp2)
         DEALLOCATE(e,f,qlmn,qlmo,qlm0)
      ENDIF
C
      RETURN
      END
