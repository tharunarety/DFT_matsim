      SUBROUTINE totale(etotal)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the valence (SOFTC=N) or the total (SOFTC=Y)      *
C   *    energy within the Spherical Cell Approximation.             *
C   *    The energy is always calculated from the input density      *
C   *    for the present iteration (chdo) and not from the output    *
C   *    density of the present iteration (chdn).                    *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts       ; USE density      ; USE dosmom ; USE fandgeq
      USE greenfunc  ; USE message      ; USE moments
      USE pota       ; USE potential    ; USE radialmesh
      USE softcore   ; USE temporary    ; USE totalenergy
      IMPLICIT NONE
      INTEGER, PARAMETER :: prnt = 1
      REAL(KIND=8), DIMENSION(mnta,nt,ns) :: qtrs
      REAL(KIND=8), DIMENSION(nq)         :: emdl
      REAL(KIND=8) :: fspin, smm, smmp, fint, fintl, potb, r
      REAL(KIND=8) :: rho, rho1, rho2, v1, v2, exc, etotal, qssm
      REAL(KIND=8) :: fintc, fintcl, twoz, vs1, vs2, concit, delq
      INTEGER :: iq, it, ita, jq, is, jws, jri, jsr, jsrp, jrn, nzc, ir
      INTEGER :: in, jt, jta, ntait, nove, lmp, lm
C
      fspin=spinfc/2.d0
C
C     1) One-electron energies
C
      eone=0.d0
      DO 10 it=1,nt
      DO 10 ita=1,nta(it)
   10 eone(ita,it)=eone(ita,it)+SUM(emom(ita,it,0:lmax,1:ns))
      IF(softc.EQ.'Y') THEN
         DO 11 it=1,nt
         DO 11 ita=1,nta(it)
   11    eone(ita,it)=eone(ita,it)+eonec(ita,it)
      ENDIF
C
C     Set up the number of electrons per site and spin
C
      qtrs=0.d0
      DO 12 is=1,ns
      DO 12 it=1,nt
      DO 12 ita=1,nta(it)
   12 qtrs(ita,it,is)=qtrs(ita,it,is)+SUM(tnos(ita,it,0:lmax,is))
C
C     2) Madelung energy
C
      emdl=0.d0
      DO 13 iq=1,nq
      smm=0.d0
      DO 14 lm=1,diml
      smmp=0.d0
      DO 15 jq=1,nq
      DO 15 lmp=1,diml
   15 smmp=smmp+vmad(iq,jq,lm,lmp)*qlmo(jq,lmp)
   14 smm=smm+qlmo(iq,lm)*smmp
      smm=0.5d0*smm/sws
   13 emdl(iq)=emdl(iq)+smm
C
      DO 16 iq=1,nq
      it=itq(iq)
      ntait=nta(it)
      IF(ntait.EQ.1) THEN
         emadl(ntait,it)=emdl(iq)
      ELSE
C
C        Screening correction by the SIM
C
         DO 17 ita=1,ntait
         delq=qcpa(ita,it)
   17    emadl(ita,it)=emdl(iq)-delq*delq*alphmd/sws
      ENDIF
   16 CONTINUE
C
C     Kinetic, Coulombic and exchange correlation terms per site
C
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jws=jwss(ita,it)
      jsr=jsrs(ita,it)
      jri=jws+2
      jrn=jrsm(ita,it)
      twoz=2.d0*nz(ita,it)
C
      rhoa(1:jri)=cor(ita,it,1:jri)
      rhoa(jri+1:jsr)=0.d0
C
C     3) Potential * density part of the kinetic energy
C
C        Vint = sum_S int_S [ v(r) - vmtz ] * n(r) dr + Q * vmtz
C
C
      vint(ita,it)=0.d0
      vintc(ita,it)=0.d0
      DO 21 is=1,ns
C
C     3a) core part for soft-core and FCD
C
      fi(1:jsr)=(v(1:jsr,ita,it,is)-vmtz(is)*ri2(1:jsr,ita,it))*
     .          rhoa(1:jsr)/ri(1:jsr,ita,it)
      CALL simpn(fi,dx,jsr,fintc)
      fintc=fintc+(nz(ita,it)-eln(ita,it))*vmtz(is)
C
C     Add int_S^Sm [ vmtzR - vmtz ] * n(r) dr ~ [ vmtzR - vmtz ] * Q_S^Sm
C
C     where Q_S^Sm is the charge between the S and S^m (>S)
C
      IF(localmt(ita,it).NE.0) THEN
         qssm=rhoa(jsr)/hsr(ita,it)/hsr(ita,it)
         fintc=fintc+(vmtzr(it,is)-vmtz(is))*
     .        qssm*(wsm(ita,it)**3.d0-hsr(ita,it)**3.d0)/3.d0
      ENDIF
      vintc(ita,it)=vintc(ita,it)+fspin*fintc
C
C     3b) valence part (For overlapping muffin-tin wells this term
C         should be calculated from the spherical part of the
C         full charge density (~chdo) and not from the renormalized
C         density (chde). However, from the spherical part we neglect
C         the higher tail components, therefore this is approximated 
C         by the renormalized density.)
C
      rhov(1:jsr)=chde(ita,it,is,1:jsr)-fspin*rhoa(1:jsr)
      fi(1:jsr)=(v(1:jsr,ita,it,is)-vmtz(is)*ri2(1:jsr,ita,it))*
     .          rhov(1:jsr)/ri(1:jsr,ita,it)
      CALL simpn(fi,dx,jsr,fint)
C
C     Add int_S^Sm [ vmtzR - vmtz ] * n(r) dr ~ [ vmtzR - vmtz ] * Q_S^Sm
C
C     where Q_S^Sm is the charge between the S and S^m (>S)
C
      IF(localmt(ita,it).NE.0) THEN
         qssm=rhov(jsr)/hsr(ita,it)/hsr(ita,it)
         fint=fint+(vmtzr(it,is)-vmtz(is))*
     .        qssm*(wsm(ita,it)**3.d0-hsr(ita,it)**3.d0)/3.d0
      ENDIF
      vint(ita,it)=vint(ita,it)+fint+qtrs(ita,it,is)*vmtz(is)
C
      IF(softc.EQ.'Y') THEN
         vint(ita,it)=vint(ita,it)+fspin*fintc
      ENDIF
   21 CONTINUE
C
      IF(ns.EQ.1) THEN
         rhov(1:jri)=chde(ita,it,1,1:jri)-rhoa(1:jri)
      ELSE
         rhov(1:jri)=chde(ita,it,1,1:jri)+chde(ita,it,2,1:jri)-
     .               rhoa(1:jri)
      ENDIF
C
C     4) Nuclear interaction
C
      enuc(ita,it)=0.d0
      CALL simpn(rhov,dx,jws,fint)
      enuc(ita,it)=enuc(ita,it)-twoz*fint
      IF(softc.EQ.'Y') THEN
         CALL simpn(rhoa,dx,jws,fint)
         enuc(ita,it)=enuc(ita,it)-twoz*fint
      ENDIF
C
C     Hartree potential for the core
C
      potb=2.d0*(nz(ita,it)-eln(ita,it))/ws(ita,it)
      nzc=0
      CALL poisson(rhoa,nzc,wc,potb,ws(ita,it),it,ita,jws,1)
C
C     5) Core-valence interaction
C
      ecor(ita,it)=0.d0
      fi(1:jri)=wc(1:jri)*rhov(1:jri)*ri(1:jri,ita,it)
      CALL simpn(fi,dx,jws,fint)
      ecor(ita,it)=ecor(ita,it)+fint
      IF(softc.EQ.'Y') THEN
         fi(1:jri)=wc(1:jri)*rhoa(1:jri)*ri(1:jri,ita,it)
         CALL simpn(fi,dx,jws,fint)
         ecor(ita,it)=ecor(ita,it)+0.5d0*fint
      ENDIF
C
C     Hartree potential from the valence
C
      potb=2.d0*(qtr(ita,it)+eln(ita,it))/ws(ita,it)
      nzc=0
      CALL poisson(rhov,nzc,wc,potb,ws(ita,it),it,ita,jws,1)
C
C     5) Valence-valence interaction
C
      eval(ita,it)=0.d0 
      fi(1:jri)=wc(1:jri)*rhov(1:jri)*ri(1:jri,ita,it)
      CALL simpn(fi,dx,jws,fint)
      eval(ita,it)=0.5d0*fint+eval(ita,it)
C
C     6) Total exchange-correlation energy
C
      exct(ita,it)=0.d0
      IF(ns.EQ.1) THEN
         chdt(1:jri)=chde(ita,it,1,1:jri)/ri2(1:jri,ita,it)
         chdt=chdt/fourpi
         IF(ixc.GE.8) THEN
            CALL diffn(chdt,rhop,rhopp,jri,dx)
            rhopp(1:jri)=(rhopp(1:jri)-rhop(1:jri))/ri2(1:jri,ita,it)
            rhop(1:jri)=rhop(1:jri)/ri(1:jri,ita,it)
         ENDIF
      ELSE
         rhosp1(1:jri)=chde(ita,it,1,1:jri)/ri2(1:jri,ita,it)
         rhosp2(1:jri)=chde(ita,it,2,1:jri)/ri2(1:jri,ita,it)
         rhosp1=rhosp1/fourpi
         rhosp2=rhosp2/fourpi
         chdt=rhosp1+rhosp2
         IF(IXC.GE.8) THEN
            CALL diffn(rhosp1,rhxcd1,rhxcdd1,jri,dx)
            CALL diffn(rhosp2,rhxcd2,rhxcdd2,jri,dx)
            rhxcdd1(1:jri)=(rhxcdd1(1:jri)-rhxcd1(1:jri))/
     .                      ri2(1:jri,ita,it)
            rhxcdd2(1:jri)=(rhxcdd2(1:jri)-rhxcd2(1:jri))/
     .                      ri2(1:jri,ita,it)
            rhxcd1(1:jri)=rhxcd1(1:jri)/ri(1:jri,ita,it)
            rhxcd2(1:jri)=rhxcd2(1:jri)/ri(1:jri,ita,it)
         ENDIF
      ENDIF
      DO 26 ir=1,jri
      r=ri(ir,ita,it)
      rho=chdt(ir)
      IF(ns.EQ.1) THEN
         rho1=0.5d0*rho
         rho2=rho1
         IF(ixc.GE.8) THEN
            rhod(1)=0.5*rhop(ir)
            rhodd(1)=0.5*rhopp(ir)
            rhod(2)=rhod(1)
            rhodd(2)=rhodd(1)
         ENDIF
      ELSE
         rho1=rhosp1(ir)
         rho2=rhosp2(ir)
         IF(ixc.GE.8) THEN
            rhod(1)=rhxcd1(ir)
            rhod(2)=rhxcd2(ir)
            rhodd(1)=rhxcdd1(ir)
            rhodd(2)=rhxcdd2(ir)
         ENDIF
      ENDIF
      CALL xcpot(ixc,rho1,rho2,rho,rhod,rhodd,r,v1,v2,exc)
      fi(ir)=exc*(rhoa(ir)+rhov(ir))*r
   26 CONTINUE
      CALL simpn(fi,dx,jws,fint)
      exct(ita,it)=exct(ita,it)+fint
C
C     7) Core exchange-correlation energy
C
      excc(ita,it)=0.d0
      IF(softc.EQ.'Y') GO TO 28
      fip(1:jri)=rhoa(1:jri)/ri2(1:jri,ita,it)
      fip=fip/fourpi
      IF(ixc.GE.8) THEN
         CALL diffn(fip,rhop,rhopp,jri,dx)
         rhopp(1:jri)=(rhopp(1:jri)-rhop(1:jri))/ri2(1:jri,ita,it)
         rhop(1:jri)=rhop(1:jri)/ri(1:jri,ita,it)
      ENDIF
      DO 27 ir=1,jri
      r=ri(ir,ita,it)
      rho=fip(ir)
      rho1=0.5d0*rho
      rho2=rho1
      IF(ixc.GE.8) THEN
         rhod(1)=0.5d0*rhop(ir)
         rhodd(1)=0.5d0*rhopp(ir)
         rhod(2)=rhod(1)
         rhodd(2)=rhodd(1)
      ENDIF
      CALL xcpot(ixc,rho1,rho2,rho,rhod,rhodd,r,v1,v2,exc)
      fi(ir)=exc*rhoa(ir)*r
   27 CONTINUE
      CALL simpn(fi,dx,jws,fint)
      excc(ita,it)=excc(ita,it)+fint
   28 CONTINUE
   20 CONTINUE
C
C     8) Kinetic energy error
C
      okae=0.d0
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      jws=jwss(ita,it)
      jsr=jsrs(ita,it)
      nove=nov(ita,it)
      DO 31 is=1,ns
      vs1=v(jsr,ita,it,is)/ri2(jsr,ita,it)-vmtz(is)
      rho2=chde(ita,it,1,jws)/ws(ita,it)/ws(ita,it)/fourpi
C
      DO 31 in=1,nove
      jt=type(ita,it,in)
      vs2=0.d0
      DO 32 jta=1,nta(jt)
      jsrp=jsrs(jta,jt)
      vs2=vs2+conc(jta,jt)*(v(jsrp,jta,jt,is)/ri2(jsrp,jta,jt)-vmtz(is))
   32 CONTINUE
   31 okae(ita,it)=okae(ita,it)+
     .             wov(ita,it,in)*okac(ita,it,in)*vs1*vs2*rho2
   30 CONTINUE
C
C     Calculate the total energy
C
      IF(softz.EQ.1) eone=eone+etotcore
      ekin=eone-vint
      etot=ekin+emadl+enuc+ecor+eval+exct-excc
      etotal=0.d0
      DO 50 it=1,nt
      DO 50 ita=1,nta(it)
      concit=conc(ita,it)
      IF(msgl.NE.0.AND.prnt.EQ.1.AND.concit.GT.1.d-6) THEN
         WRITE(msgio,100) it,eone(ita,it),vint(ita,it),ekin(ita,it),
     .         ita,ecor(ita,it),enuc(ita,it),emadl(ita,it),eval(ita,it),
     .             exct(ita,it),excc(ita,it),etot(ita,it),
     .             okae(ita,it),etot(ita,it)+okae(ita,it)
      ENDIF
   50 etotal=etotal+concit*mmt(it)*etot(ita,it)
      etotal=etotal/nq
C
  100 FORMAT(/,' IT =',i2,' EONE =',f14.6,' VINT =',f14.6,' EKIN =',
     .     f14.6,/,' ITA=',i2,' ECOR =',f14.6,' ENUC =',f14.6,' EMADL=',
     .     f14.6,/,8x,'EVAL =',f14.6,' EXCT =',f14.6,' EXCC =',
     .     f14.6,/,8x,'ETOT =',f14.6,' OKAC =',f14.6,' E[n] =',f14.6)
      RETURN
      END
