      SUBROUTINE start(ef)
C   ******************************************************************
C   *                                                                *
C   *    Construct charge density and the ASA potential for first    *
C   *    loop from input charge densities.                           *
C   *                                                                *
C   *   *ATMC = Charge density from atomic calculation,renormalized  *
C   *   *CHDA = Charge Density used generate the potential           *
C   *                                                                *
C   *   *V    = Potential*Radius**2                                  *
C   *                                                                *
C   *    Variables marked with an * are set on exit.                 *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts       ; USE density      ; USE dosmom ; USE message
      USE moments    ; USE pota         ; USE potential
      USE radialmesh ; USE temporary    ; USE text
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(nq) :: potq
      REAL(KIND=8) :: v1, v2, exc, rho, rho1, rho2, fac1, fac2, qsit
      REAL(KIND=8) :: r, rce, smms, potb, wcws, s, addz, ef, qtrit
      REAL(KIND=8), PARAMETER :: tol=1.d-06
      INTEGER :: iq, it, ntait, ita, jq, is, ir, iv
      INTEGER :: jrn, jsr, jws, jri, lmp
C
      ALLOCATE(zn0(mnta,nt))
C
      DO 10 it=1,nt
      DO 10 ita=1,nta(it)
      jrn=jrsm(ita,it)
      jri=jwss(ita,it)+2
      IF(ns.EQ.1) THEN
         spinfc=2.d0
         IF(tpot.EQ.'N') chdo(ita,it,1,jri+1:jrn)=chdo(ita,it,1,jri)
         atmc(ita,it,1:jrn)=chdo(ita,it,1,1:jrn)
      ELSE
         spinfc=1.d0
         IF(tpot.EQ.'N') THEN
            chdo(ita,it,1,jri+1:jrn)=chdo(ita,it,1,jri)
            chdo(ita,it,2,jri+1:jrn)=chdo(ita,it,2,jri)
         ENDIF
         atmc(ita,it,1:jrn)=chdo(ita,it,1,1:jrn)+chdo(ita,it,2,1:jrn)
      ENDIF
   10 CONTINUE
C
      CALL renorm
C
C     Initialize moments and Madelung potential
C
      IF(tpot.EQ.'N') THEN
         DO 11 iq=1,nq
         it=itq(iq)
         qlmo(iq,1)=SUM(conc(1:nta(it),it)*qtr(1:nta(it),it))
   11    qlmo(iq,2:diml)=0.d0
         chde=chdo
C
         potq=0.d0
         DO 12 iq=1,nq
         smms=0.d0
         DO 13 jq=1,nq
         DO 13 lmp=1,diml
   13    smms=smms+vmad(iq,jq,1,lmp)*qlmo(jq,lmp)
   12    potq(iq)=smms/sws
C
         DO 15 it=1,nt
         ntait=nta(it)
         qsit=SUM(conc(1:ntait,it)*qs(1:ntait,it))   !! Qsca
         qtrit=SUM(conc(1:ntait,it)*qtr(1:ntait,it)) !! Qlm
         qsca(it)=qsit-qtrit                         !! Qsca - Qlm
         qcpa(1:ntait,it)=qs(1:ntait,it)-qsit        !! Qsim
   15    CONTINUE
C
C        SCA and screening correction by the SIM
C
         DO 16 iq=1,nq
         it=itq(iq)
         ntait=nta(it)
         DO 17 ita=1,ntait
   17    potmc(ita,it)=-madc(ita,iq)*qsca(it)
     .                 -2.d0*alphmd*qcpa(ita,it)/sws
   16    CONTINUE
C
         DO 18 iq=1,nq
         it=itq(iq)
         DO 18 ita=1,nta(it)
   18    potm(ita,it)=potq(iq)+potmc(ita,it)
C
         IF(afm.EQ.'M') THEN
            tmagoo=mmom
            tmago=mmom
            tmag=0.d0
            DO it=1,nt
            DO ita=1,nta(it)
            tmag=tmag+conc(ita,it)*mmt(it)*split(ita,it)
            ENDDO
            ENDDO
            DO it=1,nt
            DO ita=1,nta(it)
            dexcho(ita,it)=0.d0
            CALL fxspin(ita,it,1)
            ENDDO
            ENDDO
         ELSEIF(afm.EQ.'m') THEN
            DO it=1,nt
            DO ita=1,nta(it)
            IF(fixst(ita,it).EQ.'Y') THEN
               mmom=split(ita,it)
               tmagoo=mmom
               tmago=mmom
               tmag=mmom
               dexcho(ita,it)=0.d0
               CALL fxspin(ita,it,1)
            ELSE
               dexcho(ita,it)=0.d0
               dexch(ita,it)=0.d0
            ENDIF
            ENDDO
            ENDDO
         ENDIF
      ELSE
         IF(afm.EQ.'M') THEN
            tmagoo=mmom
            tmago=mmom
            dexcho=dexch
         ELSEIF(afm.EQ.'m') THEN
            DO it=1,nt
            DO ita=1,nta(it)
            IF(fixst(ita,it).EQ.'Y') THEN
               tmagoo=split(ita,it)
               tmago=split(ita,it)
               tmag=split(ita,it)
               dexcho(ita,it)=dexch(ita,it)
            ELSE
               dexcho(ita,it)=0.d0
               dexch(ita,it)=0.d0
            ENDIF
            ENDDO
            ENDDO
         ELSEIF(afm.EQ.'F') THEN
            IF(ABS(mmom).LT.1.d-6) dexch=0.d0
         ENDIF
      ENDIF
C
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      jws=jwss(ita,it)
      jri=jws+2
      jsr=jsrs(ita,it)
      jrn=jrsm(ita,it)
      IF(tpot.EQ.'N'.AND.ns.EQ.1) THEN
C
C        Non-polarized case
C
         rhoa(1:jrn)=chdo(ita,it,1,1:jrn)
         s=ws(ita,it)
         potb=2.d0*qs(ita,it)/s
         CALL poisson(rhoa,nz(ita,it),wc,potb,s,it,ita,jws,1)
         IF(jsr.GE.jws) THEN
            wcws=wc(jws)*ws(ita,it)
            wc(jws+1:jrn)=wcws/ri(jws+1:jrn,ita,it)
         ENDIF
         wc=wc+potm(ita,it)
C
         rhoa(1:jrn)=rhoa(1:jrn)/ri2(1:jrn,ita,it)
         rhoa=rhoa/fourpi
         IF(ixc.GE.8) THEN
            CALL diffn(rhoa,rhop,rhopp,jrn,dx)
            rhopp(1:jrn)=(rhopp(1:jrn)-rhop(1:jrn))/ri2(1:jrn,ita,it)
            rhop(1:jrn)=rhop(1:jrn)/ri(1:jrn,ita,it)
         ENDIF
         DO 22 ir=1,jrn
         r=ri(ir,ita,it)
         rce=r*r
         rho=rhoa(ir)
         IF(rho.LT.0.d0) rho=1.d-10
         rho1=0.5d0*rho
         rho2=rho1
         IF(ixc.GE.8) THEN
            rhod(1:2)=0.5d0*rhop(ir)
            rhodd(1:2)=0.5d0*rhopp(ir)
         ENDIF
         CALL xcpot(ixc,rho1,rho2,rho,rhod,rhodd,r,v1,v2,exc)
   22    v(ir,ita,it,1)=rce*(wc(ir)+v1)
      ELSEIF(tpot.EQ.'N'.AND.ns.EQ.2) THEN
C
C        Spin polarized case
C
         IF(ABS(eln(ita,it)).LT.tol) THEN
            fac1=0.d0
         ELSE
            fac1=split(ita,it)/eln(ita,it)/2.d0
         ENDIF
         fac2=.5d0-fac1
         rhoa(1:jrn)=chdo(ita,it,1,1:jrn)+chdo(ita,it,2,1:jrn)
         s=ws(ita,it)
         potb=2.d0*qs(ita,it)/s
         CALL poisson(rhoa,nz(ita,it),wc,potb,s,it,ita,jws,1)
         IF(jsr.GE.jws) THEN
            wcws=wc(jws)*ws(ita,it)
            wc(jws+1:jrn)=wcws/ri(jws+1:jrn,ita,it)
         ENDIF
         wc=wc+potm(ita,it)
C
         chdo(ita,it,1,1:jrn)=fac1*cor(ita,it,1:jrn)+fac2*rhoa(1:jrn)
         chdo(ita,it,2,1:jrn)=rhoa(1:jrn)-chdo(ita,it,1,1:jrn)
         rhosp1(1:jrn)=chdo(ita,it,1,1:jrn)/ri2(1:jrn,ita,it)
         rhosp2(1:jrn)=chdo(ita,it,2,1:jrn)/ri2(1:jrn,ita,it)
         rhosp1=rhosp1/fourpi
         rhosp2=rhosp2/fourpi
         IF(ixc.GE.8) THEN
            CALL diffn(rhosp1,rhxcd1,rhxcdd1,jrn,dx)
            CALL diffn(rhosp2,rhxcd2,rhxcdd2,jrn,dx)
            rhxcdd1(1:jrn)=(rhxcdd1(1:jrn)-rhxcd1(1:jrn))/
     .                      ri2(1:jrn,ita,it)
            rhxcdd2(1:jrn)=(rhxcdd2(1:jrn)-rhxcd2(1:jrn))/
     .                      ri2(1:jrn,ita,it)
            rhxcd1(1:jrn)=rhxcd1(1:jrn)/ri(1:jrn,ita,it)
            rhxcd2(1:jrn)=rhxcd2(1:jrn)/ri(1:jrn,ita,it)
         ENDIF
         DO 24 ir=1,jrn
         r=ri(ir,ita,it)
         rce=r*r
         rho1=rhosp1(ir)
         rho2=rhosp2(ir)
         rho=rho1+rho2
         IF(rho.LT.0.d0) THEN
            rho=1.d-10
            rho1=1.d-10
            rho2=1.d-10
         ENDIF
         IF(ixc.GE.8) THEN
            rhod(1)=rhxcd1(ir)
            rhod(2)=rhxcd2(ir)
            rhodd(1)=rhxcdd1(ir)
            rhodd(2)=rhxcdd2(ir)
         ENDIF
         CALL xcpot(ixc,rho1,rho2,rho,rhod,rhodd,r,v1,v2,exc)
         v(ir,ita,it,1)=rce*(wc(ir)+v1+dexch(ita,it))
   24    v(ir,ita,it,2)=rce*(wc(ir)+v2-dexch(ita,it))
      ELSE
         CONTINUE
      ENDIF
C
      DO 25 is=1,ns
      CALL diff(dx,jrn,it,ita,is)
      potw(ita,it,is)=v(jws,ita,it,is)/ri2(jws,ita,it)
      pots(ita,it,is)=v(jsr,ita,it,is)/ri2(jsr,ita,it)
   25 CONTINUE
C
      fullp(1:dimr,ita,it,1:ns)=v(1:dimr,ita,it,1:ns)
C
   21 CONTINUE
C
C     Calculate muffin-tin zero
C
      IF(tpot.EQ.'N') CALL optpot(fixvmtz,1,ef,1.d0)
C
C     Set up the normalization function for the ASA Madelung
C     and for the ASA total energy
C
      zn0=0.d0
      znt=0.d0
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      znt=znt+conc(ita,it)*mmt(it)*nz(ita,it)
      jrn=jrsm(ita,it)
      fi(1:jrn)=ri(1:jrn,ita,it)*ri2(1:jrn,ita,it)
      jws=jwss(ita,it)
      CALL simpn(fi,dx,jws,addz)
      zn0(ita,it)=addz
      chdr(ita,it,1:jrn)=fi(1:jrn)/ri(1:jrn,ita,it)
   30 CONTINUE
C
      RETURN
      END
