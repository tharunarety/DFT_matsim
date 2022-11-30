      SUBROUTINE newpot(prnt,lin)
C   ******************************************************************
C   *                                                                *
C   *    Calculation of the one electron potential within the ASA.   *
C   *    The potential is calculated from the input density (chdo).  *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text ; USE csts
      USE density    ; USE dosmom       ; USE message
      USE moments    ; USE pota         ; USE potential
      USE radialmesh ; USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(nq) :: potq
      REAL(KIND=8) :: v1, v2, exc, rho, rho1, rho2
      REAL(KIND=8) :: r, rce, smms, potb, wcws, s
      INTEGER :: prnt, lin, iq, it, ita, jta, jq, iv, is, ir, lmp
      INTEGER :: jrn, jsr, jws
C
C     Madelung potential
C
      potq=0.d0
      DO 10 iq=1,nq
      smms=0.d0
      DO 11 jq=1,nq
      DO 11 lmp=1,diml
   11 smms=smms+vmad(iq,jq,1,lmp)*qlmo(jq,lmp)
   10 potq(iq)=smms/sws
C
C     SCA and screening correction by the SIM
C
      DO 12 iq=1,nq
      it=itq(iq)
      DO 12 ita=1,nta(it)
      potmc(ita,it)=-madc(ita,iq)*qsca(it)
     .              -2.d0*alphmd*qcpa(ita,it)/sws
   12 CONTINUE
C
      DO 13 iq=1,nq
      it=itq(iq)
      DO 13 ita=1,nta(it)
   13 potm(ita,it)=potq(iq)+potmc(ita,it)
C
C     Establish exchange splitting based on fixed input or DLM moments
C
      IF(afm.EQ.'M') THEN
         DO it=1,nt
         DO ita=1,nta(it)
         CALL fxspin(ita,it,lin)
         ENDDO
         ENDDO
      ELSEIF(afm.EQ.'m') THEN
         DO it=1,nt
         DO ita=1,nta(it)
         IF(fixst(ita,it).EQ.'Y') THEN
            mmom=split(ita,it)
            tmag=amag(ita,it)
            tmago=splito(ita,it)
            tmagoo=splitoo(ita,it)
            CALL fxspin(ita,it,lin)
            splitoo(ita,it)=tmagoo
            splito(ita,it)=tmago
         ELSE
            jta=fixsjta(ita,it)
            IF(jta.NE.0) THEN
               mmom=-amag(jta,it)
               tmag=amag(ita,it)
               IF(ABS(mmom-tmag).LT.1.d-6) THEN
                  tmago=mmom
                  tmagoo=mmom
               ELSE
                  tmago=splito(ita,it)
                  tmagoo=splitoo(ita,it)
               ENDIF
               CALL fxspin(ita,it,lin)
               splitoo(ita,it)=tmagoo
               splito(ita,it)=tmago
            ENDIF
         ENDIF
         ENDDO
         ENDDO
      ENDIF
C
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jws=jwss(ita,it)
      jsr=jsrs(ita,it)
      jrn=jrsm(ita,it)
      IF(ns.EQ.1) THEN
         rhoa(1:jrn)=chdo(ita,it,1,1:jrn)
      ELSE
         rhoa(1:jrn)=chdo(ita,it,1,1:jrn)+chdo(ita,it,2,1:jrn)
      ENDIF
C
C     Boundary condition at W
C
      s=ws(ita,it)
      potb=2.d0*qs(ita,it)/s
      CALL poisson(rhoa,nz(ita,it),wc,potb,s,it,ita,jws,1)
      IF(jsr.GE.jws) THEN
         wcws=wc(jws)*ws(ita,it)
         wc(jws+1:jrn)=wcws/ri(jws+1:jrn,ita,it)
      ENDIF
      wc=wc+potm(ita,it)
C
C     Exchange-correlation potential
C
      IF(ns.EQ.1) THEN
         rhoa(1:jrn)=rhoa(1:jrn)/ri2(1:jrn,ita,it)
         rhoa=rhoa/fourpi
         IF(ixc.GE.8) THEN
            CALL diffn(rhoa,rhop,rhopp,jrn,dx)
            rhopp(1:jrn)=(rhopp(1:jrn)-rhop(1:jrn))/ri2(1:jrn,ita,it)
            rhop(1:jrn)=rhop(1:jrn)/ri(1:jrn,ita,it)
         ENDIF
      ELSE
         rhosp1(1:jrn)=chdo(ita,it,1,1:jrn)/ri2(1:jrn,ita,it)
         rhosp2(1:jrn)=chdo(ita,it,2,1:jrn)/ri2(1:jrn,ita,it)
         rhosp1=rhosp1/fourpi
         rhosp2=rhosp2/fourpi
         rhoa=rhosp1+rhosp2
         IF(ixc.GE.8) THEN
            CALL diffn(rhosp1,rhxcd1,rhxcdd1,jrn,dx)
            CALL diffn(rhosp2,rhxcd2,rhxcdd2,jrn,dx)
            rhxcdd1(1:jrn)=
     .      (rhxcdd1(1:jrn)-rhxcd1(1:jrn))/ri2(1:jrn,ita,it)
            rhxcdd2(1:jrn)=
     .      (rhxcdd2(1:jrn)-rhxcd2(1:jrn))/ri2(1:jrn,ita,it)
            rhxcd1(1:jrn)=rhxcd1(1:jrn)/ri(1:jrn,ita,it)
            rhxcd2(1:jrn)=rhxcd2(1:jrn)/ri(1:jrn,ita,it)
         ENDIF
      ENDIF
      DO 21 ir=1,jrn
      r=ri(ir,ita,it)
      rce=r*r
      rho=rhoa(ir)
      IF(ns.EQ.1) THEN
         rho1=0.5d0*rho
         rho2=rho1
         IF(ixc.GE.8) THEN
            rhod(1:2)=0.5d0*rhop(ir)
            rhodd(1:2)=0.5d0*rhopp(ir)
         ENDIF
         CALL xcpot(ixc,rho1,rho2,rho,rhod,rhodd,r,v1,v2,exc)
         v(ir,ita,it,1)=rce*(wc(ir)+v1)
      ELSE
         rho1=rhosp1(ir)
         rho2=rhosp2(ir)
         IF(ixc.GE.8) THEN
            rhod(1)=rhxcd1(ir)
            rhod(2)=rhxcd2(ir)
            rhodd(1)=rhxcdd1(ir)
            rhodd(2)=rhxcdd2(ir)
         ENDIF
         CALL xcpot(ixc,rho1,rho2,rho,rhod,rhodd,r,v1,v2,exc)
         v(ir,ita,it,1)=rce*(wc(ir)+v1+dexch(ita,it))
         v(ir,ita,it,2)=rce*(wc(ir)+v2-dexch(ita,it))
      ENDIF
   21 CONTINUE
C
      DO 22 is=1,ns
      CALL diff(dx,jrn,it,ita,is)
      potw(ita,it,is)=v(jws,ita,it,is)/ri2(jws,ita,it)
      pots(ita,it,is)=v(jsr,ita,it,is)/ri2(jsr,ita,it)
   22 CONTINUE
C
      fullp(1:dimr,ita,it,1:ns)=v(1:dimr,ita,it,1:ns)
C
   20 CONTINUE
C
      IF(prnt.EQ.1) THEN
         DO 30 it=1,nt
         DO 30 ita=1,nta(it)
         jrn=jrsm(ita,it)
         DO 30 is=1,ns
         WRITE(m6,100) it,ita,is
         DO 30 ir=1,jrn
         r=ri(ir,ita,it)
   30    WRITE(m6,110) r,v(ir,ita,it,is)/r/r,vp(ir,ita,it,is)/r/r
      ENDIF
C
  100 FORMAT(/,' NEWPOT:  IT = ',i3,' ITA =',i3,' IS = ',i3,
     .      //,15x,'r',15x,'v(r)',12x,'vp(r)',/)
  110 FORMAT(11x,f10.6,2f16.6)
      RETURN
      END
