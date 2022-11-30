      SUBROUTINE fcdupd(ef)
C   ******************************************************************
C   *                                                                *
C   *    Save charge densities and total energies.                   *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE density    ; USE force        ; USE lattice
      USE message    ; USE moments      ; USE pota
      USE potential  ; USE radialmesh   ; USE softcore
      USE text       ; USE totalenergy
      IMPLICIT NONE
      REAL(KIND=8) :: ef, elnit, eonets
      INTEGER :: iq, it, ita, is, lm, jri, jrc, jsr, ir
C
      WRITE(10) 'KGRN',nq,lmaxh+1,lmaxf+1,dimr,mnta,alphmd,ef
      WRITE(10) sws,ns,bsx,bsy,bsz,qx,qy,qz
      WRITE(10) dato
C
      DO 20 iq=1,nq
      it=itq(iq)
      WRITE(10) nta(it)
      DO 21 ita=1,nta(it)
C
      txtp(ita,it)='Potential parameters from KGRN for'//ttxt(ita,it)
      txtp(ita,it)(59:69)=dato//'  '
      WRITE(10) txtp(ita,it),ttxt(ita,it)
      jri=jris(ita,it)
      jrc=jwsc(ita,it)
      elnit=eln(ita,it)
      IF(softc.EQ.'Y') elnit=nz(ita,it)
      WRITE(10) r1(ita,it),dx,elnit,qtr(ita,it),qcpa(ita,it),ws(ita,it),
     .          nz(ita,it),ion(ita,it),jrc,jri,jwss(ita,it),ixc,
     .          conc(ita,it),dexch(ita,it)
      eonets=eone(ita,it)-ents(ita,it)
      WRITE(10) eonets,emadl(ita,it),vint(ita,it),vintc(ita,it),
     .          enuc(ita,it),ecor(ita,it),eval(ita,it),exct(ita,it),
     .          excc(ita,it),ekin(ita,it),etot(ita,it),okae(ita,it)
      IF(softc.EQ.'Y') THEN
         WRITE(10) (0.d0,ir=1,jri)
      ELSE
         WRITE(10) cor(ita,it,1:jri)
      ENDIF
C
      DO 22 is=1,ns
      DO 23 lm=1,nlmh
      WRITE(10) chd0(ita,iq,is,lm,1:jrc)
   23 CONTINUE
      DO 24 lm=1,nlmf
      WRITE(10) chdl(ita,iq,is,lm,1:jrc)
   24 CONTINUE
   22 CONTINUE
C
   21 CONTINUE
   20 CONTINUE
C
      DO 30 iq=1,nq
      it=itq(iq)
      DO 30 ita=1,nta(it)
      jri=jris(ita,it)
      DO 30 is=1,ns
   30 WRITE(10) chde(ita,it,is,1:jri)
C
      DO 40 iq=1,nq
      it=itq(iq)
      DO 40 ita=1,nta(it)
   40 WRITE(10) jwsi(ita,it),fsx(iq),fsy(iq),fsz(iq),shf
C
      DO 50 is=1,ns
      WRITE(10) vmtz(is)
      DO 50 iq=1,nq
      it=itq(iq)
      DO 50 ita=1,nta(it)
      jsr=jsrs(ita,it)
      WRITE(10) jsr,hsr(ita,it),wsm(ita,it),vmtzr(it,is)
   50 WRITE(10) v(1:jsr,ita,it,is)
C
      CALL htimer(dato,clock)
      WRITE(m6,'(/,2a,4x,a)')
     .      ' Full charge density stored on:',dato
      IF(msgl.EQ.1) WRITE(msgio,'(/,a)')
     .         ' FCDUPD: Full charge density stored on FOR010'
C
      RETURN
      END
