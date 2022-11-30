      SUBROUTINE nxtsws(isws,efgs,hx)
C   ******************************************************************
C   *                                                                *
C   *    Prepare volume (Sws) dependent input.                       *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE control_data
      USE control_text
      USE csts
      USE density
      USE lattice
      USE message
      USE moments
      USE pota
      USE potential
      USE radialmesh
      USE text
      USE volume
      IMPLICIT NONE
      INTEGER :: isws, in, ixch, ir, it, ita, iq, is
      REAL(KIND=8) :: efgs, hx, hx0, bohr, alata
      REAL(KIND=8) :: dxa, dr1a
      SAVE hx0
C
      sws=swsl(isws)
C
C     Guess next Fermi level
C
      IF(isws.EQ.1) THEN
         hx0=hx
      ELSEIF(isws.EQ.2) THEN
         efgs=efl(1)
         hx=hx0
      ELSE
         efgs=2.d0*efl(isws-1)-efl(isws-2)
         hx=hx0
      ENDIF
C
C     Determine the individual atomic-sphere and potential-sphere radii
C
      CALL scaws
C
      IF(strt.EQ.'A') THEN
C
C        Prepare and perform atomic calculations
C
         CALL atomp(efgs)
C
      ELSE
         CALL atomd(txch,dxa,dr1a)
      ENDIF
C
      IF(strt.EQ.'A'.OR.strt.EQ.'C') THEN
         in=3
         CALL oldin(efgs,ixch,in)
      ELSEIF(strt.EQ.'B'.OR.strt.EQ.'N') THEN
         in=2
         CALL oldin(efgs,ixch,in)
      ENDIF
C
      ixc=ixch+1
C
C     Set charge transfer and redefine ELN as the number of conduction
C     electrons for the neutral atom
C
      qtr=qtro
      eln=eln+ion
      elt=0.d0
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
   20 elt=elt+conc(ita,it)*mmt(it)*eln(ita,it)
C
C     Get average ASA radius SWS.
C
      vol=0.d0
      DO 25 it=1,nt
      DO 25 ita=1,nta(it)
   25 vol=vol+conc(ita,it)*mmt(it)*ws(ita,it)**3*fourpi/3.d0
      sws=(3.d0*vol/nq/fourpi)**(1.d0/3.d0)
      bohr=0.529177d0
      alata=bohr*alat
C
      CALL setngb
C
      WRITE(m6,100)
      DO 26 iq=1,nq
      it=itq(iq)
      DO 26 ita=1,nta(it)
   26 WRITE(m6,101) iq,it,ita,mmt(it),ttxt(ita,it),nz(ita,it),
     .              ion(ita,it),eln(ita,it),qtro(ita,it),
     .              split(ita,it),fixst(ita,it),conc(ita,it)
      WRITE(m6,102) sws,alat,alata,vol,elt
C
C     Generate r-mesh
C
      ALLOCATE(ri(dimr,mnta,nt),ri2(dimr,mnta,nt))
      DO 27 iq=1,nq
      it=itq(iq)
      DO 27 ita=1,nta(it)
      DO 28 ir=1,dimr
   28 ri(ir,ita,it)=r1(ita,it)*exp(dx*(ir-1))
      IF(ri(jwsc(ita,it),ita,it).LT.wsc(iq)) THEN
         WRITE(m6,110) ri(jwsc(ita,it),ita,it),wsc(iq)
         STOP
      ENDIF
      IF(ri(jwsi(ita,it),ita,it).GT.wsi(iq)) THEN
         WRITE(m6,111) ri(jwsi(ita,it),ita,it),wsi(iq)
         STOP
      ENDIF
      ri2(1:dimr,ita,it)=ri(1:dimr,ita,it)*ri(1:dimr,ita,it)
   27 CONTINUE
C
C     Set up the operlapping neighbourhood
C
      CALL ovrlps(1)
C
  100 FORMAT(/,' NXTSWS:',3x,'IQ IT ITA MMT Type  NZ ION ',
     .       ' ELN   QTR   SPLIT  FIX  CONC',/)
  101 FORMAT(10x,3i3,i4,3x,a,2i3,2x,f4.1,2f7.3,3x,a,f7.3)
  102 FORMAT(/,11x,'SWS =',f10.4,' Alat =',f8.4,' Bohr ',f8.4,' AA',
     .       /,11x,'VOL =',f10.4,' ELT  =',f8.4)
  110 FORMAT(/,' NXTSWS:***  RI(JSC) = ',f10.6,' SC = ',f10.6)
  111 FORMAT(/,' NXTSWS:***  RI(JSI) = ',f10.6,' SI = ',f10.6)
      RETURN
      END
