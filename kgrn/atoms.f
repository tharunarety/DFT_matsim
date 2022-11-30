      SUBROUTINE ATOMS(IT,ita,nt,mnta,JWSS,JRIS,R1)
C   ******************************************************************
C   *                                                                *
C   *    Save atomic charge densities and total energies to be used  *
C   *    as a start of the selfconsistent band calculations.         *
C   *                                                                *
C   ******************************************************************
      USE atomb
      USE message
      USE text
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30,ML=4)
      CHARACTER TEXTA*30,TITLE*4,TXCH1*20
      COMMON/ATMM/DEN(MO),DQ1(MO),DFL(MO),DGC(MW,MO),DPC(MW,MO),RMAX,
     2            DR1,FOURPI,RWAT,CHDT(MW),CHDC(MW),DEXV,DEXE,VMIX,
     3            TESTE,TESTY,TESTV,TETS,NODEA(ML),NCORB(MO),NQN(MO),
     4            NQL(MO),NK(MO),NMAX(MO),NEL(MO),NORB,IWAT,IEX,JWSA,
     5            JWOA,JRIA,NITER,IZ,ION,ICEL,NPRNA
      COMMON/DIRA/DV(MW),DR(MW),DP(MW),DQ(MW),DX,Z,NSTOP,NES,TEST,NP
      COMMON/TXTA/TEXTA,TITLE(MO),TXCH1
      DIMENSION JWSS(mnta,nt),JRIS(mnta,nt),R1(mnta,nt)
C
C     Store potential parameters
C
      HELN=ICEL
      IXCH=IEX
      IF(it.EQ.1.AND.ita.EQ.1) THEN
         DUM=0.D0
         REWIND 3
         WRITE(3) 'ATOM',KNQ,KNL,'N','ASA'
         WRITE(3) HSWS,HWST(1:KNQ),HWSI(1:KNQ),HWSC(1:KNQ),DX
         WRITE(3) DUM,HEFGS,KNS,KNT,KITQ(1:KNQ),knta(1:knt)
         WRITE(3) DUM,DUM,DUM,HDEXCH
         DO 24 ITT=1,KNT
         DO 24 itta=1,knta(itt)
         WRITE(3) HWS(itta,ITT),HHSR(itta,ITT),R1(itta,ITT),
     .            JWSS(itta,ITT),JRIS(itta,ITT),TTXT(itta,ITT)
   24    CONTINUE
      ENDIF
      TXTP(ita,IT)=TEXTA
      WRITE(3) TXTP(ita,IT)
      WRITE(3) 'CHRD:ATOM ',DATO
      WRITE(3) HELN,HQTR(ita,IT),DUM,DUM,DUM,IZ,ION,IXCH,JRIA
      WRITE(3) DUM,DUM,DUM,DUM,DUM,DUM,DUM,DUM,DUM,DUM
      WRITE(3) CHDC(1:JRIA)
      IF(KNS.EQ.2) THEN
         DO 27 IR=1,JRIA
   27    CHDT(IR)=0.5D0*CHDT(IR)
      ENDIF
      DO 28 IS=1,KNS
      WRITE(3) CHDT(1:JRIA)
   28 WRITE(3) DV(1:JRIA)
      WRITE(M6,'(/,2A)') ' Atomic charge densities,',
     1                   ' and potential stored on FOR003'
      IF(MSGL.EQ.1) WRITE(MSGIO,'(/,A)')
     1         ' ATOMS: Potential etc. stored on FOR003'
      RETURN
      END
