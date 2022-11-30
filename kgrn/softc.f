      SUBROUTINE softc(WSATM,corn,icore,IT,ita,is,ns)
C   ******************************************************************
C   *                                                                *
C   *   Recalculate the core density for soft core.                  *
C   *                                                                *
C   ******************************************************************
      USE softcore
      USE message
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
      DIMENSION FINT(MW),D(MW),corn(mw)
  101 FORMAT(/,' SOFTC:    ',2A,/)
  104 FORMAT(/,' SOFTC:**  Too few iterations. Increase NITER =',I4)
  105 FORMAT('  NSTOP=',I4,' for orbital ',I3,A2)
  106 FORMAT(/,11X,'Orbital',3X,'Occup.',3X,'Valence',7X,'Final energ',
     1       'y (Rydberg)',/)
  107 FORMAT(11X,I2,A,5X,I2,9X,I1,8X,1PE15.7)
C
      IF(IZ.EQ.0) RETURN
C
      DO 19 I=1,NP
      DO 19 J=1,NORB
      DGC(I,J)=0.D0
  19  DPC(I,J)=0.D0
      N=-(ION+1)
C
      TETS=TEST
      EMAX=0.D0
      DO 21 I=1,NP
   21 D(I)=0.D0
C
C     Solve the Dirac equation for core orbitals
C
      DO 25  J=1,NORB
      IF(NCORB(J).NE.0) GO TO 25
      DE=DEN(J)
   22 CALL DIRACA(NQN(J),NQL(J),NK(J),IMAX,DEN(J),DFL(J),DQ1(J),J)
      IF(NSTOP.NE.0) THEN
         IF(NSTOP.NE.1.OR.TETS.GT.TEST) THEN
            WRITE(M6,105) NSTOP,NQN(J),TITLE(J)
            STOP
         ENDIF
         TETS=TESTV
         GO TO 22
      ENDIF
      VAL=ABS((DEN(J)-DE)/DE)
      IF(VAL.GT.EMAX) EMAX=VAL
      DO 24 I=1,NP
      DGC(I,J)=DP(I)
   24 DPC(I,J)=DQ(I)
   25 CONTINUE
C
C     Convergence test
C
      IF(TETS.LE.TEST.AND.EMAX.LE.TESTE) THEN
         icore=1
      ELSE
         icore=0
      ENDIF
C
C     Find new core charge density
C
      IF(ns.EQ.2) THEN
         facns=0.5d0 
      ELSE
         facns=1.0d0 
      ENDIF
      DO 33 J=1,NORB
      IF(NCORB(J).EQ.0) THEN
         DO 34 I=1,NP
   34    corn(I)=corn(I)+facns*NEL(J)*
     .           (DGC(I,J)*DGC(I,J)+DPC(I,J)*DPC(I,J))
      ENDIF
   33 CONTINUE
C
C     Renormalize the core
C
      IF(is.EQ.ns) THEN
         DO 35 IR=1,NP
   35    FINT(IR)=corn(IR)*DR(IR)
         CALL SIMPN(FINT,DX,JWSA,CWS)
         DCORE=Z-ICEL-CWS
         FACC=DCORE*3.D0/WSATM**3
         DO 36 IR=1,NP
   36    corn(IR)=corn(IR)+FACC*DR(IR)**2
      ENDIF
C
C     Calculate the new core one-electron energies
C
      DO 37 J=1,NORB
      IF(NCORB(J).EQ.0) THEN
         eonec(ita,it)=eonec(ita,it)+2.d0*facns*NEL(J)*DEN(J)
      ENDIF
   37 CONTINUE
C
      dens(ita,it,is,1:norb)=den(1:norb)
      dq1s(ita,it,is,1:norb)=dq1(1:norb)
      RETURN
      END
