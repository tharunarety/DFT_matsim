      SUBROUTINE ATOMC(WSATM,IT,ita,NS,etotc)
C   ******************************************************************
C   *                                                                *
C   *   Generate atomic charge densities for cold start              *
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
      DIMENSION FINT(MW),D(MW),DC(MW)
  101 FORMAT(/,' ATOMC:    ',2A,/)
  102 FORMAT(' ITER',4X,'DVMAX',10X,'R',14X,'VN-1',13X,'VN',10X,
     1'DEMAX',6X,'DPMAX',9X,'PN-1',13X,'PN')
  103 FORMAT(I5,1PE11.2,3(1PE16.6),2(1PE11.2),2(1PE16.6))
  104 FORMAT(/,' ATOMC:**  Too few iterations. Increase NITER =',I4)
  105 FORMAT('  NSTOP=',I4,' for orbital ',I3,A2)
  106 FORMAT(/,11X,'Orbital',3X,'Occup.',3X,'Valence',3X,'deltaQ',3x,
     .        'Final energy (Rydberg)',/)
  107 FORMAT(11X,I2,A,5X,I2,9X,I1,5X,f6.2,4x,1PE15.7)
  108 FORMAT(/,11x,'ETOTC = ',f16.6)
C
      IF(IZ.EQ.0) THEN
            IF(MSGL.EQ.1) WRITE(MSGIO,'(2A,/)')
     1            ' ATOMC: No atomic calculations for: ',TEXTA(1:4)
      RETURN
      ENDIF
      IF(MSGL.EQ.1) WRITE(MSGIO,'(2A,/)')
     1            ' ATOMC: Atomic calculations for: ',TEXTA(1:4)
      DO 19 I=1,NP
      DO 19 J=1,NORB
      DGC(I,J)=0.D0
  19  DPC(I,J)=0.D0
      WRITE(M6,101) TEXTA,TXCH1
      IF(NPRNA.EQ.1) WRITE(M6,102)
      N=-(ION+1)
C
C     Iteration loop
C*
      DO 20 ITER=1,NITER
      IF(ITER.EQ.1.OR.(ITER/5)*5.EQ.ITER) THEN
         IF(MSGL.EQ.1) WRITE(MSGIO,'(A,2I3)')
     1                ' NITER,ITER =',NITER,ITER
      ENDIF
      TETS=TEST
      YMAX=0.D0
      VMAX=0.D0
      EMAX=0.D0
      DO 21 I=1,NP
   21 D(I)=0.D0
C
C     Solve the Dirac equation for each orbital
C**
      DO 25  J=1,NORB
      DE=DEN(J)
   22 CALL DIRACA(NQN(J),NQL(J),NK(J),IMAX,DEN(J),DFL(J),DQ1(J),J)
      IF(NSTOP.NE.0) THEN
         IF(NSTOP.NE.1.OR.ITER.GE.10.OR.TETS.GT.TEST) THEN
            WRITE(M6,105) NSTOP,NQN(J),TITLE(J)
            STOP
         ENDIF
         TETS=TESTV
         GO TO 22
      ENDIF
      VAL=ABS((DEN(J)-DE)/DE)
      IF(VAL.GT.EMAX) EMAX=VAL
      NMAX(J)=IMAX
      DO 24 I=1,NP
      VAL=DGC(I,J)-DP(I)
      IF(ABS(DP(I)).GT.1.D0) VAL=VAL/DP(I)
      IF(ABS(VAL).GE.ABS(YMAX)) THEN
         YMAX=VAL
         Y=DP(I)
         YN=DGC(I,J)
      ENDIF
      VAL=DPC(I,J)-DQ(I)
      IF(ABS(DQ(I)).GT.1.D0) VAL=VAL/DQ(I)
      IF(ABS(VAL).GE.ABS(YMAX)) THEN
         YMAX=VAL
         Y=DQ(I)
         YN=DPC(I,J)
      ENDIF
      DGC(I,J)=DP(I)
      DPC(I,J)=DQ(I)
   24 D(I)=D(I)+NEL(J)*(DP(I)*DP(I)+DQ(I)*DQ(I))
   25 CONTINUE
C**
      CALL ATOMV(DC,D,DR,DX,Z,FOURPI,RWAT,NP,ION,IEX,IWAT)
C
      DO 23 I=1,NP
      DVAL=ABS(DC(I)-DV(I))
      IF((DR(I)*DC(I)).LE.N) DVAL=-DVAL/DC(I)
      IF(DVAL.LE.VMAX) GO TO 23
      VMAX=DVAL
      J=I
   23 CONTINUE
      IF(NPRNA.EQ.1) WRITE(M6,103) ITER,VMAX,DR(J),DV(J),DC(J),
     1                           EMAX,YMAX,YN,Y
C
C     Convergence test
C
      IF(TETS.LE.TEST.AND.
     1   EMAX.LE.TESTE.AND.
     2   VMAX.LE.TESTV.AND.
     3   YMAX.LE.TESTY) GO TO 28
C
C     Mix potentials for next iteration
C
      DMIX=1.D0-VMIX
      DO 26 I=1,NP
   26 DV(I)=DMIX*DV(I)+VMIX*DC(I)
   20 CONTINUE
C*
      WRITE(M6,104) NITER
      STOP
C
C     Calculation converged
C
   28 CONTINUE
      IF(MSGL.EQ.1) WRITE(MSGIO,'(A,/)')
     1            ' ATOMC: Atomic calculations converged'
      WRITE(M6,106)
      DO 29 I=1,NORB
C
C     Find the missing one-electron charge
C     
      fint(1:np)=dr(1:np)*(dgc(1:np,i)*dgc(1:np,i)+
     .           dpc(1:np,i)*dpc(1:np,i))
      CALL simpn(fint,dx,jwsa,cws)
      dcore=1.d0-cws
C
   29 WRITE(M6,107) NQN(I),TITLE(I),NEL(I),NCORB(I),dcore,2.d0*DEN(I)
C
C     Find total charge density
C
      DO 30 I=1,MW
      CHDT(I)=0.D0
   30 CHDC(I)=0.D0
      DO 31 J=1,NORB
      DO 31 I=1,NP
   31 CHDT(I)=CHDT(I)+NEL(J)*(DGC(I,J)*DGC(I,J)+DPC(I,J)*DPC(I,J))
C
C     Find core charge density
C
      DO 33 J=1,NORB
      IF(NCORB(J).EQ.0) THEN
         DO 34 I=1,NP
   34    CHDC(I)=CHDC(I)+NEL(J)*(DGC(I,J)*DGC(I,J)+DPC(I,J)*DPC(I,J))
      ENDIF
   33 CONTINUE
C
C     Calculate valence (softz=0), total (softz=2) or core (softz=1) part 
C     of the total energy
C
      IF(softz.EQ.0) THEN
         CALL ATOME(IT,ita,etotv)
      ELSEIF(softz.EQ.1) THEN
         CALL ATOMT(IT,ita,etota)
         CALL ATOME(IT,ita,etotv)
         etotc=etota-etotv
         WRITE(m6,108) etotc
      ELSEIF(softz.EQ.2) THEN
         CALL ATOMT(IT,ita,etota)
      ENDIF
C
C     Renormalize the core
C
      DO 35 IR=1,NP
   35 FINT(IR)=CHDC(IR)*DR(IR)
      CALL SIMPN(FINT,DX,JWSA,CWS)
      DCORE=Z-ICEL-CWS
      FACC=DCORE*3.D0/WSATM**3
      DO 36 IR=1,NP
   36 CHDC(IR)=CHDC(IR)+FACC*DR(IR)**2
C
      dens(ita,it,1,1:norb)=den(1:norb)
      dq1s(ita,it,1,1:norb)=dq1(1:norb)
      IF(ns.EQ.2) THEN
         dens(ita,it,2,1:norb)=den(1:norb)
         dq1s(ita,it,2,1:norb)=dq1(1:norb)
      ENDIF
C
      RETURN
      END
