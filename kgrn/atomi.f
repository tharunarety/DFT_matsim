      SUBROUTINE ATOMI(WSATM,R1IT,JWSAIT,NZIT,IT,ita)
C   ******************************************************************
C   *                                                                *
C   *   Read atomic data                                             *
C   *                                                                *
C   ******************************************************************
      USE softcore
      USE message
      USE text
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30,ML=4)
      CHARACTER SYMBOL*4,CONFIG*24,TORBT*4,TXAL*1,TEXTA*30,TITLE*4,
     1          TXCH1*20
      COMMON/ADAM/DEP(5),DEQ(5),BE,VC,TVC,DK,DM
      COMMON/ATMM/DEN(MO),DQ1(MO),DFL(MO),DGC(MW,MO),DPC(MW,MO),RMAX,
     2            DR1,FOURPI,RWAT,CHDT(MW),CHDC(MW),DEXV,DEXE,VMIX,
     3            TESTE,TESTY,TESTV,TETS,NODEA(ML),NCORB(MO),NQN(MO),
     4            NQL(MO),NK(MO),NMAX(MO),NEL(MO),NORB,IWAT,IEX,JWSA,
     5            JWOA,JRIA,NITER,IZ,ION,ICEL,NPRNA
      COMMON/DIRA/DV(MW),DR(MW),DP(MW),DQ(MW),DX,Z,NSTOP,NES,TEST,NP
      COMMON/TXTA/TEXTA,TITLE(MO),TXCH1
      DIMENSION TORBT(9),TXAL(5),NQNBL(5)
      DATA NQNBL/1,2,3,4,5/
      DATA TORBT/'s1/2','p3/2','p1/2','d5/2','d3/2','f7/2','f5/2',
     1           'g9/2','g7/2'/,TXAL/'s','p','d','f','g'/
  101 FORMAT(5X,30I3)
  102 FORMAT(/,11X,'Valence orbital Nodes',/)
  103 FORMAT(/,' ATOMI:    Atomic calculation for ',A,' WSATM =',F10.6,
     1       //,11X,'NITER =',I3,' IEX   =',I3,' IWAT  =',I3,
     2       ' IZ    =',I3,' ION   =',I3)
  104 FORMAT(/,11X,'Orbital',3X,'Occup.',3X,'Valence',7X,'Start energ',
     1       'y (Hartree)',/)
  105 FORMAT (11X,I2,A,5X,I2,9X,I1,8X,1PE15.7)
  106 FORMAT(/,11X,'Precision on energy, wavefunction, and potential:',
     1       //,11X,'TESTE =',1PE10.2,' TESTY =',E10.2,' TESTV =',
     2       E10.2,/,11X,'VMIX  =',0PF10.6)
  108 FORMAT(/,11X,'Radial mesh given by (R(1)=DR1/IZ) : ',//,11X,
     1       'JWSA  =',I4,'   JRIA  =',I4,'   JWOA  =',I4,'   NP    =',
     2       I4,/,11X,'DR1   =',F10.6,' DX    =',F10.6,' RMAX  =',
     3       F10.6,/,11X,'R(1)  =',E15.7,9X,'R(NP) =',E15.7)
  109 FORMAT(/,11X,'Precision on energy and number of iterations in',
     1       ' DIRACA:',//,11X,'TEST  =',1PE10.2,' NES   =',I3)
  110 FORMAT(/,11X,'Exchange-correlation: ',A)
  111 FORMAT(11X,'Slater X-alpha: DEXV  = ',F10.6,' DEXE  =',F10.6)
  112 FORMAT(/,11X,'Radius of Watson sphere, RWAT =',F10.6)
  113 FORMAT(/,' ATOMI:**  ',A,'. Probable cause: NORB too large.')
  114 FORMAT(/,' ATOMI:**   El. count',I3,' does not add up to',
     1       ' IZ-ION =',I3)
  115 FORMAT(/,' ATOMI:**  NORB =',I3,' is too large.')
  116 FORMAT(/,' ATOMI:**   NP =',I3,' is too large. Increase MW')
  117 FORMAT(/,' ATOMI:**  Wrong configuration')
  118 FORMAT(/,' ATOMI:**  Input error on the line:',E15.8,I1,2I2)
  119 FORMAT(/,' ATOMI:**  ',A,' =',F10.6,' too small')
  120 FORMAT(/,' ATOMI:**  Radius of circumscribed sphere exceeds',
     1       ' the radial mesh. ',A,I5,' MW =',I5,
     2       '. Increase MW')
C
C     Read orbital data
C
      READ(8,'(A)') SYMBOL
      symbols(ita,it)=symbol
      READ(8,'(3X,I4,6X,I3,5X,I3,8X,A)') IZ,NORB,ION,CONFIG
      izs(ita,it)=iz
      norbs(ita,it)=norb
      ions(ita,it)=ion
      configs(ita,it)=config
      IF(NZIT.NE.IZ) THEN
         WRITE(M6,'(/,2(A,I5))') 'ATOMI:** NZ =',NZIT,' IZ =',IZ
         STOP
      ENDIF
      WRITE(M6,103) SYMBOL,WSATM,NITER,IEX,IWAT,IZ,ION
      TEXTA(1:30)=SYMBOL//'  '//CONFIG
      IF(NORB.GT.MO) THEN
         WRITE(M6,115) NORB
         STOP
      ENDIF
      Z=IZ
      READ(8,101) (NQN(I),I=1,NORB)
      READ(8,101) (NK(I),I=1,NORB)
      READ(8,101) (NEL(I),I=1,NORB)
      READ(8,101) (NCORB(I),I=1,NORB)
      nqns(ita,it,1:norb)=nqn(1:norb)
      nks(ita,it,1:norb)=nk(1:norb)
      nels(ita,it,1:norb)=nel(1:norb)
      ncorbs(ita,it,1:norb)=ncorb(1:norb)
C
C     Find number of conduction electrons
C
      ICEL=0
      DO 31 J=1,NORB
      IF(NCORB(J).EQ.1) THEN
         ICEL=ICEL+NEL(J)
      ENDIF
   31 CONTINUE
C
C     Generate radial mesh
C
      IF(WSATM.LT.0.1) THEN
         WRITE(M6,119) 'WSATM',WSATM
         STOP
      ENDIF
      IF(RMAX.LT.0.1) THEN
         WRITE(M6,119) 'RMAX',RMAX
         STOP
      ENDIF
      R1=R1IT
      JWSA=JWSAIT
      NP=1.D0+LOG(RMAX/R1)/DX
      NP=2*(NP/2)+1
      IF(NP.GT.MW) THEN
         WRITE(M6,116) NP
         STOP
      ENDIF
      JRIA=JWSA+2
C
      IF(JRIA.GT.NP) THEN
         WRITE(M6,120) ' JRIA =',JRIA,NP
         STOP
      ENDIF
      DR(1)=R1
      DO 20 I=2,NP
  20  DR(I)=DR(1)*EXP((I-1)*DX)
C
      WRITE(M6,106) TESTE,TESTY,TESTV,VMIX
      WRITE(M6,109) TEST,NES
      WRITE(M6,108) JWSA,JRIA,JWOA,NP,DR1,DX,RMAX,DR(1),DR(NP)
      WRITE(M6,110) TXCH1
      IF(IEX.EQ.1) WRITE(M6,111) DEXV,DEXE
      IF(IWAT.EQ.1) WRITE(M6,112) RWAT
      DO 21 I=1,NORB
      IF(NQN(I).EQ.0) THEN
         WRITE(M6,113) 'NQN = 0'
         STOP
      ELSEIF(NK(I).EQ.0) THEN
         WRITE(M6,113) 'NK = 0'
         STOP
      ELSEIF(NEL(I).EQ.0) THEN
         WRITE(M6,113) 'NEL = 0'
         STOP
      ENDIF
   21 CONTINUE
C
      IF(NPRNA.EQ.1) WRITE(M6,104)
      K=0
      DVAL=Z*Z/(VC*VC)
      DO 22 I=1,NORB
      K=K+NEL(I)
      DEN(I)=-Z*Z/(4.D0*NQN(I)*NQN(I))
      NQL(I)=IABS(NK(I))
      IF(NK(I).LT.0) NQL(I)=NQL(I)-1
      DFL(I)=NK(I)*NK(I)
      DFL(I)=SQRT(DFL(I)-DVAL)
      L=2*IABS(NK(I))
      IF(NQL(I).GE.NQN(I).AND.NEL(I).GT.L.AND.NQN(I).LE.0.AND.NQL(I)
     1   .GT.4) THEN
         WRITE(M6,118) DEN(I),NQN(I),NQL(I),J,NEL(I)
         STOP
      ENDIF
      J=NQL(I)+IABS(NK(I))
      TITLE(I)=TORBT(J)
      IF(NPRNA.EQ.1) WRITE(M6,105) NQN(I),TITLE(I),NEL(I),NCORB(I),
     1              DEN(I)
   22 CONTINUE
C
C     Find the number of nodes in the valence wavefunctions
C
      DO 19 IL=1,ML
   19 NODEA(IL)=0
      DO 23 I=1,NORB
      DO 23 LP=1,ML
      L=LP-1
      IF(NQL(I).EQ.L) THEN
         IF(NCORB(I).EQ.1) THEN
            NODEA(LP)=NQN(I)-NQL(I)-1
            NQNBL(LP)=NQN(I)
         ELSE
            NODEA(LP)=NQN(I)-NQL(I)
            NQNBL(LP)=NQN(I)+1
         ENDIF
      ENDIF
   23 CONTINUE
      WRITE(M6,102)
      DO 24 LP=1,ML
   24 WRITE(M6,'(17X,I1,A,8X,I3)') NQNBL(LP),TXAL(LP),NODEA(LP)
C
      ITOTEL=IZ-ION
      IF (K.NE.ITOTEL) THEN
         WRITE(M6,114) K,ITOTEL
         STOP
      ENDIF
C
C     Initialize charge densities and potential for empty sphere
C
      IF(IZ.EQ.0) THEN
         DO 30 IP=1,NP
         DV(IP)=0.D0
         CHDT(IP)=0.D0
   30    CHDC(IP)=0.D0
         RETURN
      ENDIF
C
C     Generate start potential
C
      VAL=-ION-1
      DO 25 I=1,NP
      R=DR(I)
   25 DV(I)=THFPOT(R,Z,VAL)
      IF(IWAT.EQ.1) THEN
         DO 26 I=1,NP
         IF(DR(I).LE.RWAT) THEN
            DV(I)=DV(I)+ION/RWAT
         ELSE
            DV(I)=DV(I)+ION/DR(I)
         ENDIF
   26    CONTINUE
      ENDIF
      IF(NORB.GT.1) THEN
         DO 27 I=2,NORB
         K=I-1
         DO 27 J=1,K
         IF(NQN(I).NE.NQN(J).OR.NK(I).NE.NK(J)) GO TO 27
         WRITE(M6,117)
         STOP
   27    CONTINUE
      ENDIF
      DO 28 I=1,NORB
      NMAX(I)=NP
      L=1
      J=NQN(I)-NQL(I)
      IF((J-2*(J/2)).EQ.0) L=-L
      NKAI=IABS(NK(I))
      DQ1(I)=(L*NK(I)/NKAI)
   28 CONTINUE
      RETURN
      END
