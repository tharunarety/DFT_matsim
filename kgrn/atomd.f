      SUBROUTINE ATOMD(TXCH,DXA,DR1A)
C   ******************************************************************
C   *                                                                *
C   *   Set default values for the atomic calculation.               *
C   *                                                                *
C   ******************************************************************
      USE message
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30,ML=4)
      PARAMETER(ONE=1)
      CHARACTER TEXTA*30,TITLE*4,TXCH*3,TXCH1*20,TXCH2(10)*20,LABEL*5
      COMMON/ADAM/DEP(5),DEQ(5),BE,VC,TVC,DK,DM
      COMMON/ATMM/DEN(MO),DQ1(MO),DFL(MO),DGC(MW,MO),DPC(MW,MO),RMAX,
     2            DR1,FOURPI,RWAT,CHDT(MW),CHDC(MW),DEXV,DEXE,VMIX,
     3            TESTE,TESTY,TESTV,TETS,NODEA(ML),NCORB(MO),NQN(MO),
     4            NQL(MO),NK(MO),NMAX(MO),NEL(MO),NORB,IWAT,IEX,JWSA,
     5            JWOA,JRIA,NITER,IZ,ION,ICEL,NPRNA
      COMMON/DIRA/DV(MW),DR(MW),DP(MW),DQ(MW),DX,Z,NSTOP,NES,TEST,NP
      COMMON/POISA/A,C,C2,F1,EDL
      COMMON/TXTA/TEXTA,TITLE(MO),TXCH1
      DATA TXCH2/'Barth-Hedin         ','X-Alpha             ',
     1           'Bath-Hedin-Janak    ','Vosko-Wilk-Nusair-CA',
     2           'Perdew-Wang-CA      ','Wigner XC           ',
     3           'Perdew-Zunger-CA    ','Perdew-B-E GGA      ',
     4           'Local Airy Gas      ','Perdew et al. PBEsol'/
  101 FORMAT(4(10X,F10.6))
  102 FORMAT(/,11X,'Exchange-correlation: ',A)
C
C     Set default parameters
C
      PI=ACOS(-ONE)
      FOURPI=4.D0*PI
      NSTOP=MO
      DEXV=1.D0
      DEXE=1.5D0
      VC=137.0373D0
      TVC=2.D0*VC
C
C     Read control parameters
C
      REWIND 8
   20 READ(8,'(A)',END=22) LABEL
      IF(LABEL.NE.'Atom:') THEN
         GO TO 20
      ELSE
         GO TO 21
      ENDIF
   22 WRITE(M6,'(A)') ' ATOM:** No header for Atom:'
      STOP
C
   21 READ(8,'(6(7X,I3))') IEX,NP,NES,NITER,IWAT,NPRNA
      READ(8,101) VMIX,RWAT,RMAX
      TXCH1=TXCH2(IEX+1)
C
      READ(8,101) DX,DR1,TEST
      DXA=DX
      DR1A=DR1
      READ(8,101) TESTE,TESTY,TESTV
C
C     Set constants for POISON
C
      DX2=DX*DX
      EDL=EXP(DX/2.D0)
      C=DX2/6.D0
      AA=0.25D0
      A=1.D0-DX2*AA/12.D0
      B=-2.D0-5.D0*DX2*AA/6.D0
      C2=-B/A
      F1=EXP(DX*SQRT(AA))
C
C     Set exchange-correlation
C
      CALL SETXCP(IEX+1,TXCH,DEXV,1)
      WRITE(M6,102) TXCH
      RETURN
      END
