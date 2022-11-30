      SUBROUTINE SMTRX
C   ******************************************************************
C   *                                                                *
C   *   Calculate the conventional unscreened ( alpha=0 ) LMTO       *
C   *   structure constant matrix.                                   *
C   *                                                                *
C   *    NPRN  : Print vectors if = 1 or 3.                          *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE control_data
      USE csts
      USE factorial
      USE gaunt_coeff
      USE madelung_lattice
      USE message
      USE s_matrix
      IMPLICIT NONE
      COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: YLM
      COMPLEX(KIND=8), DIMENSION(NLM2) :: SR, SUMR, SUMG
      COMPLEX(KIND=8) :: CIM, CCIM, CIL, SRC
      REAL(KIND=8), DIMENSION(NUMVR,NL2) :: GLHS
      REAL(KIND=8), DIMENSION(NUMVG,NL2) :: FBLN
      REAL(KIND=8), DIMENSION(NUMVR) :: USX, USY, USZ, DRQ
      REAL(KIND=8), DIMENSION(NUMVG) :: UKX, UKY, UKZ, BDOTQ
      REAL(KIND=8), DIMENSION(NL2) :: GLH
      REAL(KIND=8) :: SWSL, FPL3OM, BACKGR, ALFA, DRI, TWOLAM, X, Y, Z,
     1                D, BETA, EXPB, AR, AG, QPX, QPY, QPZ, QPPX, QPPY,
     2                QPPZ
      INTEGER :: I, LP, KLM, L, M, NM, JLM, ILM, LPP, MPP, JQ, JQQ,
     1           JLMQ, ILMQ, NQM, JQP, IQ, IQQ, J, JP
  100 FORMAT(16F6.2)
  102 FORMAT(/,11X,'QQ diagonal block')
  103 FORMAT(10(1PE11.3))
  104 FORMAT(/,11X,'QP Q off diagonal block , QP,Q =',2I5,/)
  105 FORMAT(/,25X,'Convergence test')
  106 FORMAT(/,'Reciprocal space sums',/)
  107 FORMAT(/,'Real space sums',/)
C
      ALLOCATE(SRLMQ(NLMQ,NLMQ),SILMQ(NLMQ,NLMQ),YLM(NLM2))
C
      SWSL=SWS*ALAMDA
      FPL3OM=PI*SQRTPI/ALAMDA**3./VOL
      BACKGR=PI*SWS/VOL/ALAMDA**2
      CIM=(0.D0,1.D0)
      CCIM=(0.D0,-1.D0)
C
C     QQ diagonal blocks
C
      DO 22 I=2,NR0
      DRI=DR(I)
      USX(I)=ASX(I)/DRI
      USY(I)=ASY(I)/DRI
      USZ(I)=ASZ(I)/DRI
      ALFA=ALAMDA*DRI
C
C     Convergence function, real space
C
      CALL GAMFC(ALFA,GLH)
      DO 22 LP=1,NL2
   22 GLHS(I,LP)=GLH(LP)
C
      TWOLAM=2.D0*ALAMDA
      DO 35 I=2,NUMVG
      X=AKX(I)
      Y=AKY(I)
      Z=AKZ(I)
      D=SQRT(X*X+Y*Y+Z*Z)
      IF(D.GT.GMAX) GO TO 35
      BETA=D/TWOLAM
      EXPB=EXP(-BETA*BETA)
      UKX(I)=X/D
      UKY(I)=Y/D
      UKZ(I)=Z/D
C
C     Convergence function, recip. space
C
      DO 33 LP=1,NL2
   33 FBLN(I,LP)=BETA**(LP-3)*EXPB
   35 DG(I)=D
C
C     Evaluation of lattice sums for QPP .EQ. 0
C
      SUMR=(0.D00,0.D00)
      SUMG=(0.D00,0.D00)
C
C     Real space
C
      DO 26 I=2,NR0
      CALL YLMRK(USX(I),USY(I),USZ(I),NL2,YLM)
      KLM=0
      DO 26 LP=1,NL2
      L=LP-1
      NM=L*2+1
      DO 26 M=1,NM
      KLM=KLM+1
   26 SUMR(KLM)=SUMR(KLM)+GLHS(I,LP)*YLM(KLM)
C
C     Reciprocal space
C
      DO 34 I=2,NUMVG
      IF(DG(I).GT.GMAX) GO TO 34
      CALL YLMRK(UKX(I),UKY(I),UKZ(I),NL2,YLM)
      KLM=0
      DO 28 LP=1,NL2
      L=LP-1
      NM=L*2+1
      DO 28 M=1,NM
      KLM=KLM+1
   28 SUMG(KLM)=SUMG(KLM)+FBLN(I,LP)*YLM(KLM)
   34 CONTINUE
C
C     Set up structure constant matrix from lattice sums
C
      CIL=CIM
      KLM=0
      DO 51 LP=1,NL2
      L=LP-1
      NM=L*2+1
      CIL=CIL*CCIM
      AR=SWSL**LP*2.**L/DAC(LP)/SQRTPI
      AG=FPL3OM*AR
      DO 45 M=1,NM
      KLM=KLM+1
   45 SR(KLM)=AR*SUMR(KLM)*CIL+AG*SUMG(KLM)
   51 CONTINUE
      SR(1)=SR(1)-2.0*SWSL/SQRTPI-BACKGR
C
C     Insert into lower triangle
C
      DO 27 JLM=1,NLM
      L=LL(JLM)
      M=MM(JLM)
      DO 27 ILM=JLM,NLM
      LPP=L+LL(ILM)
      MPP=MM(ILM)-M
      KLM=LPP*LPP+LPP+MPP+1
      SRC=SR(KLM)*C(ILM,JLM)
      SRLMQ(ILM,JLM)=REAL(SRC)
      SILMQ(ILM,JLM)=AIMAG(SRC)
   27 CONTINUE
C
      IF(NPRN.EQ.1) THEN
C
C        Print structure constant matrix if requested
C
         WRITE(M6,102)
         DO 42 ILM=1,NLM
   42    WRITE(M6,100) (SRLMQ(ILM,JLM),JLM=1,ILM)
         WRITE(M6,*)
         DO 43 ILM=1,NLM
   43    WRITE(M6,100) (SILMQ(ILM,JLM),JLM=1,ILM)
         WRITE(M6,*)
      ENDIF
      IF(NQ.EQ.1) GO TO 58
C
C     Repeat first diagonal block
C
      DO 24 JQ=2,NQ
      JQQ=(JQ-1)*NLM
      DO 24 JLM=1,NLM
      JLMQ=JQQ+JLM
      DO 24 ILM=JLM,NLM
      ILMQ=JQQ+ILM
      SRLMQ(ILMQ,JLMQ)=SRLMQ(ILM,JLM)
   24 SILMQ(ILMQ,JLMQ)=SILMQ(ILM,JLM)
C
C     QQP off diagonal blocks
C
      NQM=NQ-1
      DO 20 JQ=1,NQM
      JQP=JQ+1
      JQQ=(JQ-1)*NLM
      QPX=QX(JQ)
      QPY=QY(JQ)
      QPZ=QZ(JQ)
      DO 20 IQ=JQP,NQ
      IQQ=(IQ-1)*NLM
      QPPX=QX(IQ)-QPX
      QPPY=QY(IQ)-QPY
      QPPZ=QZ(IQ)-QPZ
      DO 36 I=1,NUMVR
      X=ASX(I)-QPPX
      Y=ASY(I)-QPPY
      Z=ASZ(I)-QPPZ
      D=SQRT(X*X+Y*Y+Z*Z)
      DRQ(I)=D
      IF(D.GT.RMAX) GO TO 36
      ALFA=ALAMDA*D
      USX(I)=X/D
      USY(I)=Y/D
      USZ(I)=Z/D
C
C     Convergence function, real space
C
      CALL GAMFC(ALFA,GLH)
      DO 37 LP=1,NL2
   37 GLHS(I,LP)=GLH(LP)
   36 CONTINUE
C
      DO 38 I=2,NUMVG
      D=DG(I)
      IF(D.GT.GMAX) GO TO 38
      BDOTQ(I)=(UKX(I)*QPPX+UKY(I)*QPPY+UKZ(I)*QPPZ)*D
      BETA=D/TWOLAM
      EXPB=EXP(-BETA*BETA)
C
C     Convergence function, reciprocal space
C
      DO 39 LP=1,NL2
   39 FBLN(I,LP)=BETA**(LP-3)*EXPB
   38 CONTINUE
C
C     Evaluation of lattice sums for QPP .NE. 0
C
      SUMR=(0.D00,0.D00)
      SUMG=(0.D00,0.D00)
C
C     Real space
C
      DO 40 I=1,NUMVR
      IF(DRQ(I).GT.RMAX) GO TO 40
      CALL YLMRK(USX(I),USY(I),USZ(I),NL2,YLM)
      KLM=0
      DO 21 LP=1,NL2
      L=LP-1
      NM=L*2+1
      DO 21 M=1,NM
      KLM=KLM+1
   21 SUMR(KLM)=SUMR(KLM)+GLHS(I,LP)*YLM(KLM)
   40 CONTINUE
C
C     Reciprocal space
C
      DO 41 I=2,NUMVG
      IF(DG(I).GT.GMAX) GO TO 41
      CALL YLMRK(UKX(I),UKY(I),UKZ(I),NL2,YLM)
      KLM=0
      DO 47 LP=1,NL2
      L=LP-1
      NM=L*2+1
      DO 47 M=1,NM
      KLM=KLM+1
   47 SUMG(KLM)=SUMG(KLM)+FBLN(I,LP)*YLM(KLM)*EXP(CIM*BDOTQ(I))
   41 CONTINUE
C
C     Set up structure constant matrix from lattice sums
C
      CIL=CIM
      KLM=0
      DO 56 LP=1,NL2
      L=LP-1
      NM=L*2+1
      CIL=CIL*CCIM
      AR=SWSL**LP*2.**L/DAC(LP)/SQRTPI
      AG=FPL3OM*AR
      DO 48 M=1,NM
      KLM=KLM+1
   48 SR(KLM)=AR*SUMR(KLM)*CIL+AG*SUMG(KLM)
   56 CONTINUE
      SR(1)=SR(1)-BACKGR
C
C     Insert into lower triangle
C
      DO 23 JLM=1,NLM
      JLMQ=JQQ+JLM
      L=LL(JLM)
      M=MM(JLM)
      DO 23 ILM=1,NLM
      ILMQ=IQQ+ILM
      LPP=L+LL(ILM)
      MPP=MM(ILM)-M
      KLM=LPP*LPP+LPP+MPP+1
      SRC=SR(KLM)*C(ILM,JLM)
      SRLMQ(ILMQ,JLMQ)=REAL(SRC)
   23 SILMQ(ILMQ,JLMQ)=AIMAG(SRC)
C
      IF(NPRN.EQ.1) THEN
C
C        Print structure constant matrix if requested
C
         WRITE(M6,104) IQ,JQ
         DO 31 ILMQ=IQQ+1,IQQ+NLM
   31    WRITE(M6,100) (SRLMQ(ILMQ,JLMQ),JLMQ=JQQ+1,JQQ+NLM)
         WRITE(M6,*)
         DO 32 ILMQ=IQQ+1,IQQ+NLM
   32    WRITE(M6,100) (SILMQ(ILMQ,JLMQ),JLMQ=JQQ+1,JQQ+NLM)
      ENDIF
   20 CONTINUE
C
C     Fill in upper triangle
C
   58 DO 57 J=1,NLMQ-1
      JP=J+1
      DO 57 I=JP,NLMQ
      SRLMQ(J,I)=SRLMQ(I,J)
   57 SILMQ(J,I)=-SILMQ(I,J)
      DEALLOCATE(YLM)
C
      RETURN
      END
