      SUBROUTINE DBAR(PH,TH,GM,L,MY,MYP,DB)
C   ******************************************************************
C   *                                                                *
C   *   Calculate the real matrix elements of finite rotation        *
C   *   for Euler angles (ph,th,gm) and quantum numbers (l,my,myp),  *
C   *   as defined in A.R. Edmons:Angular Momentum in Quantum        *
C   *                             Mechanics(Princeton 1957).         *
C   ******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(SQ2=1.414213562D0)
      DIMENSION DR(2),DI(2)
C
      MZP=ABS(MYP)
      MZ=ABS(MY)
      CALL GEND(PH,TH,GM,L,MZ,MZP,DR,DI)
C
      ONE=(-1.D0)**(MY+MYP)
      ONEP=(-1.D0)**MYP
      IF(MY.GT.0.AND.MYP.GT.0) THEN
         DB=ONE*(DR(2)+ONEP*DR(1))
      ELSEIF(MY.GT.0.AND.MYP.LT.0) THEN
         DB=ONE*(DI(2)-ONEP*DI(1))
      ELSEIF(MY.LT.0.AND.MYP.GT.0) THEN
         DB=ONE*(-DI(2)-ONEP*DI(1))
      ELSEIF(MY.LT.0.AND.MYP.LT.0) THEN
         DB=ONE*(DR(2)-ONEP*DR(1))
      ELSEIF(MY.GT.0.AND.MYP.EQ.0) THEN
         DB=SQ2*ONE*DR(2)
      ELSEIF(MY.LT.0.AND.MYP.EQ.0) THEN
         DB=-SQ2*ONE*DI(2)
      ELSEIF(MY.EQ.0.AND.MYP.GT.0) THEN
         DB=SQ2*ONE*DR(2)
      ELSEIF(MY.EQ.0.AND.MYP.LT.0) THEN
         DB=SQ2*ONE*DI(2)
      ELSEIF(MY.EQ.0.AND.MYP.EQ.0) THEN
         DB=DR(2)
      ENDIF
C
      RETURN
      END
      SUBROUTINE GEND(PH,TH,GM,L,MP,M,DR,DI)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION DR(2),DI(2)
C
      DO 30 IMP=1,2
      IF(IMP.EQ.1) MZ=-M
      IF(IMP.EQ.2) MZ=M
C
      DL=FDL(L,MP,MZ,TH)
      DR(IMP)=DL*DCOS(MP*PH+MZ*GM)
      DI(IMP)=DL*DSIN(MP*PH+MZ*GM)
   30 CONTINUE
      RETURN
      END
C
      REAL*8 FUNCTION FDL(J,MP,M,TH)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MFC=170,PI=3.141592653D0)
      COMMON/RFAC/FAC(MFC)
C
      DTH=DABS(TH)
      IF(DTH.LE.1.D-8) THEN
        IF(MP.EQ.M) THEN
          FDL=1.D0
        ELSE
          FDL=0.D0
        ENDIF
      ELSEIF(DTH.GE.PI) THEN
        IF((MP+M).EQ.0) THEN
          FDL=(-1.D0)**(J+M)
        ELSE
          FDL=0.D0
        ENDIF
      ELSE
        TH2=TH/2.D0
        SUM=0.D0
        DO 10 ISG=0,J-MP
        CTH=(DCOS(TH2))**(2*ISG+MP+M)
        STH=(DSIN(TH2))**(2*J-2*ISG-MP-M)
        ONE=(-1.D0)**(J-MP-ISG)
        C1=COMB(J+M,J-MP-ISG)
        C2=COMB(J-M,ISG)
   10   SUM=SUM+C1*C2*ONE*CTH*STH
C
        COEF=FAC(J+MP+1)*FAC(J-MP+1)/FAC(J+M+1)/
     1       FAC(J-M+1)
        FDL=DSQRT(COEF)*SUM
      ENDIF
C
      RETURN
      END
C
      REAL*8 FUNCTION COMB(N,K)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MFC=170)
      COMMON/RFAC/FAC(MFC)
C
      IF(K.LT.0.OR.K.GT.N) THEN
        COMB=0.D0
        RETURN
      ELSEIF(K.EQ.0.OR.K.EQ.N) THEN
        COMB=1.D0
        RETURN
      ENDIF
      COMB=FAC(N+1)/FAC(K+1)/FAC(N-K+1)
C
      RETURN
      END
      SUBROUTINE ROTFACT
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MFC=170)
      COMMON/RFAC/FAC(MFC)
C
      FAC(1)=1.0D00
      DO 10 I=1,MFC-1
   10 FAC(I+1)=I*FAC(I)
C
      RETURN
      END
