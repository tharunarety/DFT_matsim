      SUBROUTINE ATOMO(DP,DQ,DR,DQ1,DFL,DV,Z,TEST,JC)
C   ******************************************************************
C   *                                                                *
C   *   Find starting values for the outward integration             *
C   *                                                                *
C   ******************************************************************
      USE message
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30)
      COMMON/ADAM/DEP(5),DEQ(5),BE,VC,TVC,DK,DM
      COMMON/INIH/DPNO(4,MO),DQNO(4,MO)
      DIMENSION DP(MW),DQ(MW),DR(MW)
      DO 20 I=1,10
      DP(I)=0.D0
   20 DQ(I)=0.D0
      DVAL=Z/VC
      DEVA1=-DVAL
      DEVA2=DV/VC+DVAL/DR(1)-BE
      DEVA3=0.D0
      IF (DK) 21,21,22
   21 DBE=(DK-DFL)/DVAL
      GO TO 23
   22 DBE=DVAL/(DK+DFL)
   23 DQ(10)=DQ1
      DP(10)=DBE*DQ1
      DO 24 I=1,5
      DP(I)=DP(10)
      DQ(I)=DQ(10)
      DEP(I)=DP(I)*DFL
   24 DEQ(I)=DQ(I)*DFL
      M=1
   25 DM=M+DFL
      DSUM=DM*DM-DK*DK+DEVA1*DEVA1
      DQR=(TVC-DEVA2)*DQ(M+9)-DEVA3*DQ(M+7)
      DPR=DEVA2*DP(M+9)+DEVA3*DP(M+7)
      DVAL=((DM-DK)*DQR-DEVA1*DPR)/DSUM
      DSUM=((DM+DK)*DPR+DEVA1*DQR)/DSUM
      J=-1
      DO 26 I=1,5
      DRT=DR(I)
      DPR=DRT**M
      DQR=DSUM*DPR
      DPR=DVAL*DPR
      IF (M.EQ.1) GO TO 27
      IF (ABS(DPR/DP(I)).LE.TEST.AND.ABS(DQR/DQ(I)).LE.TEST) J=1
   27 DP(I)=DP(I)+DPR
      DQ(I)=DQ(I)+DQR
      DEP(I)=DEP(I)+DPR*DM
   26 DEQ(I)=DEQ(I)+DQR*DM
      IF (J.EQ.1) GO TO 28
      DP(M+10)=DVAL
      DQ(M+10)=DSUM
      M=M+1
      IF (M.LE.20) GO TO 25
      WRITE(M6,'(A)')
     1 ' ATOMO:  The expansion at the origin does no converge'
      STOP
   28 DO 29 I=1,4
      DPNO(I,JC)=DP(I+9)
   29 DQNO(I,JC)=DQ(I+9)
      RETURN
      END
