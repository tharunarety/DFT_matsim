      SUBROUTINE DIRACA(NQN,NQL,NK,IMAX,DE,DFL,DQ1,JC)
C   ******************************************************************
C   *                                                                *
C   *    Solve the Dirac equation for each orbital                   *
C   *                                                                *
C   *    NQN : Principal quantum number                              *
C   *    NQL : Orbital quantum number                                *
C   *    NK  : Combined spin 1/2 and orbital quantum number          *
C   *    IMAX: Last mesh point                                       *
C   *    DE  : Energy                                                *
C   *    DFL : Exponent of the first term in the expansion           *
C   *    DQ1 : Slope at the origin                                   *
C   *    NES : Maximun number of adjustments of the energy           *
C   *    TEST: Precission on the energy                              *
C   *                                                                *
C   ******************************************************************
      USE message
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30)
      PARAMETER(ONE=1,DKOEF=ONE/720)
      COMMON/ADAM/DEP(5),DEQ(5),BE,VC,TVC,DK,DM
      COMMON/DIRA/DV(MW),DR(MW),DP(MW),DQ(MW),DX,Z,NSTOP,NES,TEST,NP
      COMMON/INIH/DPNO(4,MO),DQNO(4,MO)
      NSTOP=0
      IMM=0
      IES=0
      DK=NK
      LLL=(NQL*(NQL+1))/2
      ND=0
      NOEUD=NQN-NQL
      IF (LLL.NE.0) GO TO 20
      ELIM=-Z*Z/(1.5D0*NQN*NQN)
      GO TO 22
   20 ELIM=DV(1)+LLL/(DR(1)*DR(1))
      DO 23 I=2,NP
      VAL=DV(I)+LLL/(DR(I)*DR(I))
      IF (VAL.LE.ELIM) ELIM=VAL
   23 CONTINUE
      IF(ELIM.GT.0.D0) THEN
         WRITE(M6,'(A)')
     1    ' DIRACA:  2*V+L*(L+1)/R**2 is positive everywhere'
         STOP
      ENDIF
C
   22 IF(DE.LE.ELIM) DE=ELIM*0.5D0
   21 IF (IMM.EQ.1) GO TO 25
      DO 24 I=7,NP,2
      IMAT=NP+1-I
      IF((DV(IMAT)+LLL/(DR(IMAT)*DR(IMAT))-DE).LE.0.D0) GO TO 26
   24 CONTINUE
   26 IF (IMAT.GT.5) GO TO 25
      DE=DE*0.5D0
      IF(DE.LT.-TEST.AND.ND.LE.NOEUD) GO TO 21
      WRITE(M6,'(A)')
     1 ' DIRACA:  2*V+L*(L+1)/R**2-2*E is positive everywhere'
      STOP
C
C     Starting values for outward integration
C
   25 BE=DE/VC
      CALL ATOMO(DP,DQ,DR,DQ1,DFL,DV(1),Z,TEST,JC)
      ND=1
      DO 27 I=1,5
      DVAL=DR(I)**DFL
      IF(I.EQ.1) GO TO 28
      IF(DP(I-1).EQ.0.D0) GO TO 28
      IF((DP(I)/DP(I-1)).GT.0.D0) GO TO 28
      ND=ND+1
   28 DP(I)=DP(I)*DVAL
      DQ(I)=DQ(I)*DVAL
      DEP(I)=DEP(I)*DVAL
   27 DEQ(I)=DEQ(I)*DVAL
      K=-1+2*(NOEUD-2*(NOEUD/2))
      IF((DP(1)*K).GT.0.D0) GO TO 29
   30 WRITE(M6,'(A)')
     1 ' DIRACA:  Error in the expansion at the origin'
      STOP
C
   29 IF((K*NK*DQ(1)).LT.0.D0) GO TO 30
      DM=DX*DKOEF
C
      DO 31 I=6,IMAT
      DP(I)=DP(I-1)
      DQ(I)=DQ(I-1)
      CALL ADAM5(DP(I),DQ(I),DV(I),DR(I))
      IF(DP(I-1).EQ.0.D0) GO TO 31
      IF((DP(I)/DP(I-1)).GT.0.D0) GO TO 31
      ND=ND+1
      IF(ND.GT.NOEUD) GO TO 32
   31 CONTINUE
      IF(ND.EQ.NOEUD) GO TO 33
      DE=0.8D0*DE
      IF(DE.LT.-TEST) GO TO 21
      WRITE(M6,'(A)')
     1 ' DIRACA:  The number of nodes too small'
      STOP
C
   32 DE=1.2D0*DE
      IF(DE.GT.ELIM) GO TO 21
      WRITE(M6,'(A)')
     1 ' DIRACA:  The number of nodes too large'
      STOP
C
C     Starting values for inward integration
C
   33 DQM=DQ(IMAT)
      DPM=DP(IMAT)
      IF (IMM.EQ.1) GO TO 34
      DO 35 I=1,NP,2
      IMAX=NP+1-I
      IF(((DV(IMAX)-DE)*DR(IMAX)*DR(IMAX)).LE.300.D0) GO TO 34
   35 CONTINUE
   34 DD=SQRT(-DE*(2.D0+BE/VC))
      DPQ=-DD/(TVC+BE)
      DM=-DM
      DO 36 I=1,5
      J=IMAX+1-I
      DP(J)=EXP(-DD*DR(J))
      DEP(I)=-DD*DP(J)*DR(J)
      DQ(J)=DPQ*DP(J)
   36 DEQ(I)=DPQ*DEP(I)
      M=IMAX-5
C
      DO 37 I=IMAT,M
      J=M+IMAT-I
      DP(J)=DP(J+1)
      DQ(J)=DQ(J+1)
   37 CALL ADAM5(DP(J),DQ(J),DV(J),DR(J))
C
C     Match the large component
C
      DVAL=DPM/DP(IMAT)
      IF(DVAL.GT.0.D0) GO TO 38
      WRITE(M6,'(A)')
     1 ' DIRACA:  Error in the sign of the large component'
      STOP
C
   38 DO 39 I=IMAT,IMAX
      DP(I)=DP(I)*DVAL
   39 DQ(I)=DQ(I)*DVAL
C
C     Calculate the normalization
C
      DSUM=3.D0*DR(1)*(DP(1)**2+DQ(1)**2)/(DX*(DFL+DFL+1.D0))
      DO 40 I=3,IMAX,2
      DSUM=DSUM+DR(I)*(DP(I)**2+DQ(I)**2)+4.D0*DR(I-1)*(DP(I-1)**2
     1    +DQ(I-1)**2)+DR(I-2)*(DP(I-2)**2+DQ(I-2)**2)
   40 CONTINUE
      DSUM=DX*(DSUM+DR(IMAT)*(DQM*DQM-DQ(IMAT)*DQ(IMAT)))
     1    *(1.D0/3.D0)
C
C     Energy correction
C
      DBE=DP(IMAT)*(DQM-DQ(IMAT))*VC/DSUM
      IMM=0
      VAL=ABS(DBE/DE)
      IF(VAL.LE.TEST) GO TO 41
   42 DVAL=DE+DBE
      IF(DVAL.LT.0.D0) GO TO 43
      DBE=DBE*0.5D0
      VAL=VAL*0.5D0
      IF(VAL.GT.TEST) GO TO 42
      WRITE(M6,'(A)')
     1 ' DIRACA:  Zero energy'
      STOP
C
   43 DE=DVAL
      IF(VAL.LE.0.1) IMM=1
      IES=IES+1
      IF(IES.LE.NES) GO TO 21
      WRITE(M6,'(A)')
     1 ' DIRACA:  Too many iterations'
      NSTOP=1
      RETURN
C
   41 DSUM=SQRT(DSUM)
      DQ1=DQ1/DSUM
      DO 44 I=1,IMAX
      DP(I)=DP(I)/DSUM
   44 DQ(I)=DQ(I)/DSUM
      DO 45 I=1,4
      DPNO(I,JC)=DPNO(I,JC)/DSUM
   45 DQNO(I,JC)=DQNO(I,JC)/DSUM
      IF(IMAX.EQ.NP) GO TO 48
      J=IMAX+1
      DO 47 I=J,NP
      DP(I)=0.D0
   47 DQ(I)=0.D0
   48 NSTOP=0
      RETURN
      END
