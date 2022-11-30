      SUBROUTINE INVERC(NM,N,AR,AI,BR,BI)
C    **************************************************************
C    *                                                            *
C    *     SOLVES X*A=E, X AND A ARE HERMITIAN  COMPLEX N*N       *
C    *                   E IS A UNIT N*N MATRIX                   *
C    *                   ON EXIT X IS STORED IN A                 *
C    *                   L(I,I) IS REAL AND                       *
C    *                   U(I,J)=CONJG(L(J,I))/L(I,I)              *
C    *                                                            *
C    **************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N)
      DATA T0,T1/0.D0,1.D0/
      DATA IWR/0/,ICHECK/0/
C >>>>>>>>>>>>>>>>>>>>  CROUT FACTORIZATION   >>>>>>>>>>>>>>>>>>>>>>>>
C >>   A=L*U, L AND U ARE STORED IN A EXCEPT THE UNIT DIAGONAL OF U >>
      IF(ICHECK.EQ.0)GO TO 61
      DO 60 I=1,N
      DO 60 J=1,N
      BR(J,I)=AR(J,I)
   60 BI(J,I)=AI(J,I)
   61 CONTINUE
      DO 1 I=1,N
      IM1=I-1
      IF(I.EQ.1)GO TO 4
      DO 2 K=1,IM1
CDIR$ IVDEP
      DO 2 J=I,N
      AR(J,I)=AR(J,I)-(AR(J,K)*AR(K,I)-AI(J,K)*AI(K,I))
    2 AI(J,I)=AI(J,I)-(AR(J,K)*AI(K,I)+AI(J,K)*AR(K,I))
    4 CONTINUE
      IF(I.EQ.N)GO TO 1
      IP1=I+1
      XR=T1/AR(I,I)
CDIR$ IVDEP
      DO 11 J=IP1,N
      AR(I,J)=XR*AR(J,I)
   11 AI(I,J)=-XR*AI(J,I)
C     IF(IWR.EQ.1)WRITE(8,6070)I,J,AR(I,J),AI(I,J)
   1  CONTINUE
C  <<<<<<<<    END OF CROUT FACTORIZATION    <<<<<<<<<<
C  >>>>  FIND Y FROM Y*U=E ; Y STORED IN A   >>>>>>>>>>>>>
      NM1=N-1
      DO 20 I=1,NM1
      IP1=I+1
      DO 21 J=IP1,N
      AR(I,J)=-AR(I,J)
   21 AI(I,J)=-AI(I,J)
      DO 20 J=IP1,NM1
      JP1=J+1
CDIR$ IVDEP
      DO 20 K=JP1,N
      AR(I,K)=AR(I,K)-(AR(I,J)*AR(J,K)-AI(I,J)*AI(J,K))
  20  AI(I,K)=AI(I,K)-(AR(I,J)*AI(J,K)+AI(I,J)*AR(J,K))
C  <<<<<<<<<<<  Y HAS BEEN FOUND  <<<<<<<<<<<<<<<<<<<
C  >>>>>>>>>>> FIND X FROM X*L=Y, X IS STORED IN A >>>>>>
      DO 30 J=1,N
      JI=N-J+1
      XR=T1/AR(JI,JI)
      AR(JI,JI)=T1
      AI(JI,JI)=T0
      IF(J.EQ.1)GO TO 32
      JP1=JI+1
      DO 31 K=JP1,N
CDIR$ IVDEP
      DO 31 I=1,JI
      AR(I,JI)=AR(I,JI)-(AR(I,K)*AR(K,JI)-AI(I,K)*AI(K,JI))
  31  AI(I,JI)=AI(I,JI)-(AR(I,K)*AI(K,JI)+AI(I,K)*AR(K,JI))
  32  CONTINUE
      DO 33 I=1,JI
      AR(I,JI)=AR(I,JI)*XR
   33 AI(I,JI)=AI(I,JI)*XR
   30 CONTINUE
C <<<<<<<<<   X HAS BEEN FOUND   <<<<<<<<<<<<<
      DO 80 J=1,N
      JP1=J+1
CDIR$ IVDEP
      DO 80 I=JP1,N
      AR(I,J)=AR(J,I)
   80 AI(I,J)=-AI(J,I)
      IF(ICHECK.EQ.0)RETURN
      DO 70 J=1,N
      DO 70 I=1,N
      XR=T0
      XI=T0
      DO 71 K=1,N
      XR=XR+AR(I,K)*BR(K,J)-AI(I,K)*BI(K,J)
   71 XI=XI+AI(I,K)*BR(K,J)+AR(I,K)*BI(K,J)
      WRITE(8,6070)I,J,XR,XI
   70 CONTINUE
      IF(IWR.EQ.0)RETURN
      WRITE(8,6020)
      DO 50 I=1,N
  50  WRITE(8,6030)(AR(I,J),J=1,N)
      WRITE(8,6040)
      DO 51 I=1,N
  51  WRITE(8,6030)(AI(I,J),J=1,N)
 6020 FORMAT(1H ,' REAL PART OF RESULT MATRIX ')
 6030 FORMAT(1H ,8D16.9)
 6040 FORMAT(1H ,' IMAGINARY PART OF RESULT MATRIX ')
 6070 FORMAT(1H ,' I=',I3,' J=',I3,' UR=',D16.8,' UI=',D16.8)
      RETURN
      END
