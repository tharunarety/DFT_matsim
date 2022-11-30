      SUBROUTINE YLMRK(XX,YY,ZZ,NLMAX,YLM)
C   ******************************************************************
C   *                                                                *
C   *    Calculates spherical harmonics.                             *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    NLMAX : Calculates Ylm for l = 0,..,NLMAX.                  *
C   *    XX                                                          *
C   *    YY    : Unit vector.                                        *
C   *    ZZ                                                          *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    YLM    := (sqrt(4*Pi)/CYLM)*c.c(Ylm(x,y,z)) where CYLM are  *
C   *              constants calculated in GAUNT.                    *
C   *                                                                *
C   *                                                                *
C   *    Ylm  = YLM(lm)*CYLM/SQRT(4*PI)                              *
C   *                                                                *
C   *    CYLM = (-1)^((m+|m|)/2)*SQRT{(2*l+1)*(l-|m|)!/(l+|m|)!}     *
C   *                                                                *
C   *    Ylm  = sum(m'){CHY(lm,lm')*RYlm'}                           *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE factorial
      USE gaunt_coeff
      IMPLICIT NONE
      COMPLEX(KIND=8), DIMENSION(*) :: YLM
      REAL(KIND=8), DIMENSION(NL22) :: PLM
      REAL(KIND=8), DIMENSION(NLM2) :: COSMP, SINMP
      REAL(KIND=8) :: XX, YY, ZZ, P, X, Y, XA, TA, A, B, COSPHI,
     1                SINPHI
      INTEGER :: NLMAX, NLP, LP, L, LA, M, K, IA, IB, IC, JC, MP, J,
     1           L2, N, KLM, NM, MA, LB
C
C     Calculate Legendre polynomials by recursion
C
      NLP=NLMAX
      P=SQRT(XX*XX+YY*YY)
      X=ZZ
      Y=P
      XA=ABS(X)
      IF(XA.GT.1.D-06) GO TO 10
C
C     ABS(X) = 0
C
      DO 11 LP=1,NLP
      L=LP-1
      LA=L*(L+1)/2+1
      TA=2.D0**L
      DO 11 MP=1,LP
      M=MP-1
      K=L+M
      IF(K-2*(K/2).EQ.0) GO TO 12
      J=LA+M
      PLM(J)=0.D0
      GO TO 11
   12 IA=K+1
      IB=K/2+1
      JC=(L-M)/2
      IC=JC+1
      J=LA+M
      PLM(J)=(((-1)**JC)*FAC(IA))/(TA*FAC(IB)*FAC(IC))
   11 CONTINUE
      GO TO 32
   10 IF(XA.LT.0.999999D0) GO TO 20
C
C     ABS(X) = 1
C
      PLM(1)=1.D0
      PLM(2)=X
      DO 13 LP=3,NLP
      L=LP-1
      J=L*(L+1)/2+1
      L2=2*L-1
      K=J-L
      M=J-L2
   13 PLM(J)=(L2*X*PLM(K)-(L-1)*PLM(M))/L
      DO 14 LP=2,NLP
      L=LP-1
      LA=L*(L+1)/2
      DO 14 MP=2,LP
      J=LA+MP
   14 PLM(J)=0.D0
      GO TO 32
C
C     0 < ABS(X) < 1
C
   20 PLM(1)=1.D0
      PLM(2)=X
      PLM(3)=Y
      PLM(5)=3.D0*Y*X
      DO 21 LP=3,NLP
      L=LP-1
      J=L*(L+1)/2+1
      L2=2*L-1
      K=J-L
      M=J-L2
   21 PLM(J)=(L2*X*PLM(K)-(L-1)*PLM(M))/L
      DO 22 LP=4,NLP
      L=LP-1
      J=L*(L+1)/2+2
      L2=2*L-1
      K=J-L
      M=J-L2
   22 PLM(J)=(L2*X*PLM(K)-L*PLM(M))/(L-1)
      DO 23 LP=3,NLP
      L=LP-1
      LA=L*(L+1)/2
      DO 23 MP=3,LP
      M=MP-1
      J=LA+MP
      K=J-1
      N=K-1
      A=(M-1)*2.D0*X/Y
      B=(L+M-1)*(L-M+2)
   23 PLM(J)=A*PLM(K)-B*PLM(N)
   32 CONTINUE
C
C     Form spherical harmonics
C
      IF(P.GT.1.D-06) GO TO 34
      COSPHI=1.D0
      SINPHI=0.D0
      GO TO 35
   34 COSPHI=XX/P
      SINPHI=YY/P
   35 COSMP(1)=1.D0
      SINMP(1)=0.D0
      DO 33 MP=2,NLP
      COSMP(MP)=COSMP(MP-1)*COSPHI-SINMP(MP-1)*SINPHI
   33 SINMP(MP)=SINMP(MP-1)*COSPHI+COSMP(MP-1)*SINPHI
      KLM=0
      DO 36 LP=1,NLP
      L=LP-1
      NM=L*2+1
      DO 36 MP=1,NM
      KLM=KLM+1
      M=MM(KLM)
      MA=IABS(M)+1
      LB=L*(L+1)/2+MA
      IF(M.LE.0) GO TO 37
      YLM(KLM)=PLM(LB)*CMPLX(COSMP(MA),-SINMP(MA),KIND=8)
      GO TO 36
   37 YLM(KLM)=PLM(LB)*CMPLX(COSMP(MA),SINMP(MA),KIND=8)
   36 CONTINUE
      RETURN
      END
