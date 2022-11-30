      FUNCTION CLBSHG(LPP,LP,L,MPP,MP,M)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the Clebsch-Gordan coefficients as given for      *
C   *    instance by Tinkham Eq. (5-47).                             *
C   *                                                                *
C   ******************************************************************
      USE clebsch_gordon
      IMPLICIT NONE
      REAL(KIND=8) :: CLBSHG, CLEBSH, X, Y, Z
      INTEGER :: LPP, LP, L, MPP, MP, M, IPARF, I, K1, K2, K3, K4, K5,
     1           K6, KEY, IJ, IJPAR, M1, M2, M3, M4, M5, M6, M7, M8,
     2           M9, M10, IX, IY, IZ, JQ, N4, N5, N5PAR, MM1, MM2, MM3,
     3           MM4, MM5
      IPARF(I)=4*(I/4)-I+1
C
      K1=LPP+LPP
      K2=LP+LP
      K3=L+L
      KEY=ABS(MPP)+ABS(MP)+ABS(M)
      IF(KEY.EQ.0) THEN
C
C        All magnetic quantum numbers are zero
C
         IF((K1-(K1/2)-(K1/2)).NE.0) GO TO 420
         IF((K2-(K2/2)-(K2/2)).NE.0) GO TO 420
         IF((K3-(K3/2)-(K3/2)).NE.0) GO TO 420
         IJ=K1+K2+K3
         IJPAR=IPARF(IJ)
         IF(IJPAR.LE.0) GO TO 420
         M1=IJ-K1-K1
         M2=IJ-K2-K2
         M3=IJ-K3-K3
         IF(M1.LT.0) GO TO 420
         IF(M2.LT.0) GO TO 420
         IF(M3.LT.0) GO TO 420
         M1=M1/2+1
         M2=M2/2+1
         M3=M3/2+1
         M4=IJ/2+2
         Y=SQRT(H(M1)*H(M2)*H(M3)/H(M4))
         IY=(J(M1)+J(M2)+J(M3)-J(M4))/2
         IJ=IJ/2
         IJPAR=IPARF(IJ)
         IJ=IJ/2+1
         M1=M1/2+1
         M2=M2/2+1
         M3=M3/2+1
         Z=H(IJ)/(H(M1)*H(M2)*H(M3))
         IZ=J(IJ)-J(M1)-J(M2)-J(M3)
         IZ=IZ+IY
         CLEBSH=IJPAR*Y*Z*10.0D0**IZ
         JQ=K2-K1
         IF(JQ.LT.0) JQ=-JQ
         IJPAR=IPARF(JQ)
         CLBSHG=CLEBSH*IJPAR*SQRT(K3+1.D0)
         RETURN
  420    CLBSHG=0.D0
         RETURN
      ELSE
C
C        Non-zero magnetic quantum numbers
C
         K4=MPP+MPP
         K5=MP+MP
         K6=M+M
         IF((K4+K5-K6).NE.0) GO TO 710
         M1=K1+K2-K3
         M2=K2+K3-K1
         M3=K3+K1-K2
         M4=K1+K4
         M5=K1-K4
         M6=K2+K5
         M7=K2-K5
         M8=K3+K6
         M9=K3-K6
         M10=K1+K2+K3+2
         IF(M1.LT.0) GO TO 710
         IF(M2.LT.0) GO TO 710
         IF(M3.LT.0) GO TO 710
         IF(M4.LT.0) GO TO 710
         IF(M5.LT.0) GO TO 710
         IF(M6.LT.0) GO TO 710
         IF(M7.LT.0) GO TO 710
         IF(M8.LT.0) GO TO 710
         IF(M9.LT.0) GO TO 710
         IF((M4-(M4/2)-(M4/2)).NE.0) GO TO 710
         IF((M10-(M10/2)-(M10/2)).NE.0) GO TO 710
         Y=K3+1
         M1=M1/2+1
         M2=M2/2+1
         M3=M3/2+1
         M4=M4/2+1
         M5=M5/2+1
         M6=M6/2+1
         M7=M7/2+1
         M8=M8/2+1
         M9=M9/2+1
         M10=M10/2+1
         Y=SQRT(Y*H(M1)*H(M2)*H(M3)*H(M4)*H(M5)
     1     *H(M6)*H(M7)*H(M8)*H(M9)/H(M10))
         IY=(J(M1)+J(M2)+J(M3)+J(M4)+J(M5)
     1     +J(M6)+J(M7)+J(M8)+J(M9)-J(M10))/2
         N4=M1
         IF(N4.GT.M5)N4=M5
         IF(N4.GT.M6)N4=M6
         N4=N4-1
         M2=K2-K3-K4
         M3=K1+K5-K3
         N5=0
         IF(N5.LT.M2) N5=M2
         IF(N5.LT.M3) N5=M3
         N5PAR=IPARF(N5)
         N5=N5/2
         Z=0.0D0
         GO TO 610
  700    MM1=M1-N5
         MM2=M5-N5
         MM3=M6-N5
         MM4=N5-(M2/2)+1
         MM5=N5-(M3/2)+1
         X=1./(H(MM1)*H(MM2)*H(MM3)*H(MM4)*H(MM5)*H(N5+1))
         IX=-J(MM1)-J(MM2)-J(MM3)-J(MM4)-J(MM5)-J(N5+1)
  800    IF(IX+IY)900,210,110
  900    X=0.1D0*X
         IX=IX+1
         GO TO 800
  110    X=10.0D0*X
         IX=IX-1
         GO TO 800
  210    IF(N5PAR.LT.0) X=-X
         Z=Z+X
  510    N5PAR=-N5PAR
         N5=N5+1
  610    IF(N5-N4)700,700,810
  710    CLEBSH=0.0D0
         GO TO 220
  810    CLEBSH=Z*Y
  220    CLBSHG=CLEBSH
         RETURN
      ENDIF
      END
