      SUBROUTINE F03AHF(N, A, IA, DETR, DETI, IDETE, RINT, IFAIL)
C   ******************************************************************
C   * COMPDET  MARK 2 RELEASE. NAG COPYRIGHT 1972                    *
C   *          MARK 11 REVISED. VECTORISATION (JAN 1984).            *
C   *          DOUBLE PRECISION VERSION WITH COMPLEX*16              *
C   * THE COMPLEX UNSYMMETRIC MATRIX, A, IS STORED IN THE ARRAY      *
C   * A(N,N). THE DECOMPOSITION A = LU, WHERE L IS A LOWER TRIANGU-  *
C   * LAR MATRIX AND U IS A UNIT UPPER TRIANGULAR MATRIX, IS PER-    *
C   * FORMED AND OVERWRITTEN ON A, OMITTING THE UNIT DIAGONAL OF U.  *
C   * A RECORD OF ANY INTERCHANGES MADE TO THE ROWS OF A IS KEPT IN  *
C   * RINT(I), I=1,N, SUCH THAT THE I-TH ROW AND THE RINT(I)-TH ROW  *
C   * WERE INTERCHANGED AT THE I-TH STEP. THE DETERMINANT,           *
C   * (DETR + I * DETI) * 2.0**IDETE, OF A IS ALSO COMPUTED. THE     *
C   * SUBROUTINE WILL FAIL IF A, MODIFIED BY THE ROUNDING ERRORS, IS *
C   * SINGULAR.                                                      *
C   * 1ST DECEMBER 1971 , HLS: 18-JAN-87.                            *
C   ******************************************************************
      INTEGER IFAIL, I, N, IA, IDETE, K, L, J, P01AAF
      DOUBLE PRECISION SRNAME
      COMPLEX*16 CZ, A(IA,N)
      DOUBLE PRECISION DIMAG,DREAL,DETR, DETI, Z, X, Y, W, RINT(N)
C     ..
      DATA SRNAME /8H F03AHF /
      DO 5 I = 1,N
         RINT(I) = 0.0D0
    5 CONTINUE
      DO 15 J=1,N
         DO 10 I=1,N
            RINT(I) = RINT(I) + DREAL(A(I,J))**2 + DIMAG(A(I,J))**2
   10    CONTINUE
   15 CONTINUE
      DO 20 I=1,N
         IF (RINT(I).LE.0.0D0) GO TO 200
   20 CONTINUE
      DETR = 1.0D0
      DETI = 0.0D0
      IDETE = 0
      DO 180 K=1,N
         L = K
         Z = 0.0D0
         DO 40 I=K,N
            X = DREAL(A(I,K))
            Y = DIMAG(A(I,K))
            X = (X*X+Y*Y)/RINT(I)
            IF (X.LE.Z) GO TO 40
            Z = X
            L = I
   40    CONTINUE
         IF (L.EQ.K) GO TO 80
         DETR = -DETR
         DETI = -DETI
         DO 60 J=1,N
            CZ = A(K,J)
            A(K,J) = A(L,J)
            A(L,J) = CZ
   60    CONTINUE
         RINT(L) = RINT(K)
   80    RINT(K) = L
         X = DREAL(A(K,K))
         Y = DIMAG(A(K,K))
         Z = X*X + Y*Y
         W = X*DETR - Y*DETI
         DETI = X*DETI + Y*DETR
         DETR = W
         IF (DABS(DETR).LE.DABS(DETI)) W = DETI
         IF (W.EQ.0.0D0) GO TO 200
  100    IF (DABS(W).LT.1.0D0) GO TO 120
         W = W*0.0625D0
         DETR = DETR*0.0625D0
         DETI = DETI*0.0625D0
         IDETE = IDETE + 4
         GO TO 100
  120    IF (DABS(W).GE.0.0625D0) GO TO 140
         W = W*16.0D0
         DETR = DETR*16.0D0
         DETI = DETI*16.0D0
         IDETE = IDETE - 4
         GO TO 120
  140    IF (K.LT.N) CALL F03AHZ(A, IA, N, K, A(1,K+1))
  180 CONTINUE
      IFAIL = 0
      RETURN
  200 IDETE = 0
      IFAIL = P01AAF(IFAIL,1,SRNAME)
      RETURN
      END
