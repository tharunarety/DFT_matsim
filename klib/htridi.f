      SUBROUTINE HTRIDI(NM,N,AR,AI,D,E,E2,TAU)
C
      INTEGER I,J,K,L,N,II,NM,JP1
      REAL*8 AR(NM,N),AI(NM,N),D(N),E(N),E2(N),TAU(2,N)
      REAL*8 F,G,H,FI,GI,HH,SI,SCALE
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE TRED1, NUM. MATH. 11, 181-195(1968)
C     BY MARTIN, REINSCH, AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
C
C     THIS SUBROUTINE REDUCES A COMPLEX HERMITIAN MATRIX
C     TO A REAL SYMMETRIC TRIDIAGONAL MATRIX USING
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRIX,
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX HERMITIAN INPUT MATRIX.
C          ONLY THE LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.
C
C     ON OUTPUT-
C
C        AR AND AI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION IN THEIR FULL LOWER
C          TRIANGLES.  THEIR STRICT UPPER TRIANGLES AND THE
C          DIAGONAL OF AR ARE UNALTERED,
C
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE THE TRIDIAGONAL MATRIX,
C
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,
C
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED,
C
C        TAU CONTAINS FURTHER INFORMATION ABOUT THE TRANSFORMATIONS.
C
C     ARITHMETIC IS REAL EXCEPT FOR THE USE OF THE SUBROUTINES
C     CABS AND CMPLX IN COMPUTING COMPLEX ABSOLUTE VALUES.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     HLS: 18-JAN-87
C     ------------------------------------------------------------------
C
      TAU(1,N) = 1.0D0
      TAU(2,N) = 0.0D0
C
      DO 100 I = 1, N
  100 D(I) = AR(I,I)
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********
      DO 300 II = 1, N
         I = N + 1 - II
         L = I - 1
         H = 0.0D0
         SCALE = 0.0D0
         IF (L .LT. 1) GO TO 130
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********
         DO 120 K = 1, L
  120    SCALE = SCALE + ABS(AR(I,K)) + ABS(AI(I,K))
C
         IF (SCALE .NE. 0.0D0) GO TO 140
         TAU(1,L) = 1.0D0
         TAU(2,L) = 0.0D0
  130    E(I) = 0.0D0
         E2(I) = 0.0D0
         GO TO 290
C
  140    DO 150 K = 1, L
            AR(I,K) = AR(I,K) / SCALE
            AI(I,K) = AI(I,K) / SCALE
            H = H + AR(I,K) * AR(I,K) + AI(I,K) * AI(I,K)
  150    CONTINUE
C
         E2(I) = SCALE * SCALE * H
         G = SQRT(H)
         E(I) = SCALE * G
         F = CABS(CMPLX(AR(I,L),AI(I,L)))
C     ********** FORM NEXT DIAGONAL ELEMENT OF MATRIX T **********
         IF (F .EQ. 0.0D0) GO TO 160
         TAU(1,L) = (AI(I,L) * TAU(2,I) - AR(I,L) * TAU(1,I)) / F
         SI = (AR(I,L) * TAU(2,I) + AI(I,L) * TAU(1,I)) / F
         H = H + F * G
         G = 1.0D0 + G / F
         AR(I,L) = G * AR(I,L)
         AI(I,L) = G * AI(I,L)
         IF (L .EQ. 1) GO TO 270
         GO TO 170
  160    TAU(1,L) = -TAU(1,I)
         SI = TAU(2,I)
         AR(I,L) = G
  170    F = 0.0D0
C
         DO 240 J = 1, L
            G = 0.0D0
            GI = 0.0D0
C     ********** FORM ELEMENT OF A*U **********
            DO 180 K = 1, J
               G = G + AR(J,K) * AR(I,K) + AI(J,K) * AI(I,K)
               GI = GI - AR(J,K) * AI(I,K) + AI(J,K) * AR(I,K)
  180       CONTINUE
C
            JP1 = J + 1
            IF (L .LT. JP1) GO TO 220
C
            DO 200 K = JP1, L
               G = G + AR(K,J) * AR(I,K) - AI(K,J) * AI(I,K)
               GI = GI - AR(K,J) * AI(I,K) - AI(K,J) * AR(I,K)
  200       CONTINUE
C     ********** FORM ELEMENT OF P **********
  220       E(J) = G / H
            TAU(2,J) = GI / H
            F = F + E(J) * AR(I,J) - TAU(2,J) * AI(I,J)
  240    CONTINUE
C
         HH = F / (H + H)
C     ********** FORM REDUCED A **********
         DO 260 J = 1, L
            F = AR(I,J)
            G = E(J) - HH * F
            E(J) = G
            FI = -AI(I,J)
            GI = TAU(2,J) - HH * FI
            TAU(2,J) = -GI
C
            DO 260 K = 1, J
               AR(J,K) = AR(J,K) - F * E(K) - G * AR(I,K)
     X                           + FI * TAU(2,K) + GI * AI(I,K)
               AI(J,K) = AI(J,K) - F * TAU(2,K) - G * AI(I,K)
     X                           - FI * E(K) - GI * AR(I,K)
  260    CONTINUE
C
  270    DO 280 K = 1, L
            AR(I,K) = SCALE * AR(I,K)
            AI(I,K) = SCALE * AI(I,K)
  280    CONTINUE
C
         TAU(2,L) = -SI
  290    HH = D(I)
         D(I) = AR(I,I)
         AR(I,I) = HH
         AI(I,I) = SCALE * SQRT(H)
  300 CONTINUE
C
      RETURN
      END
