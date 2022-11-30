      SUBROUTINE GAMFC(ALPHA,GLH)
C   ******************************************************************
C   *                                                                *
C   *    Calculation of convergence functions used in the Ewald      *
C   *    method.                                                     *
C   *                                                                *
C   *    GLH = GAMMA(ALPHA**2,L+1/2)/ALPHA**(L+1)                    *
C   *                                                                *
C   ******************************************************************
      USE csts
      USE message
      USE control_data
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(*) :: GLH
      REAL(KIND=8) :: ALPHA, ARG, FACL, FEX, DEL, AINV, DERFC
      INTEGER :: LP
C
      DEL=ALPHA
      ARG=ALPHA*ALPHA
      FEX=EXP(ARG)
      GLH(1)=SQRTPI*FEX*DERFC(ALPHA)
      FACL=0.5D0
C
C     Incomplete gamma function by recursion
C
      DO 20 LP=2,NL2
      GLH(LP)=DEL+FACL*GLH(LP-1)
      DEL=DEL*ARG
   20 FACL=FACL+1.0D0
      FEX=1.0D0/FEX
      AINV=1.0D0/ALPHA
      DO 21 LP=1,NL2
      FEX=FEX*AINV
   21 GLH(LP)=GLH(LP)*FEX
      RETURN
      END
