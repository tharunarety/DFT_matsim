      SUBROUTINE GAUNT
C   ******************************************************************
C   *                                                                *
C   *   Calculation of Gaunt coefficients.                           *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE factorial
      USE gaunt_coeff
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(NLM2) :: CYLM
      REAL(KIND=8) :: ARG, SETN, SETD, SET3, SET4, SETF
      INTEGER :: LP, L, LA, MP, NM, M, L1, L2, KLM, ILM, JLM, MABS,
     1           LPP, MPP, MPPA, ISIGN, LPM, LIM, LM1P, LM1M, LM2P,
     2           LM2M
      ALLOCATE(C(NLM,NLM))
C
C     CYLM = (-1)**M*SQRT((2L+1)(L-M)!/(L+M)!)
C
      DO 20 LP=1,NL2
      L=LP-1
      LA=L*(L+1)/2
      NM=L*2+1
      DO 20 MP=1,LP
      M=MP-1
      L1=L-M+1
      L2=L+MP
      KLM=LA+MP
      ARG=NM*FAC(L1)/FAC(L2)
   20 CYLM(KLM)=SQRT(ARG)*(-1)**M
C
C     Calculate CYLM*GAUNT-coefficients
C
      DO 21 ILM=1,NLM
      LP=LL(ILM)
      MP=MM(ILM)
      DO 21 JLM=1,NLM
      L=LL(JLM)
      M=MM(JLM)
      MABS=ABS(M+1)
      LPP=LP+L
      MPP=MP-M
      MPPA=IABS(MPP)
      KLM=LPP*(LPP+1)/2+MPPA+1
C
C     Sign of YLM for M .LT. 0
C
      ISIGN=1
      IF(MPP.LT.0) ISIGN=(-1)**MPPA
      LPM=LPP+MPP+1
      LIM=LPP-MPP+1
      LM1P=LP+MP+1
      LM1M=LP-MP+1
      LM2P=L+M+1
      LM2M=L-M+1
      SETN=FAC(LPM)*FAC(LIM)
      SETD=FAC(LM1P)*FAC(LM1M)*FAC(LM2P)*FAC(LM2M)
      SET3=(2*LP+1)*(2*L+1)
      SET4=(2*LPP+1)
      SETF=SET3*SETN/SET4/SETD
   21 C(ILM,JLM)=((-1)**MABS)*ISIGN*CYLM(KLM)*SQRT(SETF)*2.0D0
      RETURN
      END
