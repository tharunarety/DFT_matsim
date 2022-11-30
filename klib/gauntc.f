      FUNCTION GAUNTC(K,L,M,LP,MP)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the Gaunt coefficients ck(lm,l'm') as defined     *
C   *    for instance by Tinkham Eq. (6-37a).                        *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      REAL(KIND=8) :: GAUNTC, FAC, AM, A0, CLBSHG
      INTEGER :: K, L, M, LP, MP, MK
      MK=M-MP
      FAC=(1.D0+LP+LP)/(1.D0+L+L)
      AM=CLBSHG(K,LP,L,MK,MP,M)
      A0=CLBSHG(K,LP,L,0,0,0)
      GAUNTC=SQRT(FAC)*A0*AM
      RETURN
      END

