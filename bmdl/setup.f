      SUBROUTINE setup
C   ******************************************************************
C   *                                                                *
C   *    Initialize various constants and factorial arrays.          *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    FAC(N)=1*2*3*...*(N-1)                                      *
C   *    DAC(N)=1*1*3*5*...*(2*(N-1)-1)                              *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE control_data
      USE factorial
      USE gaunt_coeff
      USE madelung_lattice
      USE madelung_matrix
      IMPLICIT NONE
      INTEGER :: I, L, M, LM
C
C     Initialize various constants
C
      NLM=NL*NL
      NLMQ=NLM*NQ
C
      ALLOCATE(VMD(NQ,NQ,NLM,NLM),VMAD(NLMQ,NLMQ))
C
      NL2=2*NL-1
      NL22=NL*NL2+1
      NLM2=NL2*NL2
      RMAX=AMAX/ALAMDA
      GMAX=2.*ALAMDA*BMAX
C
C     Calculate and store factorials
C
      NFCTRL=4*NL+4
      ALLOCATE(FAC(NFCTRL),DAC(NFCTRL))
C
      FAC(1)=1.D0
      FAC(2)=1.D0
      DO I=3,NFCTRL
         FAC(I)=FAC(I-1)*(I-1)
      ENDDO
      DAC(1)=1.D0
      DAC(2)=1.D0
      DO I=3,NFCTRL
         DAC(I)=DAC(I-1)*(2*I-3)
      ENDDO
C
      ALLOCATE(LL(NLM2),MM(NLM2))
      DO L=0,NL2-1
      DO M=-L,L
         LM=L*L+L+M+1
         LL(LM)=L
         MM(LM)=M
      END DO 
      END DO 
C
      RETURN
      END
