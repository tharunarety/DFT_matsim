C   ******************************************************************
C   *                                                                *
C   *                         *** BMDL ***                           *
C   *                                                                *
C   *          Calculation of Madelung potentials for bulk.          *
C   *                                                                *
C   *                   H.L.Skriver and L. Vitos                     *
C   *                      Physics Department                        *
C   *                Technical University of Denmark                 *
C   *                        DK-2800 Lyngby                          *
C   *                                                                *
C   *                          J. Kollar                             *
C   *             Central Research Institute for Physics             *
C   *                       H-1525 Budapest                          *
C   *                                                                *
C   ******************************************************************
      PROGRAM BMDL
      USE basis              ; USE control_data    ; USE control_text
      USE madelung_lattice   ; USE madelung_matrix ; USE message
      IMPLICIT NONE
      CHARACTER(LEN=9) :: DATO
      CHARACTER(LEN=5) :: CLOCK
      REAL(KIND=8) :: CMDL
      INTEGER :: IQ, JQ, LM, LMP, I, J, LP, L, NLMAX
C
      VERSION='LV  05 Mar 1999'
C
  101 FORMAT(/,' BMDL:',4X,' NPRN  =',I3,/,11X,'NL    =',I3,
     1     ' NQ    =',I3,' NLM   =',I3,' NLMQ  =',I3,/,11X,
     2     'MSGL  =',I3,//,11X,'AMAX  =',F10.6,' BMAX  =',F10.6,
     3      ' ALAMDA=',F10.6,/,11X,'RMAX  =',F10.6,' GMAX  =',F10.6)
  102 FORMAT(3F10.6,A4)
  103 FORMAT(/,11X,'Madelung-potential matrices stored on FOR001')
  110 FORMAT(16F6.2)
  111 FORMAT(/,11X,'QP Q block , QP,Q =',2I5,/)
C
      CALL FHNDLR(DATO,CLOCK)
      CALL SETCST
C
      IF(MSGL.EQ.1) WRITE(MSGIO,'(A,I5,/)') ' BMDL:  Message level =',
     1              MSGL
C
      READ(5,'(8X,I2,10X,F10.6)') NL
      READ(5,'(3(10X,F10.6))') ALAMDA,AMAX,BMAX
C
      CALL PRIMV
      CALL PRIMKR
      CALL SET3D
      CALL SETUP
      CALL GAUNT
C
      WRITE(M6,101) NPRN,NL,NQ,NLM,NLMQ,MSGL,AMAX,BMAX,ALAMDA,RMAX,
     1              GMAX
      CLOSE(5)
C
      IF(MSGL.EQ.1) WRITE(MSGIO,'(A,I3,/)')
     1             '  Ewald technique for l <=',NL-1
C
      CALL LATT3M
      CALL SMTRX
      CALL MADL3M
C
      WRITE(1) TXT
      WRITE(1) NQ,NL
      DO 20 IQ=1,NQ
      DO 20 JQ=1,NQ
      DO 30 LM=1,NLM
      I=(IQ-1)*NLM+LM
      DO 30 LMP=1,NLM
      J=(JQ-1)*NLM+LMP
   30 VMD(IQ,JQ,LM,LMP)=VMAD(I,J)
   20 WRITE(1) ((VMD(IQ,JQ,L,LP),L=1,NLM),LP=1,NLM)
C
      CMDL=0.D0
      DO 21 IQ=1,NQ
      DO 21 JQ=1,NQ
   21 CMDL=CMDL-VMD(IQ,JQ,1,1)
      CMDL=0.5D0*CMDL/NQ
C
      WRITE(M6,'(/,A,F16.10)') ' BMDL:    CMDL=',CMDL
      IF(NPRN.EQ.3) THEN
        NLMAX=NLM
        WRITE(M6,'(/,A,/)') ' BMDL:    Madelung matrice'
        DO 40 IQ=1,NQ
        DO 40 JQ=1,NQ
        WRITE(M6,111) IQ,JQ
        DO 40 LM=1,NLMAX
   40   WRITE(M6,110) (VMD(IQ,JQ,LM,LMP),LMP=1,LM)
      ENDIF
      WRITE(M6,103)
C
      STOP
      END

