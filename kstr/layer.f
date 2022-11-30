      SUBROUTINE LAYER(NPA)
C   ******************************************************************
C   *                                                                *
C   *    Sort vectors into layers and find the number of layers.     *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE lattice
      USE message
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: TOL=.5D-6
      REAL(KIND=8), DIMENSION(NR) :: P
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: PTZ, PTL
      REAL(KIND=8) :: DBK3, UKX, UKY, UKZ, UKMAG, TOLP, RDOTN, PA, PMAX
      INTEGER, DIMENSION(NR) :: INDEX, LAMBDA, LIQA
      INTEGER, DIMENSION(:), ALLOCATABLE :: LAMQ, INDX, NIQA
      INTEGER :: NPA, IR, IQA, NQA, ISRF, IQAS, INX, LAM, NLAM, IR1, IR2
  101 FORMAT(/,' LAYER:    Surface normal :  (',F10.5,',',F10.5,',',
     1       F10.5,' )')
  102 FORMAT(' ',9X,I3,2I4,I3,2I4,F10.6,2X,3F10.6)
  103 FORMAT(//,11X,'IQA LAM IQ JQ  IR INX',6X,'P',10X,'TX',8X,'TY',
     1       8X,'TZ')
  104 FORMAT(/,11X,'Surface layer at  IQA =',I4,'   NQA =',I4,
     1      '  NLAM =',I4)
  105 FORMAT(/,12X,'----------------------', ' Surface ',
     1       '-----------------------')
  106 FORMAT(/,' LAYER:**  NQA = 0. Increase DMAX.')
C
      DBK3=SQRT(BKX(3)**2+BKY(3)**2+BKZ(3)**2)
      UKX=BKX(3)/DBK3
      UKY=BKY(3)/DBK3
      UKZ=BKZ(3)/DBK3
      UKMAG=UKX*UKX+UKY*UKY
      WRITE(M6,101) UKX,UKY,UKZ
      IF(UKMAG.GT.TOL) THEN
         WRITE(M6,'(/,2A)') ' LAYER:**  The two first primitive',
     1                    ' translations are not in the surface'
         STOP
      ENDIF
C
C     Sort vectors into layers
C
      DO 20 IR=1,NR
      TOLP=(JQBAS(IR)-IQBAS(IR))*TOL
      RDOTN=BKX(3)*RX(IR)+BKY(3)*RY(IR)+BKZ(3)*RZ(IR)
      P(IR)=RDOTN/DBK3+SIGN(2.D0*TOLP,RDOTN)
   20 CONTINUE
      CALL QSORT(P,INDEX,NR)
      PA=-100.
      IQA=0
      DO 21 IR=1,NR
      IF(P(IR)-PA.GT.TOL) IQA=IQA+1
      LIQA(IR)=IQA
      PA=P(IR)
   21 CONTINUE
      PMAX=(1.+TOL)*PA
      NQA=IQA
      IF(NQA.EQ.0) THEN
         WRITE(M6,106)
         STOP
      ENDIF
      ALLOCATE(PTZ(NQA),PTL(NQA),LAMQ(NQA),INDX(NQA),NIQA(0:NQA))
C
C     Find surface layer
C
      DO 22 IR=1,NR
   22 IF(P(IR).LT.TOL) ISRF=IR
      IQAS=LIQA(ISRF)
      NIQA(0)=0
      DO 23 IR=1,NR
      IQA=LIQA(IR)
   23 NIQA(IQA)=IR
C
C     Find all layers PTL and primitive layers PTZ
C
      DO 24 IQA=1,NQA
      IR=NIQA(IQA)
      PTL(IQA)=P(IR)
      INX=INDEX(IR)
   24 PTZ(IQA)=TZ(INX)
      CALL QSORT(PTZ,INDX,NQA)
      PA=-100.
      LAM=0
      DO 25 IQA=1,NQA
      IF(PTZ(IQA)-PA.GT.TOL) THEN
         LAM=LAM+1
      ENDIF
      LAMQ(IQA)=LAM
      PA=PTZ(IQA)
   25 CONTINUE
      NLAM=LAM
C
      DO 26 IR=1,NR
      IQA=LIQA(IR)
      INX=INDX(IQA)
   26 LAMBDA(IR)=LAMQ(INX)
C
C     Print the layers generated
C
      WRITE(M6,103)
      DO 30 IQA=1,NQA
      IF(IQA.EQ.IQAS+1) WRITE(M6,105)
      IR1=NIQA(IQA-1)+1
      IR2=NIQA(IQA)
      WRITE(M6,'(A)') ' '
      DO 30 IR=IR1,IR2
      INX=INDEX(IR)
   30 WRITE(M6,102) LIQA(IR),LAMBDA(IR),IQBAS(INX),JQBAS(INX),IR,INX,
     1             P(IR),TX(INX),TY(INX),TZ(INX)
C
C     Establish the width of a principal layer
C
      WRITE(M6,104) IQAS,NQA,NLAM
      CALL SETPLW(PTL,PMAX,NLAM,NPA,NQA)
C
      DEALLOCATE(PTZ,PTL,LAMQ,INDX,NIQA)
      RETURN
      END

