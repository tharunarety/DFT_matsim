      SUBROUTINE SETPLW(PTL,PMAX,NLAMP,NPA,NQA)
C   ******************************************************************
C   *                                                                *
C   *    Establish the width of a principal layer from the data      *
C   *    generated in LAYER and test against NPA from LAYER.         *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    BSX                                                         *
C   *    BSY   : 3D translation vectors.                             *
C   *    BSZ                                                         *
C   *    PMAX  : z-coordinate of last layer.                         *
C   *    NLAMP : Number of layers in a primitive principal layer.    *
C   *    QX                                                          *
C   *    QY    : 3D basis vectors.                                   *
C   *    QZ                                                          *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    NPA   : Number of atomic layers in a principal layer.       *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE message
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: TOL=.5D-6
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::QX2, QY2, QZ2, PZ
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE ::UQX, UQY, UQZ
      INTEGER, DIMENSION(:), ALLOCATABLE :: IQQ, IQB, LAMQX, LAMQB
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDX
      REAL(KIND=8), DIMENSION(*) :: PTL
      REAL(KIND=8) :: PMAX, BX, BY, BZ, UQXI, UQYI, UQZI, DQI, TOLSRF,
     1                TOLI, TQZ
      INTEGER :: NLAMP, NQA, NPA, I, LAM, IQ, N, ISRF, IQSRF, INX,
     1           NQAZ, IPW, IQA, LAMSRF, ILM, NPAX
  101 FORMAT(/,' SETPLW:   Dimension of principal layer =',I3)
  102 FORMAT(/,'           Basis vectors in 2D cell',//,11X,
     1       'IQA IQ  LAM INCL',5X,'QX',8X,'QY',8X,'QZ',/)
  103 FORMAT(11X,I2,3I4,1X,3F10.6)
  104 FORMAT(/,' SETPLW:** Some interface layers neglected in the',
     1       ' cluster inversion.',/,11X,'Cannot proceed.')
C
C     Establish the interface layers
C
      I=0
      DO 19 LAM=-NLAMP-2,NLAMP+2
      BX=LAM*BSX(3)
      BY=LAM*BSY(3)
      BZ=LAM*BSZ(3)
      DO 19 IQ=1,NQ
      UQXI=BX+QX(IQ)
      UQYI=BY+QY(IQ)
      UQZI=BZ+QZ(IQ)
      DQI=ABS(UQZI)
      IF(DQI.LT.PMAX) THEN
         I=I+1
      ENDIF
   19 CONTINUE
      N=I
      ALLOCATE(QX2(N),QY2(N),QZ2(N),IQQ(N),LAMQX(N),PZ(N),
     1         UQX(N),UQY(N),UQZ(N),IQB(N),LAMQB(N),INDX(N))

      I=0
      DO 20 LAM=-NLAMP-2,NLAMP+2
      BX=LAM*BSX(3)
      BY=LAM*BSY(3)
      BZ=LAM*BSZ(3)
      DO 20 IQ=1,NQ
      UQXI=BX+QX(IQ)
      UQYI=BY+QY(IQ)
      UQZI=BZ+QZ(IQ)
      DQI=ABS(UQZI)
      IF(DQI.LT.PMAX) THEN
         I=I+1
         IQB(I)=IQ
         LAMQB(I)=LAM
         UQX(I)=UQXI
         UQY(I)=UQYI
         UQZ(I)=UQZI
         PZ(I)=UQZI
      ENDIF
   20 CONTINUE
C
      CALL QSORT(PZ,INDX,N)
C
C     Find the surface layer defined as the first layer with zero PZ
C
      DO 21 I=1,N
      IF(ABS(PZ(I)).LT.TOL) THEN
         ISRF=I
      ENDIF
   21 CONTINUE
C
C     Find the dimension of the principal layer
C
      IQSRF=IQB(INDX(ISRF))
      DO 28 I=ISRF+1,N
      INX=INDX(I)
      IF(IQB(INX).EQ.IQSRF) NQAZ=I
   28 CONTINUE
      NPA=NQAZ-ISRF
C
C     Now NPA is the number of atomic layers in a principal layer
C
      WRITE(M6,101) NPA
C
      IPW=0
      DO 22 I=ISRF-NPA+1,ISRF+NPA
      IPW=IPW+1
      INX=INDX(I)
      IQQ(IPW)=IQB(INX)
      LAMQX(IPW)=LAMQB(INX)
      QX2(IPW)=UQX(INX)
      QY2(IPW)=UQY(INX)
      QZ2(IPW)=UQZ(INX)
   22 CONTINUE
      NPW=IPW
C
C     Find the surface layer and shift LAMQX
C
      TOLSRF=10.*TOL
      DO 23 IQA=1,NPW
      IF(QZ2(IQA).LT.TOLSRF) THEN
         LAMSRF=LAMQX(IQA)
      ENDIF
   23 CONTINUE
      DO 24 IQA=1,NPW
   24 LAMQX(IQA)=LAMQX(IQA)-LAMSRF
C
C     Are the all the NPA layers used in the matrix inversion?
C
      ALLOCATE(incl(npw))
      TOLI=10.*TOL
      DO 25 IQA=1,NPW
      TQZ=QZ2(IQA)
      INCL(IQA)=0
      DO 25 ILM=1,NQA
      IF(ABS(TQZ-PTL(ILM)).LT.TOLI) THEN
         INCL(IQA)=1
      ENDIF
   25 CONTINUE
C
C     Print the layers generated
C
      WRITE(M6,102)
      DO 26 IQA=1,NPW
      IQ=IQQ(IQA)
      WRITE(M6,103) IQA,IQ,LAMQX(IQA),INCL(IQA),QX2(IQA),QY2(IQA),
     1              QZ2(IQA)
      IF(IQA.EQ.NPW/2) WRITE(M6,'(11X,A)')
     1       '-----------------Surface---------------------'
   26 CONTINUE
      NPAX=0
      DO IQA=1,NPA
      NPAX=NPAX+INCL(IQA)
      ENDDO
      IF(NPAX.NE .NPA) THEN
         WRITE(M6,104)
         STOP
      ENDIF
      RETURN
      END
