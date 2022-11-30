      SUBROUTINE LATT3D
C   ******************************************************************
C   *                                                                *
C   *    Generate a lattice in real space from the primitive vectors *
C   *    of the 3D cell. For each Q' labelled by IQ the routine      *
C   *    determines the set of vectors T+Q-Q' shorter than DMAX and  *
C   *    sort these into shells of neighbours.                       *
C   *                                                                *
C   *   *On entry (BSX,BSY,BSZ) contain the primitive translation    *
C   *    vectors and (QX,QY,QZ) the NQ basis vectors.                *
C   *                                                                *
C   *   *On exit (RPX,RPY,RPZ) = T+Q-Q', where T are the translation *
C   *    vectors, sorted into shells of increasing length and given  *
C   *    a single index IR which combines IQ with the shell index IS *
C   *    and the vector index in each shell. (RX,RY,RZ) contain the  *
C   *    corresponding R=T+Q and (TX,TY,TZ) the bare translations T. *
C   *    Vectors belonging to a given Q' have indices between        *
C   *    NRIQ(IQ-1)+1 and NRIQ(IQ).                                  *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE control_data
      USE lattice
      USE message
      IMPLICIT NONE
      CHARACTER(LEN=10) :: TGENV = ' LATT3D:**'
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: CSX, CSY, CSZ, DC
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: D
      REAL(KIND=8) :: QIX, QIY, QIZ, QMQPX, QMQPY, QMQPZ, SX, SY, SZ,
     1                DX, DA
      INTEGER, DIMENSION(:), ALLOCATABLE :: INDEX, JQT
      INTEGER :: NCR, IQ, JQ, NCRQ, L, M, N, IR, IS, I, INX, NSMAX, IN
  101 FORMAT(/,' LATT3D:   Generate lattice vectors, NGHBP =',I3)
  102 FORMAT(' ',9X,2I3,I4,I3,F10.6,2X,3F10.6)
  103 FORMAT(' ',10X,'(',F10.5,',',F10.5,',',F10.5,'  )',/)
  104 FORMAT(//,11X,'IS IN  IR JQ',5X,'D',10X,'RPX',7X,'RPY',7X,'RPZ')
  107 FORMAT(' ',/,11X,'Shells of neighbours around',//,10X,' IQ =',I3,
     1       22X,'QP =',3F7.3)
  108 FORMAT(A,' DMAX =',F10.5,' is too small',/,11X,'last DX =',F10.5)
  110 FORMAT(/,11X,'Number of vectors:',I4,'. Number of shells:',I3)
C
C     Establish the dimension NCR of temporary arrays
C
      NCR=0
      DO 20 IQ=1,NQ
      QIX=QX(IQ)
      QIY=QY(IQ)
      QIZ=QZ(IQ)
      NCRQ=0
      DO 21 JQ=1,NQ
      QMQPX=QX(JQ)-QIX
      QMQPY=QY(JQ)-QIY
      QMQPZ=QZ(JQ)-QIZ
      DO 21 L=-NGHBP,NGHBP
      DO 21 M=-NGHBP,NGHBP
      DO 21 N=-NGHBP,NGHBP
      SX=QMQPX+L*BSX(1)+M*BSX(2)+N*BSX(3)
      SY=QMQPY+L*BSY(1)+M*BSY(2)+N*BSY(3)
      SZ=QMQPZ+L*BSZ(1)+M*BSZ(2)+N*BSZ(3)
      DX=SQRT(SX*SX+SY*SY+SZ*SZ)
      IF(DX.LE.DMAX) THEN
         NCRQ=NCRQ+1
      ENDIF
   21 CONTINUE
      IF(NCRQ.LE.1) THEN
         WRITE(M6,108) TGENV,DMAX,DX
         STOP
      ENDIF
      NCR=MAX(NCR,NCRQ)
   20 CONTINUE
      ALLOCATE(CSX(NCR),CSY(NCR),CSZ(NCR),DC(NCR),INDEX(NCR),JQT(NCR))
C
C     Establish the dimension NR and NSMAX of lattice arrays
C
      IR=0
      NSMAX=0
      DO 23 IQ=1,NQ
      QIX=QX(IQ)
      QIY=QY(IQ)
      QIZ=QZ(IQ)
      NCR=0
      DO 22 JQ=1,NQ
      QMQPX=QX(JQ)-QIX
      QMQPY=QY(JQ)-QIY
      QMQPZ=QZ(JQ)-QIZ
      DO 22 L=-NGHBP,NGHBP
      DO 22 M=-NGHBP,NGHBP
      DO 22 N=-NGHBP,NGHBP
      SX=QMQPX+L*BSX(1)+M*BSX(2)+N*BSX(3)
      SY=QMQPY+L*BSY(1)+M*BSY(2)+N*BSY(3)
      SZ=QMQPZ+L*BSZ(1)+M*BSZ(2)+N*BSZ(3)
      DX=SQRT(SX*SX+SY*SY+SZ*SZ)
      IF(DX.LE.DMAX) THEN
         NCR=NCR+1
         DC(NCR)=DX
         CSX(NCR)=SX
         CSY(NCR)=SY
         CSZ(NCR)=SZ
         JQT(NCR)=JQ
      ENDIF
   22 CONTINUE
C
C     Sort the R-vectors in order of increasing length
C
      CALL QSORT(DC,INDEX,NCR)
C
C     Sort  R-vectors in shells of neighbours
C
      DA=-1.
      IS=0
      DO 24 I=1,NCR
      INX=INDEX(I)
      IF(DC(I)-DA.GT.1.E-5) THEN
         IF(DC(I).GT.DMAX) GO TO 24
         IS=IS+1
      ENDIF
      IR=IR+1
      DA=DC(I)
   24 CONTINUE
      IF(IS.GT.NSMAX) NSMAX=IS
   23 CONTINUE
      NR=IR
      ALLOCATE(RX(0:NR),RY(0:NR),RZ(0:NR),TX(0:NR),TY(0:NR),TZ(0:NR))
      ALLOCATE(IQBAS(NR),JQBAS(NR),NRIQ(0:NQ),NSH(NQ))
      ALLOCATE(D(NQ,NSMAX),RWATS(NQ))
C
      WRITE(M6,101) NGHBP
C
C     Generate RP = T+Q-QP
C
      NRIQ(0)=0
      IR=0
      DO 30 IQ=1,NQ
      QIX=QX(IQ)
      QIY=QY(IQ)
      QIZ=QZ(IQ)
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) THEN
         WRITE(M6,107) IQ,QIX,QIY,QIZ
         WRITE(M6,104)
      ENDIF
      NCR=0
      DO 31 JQ=1,NQ
      QMQPX=QX(JQ)-QIX
      QMQPY=QY(JQ)-QIY
      QMQPZ=QZ(JQ)-QIZ
      DO 31 L=-NGHBP,NGHBP
      DO 31 M=-NGHBP,NGHBP
      DO 31 N=-NGHBP,NGHBP
      SX=QMQPX+L*BSX(1)+M*BSX(2)+N*BSX(3)
      SY=QMQPY+L*BSY(1)+M*BSY(2)+N*BSY(3)
      SZ=QMQPZ+L*BSZ(1)+M*BSZ(2)+N*BSZ(3)
      DX=SQRT(SX*SX+SY*SY+SZ*SZ)
      IF(DX.LE.DMAX) THEN
         NCR=NCR+1
         DC(NCR)=DX
         CSX(NCR)=SX
         CSY(NCR)=SY
         CSZ(NCR)=SZ
         JQT(NCR)=JQ
      ENDIF
   31 CONTINUE
C
C     Sort the R-vectors in order of increasing length
C
      CALL QSORT(DC,INDEX,NCR)
      RWATS(IQ)=DC(NCR)
C
C     Sort  R-vectors in shells of neighbours
C
      DA=-1.
      IS=0
      DO 32 I=1,NCR
      INX=INDEX(I)
      IF(DC(I)-DA.GT.1.E-5) THEN
         IF(DC(I).GT.DMAX) GO TO 32
         IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,'(A)') ' '
         IS=IS+1
         IN=0
         D(IQ,IS)=DC(I)
      ENDIF
      IN=IN+1
      IR=IR+1
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,102) IS,IN,IR,JQT(INX),
     1   D(IQ,IS),CSX(INX),CSY(INX),CSZ(INX)
      RX(IR)=CSX(INX)+QIX
      RY(IR)=CSY(INX)+QIY
      RZ(IR)=CSZ(INX)+QIZ
      JQ=JQT(INX)
      TX(IR)=RX(IR)-QX(JQ)
      TY(IR)=RY(IR)-QY(JQ)
      TZ(IR)=RZ(IR)-QZ(JQ)
      IQBAS(IR)=IQ
      JQBAS(IR)=JQ
      DA=DC(I)
   32 CONTINUE
      NRIQ(IQ)=IR
      NSH(IQ)=IS
      WRITE(M6,110) NRIQ(IQ),NSH(IQ)
   30 CONTINUE
C
      NR=IR
      DEALLOCATE(CSX,CSY,CSZ,DC,INDEX)
      RETURN
      END
