      SUBROUTINE PRIMV
C   ******************************************************************
C   *                                                                *
C   *    Obtain primitive translational vectors of real space and    *
C   *    read atomic positions (basis).                              *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    LAT   : Lattice according to the list in MESH.              *
C   *                                                                *
C   *   *Input:                                                      *
C   *                                                                *
C   *    A     : Lattice parameter a.                                *
C   *    ALPHA : Crystallographic angle between a and b.             *
C   *    B     : Lattice parameter b.                                *
C   *    BETA  : Crystallographic angle between b and c.             *
C   *    C     : Lattice parameter c.                                *
C   *    GAMMA : Crystallographic angle between c and a.             *
C   *    IPRIM := 0: Read non-standard vectors.                      *
C   *           = 1: Generate standard vectors.                      *
C   *           = 2: Generate supercell.                             *
C   *    NQR2  . Repeat for 2D supercell.                            *
C   *    NQ3   : Number of atoms in the 3D cell.                     *
C   *    QX3                                                         *
C   *    QY3   : Position of atoms in the 3D cell.                   *
C   *    QZ3                                                         *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    BOA   : b/a.                                                *
C   *    BSX(I)                                                      *
C   *    BSY(I): Three primitive translation vectors.                *
C   *    BSZ(I)                                                      *
C   *    COA   : c/a.                                                *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE control_data
      USE csts
      USE message
      IMPLICIT NONE
      CHARACTER(LEN=4), DIMENSION(14) :: TLAT
      REAL(KIND=8) :: A, B, C, RADF, ALPHA, BETA, GAMMA
      INTEGER :: LAT, IPRIM, I, NQR2, K, J, IQ
      DATA TLAT/'  Sc',' Fcc',' Bcc',' Hex','  St',' Bct','Trig',
     1          '  So','Baco',' Bco',' Fco','  Sm',' Bcm','Tric'/
  101 FORMAT(/,A,'Standard choice not implemented for LAT =',I3)
  102 FORMAT(3(10X,F10.6))
  103 FORMAT(/,11X,'A     =',F10.6,' B     =',F10.6,' C     =',F10.6)
  104 FORMAT(/,'IPRIM =',I3,'. Must be 0, 1, or 2.')
  105 FORMAT(11X,'(',F10.5,',',F10.5,',',F10.5,' )')
  106 FORMAT(/,11X,'Primitive vectors for',A,' lattice in',/,11X,
     1       'units of the lattice spacing a:',/)
  107 FORMAT(/,11X,'A     =',F10.6,' B     =',F10.6,' C     =',F10.6,
     1       /,11X,'ALPHA =',F10.6,' BETA  =',F10.6,' GAMMA =',F10.6)
  108 FORMAT(/,11X,'Basis vectors:',12X,'NQ3 =',I4,/)
  110 FORMAT(/,' PRIMV:**  Number of supercell basis vectors must',
     1       ' be equal to ',A,' =',I3)
  111 FORMAT(/,' PRIMV:**  Supercell calculation: ',A,I3,/,11X,
     1       'Must be smaller than or equal to',A,I3)
C
      IF(MSGL.EQ.1) WRITE(MSGIO,'(A,I5)')
     1          ' Construct lattice'
      READ(5,'(7x,i3,4(8X,I2))') NQ3,LAT,IPRIM,NQR2
C
      ALLOCATE(QX3(NQ3),QY3(NQ3),QZ3(NQ3))
C
      READ(5,102) A,B,C
      BOA=B/A
      COA=C/A
      IF(IPRIM.EQ.0) THEN
         WRITE(M6,'(/,A,4X,2A)') ' PRIMV:','Special choice of',
     1                          ' primitive vectors.'
         WRITE(M6,103) A,B,C
C
C        Read the primitive vectors on FOR005
C
         DO 20 I=1,3
   20    READ(5,102) BSX(I),BSY(I),BSZ(I)
         IF(LAT.EQ.7.AND.BSZ(1).LT.1.D-07) BOA=-1.
      ELSEIF(IPRIM.EQ.1) THEN
         WRITE(M6,'(/,A,4X,2A)') ' PRIMV:','Default choice of',
     1                          ' primitive vectors.'
C
C        Construct primitive vectors from A,B,C,ALPHA,BETA,GAMMA
C
         RADF=PI/180.D0
         READ(5,102) ALPHA,BETA,GAMMA
         WRITE(M6,107) A,B,C,ALPHA,BETA,GAMMA
         ALF=ALPHA*RADF
         BET=BETA*RADF
         GAM=GAMMA*RADF
         GO TO (1,2,3,4,5,6,7,8,9,10,11,12,13,14) LAT
C
C        Simple cubic
C
    1    BSX(1)=1.D0
         BSY(1)=0.D0
         BSZ(1)=0.D0
         BSX(2)=0.D0
         BSY(2)=1.D0
         BSZ(2)=0.D0
         BSX(3)=0.D0
         BSY(3)=0.D0
         BSZ(3)=1.D0
         GO TO 21
C
C        Face centreed cubic
C
    2    BSX(1)=0.5D0
         BSY(1)=0.5D0
         BSZ(1)=0.0D0
         BSX(2)=0.0D0
         BSY(2)=0.5D0
         BSZ(2)=0.5D0
         BSX(3)=0.5D0
         BSY(3)=0.0D0
         BSZ(3)=0.5D0
         GO TO 21
C
C        Body centred cubic
C
    3    BSX(1)=0.5D0
         BSY(1)=0.5D0
         BSZ(1)=-0.5D0
         BSX(2)=-0.5D0
         BSY(2)=0.5D0
         BSZ(2)=0.5D0
         BSX(3)=0.5D0
         BSY(3)=-0.5D0
         BSZ(3)=0.5D0
         GO TO 21
C
C        Hexagonal
C
    4    BSX(1)=1.D0
         BSY(1)=0.D0
         BSZ(1)=0.D0
         BSX(2)=-0.5D0
         BSY(2)=SQRT(3.D0)/2.D0
         BSZ(2)=0.D0
         BSX(3)=0.D0
         BSY(3)=0.D0
         BSZ(3)=COA
         GO TO 21
C
C        Simple tetragonal
C
    5    BSX(1)=1.D0
         BSY(1)=0.D0
         BSZ(1)=0.D0
         BSX(2)=0.D0
         BSY(2)=1.D0
         BSZ(2)=0.D0
         BSX(3)=0.D0
         BSY(3)=0.D0
         BSZ(3)=COA
         GO TO 21
C
C        Body centred tetragonal
C
    6    BSX(1)=0.5D0
         BSY(1)=-0.5D0
         BSZ(1)=COA/2.D0
         BSX(2)=0.5D0
         BSY(2)=0.5D0
         BSZ(2)=-COA/2.D0
         BSX(3)=-0.5D0
         BSY(3)=0.5D0
         BSZ(3)=COA/2.D0
         GO TO 21
C
C        Trigonal
C
    7    BSX(1)=0.D0
         BSY(1)=1.D0
         BSZ(1)=COA
         BSX(2)=-SQRT(3.D0)/2.D0
         BSY(2)=-0.5D0
         BSZ(2)=COA
         BSX(3)=SQRT(3.D0)/2.D0
         BSY(3)=-0.5D0
         BSZ(3)=COA
         GO TO 21
C
C        Simple orthorombic
C
    8    BSX(1)=1.D0
         BSY(1)=0.D0
         BSZ(1)=0.D0
         BSX(2)=0.D0
         BSY(2)=BOA
         BSZ(2)=0.D0
         BSX(3)=0.D0
         BSY(3)=0.D0
         BSZ(3)=COA
         GO TO 21
C
C        Basis centered orthorombic
C
    9    BSX(1)=1.D0/2.D0
         BSY(1)=-BOA/2.D0
         BSZ(1)=0.D0
         BSX(2)=1.D0/2.D0
         BSY(2)=BOA/2.D0
         BSZ(2)=0.D0
         BSX(3)=0.D0
         BSY(3)=0.D0
         BSZ(3)=COA
         GO TO 21
C
C        Body centred orthorombic
C
   10    BSX(1)=1.D0/2.D0
         BSY(1)=-BOA/2.D0
         BSZ(1)=COA/2.D0
         BSX(2)=1.D0/2.D0
         BSY(2)=BOA/2.D0
         BSZ(2)=-COA/2.D0
         BSX(3)=-1.D0/2.D0
         BSY(3)=BOA/2.D0
         BSZ(3)=COA/2.D0
         GO TO 21

C
C        Face centred orthorombic
C
   11    BSX(1)=1.D0/2.D0
         BSY(1)=0.D0
         BSZ(1)=COA/2.D0
         BSX(2)=1.D0/2.D0
         BSY(2)=BOA/2.D0
         BSZ(2)=0.D0
         BSX(3)=0.D0
         BSY(3)=BOA/2.D0
         BSZ(3)=COA/2.D0
         GO TO 21
C
C        Simple monoclinic
C
   12    BSX(1)=1.D0
         BSY(1)=0.D0
         BSZ(1)=0.D0
         BSX(2)=BOA*COS(gam)
         BSY(2)=BOA*SIN(gam)
         BSZ(2)=0.D0
         BSX(3)=0.D0
         BSY(3)=0.D0
         BSZ(3)=COA
         GO TO 21
C
C        Base centred monoclinic
C
   13    BSX(1)=0.D0
         BSY(1)=-BOA
         BSZ(1)=0.D0
         BSX(2)=0.5D0*SIN(GAM)
         BSY(2)=-0.5D0*COS(GAM)
         BSZ(2)=-0.5D0*COA
         BSX(3)=0.5D0*SIN(GAM)
         BSY(3)=-0.5D0*COS(GAM)
         BSZ(3)=0.5D0*COA
         GO TO 21
C
C        Simple triclinic
C
   14    BSX(1)=1.d0
         BSY(1)=0.d0
         BSZ(1)=0.d0
         BSX(2)=boa*COS(gam)
         BSY(2)=boa*SIN(gam)
         BSZ(2)=0.d0
         BSX(3)=coa*COS(bet)
         BSY(3)=coa*(COS(alf)-COS(bet)*COS(gam))/SIN(gam)
         BSZ(3)=coa*SQRT((1.d0-COS(gam)*COS(gam)-COS(alf)*COS(alf)-
     .   COS(bet)*COS(bet)+2.d0*COS(alf)*COS(bet)*COS(gam)))/SIN(gam)
         GO TO 21
      ELSEIF(IPRIM.EQ.2) THEN
         WRITE(M6,'(/,A,4X,2A)') ' PRIMV:','Supercell',
     1                          ' primitive vectors.'
         WRITE(M6,103) A,B,C
C
C        Read the primitive vectors on FOR005
C
         DO 22 I=1,3
   22    READ(5,102) BSX(I),BSY(I),BSZ(I)
         K=0
         DO 23 J=1,NQR2
         DO 23 I=1,NQR2
         K=K+1
         IF(K.GT.NQ3) THEN
            WRITE(M6,111) 'K =',K,' NQ3 =',NQ3
            STOP
         ENDIF
         QX3(K)=(I-1)*BSX(1)+(J-1)*BSX(2)
         QY3(K)=(I-1)*BSY(1)+(J-1)*BSY(2)
   23    QZ3(K)=(I-1)*BSZ(1)+(J-1)*BSZ(2)
         IF(K.NE.NQ3) THEN
            WRITE(M6,110) 'NQ3',NQ3
            STOP
         ENDIF
         DO 24 I=1,2
         BSX(I)=NQR2*BSX(I)
   24    BSY(I)=NQR2*BSY(I)
      ELSE
         WRITE(M6,104) 'PRIMV: **',IPRIM
         STOP
      ENDIF
   21 IF(IPRIM.NE.2) THEN
C
C        Read the NQ3 basis vectors in the primitive cell
C
         DO 25 IQ=1,NQ3
   25    READ(5,102) QX3(IQ),QY3(IQ),QZ3(IQ)
      ENDIF
C
C     Print primitive and basis vectors
C
      IF(LAT.GE.1.AND.LAT.LE.14) THEN
         WRITE(M6,106) TLAT(LAT)
      ELSE
         WRITE(M6,106) ' unspecified'
      ENDIF
      DO 26 I=1,3
   26 WRITE(M6,105) BSX(I),BSY(I),BSZ(I)
      WRITE(M6,108) NQ3
      DO 27 I=1,NQ3
   27 WRITE(M6,105) QX3(I),QY3(I),QZ3(I)
      RETURN
C
   99 WRITE(M6,101) 'PRIMV: **',LAT
      STOP
      END
