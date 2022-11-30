      SUBROUTINE rotm3d(ibz,lmax)
C  ********************************************************************
C  *                                                                  *
C  *     Generate the 3D point group rotation matrices for the        *
C  *     cubic (48 symmetry elements), tetragonal (16),               *
C  *     orthorhombic (8), monoclinic (4) or triclinic (2)            *
C  *     crystal systems.                                             *
C  *                                                                  *
C  *     Generate the 3D point group rotation matrices for the        *
C  *     hexagonal close packed (24) crystal system.                  *
C  *                                                                  *
C  *     Generate the 3D point group rotation matrices for the        *
C  *     trigonal (12) crystal system.                                *
C  *                                                                  *
C  *     Author:   Magnus Alden  921019                               *
C  *   1 Modif.:   Igor Abrikosov 931108                              *
C  *   2 Modif.:   Andrei Ruban 940503                                *
C  *   3 Modif.:   Igor Abrikosov 950503                              *
C  *   for general lmax : Levente Vitos 970411                        *
C  *   trigonal system  : Levente Vitos 000228                        *
C  *                                                                  *
C  ********************************************************************
      USE symmetry
      USE message
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: nugam
      REAL(KIND=8), PARAMETER :: err = 1.d-8
      REAL(KIND=8) :: pi, smm, dd, db, a, b, g, diffg
      INTEGER, PARAMETER :: dotest = 0
      INTEGER :: ibz, lmax, irot, ir1, ir2, one, onep
      INTEGER :: irinv, l, lm, lmp, my, myp, mypp, fnd
      pi=dacos(-1.d0)
C
      ALLOCATE(ugam(48,0:lmax,-lmax:lmax,-lmax:lmax),ibzrot(48))
C
C     Initialize rotation matrices
C
      ugam=0.d0
C
      WRITE(m6,100) ibz
C
C     Set number of proper symmetry operations for particular 
C     IBZ and index IBZROT which is 1 if IROT is a proper
C     symmetry operation for this IBZ and 0 otherwise
C
      nrot=48
      ibzrot=0
      IF(ibz.LE.3) THEN
C
C        Cubic system
C
         nprprt=48
         ibzrot=1
      ELSEIF(ibz.EQ.4) THEN
C
C        Hexagonal system
C
         nprprt=24
         ibzrot(1:24)=1
         GO TO 30
      ELSEIF(ibz.EQ.5.OR.ibz.EQ.6) THEN
C
C        Tetragonal system
C
         nprprt=16
         ibzrot(1)=1
         ibzrot(2)=1
         ibzrot(7)=1
         ibzrot(8)=1
         ibzrot(9)=1
         ibzrot(10)=1
         ibzrot(15)=1
         ibzrot(20)=1
         ibzrot(25)=1
         ibzrot(26)=1
         ibzrot(31)=1
         ibzrot(32)=1
         ibzrot(33)=1
         ibzrot(34)=1
         ibzrot(39)=1
         ibzrot(44)=1
      ELSEIF(ibz.EQ.7) THEN
C
C        Trigonal system
C
         nprprt=12
         ibzrot(1:12)=1
         GO TO 40
      ELSEIF(ibz.GE.8.AND.ibz.LE.11) THEN
C
C        Orthorhombic system
C
         nprprt=8
         ibzrot(1)=1
         ibzrot(9)=1
         ibzrot(10)=1
         ibzrot(15)=1
         ibzrot(25)=1
         ibzrot(33)=1
         ibzrot(34)=1
         ibzrot(39)=1
      ELSEIF(ibz.EQ.12.OR.ibz.EQ.13) THEN
C
C        Monoclinic system
C
         nprprt=4
         ibzrot(1)=1
         ibzrot(9)=1
         ibzrot(25)=1
         ibzrot(33)=1
      ELSEIF(ibz.EQ.14) THEN
C
C        Triclinic system
C
         nprprt=2
         ibzrot(1)=1
         ibzrot(25)=1
      ENDIF
C
C     Generate cubic (48) symmetry elements
C
      DO 20 l=0,lmax
      DO 20 my=-l,l
      one=(-1.d0)**(l+my)
      onep=((-1.d0)**my)*ISIGN(1,my)
      DO 20 myp=-l,l
C
C     1. Identity element
C
         db=0.d0
         IF(my.EQ.myp) db=1.d0
         ugam(1,l,my,myp)=db
C
C     2. Mirror in x=y
C
         a=pi/2.d0
         b=0.d0
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(2,l,my,myp)=db*onep
C
C     3. Rotation +120 around (1,1,1)   
C
         a=pi/2.d0
         b=pi/2.d0
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(3,l,my,myp)=db
C
C     4. Rotation -120 around (1,1,1)   
C
         a=pi
         b=pi/2.d0
         g=pi/2.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(4,l,my,myp)=db
C
C     5. Mirror in y=z   
C
         a=-pi/2.d0
         b=pi/2.d0
         g=pi/2.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(5,l,my,myp)=db*one
C
C     6. Mirror in x=z   
C
         a=0.d0
         b=pi/2.d0
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(6,l,my,myp)=db*one 
C
C     7. Rotation +90 around (0,0,1)  
C
         a=0.d0
         b=0.d0
         g=pi/2.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(7,l,my,myp)=db
C
C     8. Rotation -90 around (0,0,1)  
C
         a=0.d0
         b=0.d0
         g=-pi/2.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(8,l,my,myp)=db
C
C     9. Rotation +-180 around (0,0,1)   
C
         a=0.d0
         b=0.d0
         g=pi
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(9,l,my,myp)=db
C
   20 CONTINUE
C
C     10-24.  2-6 above + 7-9 above :   
C   
      irot=9
      DO 21 ir1=7,9
      DO 21 ir2=2,6
      irot=irot+1
      DO 22 l=0,lmax
      DO 22 my=-l,l
      DO 22 myp=-l,l
      smm=0.d0
      DO 23 mypp=-l,l
   23 smm=smm+ugam(ir1,l,my,mypp)*ugam(ir2,l,mypp,myp)
   22 ugam(irot,l,my,myp)=smm
   21 CONTINUE
C
C     25. Inversion through origin.   
C
      DO 24 l=0,lmax
      one=(-1.d0)**l
      DO 24 my=-l,l
   24 ugam(25,l,my,my)=one
C
C     26-48.  2-24 above + 25 (inversion through origin) 
C
      DO 25 irinv=26,48
      irot=irinv-24
      DO 26 l=0,lmax
      DO 26 my=-l,l
      DO 26 myp=-l,l
      smm=0.d0
      DO 27 mypp=-l,l
   27 smm=smm+ugam(25,l,my,mypp)*ugam(irot,l,mypp,myp)
   26 ugam(irinv,l,my,myp)=smm
   25 CONTINUE
      GO TO 50
C
C     Hexagonal close packed crystal
C
   30 CONTINUE
      DO 31 l=0,lmax
      DO 31 my=-l,l
      one=(-1.d0)**(l+my)
      onep=((-1.d0)**my)*ISIGN(1,my)
      DO 31 myp=-l,l
C
C     1. Identity element
C
         db=0.d0
         IF(my.EQ.myp) db=1.d0
         ugam(1,l,my,myp)=db
C
C     2. Rotation +120 around [0,0,0,1]_h
C
         a=2.d0*pi/3.d0
         b=0.d0
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(2,l,my,myp)=db
C
C     3. Rotation +240 around [0,0,0,1]_h
C
         a=4.d0*pi/3.d0
         b=0.d0
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(3,l,my,myp)=db
C
C     4. Rotation +180 around [-1,1,0,0]_h
C
         a=0.d0
         b=pi
         g=2.d0*pi/3.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(4,l,my,myp)=db
C
C     5. Rotation +180 around [2,1,0,0]_h  
C
         a=0.d0
         b=pi
         g=-2.d0*pi/3.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(5,l,my,myp)=db
C
C     6. Rotation +180 around [1,2,0,0]_h  
C
         a=0.d0
         b=pi
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(6,l,my,myp)=db
C
C     7. Mirror in [0,0,0,1]_h
C
         db=0.d0
         IF(my.EQ.myp) db=one
         ugam(7,l,my,myp)=db
   31 CONTINUE
C
C     8-12.  2-6 above + 7 (mirror in [0,0,0,1]_h)
C
      DO 32 irinv=8,12
      irot=irinv-6
      DO 33 l=0,lmax
      DO 33 my=-l,l
      DO 33 myp=-l,l
      smm=0.d0
      DO 34 mypp=-l,l
   34 smm=smm+ugam(7,l,my,mypp)*ugam(irot,l,mypp,myp)
   33 ugam(irinv,l,my,myp)=smm
   32 CONTINUE
C
C     13. Inversion through origin.   
C
      DO 35 l=0,lmax
      one=(-1.d0)**l
      DO 35 my=-l,l
   35 ugam(13,l,my,my)=one
C
C     14-24.  2-12 above + 13 (inversion through the origin)
C
      DO 36 irinv=14,24
      irot=irinv-12
      DO 37 l=0,lmax
      DO 37 my=-l,l
      DO 37 myp=-l,l
      smm=0.d0
      DO 38 mypp=-l,l
   38 smm=smm+ugam(13,l,my,mypp)*ugam(irot,l,mypp,myp)
   37 ugam(irinv,l,my,myp)=smm
   36 CONTINUE
      Go TO 50
C
C     Trigonal crystal
C
   40 CONTINUE
      DO 41 l=0,lmax
      DO 41 my=-l,l
      one=(-1.d0)**(l+my)
      onep=((-1.d0)**my)*ISIGN(1,my)
      DO 41 myp=-l,l
C
C     1. Identity element
C
         db=0.d0
         IF(my.EQ.myp) db=1.d0
         ugam(1,l,my,myp)=db
C
C     2. Rotation +120 around [0,0,0,1]_h
C
         a=2.d0*pi/3.d0
         b=0.d0
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(2,l,my,myp)=db
C
C     3. Rotation +240 around [0,0,0,1]_h  
C
         a=4.d0*pi/3.d0
         b=0.d0
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(3,l,my,myp)=db
C
C     4. Rotation +180 around [1,0,0,0]_h
C
         a=pi
         b=pi
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(4,l,my,myp)=db
C
C     5. Rotation +180 around [0,1,0,0]_h
C
         a=-pi/3.d0
         b=pi
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(5,l,my,myp)=db
C
C     6. Rotation +180 around [0,0,1,0]_h
C
         a=pi/3.d0
         b=pi
         g=0.d0
         CALL dbar(a,b,g,l,myp,my,db)
         ugam(6,l,my,myp)=db
C
   41 CONTINUE
C
C     7. Inversion through origin.   
C
      DO 45 l=0,lmax
      one=(-1.d0)**l
      DO 45 my=-l,l
   45 ugam(7,l,my,my)=one
C
C     8-12.  2-6 above + 7 (inversion through origin) 
C
      DO 46 irinv=8,12
      irot=irinv-6
      DO 47 l=0,lmax
      DO 47 my=-l,l
      DO 47 myp=-l,l
      smm=0.d0
      DO 48 mypp=-l,l
   48 smm=smm+ugam(7,l,my,mypp)*ugam(irot,l,mypp,myp)
   47 ugam(irinv,l,my,myp)=smm
   46 CONTINUE
C
C     Print rotation matrices
C
   50 WRITE(m6,110) nprprt
      IF(NPRN.EQ.10) THEN
         DO 51 irot=1,nrot
         IF(ibzrot(irot).EQ.0) GO TO 51
         WRITE(m6,*)
         DO 52 l=0,lmax
         WRITE(m6,*)
         DO 52 my=-l,l
         lm=l*l+l+my+1
         DO 52 myp=-l,l
         lmp=l*l+l+myp+1
         dd=ugam(irot,l,my,myp)
         IF(ABS(dd).GT.1.d-10) THEN
            WRITE(m6,120) irot,lm,lmp,dd
         ENDIF
   52    CONTINUE
   51    CONTINUE
      ENDIF
C
C     Test for the point groups
C
      IF(dotest.EQ.1) THEN
         ALLOCATE(nugam(0:lmax,-lmax:lmax,-lmax:lmax))
C
         DO 60 ir1=1,nrot
         IF(ibzrot(ir1).EQ.0) GO TO 60
         DO 61 ir2=1,nrot
         IF(ibzrot(ir2).EQ.0) GO TO 61
C
         DO 62 l=0,lmax
         DO 62 my=-l,l
         DO 62 myp=-l,l
         smm=0.d0
         DO 63 mypp=-l,l
   63    smm=smm+ugam(ir1,l,my,mypp)*ugam(ir2,l,mypp,myp) 
   62    nugam(l,my,myp)=smm
C
         irot=0
   64    irot=irot+1
         IF(ibzrot(irot).EQ.0) GO TO 65
         fnd=1
         DO 66 l=0,lmax
         DO 66 my=-l,l
         DO 66 myp=-l,l
         diffg=ABS(ugam(irot,l,my,myp)-nugam(l,my,myp))
         IF(diffg.GT.err) fnd=0
   66    CONTINUE
         IF(fnd.EQ.1) GO TO 61
   65    CONTINUE
         IF(irot.LT.nrot) GO TO 64
         WRITE(m6,130) ir1,ir2
         IF(msgl.NE.0) WRITE(msgio,130) ir1,ir2
         STOP
   61    CONTINUE
   60    CONTINUE
C
         DEALLOCATE(nugam)
      ENDIF
C
  100 FORMAT(/,' ROTM3D:   Set up the 3D rotation matrices for IBZ =',
     .       i3)
  110 FORMAT(/,11x,'Rotation matrices generated, NPRPRT =',i3)
  120 FORMAT(11x,'UGAM(',i3,',',i3,',',i3,')=',f10.6)
  130 FORMAT(/,' ROTM3D:   Error !!! UGAM(',i3,')*UGAM(',i3,') is',
     .         ' not element of the group')
      RETURN
      END
