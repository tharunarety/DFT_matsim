      SUBROUTINE kmesh(lat,npx,npy,npz,npt,dkx,dky,dkz,dhx,
     .           boa,coa,alf,bet,gam,m6,mkm,tkx,tky,tkz,
     .           kx,ky,kz,ww)
C   ******************************************************************
C   *                                                                *
C   *    Construction of a mesh in k-space to be used with state     *
C   *    density programs etc. K-vectors are in units of pi/a        *
C   *    where a is the lattice parameter.                           *
C   *                                                                *
C   ******************************************************************
C   *                                                                *
C   *    LAT Implemented Brillouin zone                              *
C   *                                                                *
C   *     1      yes     Simple cubic                                *
C   *     2      yes     Face centred cubic                          *
C   *     3      yes     Body centred cubic                          *
C   *     4      yes     Hexagonal close packed(24th part)           *
C   *     5      yes     Simple tetragonal                           *
C   *     6      yes     Body centred tetragonal                     *
C   *     7      yes     Trigonal                                    *
C   *     8      yes     Simple orthorombic(1994.06.22)              *
C   *     9      yes     Base centred orthorombic                    *
C   *    10      yes     Body centred orthorombic(1995.10.25)        *
C   *    11      yes     Face centred orthorombic                    *
C   *    12      yes     Simple monoclinic(1994.02.18)               *
C   *    13      yes     Base centred monoclinic (Check KTRNSFM)     *
C   *    14      yes     Simple triclinic(2013.10.24)                *
C   *                                                                *
C   *    HLS: 02-Nov-96                                              *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER TBZ(14)*4
      INTEGER X,Y,Z
      DIMENSION TKX(3),TKY(3),TKZ(3)
      DIMENSION KX(MKM),KY(MKM),KZ(MKM),WW(MKM)
C
      DATA TBZ/'  sc',' fcc',' bcc',' hex','  st',' bct','trig',
     1         '  so','baco',' bco',' fco','  sm',' bcm','tric'/
  101 FORMAT(/,' KMESH:**  Number of points on k-mesh exeeds MKM =',
     1       I5,/,11X,'Decrease the number of points or increase MKM')
  102 FORMAT(/,' KMESH:**  Wrong NPY, NPY =',I4)
  103 FORMAT(/,' KMESH:**  Wrong NPZ, NPZ =',I4)
  104 FORMAT(/,' KMESH:**  Wrong NPX, NPX =',I4)
  105 FORMAT(/,' KMESH:**  LAT =',I3,' not implemented')
  106 FORMAT(/,' KMESH:',4X,'Brillouin zone:',A,'. NPT =',I5)
  107 FORMAT(/,' KMESH:**  None of the parameters can be zero',
     .       ' NPX,NPY,NPZ =',3i4)
C
      NPXM=NPX-1
      NPYM=NPY-1
      NPZM=NPZ-1
      NP=0
      IF(LAT.EQ.1) THEN
C
C        Simple cubic:
C
C        0 .LE. Kz .LE. Kx .LE. Ky .LE. Pi/a
C
         DKX=1.D0/NPYM
         DKY=DKX
         DKZ=DKX
         DHX=0.D0
         DO 30 J=1,NPY
         Y=J-1
         DO 30 I=1,J
         X=I-1
         DO 30 K=1,I
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=48.D0
         IF(X.EQ.Y.OR.X.EQ.Z.OR.Y.EQ.Z) W=W/2.D0
         IF(X.EQ.Y.AND.Y.EQ.Z) W=W/3.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.0) W=W/2.D0
         IF(X.EQ.NPYM) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(X+Y+Z.EQ.3*NPYM) W=1.D0
         IF(X+Y+Z.EQ.0) W=1.D0
   30    WW(NP)=W
      ELSEIF(LAT.EQ.2) THEN
C
C        Face centred cubic:
C
C        0 .LE. Kz .LE. Kx .LE. Ky .LE. 2Pi/a
C        Kx + Ky + Kz .LE. 3Pi/a
C
         IF(NPYM.NE.4*(NPYM/4)) GO TO 998
         DKX=2.D0/NPYM
         DKY=DKX
         DKZ=DKX
         DHX=0.D0
         NPPX=NPY
         NPH=NPY/2+1
         NPTH=(NPYM/2)*3
         DO 31 I=1,NPPX
         Y=I-1
         M1=NPPX-I+NPH
         NQY=MIN0(I,M1)
         DO 31 J=1,NQY
         X=J-1
         M5=M1-J+1
         NPPZ=MIN0(J,M5)
         DO 31 K=1,NPPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=48.D0
         IF(X.EQ.Y.OR.X.EQ.Z.OR.Y.EQ.Z) W=W/2.D0
         IF(X.EQ.Y.AND.Y.EQ.Z) W=W/3.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(X+Y+Z.EQ.NPTH) W=W/2.D0
         IF(X+Y+Z.EQ.0) W=1.D0
   31    WW(NP)=W
      ELSEIF(LAT.EQ.3) THEN
C
C        Body centred cubic:
C
C        0 .LE. Kz .LE. Ky .LE. Kx .LE. 2Pi/a
C        Kx + Ky .LE. 2Pi/a
C
         IF(NPYM.NE.2*(NPYM/2)) GO TO 998
         DKX=2.D0/NPYM
         DKY=DKX
         DKZ=DKX
         DHX=0.D0
         NPYH=NPYM/2
         DO 32 I=1,NPY
         X=I-1
         JM=MIN0(I,NPY-I+1)
         DO 32 J=1,JM
         Y=J-1
         DO 32 K=1,J
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=48.D0
         IF(X.EQ.Y.OR.X.EQ.Z.OR.Y.EQ.Z) W=W/2.D0
         IF(X.EQ.Y.AND.Y.EQ.Z) W=W/3.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(X+Y.EQ.NPYM) W=W/2.D0
         IF(Y.EQ.Z.AND.X+Y.EQ.NPYM) W=8.D0
         IF(Z.EQ.NPYH) W=2.D0
         IF(X.EQ.NPYM) W=1.D0
         IF(X+Y+Z.EQ.0) W=1.D0
   32    WW(NP)=W
      ELSEIF(LAT.EQ.4) THEN
C
C        Hexagonal close-packed:
C
C        0 .LE. 2Ky .LE. Kx .LE. 4/3 Pi/a
C        0 .LE. Kz .LE. a/c Pi/a
C
         IF(NPYM.NE.2*(NPYM/2)) GO TO 998
         NPH=NPY/2+1
         NPZM=NPZ-1
         AOC=1.D0/COA
         DKX=4.D0/3.D0/NPYM
         DKY=DSQRT(3.D0)/2.D0*DKX
         IF(NPZ.GT.1) DKZ=AOC/NPZM
         DHX=-0.5D0*DKX
         DO 33 K=1,NPZ
         Z=K-1
         DO 33 J=1,NPH
         Y=J-1
         IM=2*Y+1
         DO 33 I=IM,NPY
         X=I-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=24.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(X.EQ.0) W=W/3.D0
         IF(X.EQ.2*Y) W=W/2.D0
         IF(Z.EQ.NPZM.AND.Z.GT.0) W=W/2.D0
         IF(X.EQ.NPYM.AND.Y.EQ.0) W=W/3.D0
         IF(X.EQ.NPYM.AND.Y.NE.0) W=W/2.D0
         IF(X+Y+Z.EQ.0) W=1.D0
         IF(X+Y.EQ.0.AND.Z.EQ.NPZM) W=1.D0
         WW(NP)=W
   33    CONTINUE
      ELSEIF(LAT.EQ.5) THEN
C
C        Simple tetragonal:
C
C        0 .LE. Ky .LE. Kx .LE. Pi/a
C        0 .LE. Kz .LE. a/c Pi/a
C
         AOC=1.D0/COA
         DKX=1.D0/NPYM
         DKY=DKX
         IF(NPZ.GT.1) THEN
            DKZ=AOC/NPZM
         ELSE
            DKZ=AOC
         ENDIF
         DHX=0.D0
         DO 34 I=1,NPY
         X=I-1
         DO 34 J=1,I
         Y=J-1
         DO 34 K=1,NPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=16.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.Y) W=W/2.D0
         IF(X.EQ.NPYM) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(Z.EQ.NPZM.AND.Z.GT.0) W=W/2.D0
   34    WW(NP)=W
      ELSEIF(LAT.EQ.6) THEN
C
C        Body centred tetragonal:
C
C        0 .LE. Ky .LE. Kx .LE. Pi/a
C        0 .LE. Kz .LE. 2*a/c Pi/a
C
         IF(NPZM.NE.2*(NPZM/2)) GO TO 997
         AOC=1.D0/COA
         DKX=1.D0/NPYM
         DKY=DKX
         IF(NPZ.GT.1) THEN
            DKZ=2.D0*AOC/NPZM
         ELSE
            DKZ=AOC
         ENDIF
         DHX=0.D0
         DO 35 I=1,NPY
         X=I-1
         DO 35 J=1,I
         Y=J-1
         DO 35 K=1,NPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=16.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.Y) W=W/2.D0
         IF(X.EQ.NPYM) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(Z.EQ.NPZM) W=W/2.D0
   35    WW(NP)=W
      ELSEIF(LAT.EQ.7) THEN
C
C        Trigonal:
C
C        1) Trigonal primitive translations
C        0 .LE. 2 Ky .LE. Kx .LE. 4/(3 SQRT(3)) Pi/a
C        -a/c Pi/a .LE. Kz .LE. a/c Pi/a
C
C        2) Hexagonal primitive translations
C        0 .LE. 2 Ky .LE. Kx .LE. 4/3 Pi/a
C       -a/c Pi/a .LE. Kz .LE. a/c Pi/a
C
         IF(NPYM.NE.2*(NPYM/2)) GO TO 998
         IF(NPZM.NE.2*(NPZM/2)) GO TO 997
         IF(BOA.GT.0.D0.AND.NPZM.NE.6*(NPZM/6)) GO TO 997
         AOC=1.D0/COA
         NPH=NPY/2+1
         NPZH=NPZM/2
         NPZM=NPZ-1
         DKX=4.D0/3.D0/NPYM
         IF(BOA.GT.0.D0) DKX=DKX/DSQRT(3.D0)
         DKY=DSQRT(3.D0)/2.D0*DKX
         DKZ=2.D0*AOC/NPZM
         DHX=-0.5D0*DKX
         DO 37 K=-NPZH,NPZH
         Z=K
         DO 37 J=1,NPH
         Y=J-1
         IM=2*Y+1
         DO 37 I=IM,NPY
         X=I-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=12.D0
         IF(X.EQ.0) W=W/3.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(Z.EQ.-NPZH) W=W/2.D0
         IF(X.EQ.2*Y) W=W/2.D0
         IF(Z.EQ.NPZH) W=W/2.D0
         IF(X.EQ.NPYM.AND.Y.EQ.0) W=W/3.D0
         IF(X.EQ.NPYM.AND.Y.NE.0) W=W/2.D0
   37    WW(NP)=W
      ELSEIF(LAT.EQ.8) THEN
C
C        Simple orthorombic
C
C        0 .LE. Ky .LE. Pi/a
C        0 .LE. Ky .LE. a/b Pi/a
C        0 .LE. Kz .LE. a/c Pi/a
C
         AOC=1.D0/COA
         AOB=1.D0/BOA
         DKX=1.D0/NPXM
         DKY=AOB/NPYM
         IF(NPZ.GT.1) THEN
            DKZ=AOC/NPZM
         ELSE
            DKZ=AOC
         ENDIF
         DHX=0.D0
         DO 46 I=1,NPX
         X=I-1
         DO 46 J=1,NPY
         Y=J-1
         DO 46 K=1,NPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=8.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.NPXM) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(Z.EQ.NPZM.AND.Z.GT.0) W=W/2.D0
   46    WW(NP)=W
      ELSEIF(LAT.EQ.9) THEN
C
C        Base centred orthorombic:
C
C        0 .LE. Kx .LE. Pi/a
C        0 .LE. Ky .LE. 2 Pi/b
C        0 .LE. Kz .LE. Pi/c
C
         AOC=1.D0/COA
         DKX=1.D0/NPXM
         DKY=2.D0/BOA/NPYM
         IF(NPZ.GT.1) DKZ=AOC/NPZM
         DHX=0.D0
         DO 36 I=1,NPX
         X=I-1
         DO 36 J=1,NPY
         Y=J-1
         DO 36 K=1,NPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=8.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.NPXM) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(Z.EQ.NPZM.AND.Z.GT.0) W=W/2.D0
   36    WW(NP)=W
      ELSEIF(LAT.EQ.10) THEN
C
C        Body centered orthorombic
C
C        0 .LE. Kx .LE. Pi/a
C        0 .LE. Ky .LE. a/b Pi/a
C        0 .LE. Kz .LE. 2a/c Pi/a
C
         AOC=1.D0/COA
         AOB=1.D0/BOA
         DKX=1.D0/NPXM
         DKY=AOB/NPYM
         IF(NPZ.GT.1) THEN
            DKZ=2.D0*AOC/NPZM
         ELSE
            DKZ=AOC
         ENDIF
         DHX=0.D0
         DO 47 I=1,NPX
         X=I-1
         DO 47 J=1,NPY
         Y=J-1
         DO 47 K=1,NPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=8.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.NPXM) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(Z.EQ.NPZM.AND.Z.GT.0) W=W/2.D0
   47    WW(NP)=W
      ELSEIF(LAT.EQ.11) THEN
C
C        Face centered orthorombic
C
C        0 .LE. Kx .LE. 2 Pi/a
C        0 .LE. Ky .LE. 2a/b Pi/a
C        0 .LE. Kz .LE. 2a/c Pi/a
C
         AOC=1.D0/COA
         IF(NPX.GT.1) THEN
            DKX=2.D0/NPXM
         ELSE
            DKX=2.d0
         ENDIF
         IF(NPY.GT.1) THEN
            DKY=2.D0/BOA/NPYM
         ELSE
            DKY=2.d0/BOA
         ENDIF
         IF(NPZ.GT.1) THEN
            DKZ=2.D0*AOC/NPZM
         ELSE
            DKZ=AOC
         ENDIF
         DHX=0.D0
         DO 49 I=1,NPX
         X=I-1
         DO 49 J=1,NPY
         Y=J-1
         DO 49 K=1,NPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=8.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.0) W=W/2.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.NPXM) W=W/2.D0
         IF(Y.EQ.NPYM) W=W/2.D0
         IF(Z.EQ.NPZM.AND.Z.GT.0) W=W/2.D0
   49    WW(NP)=W
      ELSEIF(LAT.EQ.12) THEN
C
C        Simple monoclinic:
C
C        0 .LE. Kx .LE. Pi/a
C        -a/b/sin(gamma) Pi/a .LE. Ky .LE. a/b/sin(gamma) Pi/a
C        0 .LE. Kz .LE. a/c Pi/a
C
         IF((NPY/2)*2.NE.NPY) GO TO 998
         NPH=NPY/2
         AOC=1.D0/COA
         AOB=1.D0/BOA
         IF(NPX.GT.1) THEN
            DKX=1.D0/NPXM
         ELSE
            DKX=1.d0
         ENDIF
         DKY=2.D0*AOB/NPY/SIN(gam)
         IF(NPZ.GT.1) THEN
            DKZ=AOC/NPZM
         ELSE
            DKZ=AOC
         ENDIF
         DHX=0.D0
         DO 45 I=1,NPX
         X=I-1
         DO 45 J=0,NPY
         Y=J-NPH
         DO 45 K=1,NPZ
         Z=K-1
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=4.D0
         IF(X.EQ.0) W=W/2.D0
         IF(Y.EQ.-NPH) W=W/2.D0
         IF(Z.EQ.0) W=W/2.D0
         IF(X.EQ.NPXM) W=W/2.D0
         IF(Y.EQ.NPH) W=W/2.D0
         IF(Z.EQ.NPZM) W=W/2.D0
   45    WW(NP)=W
      ELSEIF(LAT.EQ.13) THEN
C
C        Base centred monoclinic:
C
         IF(NPXM.NE.2*(NPXM/2)) GO TO 996
         IF(NPYM.NE.2*(NPYM/2)) GO TO 998
         NPYH=NPYM/2
         TKX(1)=2./SIN(GAM)/NPXM
         TKY(1)=0.
         TKZ(1)=0.
         TKX(2)=2./TAN(GAM)/BOA/NPYM
         TKY(2)=2.0/BOA/NPYM
         TKZ(2)=0.
         TKX(3)=0.
         TKY(3)=0.
         TKZ(3)=2./COA/NPXM
         DO 38 I=1,NPX
         X=I-1
         KM=1+NPX-I
         DO 38 K=1,KM
         Z=K-1
         DO 38 J=-NPYH,NPYH
         Y=J
         NP=NP+1
         IF(NP.GT.MKM) GO TO 999
         KX(NP)=X
         KY(NP)=Y
         KZ(NP)=Z
         W=4.D0
         IF(X.EQ.0) W=W/2.
         IF(Z.EQ.0) W=W/2.
         IF(X+Z.EQ.NPXM) W=W/2.
         IF(Y.EQ.-NPYH) W=W/2.
         IF(Y.EQ.NPYH) W=W/2.
   38    WW(NP)=W
      ELSEIF(LAT.EQ.14) THEN
C
C        Simple triclinic:
C
         ca=COS(alf)
         cb=COS(bet)
         cg=COS(gam)
         sg=SIN(gam)
         aoc=1.d0/coa
         aob=1.d0/boa
         cx=cb
         cy=(ca-cb*cg)/sg
         cz=SQRT(1.d0-ca*ca-cb*cb-cg*cg+2.d0*ca*cb*cg)/sg
C
         dkx=1.d0/npxm
         dky=2.d0*aob/sg/npy
         dkz=2.d0*aoc/cz/npz
         tkx(1)=-dkx*cg/sg
         tkx(2)=dkx*(cy*cg-cx*sg)/cz/sg
         tkx(3)=-dky*cy/cz
         DO 48 i=1,npx
         x=i-1
         DO 48 j=1,npy
         y=j-1
         DO 48 k=1,npz
         z=k-1
         np=np+1
         IF(np.GT.mkm) GO TO 999
         kx(np)=x
         ky(np)=y
         kz(np)=z
         w=2.d0
         IF(x.EQ.0) w=w/2.d0
         IF(x.EQ.npxm) w=w/2.d0
   48    ww(np)=w
      ENDIF
      NPT=NP
      WRITE(M6,106) TBZ(LAT),NPT
      RETURN
  994 WRITE(M6,107) NPX,NPY,NPZ
      STOP
  995 WRITE(M6,105) LAT
      STOP
  996 WRITE(M6,104) NPX
      STOP
  997 WRITE(M6,103) NPZ
      STOP
  998 WRITE(M6,102) NPY
      STOP
  999 WRITE(M6,101) NP,MKM
      STOP
      END
