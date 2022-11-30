      SUBROUTINE LATT3M
C   ******************************************************************
C   *                                                                *
C   *    Generate vectors of real and reciprocal space from the      *
C   *    primitive translation vectors (BSX,BSY,BSZ) to be used      *
C   *    the Madelung calculation.                                   *
C   *                                                                *
C   *    NPRN  : Print vectors if = 1 or 2.                          *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE csts
      USE message
      USE madelung_lattice
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: CSX, CSY, CSZ, DC
      REAL(KIND=8), DIMENSION(3) :: DD, DK
      REAL(KIND=8) :: R1, PQX, PQY, PQZ, X, Y, Z, DQ, RA, G1, GA, DDM,
     1                DKM, SX, SY, SZ, DX, A, B, C, DA, AMIN, DB, GX,
     2                GY, GZ
      INTEGER :: IQ, JQ, I, NUMR, NUMG, NUMRH, NUMGH, NCR, L, M,
     1           N, NR, NSH, NSHL, K, N1, NGR, NG
  101 FORMAT(/,15X,'Shell number',I5,' with',I5,' vectors',/)
  102 FORMAT(9X,I4,4F10.6)
  103 FORMAT(11X,'(',F10.5,',',F10.5,',',F10.5,')')
  104 FORMAT(/,11X,'Real space vectors:',//,11X,'NO',5X,'SX',8X,'SY',8X,
     1       'SZ',8X,'D',/)
  105 FORMAT(/,11X,'Reciprocal space vectors',//,11X,'NO',5X,'KX',8X,
     1       'KY',8X,'KZ',8X,'D',/)
  107 FORMAT(/,' LATT3M:',3X,'R1    =',F10.4,' RA    =',F10.4,
     1       ' G1    =',F10.4,/,11X,'GA    =',F10.4,/)
  108 FORMAT(/,11X,'NUMR  =',I3,' NUMG  =',I3,' NUMVR =',I3,' NUMVG =',
     1       I3)
  111 FORMAT(/,' LATT3M:** The basis vectors must be distinct:',//,11X,
     1      'Q(',I2,') =',3F10.6,/,11X,'Q(',I2,') =',3F10.6)
C
C     Calculate the radius RA of the spheres holding all the vectors
C     used in the lattice summations. R1 is the longest basis vector.
C
      R1=1.D-06
      DO 35 IQ=1,NQ
      PQX=QX(IQ)
      PQY=QY(IQ)
      PQZ=QZ(IQ)
      DO 35 JQ=IQ,NQ
      IF(IQ.NE.JQ) THEN
         X=PQX-QX(JQ)
         Y=PQY-QY(JQ)
         Z=PQZ-QZ(JQ)
         DQ=SQRT(X*X+Y*Y+Z*Z)
         IF(DQ.GE.R1) THEN
            R1=DQ
         ELSEIF(DQ.LT.1.D-8) THEN
            WRITE(M6,111) IQ,PQX,PQY,PQZ,JQ,QX(JQ),QY(JQ),QZ(JQ)
            STOP
         ENDIF
      ENDIF
   35 CONTINUE
      R1=R1*1.001
      RA=RMAX+R1
      G1=0.D0
      GA=GMAX+G1
      WRITE(M6,107) R1,RA,G1,GA
C
      DO 36 I=1,3
      DD(I)=SQRT(BSX(I)**2+BSY(I)**2+BSZ(I)**2)
   36 DK(I)=SQRT(BKX(I)**2+BKY(I)**2+BKZ(I)**2)
      DDM=MAX(DD(1),DD(2),DD(3))
      DKM=MAX(DK(1),DK(2),DK(3))
      DDM=TWOPI/DDM
      DKM=TWOPI/DKM
      NUMR=2*(INT(RA/DKM)+1)+1
      NUMG=2*(INT(GA/DDM)+1)+1
      NUMRH=NUMR/2+1
      NUMGH=NUMG/2+1
C
C     Real space
C
      NCR=0
      DO 42 L=1,NUMR
      A=L-NUMRH
      DO 42 M=1,NUMR
      B=M-NUMRH
      DO 42 N=1,NUMR
      C=N-NUMRH
      SX=A*BSX(1)+B*BSX(2)+C*BSX(3)
      SY=A*BSY(1)+B*BSY(2)+C*BSY(3)
      SZ=A*BSZ(1)+B*BSZ(2)+C*BSZ(3)
      DX=SQRT(SX*SX+SY*SY+SZ*SZ)
      IF(DX.GT.RA) GO TO 42
      NCR=NCR+1
   42 CONTINUE
      ALLOCATE(CSX(NCR),CSY(NCR),CSZ(NCR),DC(NCR))
      ALLOCATE(ASX(NCR),ASY(NCR),ASZ(NCR),DR(NCR))

      IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,104)
      NR=0
      NR0=0
      DO 22 L=1,NUMR
      A=L-NUMRH
      DO 22 M=1,NUMR
      B=M-NUMRH
      DO 22 N=1,NUMR
      C=N-NUMRH
      SX=A*BSX(1)+B*BSX(2)+C*BSX(3)
      SY=A*BSY(1)+B*BSY(2)+C*BSY(3)
      SZ=A*BSZ(1)+B*BSZ(2)+C*BSZ(3)
      DX=SQRT(SX*SX+SY*SY+SZ*SZ)
      IF(DX.GT.RA) GO TO 22
      IF(DX.LE.RMAX) NR0=NR0+1
      NR=NR+1
      DC(NR)=DX
      CSX(NR)=SX
      CSY(NR)=SY
      CSZ(NR)=SZ
   22 CONTINUE
C
C     Sort vectors in order of increasing length
C
      DA=1.D-06
      NSH=0
      NSHL=-1
      DO 23 K=1,NR
      AMIN=1000.
      DO 24 N=1,NR
      IF(DC(N)-AMIN)25,24,24
   25 AMIN=DC(N)
      N1=N
   24 CONTINUE
      NSHL=NSHL+1
      ASX(K)=CSX(N1)
      ASY(K)=CSY(N1)
      ASZ(K)=CSZ(N1)
      DB=DC(N1)
      DR(K)=DB
      IF(DB.GT.DA+1.D-06) GO TO 26
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,102)K, ASX(K),ASY(K),ASZ(K),
     1                                        DB
      GO TO 23
   26 NSH=NSH+1
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) THEN
         WRITE(M6,101) NSH,NSHL
         WRITE(M6,102)K, ASX(K),ASY(K),ASZ(K),DB
      ENDIF
      NSHL=0
      DA=DB
   23 DC(N1)=1000.
      NSH=NSH+1
      NSHL=NSHL+1
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,101) NSH,NSHL
      NUMVR=NR
      DEALLOCATE(csx,csy,csz,dc)
C
C     Reciprocal space
C
      NGR=0
      DO 47 L=1,NUMG
      A=L-NUMGH
      DO 47 M=1,NUMG
      B=M-NUMGH
      DO 47 N=1,NUMG
      C=N-NUMGH
      GX=A*BKX(1)+B*BKX(2)+C*BKX(3)
      GY=A*BKY(1)+B*BKY(2)+C*BKY(3)
      GZ=A*BKZ(1)+B*BKZ(2)+C*BKZ(3)
      DX=SQRT(GX*GX+GY*GY+GZ*GZ)
      IF(DX.GT.GA) GO TO 47
      NGR=NGR+1
   47 CONTINUE
      allocate(csx(ngr),csy(ngr),csz(ngr),dc(ngr))
      allocate(akx(ngr),aky(ngr),akz(ngr),dg(ngr))
C
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,105)
      NG=0
      DO 27 L=1,NUMG
      A=L-NUMGH
      DO 27 M=1,NUMG
      B=M-NUMGH
      DO 27 N=1,NUMG
      C=N-NUMGH
      GX=A*BKX(1)+B*BKX(2)+C*BKX(3)
      GY=A*BKY(1)+B*BKY(2)+C*BKY(3)
      GZ=A*BKZ(1)+B*BKZ(2)+C*BKZ(3)
      DX=SQRT(GX*GX+GY*GY+GZ*GZ)
      IF(DX.GT.GA) GO TO 27
      NG=NG+1
      DC(NG)=DX
      CSX(NG)=GX
      CSY(NG)=GY
      CSZ(NG)=GZ
   27 CONTINUE
C
C     Sort vectors in order of increasing length
C
      DA=1.E-06
      NSH=0
      NSHL=-1
      DO 28 K=1,NG
      AMIN=1000.
      DO 29 N=1,NG
      IF(DC(N)-AMIN)30,29,29
   30 AMIN=DC(N)
      N1=N
   29 CONTINUE
      NSHL=NSHL+1
      AKX(K)=CSX(N1)
      AKY(K)=CSY(N1)
      AKZ(K)=CSZ(N1)
      DB=DC(N1)
      DG(K)=DB
      IF(DB.GT.DA*1.000001) GO TO 31
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,102) K,AKX(K),AKY(K),AKZ(K),
     1                                         DB
      GO TO 28
   31 NSH=NSH+1
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) THEN
        WRITE(M6,101) NSH,NSHL
        WRITE(M6,102) K,AKX(K),AKY(K),AKZ(K),DB
      ENDIF
      NSHL=0
      DA=DB
   28 DC(N1)=1000.
      NSH=NSH+1
      NSHL=NSHL+1
      IF(NPRN.EQ.1.OR.NPRN.EQ.2) WRITE(M6,101) NSH,NSHL
      NUMVG=NG
      WRITE(M6,108) NUMR,NUMG,NUMVR,NUMVG
      DEALLOCATE(CSX,CSY,CSZ,DC)
      RETURN
      END
