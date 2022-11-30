      SUBROUTINE KMESH2(IBZ2,NKX2,NKY2,NKZ2,NKVEC,MPAR,MPER)
C   ******************************************************************
C   *                                                                *
C   *    Construct a 3D k-mesh suitable for the calculation of the   *
C   *    ideal Green's function for a crystal surface. The mesh is   *
C   *    combined of a grid in the 2D Brillouin zone specified by    *
C   *    IBZ2 and a subdivision of the primitive reciprocal trans-   *
C   *    lation perpendicular to the surface.                        *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    IBZ2  :=1 Hexagonal, 1/12 BZ                                *
C   *          :=201 Hexagonal, 1/6 BZ                               *
C   *          :=2 Square                                            *
C   *          :=3 Rectangular P                                     *
C   *          :=4 Rectangular C                                     *
C   *          :=5 Oblique                                           *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    AKX   : x component of 2D grid                              *
C   *    AKY   : y component of 2D grid                              *
C   *    AKZ   : Component perpendicular to the surface plane        *
C   *                                                                *
C   ******************************************************************
      USE bzmesh2
      USE lattice
      USE message
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ONE=1)
      CHARACTER TXTM*13
      COMMON/KSPC/BKX(3),BKY(3),BKZ(3)
      DIMENSION CKX(3),CKY(3),CKZ(3),DOT(3),TXTM(5)
      DATA TXTM/'Hexagonal    ','Square       ','Rectangular P',
     1          'Rectangular C','Oblique      '/
  101 FORMAT(//,'KMESH2:',4X,'Distance between surface planes,',//,11X,
     1       'D-perp.=',F10.5,' in units of a.')
  102 FORMAT(10X,'(',F10.5,',',F10.5,',',F10.5,'  )')
  105 FORMAT(/,11X,'Primitive 2D reciprocal translations',/,11X,
     1       'in units of 2*pi/a:',/)
  106 FORMAT(/,' KMESH2:** NKZ2 =',I4,' too large. Increase MPAR =',
     1       I4,'.')
  107 FORMAT(/,11X,'2D Brillouin zone: ',A,//,11X,I4,' points in',
     1       ' the 2D zone',/,11X,I4,' points along K-perpendicular')
  108 FORMAT(/,' KMESH2:** IBZ2 =',I2,' not implemented')
  109 FORMAT(/,' KMESH2:** ',A,I3,2A,I4)
  110 FORMAT(/,' KMESH2::** Angle between primitive vectors .GT. Pi/4',
     1       /,11X,'Change input to BSTR.')
C
      PI=ACOS(-ONE)
      TWOPI=2.D0*PI
      IF(IBZ2.GT.200) THEN
C
C        IBZ(hex) = 1/6
C
         IBZ2=IBZ2-200
         IHEX=2
      ELSE
C
C        IBZ(hex) = 1/12
C
         IHEX=1
      ENDIF
C
C     Check IBZ2 and NKX2
C
      IF(IBZ2.LT.1.OR.IBZ2.GT.5) THEN
         WRITE(M6,109) 'IBZ2 =',IBZ2,
     1                '. Must be a number between 1 and 5'
         STOP
      ENDIF
      IF(IBZ2.NE.1) THEN
         ANKX=NKX2
         AN=LOG(ANKX)/LOG(2.D0)
         N=AN
         DIFF=ABS(AN-N)
         IF(DIFF.GT.1E-4) THEN
            WRITE(M6,109) 'NKX2 =',NKX2,' not properly chosen.',
     1                   ' Must be a power of 2'
            STOP
         ENDIF
      ENDIF
C
C     Determine surface reciprocal translations
C
      DO 123 I=1,3
  123 DOT(I)=BKX(I)*BKX(3)+BKY(I)*BKY(3)+BKZ(I)*BKZ(3)
      DSURF=TWOPI/SQRT(DOT(3))
      DO 124 I=1,2
      CKX(I)=BKX(I)-DOT(I)*BKX(3)/DOT(3)
      CKY(I)=BKY(I)-DOT(I)*BKY(3)/DOT(3)
  124 CKZ(I)=BKZ(I)-DOT(I)*BKZ(3)/DOT(3)
      CKX(3)=BKX(3)
      CKY(3)=BKY(3)
      CKZ(3)=BKZ(3)
      WRITE(M6,101) DSURF
      WRITE(M6,105)
      DO 125 I=1,3
  125 WRITE(M6,102) CKX(I)/TWOPI,CKY(I)/TWOPI,CKZ(I)/TWOPI
C
C     Select the 2D Brillouin zone
C
      GO TO (1,2,3,4,5),IBZ2
C
C     Hexagonal
C
    1 LPAR=0
      IF(IHEX.EQ.2) THEN
         NPH1=-NPH
      ELSE
         NPH1=0
      ENDIF
      FAH=ABS(1.D0/BSX(1))
      WEIGHT=0.D0
      NDX=3**(NKX2/2+1)
      NPH=NDX/2
      DKX=TWOPI*FAH*(2.D0/3.D0)/NDX
      DKY=DSQRT(3.D0)/2.D0*DKX
      DHX=-0.5D0*DKX
      IF(2*(NKX2/2).NE.NKX2) THEN
         DO 10 LYS=NPH1,NPH
         LY=IABS(LYS)
         IXB=2*LY
         DO 10 LX=IXB,NDX
         LXY=LX+LY
         IF(3*(LXY/3).NE.LXY) THEN
            LPAR=LPAR+1
            IF(LPAR.GT.MPAR) THEN
               WRITE(M6,109) 'LPAR =',LPAR,' too large.',
     1                      'Increase MPAR =',MPAR
               STOP
            ENDIF
            AKX(LPAR)=LX*DKX+LY*DHX
            AKY(LPAR)=LYS*DKY
            IF(LY.EQ.0.AND.NPH1.EQ.0) THEN
               WKPAR=1.0D0
            ELSEIF(LX.EQ.NDX) THEN
               WKPAR=1.0D0
            ELSE
               WKPAR=2.0D0
            ENDIF
            WK(LPAR)=WKPAR
            WEIGHT=WEIGHT+WKPAR
         ENDIF
   10    CONTINUE
      ELSE
         DO 11 LYS=NPH1,NPH
         LY=IABS(LYS)
         IF(3*(LY/3).NE.LY) THEN
            IXB=2*LY
            DO 12 LX=IXB,NDX
            LXY=LX+LY
            IF(3*(LXY/3).EQ.LXY) THEN
               LPAR=LPAR+1
               IF(LPAR.GT.MPAR) THEN
                  WRITE(M6,109) 'LPAR =',LPAR,' too large.',
     1                         'Increase MPAR =',MPAR
                  STOP
               ENDIF
               AKX(LPAR)=LX*DKX+LY*DHX
               AKY(LPAR)=LYS*DKY
               IF(LX.EQ.2*LY) THEN
                  WKPAR=1.0D0
               ELSE
                  WKPAR=2.0D0
               ENDIF
               WK(LPAR)=WKPAR
               WEIGHT=WEIGHT+WKPAR
            ENDIF
   12       CONTINUE
         ENDIF
   11    CONTINUE
      ENDIF
      NPAR=LPAR
      DO 13 LPAR=1,NPAR
   13 WK(LPAR)=WK(LPAR)/WEIGHT
      GO TO 128
C
C     Square lattice
C
    2 LPAR=0
      DX=.5D0/NKX2
      XOFFST=0.5D0*DX*(CKX(1)+CKX(2))
      YOFFST=0.5D0*DX*(CKY(1)+CKY(2))
      WEIGHT=0.
      DO 26 J=1,NKX2
      B=(J-1)*DX
      X2=B*CKX(2)+XOFFST
      Y2=B*CKY(2)+YOFFST
      DO 26 I=J,NKX2
      LPAR=LPAR+1
      IF(LPAR.GT.MPAR) THEN
         WRITE(M6,109) 'LPAR =',LPAR,' too large.',
     1                'Increase MPAR =',MPAR
         STOP
      ENDIF
      A=(I-1)*DX
      AKX(LPAR)=A*CKX(1)+X2
      AKY(LPAR)=A*CKY(1)+Y2
      IF(I.EQ.J) THEN
         WKPAR=1.D0
      ELSE
         WKPAR=2.D0
      ENDIF
      WK(LPAR)=WKPAR
      WEIGHT=WEIGHT+WKPAR
   26 CONTINUE
      NPAR=LPAR
      DO 24 LPAR=1,NPAR
   24 WK(LPAR)=WK(LPAR)/WEIGHT
      GO TO 128
C
C     Rectangular p
C

    3 LPAR=0
      DX=.5D0/NKX2
      XOFFST=.5D0*DX*(CKX(1)+CKX(2))
      YOFFST=.5D0*DX*(CKY(1)+CKY(2))
      DO 30 J=1,NKX2
      B=(J-1)*DX
      X2=B*CKX(2)+XOFFST
      Y2=B*CKY(2)+YOFFST
      DO 30 I=1,NKX2
      LPAR=LPAR+1
      IF(LPAR.GT.MPAR) THEN
         WRITE(M6,109) 'LPAR =',LPAR,' too large.',
     1                'Increase MPAR =',MPAR
         STOP
      ENDIF
      A=(I-1)*DX
      AKX(LPAR)=A*CKX(1)+X2
   30 AKY(LPAR)=A*CKY(1)+Y2
      NPAR=LPAR
      DO 31 LPAR=1,NPAR
   31 WK(LPAR)=1.D0/NPAR
      GO TO 128
C
C     Rectangular C
C
    4 LPAR=0
      BDOT=(BSX(1)*BSX(2)+BSY(1)*BSY(2))
      BS1=SQRT(BSX(1)*BSX(1)+BSY(1)*BSY(1))
      BS2=SQRT(BSX(2)*BSX(2)+BSY(2)*BSY(2))
      BDOT=BDOT/BS1/BS2
      IF(BDOT.LT.COS(PI/4.D0)) THEN
         WRITE(M6,110)
         STOP
      ENDIF
C
C     Interchange primitive vectors if needed
C
      IF(BS1.LT.BS2) THEN
         DUMMY=CKX(1)
         CKX(1)=CKX(2)
         CKX(2)=DUMMY
         DUMMY=CKY(1)
         CKY(1)=CKY(2)
         CKY(2)=DUMMY
      ENDIF
C
C     Rotate coordinate system
C
      RCKY2=SQRT(CKX(2)*CKX(2)+CKY(2)*CKY(2))
      SINT=-CKX(2)/RCKY2
      COST=CKY(2)/RCKY2
      RCKX1=CKX(1)*COST+CKY(1)*SINT
      RCKY1=-CKX(1)*SINT+CKY(1)*COST
      DELTX=RCKX1/NKX2
      DELTY=RCKY2/4/NKX2
      YOFFST=0.5D0*DELTY
      XOFFST=0.5D0*DELTX
      SUMCY=RCKY1+RCKY2
      YMAX=0.5D0*(SUMCY+RCKX1*RCKX1/SUMCY)
      DO 33 I=1,NKX2
      RAKX=XOFFST+(I-1)*DELTX
      YLIN=YMAX-RCKX1*RAKX/SUMCY
      NYMAX=(YLIN-YOFFST)/DELTY
      NYMAX=NYMAX+1
      DO 34 J=1,NYMAX
      LPAR=LPAR+1
      IF(LPAR.GT.MPAR) THEN
         WRITE(M6,109) 'LPAR =',LPAR,' too large.',
     1                'Increase MPAR =',MPAR
         STOP
      ENDIF
      RAKY=YOFFST+(J-1)*DELTY
      AKX(LPAR)=COST*RAKX-SINT*RAKY
   34 AKY(LPAR)=SINT*RAKX+COST*RAKY
   33 CONTINUE
      NPAR=LPAR
C
C     Restore original primitive k-vectors
C
      IF(BS1.LT.BS2) THEN
         DUMMY=CKX(1)
         CKX(1)=CKX(2)
         CKX(2)=DUMMY
         DUMMY=CKY(1)
         CKY(1)=CKY(2)
         CKY(2)=DUMMY
      ENDIF
      DO 35 LPAR=1,NPAR
   35 WK(LPAR)=1.D0/NPAR
      GO TO 128
C
    5 WRITE(M6,108) IBZ2
      STOP
C
C     Subdivision perpendicular to the surface
C
  128 NPER=NKZ2-1
      IF(NPER.GT.MPER) THEN
         WRITE(M6,106) NPER,MPER
         STOP
      ENDIF
      NKVEC=NPAR*NPER
      if(nper.EQ.0) nkvec=npar
      DKP=1.D0/CKZ(3)
      IF(NKZ2.GT.1) THEN
         DKZ2=CKZ(3)/(NKZ2-1)
      ELSE
         DKZ2=0.D0
      ENDIF
      DO 127 K=1,NKZ2
  127 AKZ(K)=(K-1)*DKZ2
      WRITE(M6,107) TXTM(IBZ2),NPAR,NPER
C
      WRITE(M6,'(/,12X,A,4X,A,7X,A,7X,A,/)') 'LPAR','AKX','AKY','WK'
      DO 129 LPAR=1,NPAR
  129 WRITE(M6,'(10X,I5,3F10.6)') LPAR,AKX(LPAR)/TWOPI,
     1                           AKY(LPAR)/TWOPI,WK(LPAR)
      WRITE(M6,'( )')
      RETURN
      END
