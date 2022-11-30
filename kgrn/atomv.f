      SUBROUTINE ATOMV(DV,D,DR,DX,Z,FOURPI,RWAT,NP,ION,IEX,IWAT)
C   ******************************************************************
C   *                                                                *
C   *    Solve Poisson's equation and add exchange-correlation       *
C   *    potential                                                   *
C   *                                                                *
C   *    DV   : Potential                                            *
C   *    D    : Charge density                                       *
C   *    DR   : Radial mesh                                          *
C   *    DX   : Step on exponential mesh                             *
C   *    Z    : Atomic number                                        *
C   *    NP   : Number of points                                     *
C   *    ION  : Z - number of electrons                              *
C   *                                                                *
C   ******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999)
      DIMENSION DRHO(MW),DRHOP(MW),DRHOPP(MW),DV(MW),D(MW),DR(MW),
     1          RHOD(2),RHODD(2)
C
      POTB=2.D0*Z/DR(NP)
      CALL POISON(D,DR,DV,POTB,NP,DX)
      DO 20 I=1,NP
   20 DV(I)=0.5D0*DR(I)*DV(I)
C
      DLO=-ION-1
      IXC=IEX+1
C
C     Exchange-correlation potential
C
      IF(IXC.EQ.5.OR.IXC.GE.8) THEN
         DO 42 I=1,NP
         DFPR=FOURPI*DR(I)**2
   42    DRHO(I)=D(I)/DFPR
         CALL DIFFN(DRHO,DRHOP,DRHOPP,NP,DX)
         DO 23 I=1,NP
         DRI=DR(I)
         DR2=DRI*DRI
         RHO=DRHO(I)
         RHO1=RHO/2.D0
         RHO2=RHO1
         RHOD(1)=DRHOP(I)/DRI
         RHODD(1)=(DRHOPP(I)-DRHOP(I))/DR2
         RHOD(1)=RHOD(1)/2.D0
         RHOD(2)=RHOD(1)
         RHODD(1)=RHODD(1)/2.D0
         RHODD(2)=RHODD(1)
         CALL XCPOT(IXC,RHO1,RHO2,RHO,RHOD,RHODD,DRI,VXC1,VXC2,EXC)
         DV(I)=DV(I)-(Z-0.5D0*DRI*VXC1)
C
C        Correct the potential by -(ION+1)/R
C
         IF (DV(I).GT.DLO) DV(I)=DLO
   23    DV(I)=DV(I)/DR(I)
      ELSE
         DO 24 I=1,NP
         DRI=DR(I)
         DR2=FOURPI*DRI*DRI
         RHO=D(I)/DR2
         RHO1=RHO/2.D0
         RHO2=RHO1
         CALL XCPOT(IXC,RHO1,RHO2,RHO,0.D0,0.D0,DRI,VXC1,VXC2,EXC)
         DV(I)=DV(I)-(Z-0.5D0*DRI*VXC1)
C
C        Correct the potential by -(ION+1)/R
C
         IF (DV(I).GT.DLO) DV(I)=DLO
   24    DV(I)=DV(I)/DR(I)
      ENDIF
      IF(IWAT.EQ.0) RETURN
      DO 26 I=1,NP
      IF(DR(I).GT.RWAT) GO TO 27
      DV(I)=DV(I)+ION/RWAT
      GO TO 26
   27 DV(I)=DV(I)+ION/DR(I)
   26 CONTINUE
      RETURN
      END
