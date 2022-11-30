      SUBROUTINE POISON(YSQ,DR,W,POTB,JRI,DX)
C   ******************************************************************
C   *                                                                *
C   *    Solve the Poisson equation as described by T. Loucks in:    *
C   *    "The Augmented Plane Wave Method", Benjamin, New York       *
C   *    (1967)  P.98                                                *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    YSQ   : 4.*Pi*Radius**2*Charge density.                     *
C   *    Z     : Atomic number.                                      *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    W     : Electrostatic potential from nucleus and electron   *
C   *            charge density including Madelung contribution.     *
C   *                                                                *
C   ******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999)
      COMMON/POISA/A,C,C2,F1,EDL
      DIMENSION YSQ(MW),DR(MW),W(MW),E(MW),F(MW)
C
      E(1)=.0D0
      F(1)=F1
C
C     Find charge inside R(JRI)
C
      DO 24 IR=1,JRI
   24 W(IR)=YSQ(IR)*DR(IR)
      CALL SIMPN(W,DX,JRI,QJRI)
C
C     Solve Poisson
C
      ITOP=JRI-1
      DO 21 J=2,ITOP
      D=C*SQRT(DR(J))*(EDL*YSQ(J+1)+10.D0*YSQ(J)+YSQ(J-1)/EDL)
      F(J)=C2-1.D0/F(J-1)
   21 E(J)=(D/A+E(J-1))/F(J)
      W(JRI)=2.D0*QJRI/SQRT(DR(JRI))
      DO 22 J=1,ITOP
      JV=JRI-J
   22 W(JV)=E(JV)+W(JV+1)/F(JV)
C
C     Impose proper boundary conditions
C
      VS=W(JRI)/SQRT(DR(JRI))
      DV=POTB-VS
      DO 23 IR=1,JRI
      R=DR(IR)
   23 W(IR)=W(IR)/SQRT(R)+DV
      RETURN
      END
