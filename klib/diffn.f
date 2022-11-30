      SUBROUTINE DIFFN(RHO,RHOP,RHOPP,N,H)
C   *******************************************************************
C   *                                                                 *
C   *   Differentiate charge density for use in XCPOT.                *
C   *                                                                 *
C   *   HLS:  24-Mar-88                                               *
C   *******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RHO(N),RHOP(N),RHOPP(N)
      NM2=N-2
      SIXH=6.D0*H
      TWLWH=12.D0*H
      TWLWH2=TWLWH*H
      H2=H*H
C
C     Forward difference at the beginning of the table
C
      RHOP(1)=((2.D0*RHO(4)+18.D0*RHO(2))-(9.D0*RHO(3)
     1       +11.D0*RHO(1)))/SIXH
      RHOP(2)=((2.D0*RHO(5)+18.D0*RHO(3))-(9.D0*RHO(4)
     1       +11.D0*RHO(2)))/SIXH
      RHOPP(1)=((2.D0*RHO(1)+4.D0*RHO(3))-(5.D0*RHO(2)+RHO(4)))/H2
      RHOPP(2)=((2.D0*RHO(2)+4.D0*RHO(4))-(5.D0*RHO(3)+RHO(5)))/H2
C
C     Central difference at the interior of the table
C
      DO 20 I=3,NM2
      RHOP(I)=((RHO(I-2)+8.D0*RHO(I+1))-
     1        (8.D0*RHO(I-1)+RHO(I+2)))/TWLWH
      RHOPP(I)=((16.D0*RHO(I+1)+16.D0*RHO(I-1))-
     1         (RHO(I+2)+RHO(I-2)+30.D0*RHO(I)))/TWLWH2
   20 CONTINUE
C
C     Backward difference at the end of the table
C
      RHOP(N)=((11.D0*RHO(N)+9.D0*RHO(N-2))-
     1        (18.D0*RHO(N-1)+2.D0*RHO(N-3)))/SIXH
      RHOP(N-1)=((11.D0*RHO(N-1)+9.D0*RHO(N-3))-
     1          (18.*RHO(N-2)+2.d0*RHO(N-4)))/SIXH
      RHOPP(N)=((2.D0*RHO(N)+4.D0*RHO(N-2))-
     1         (5.D0*RHO(N-1)+RHO(N-3)))/H2
      RHOPP(N-1)=((2.D0*RHO(N-1)+4.*RHO(N-3))-
     1           (5.D0*RHO(N-2)+RHO(N-4)))/H2
      RETURN
      END
