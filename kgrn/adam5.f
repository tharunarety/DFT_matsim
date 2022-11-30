      SUBROUTINE ADAM5(P,Q,V,R)
C   ******************************************************************
C   *                                                                *
C   *   Adams 5 point integration.                                   *
C   *                                                                *
C   *   P:   The large component                                     *
C   *   Q:   The small component                                     *
C   *   R:   Mesh point                                              *
C   *   V:   Potential                                               *
C   *   DEP: Derivative of the large component                       *
C   *   DEQ: Derivative of the small component                       *
C   *   BE:  Energy/Speed of light                                   *
C   *   VC:  Speed of light                                          *
C   *   TVC: 2.* Speed of light                                      *
C   *   DK:  Combined spin and orbital quantum number                *
C   *   DM:  Exponential step/720.                                   *
C   *                                                                *
C   ******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(AA=475.D0/502,AB=27.D0/502)
      COMMON/ADAM/DEP(5),DEQ(5),BE,VC,TVC,DK,DM
C
      DPR=P+DM*((251.D0*DEP(1)+2616.D0*DEP(3)+1901.D0*DEP(5))
     1   -(1274.D0*DEP(2)+2774.D0*DEP(4)))
      DQR=Q+DM*((251.D0*DEQ(1)+2616.D0*DEQ(3)+1901.D0*DEQ(5))
     1   -(1274.D0*DEQ(2)+2774.D0*DEQ(4)))
      DO 20 I=2,5
      DEP(I-1)=DEP(I)
   20 DEQ(I-1)=DEQ(I)
      DSUM=(BE-V/VC)*R
      DEP(5)=-DK*DPR+(TVC*R+DSUM)*DQR
      DEQ(5)=DK*DQR-DSUM*DPR
      P=P+DM*((106.D0*DEP(2)+646.D0*DEP(4)+251.D0*DEP(5))
     1  -(19.D0*DEP(1)+264.D0*DEP(3)))
      Q=Q+DM*((106.D0*DEQ(2)+646.D0*DEQ(4)+251.D0*DEQ(5))
     1  -(19.D0*DEQ(1)+264.D0*DEQ(3)))
      P=AA*P+AB*DPR
      Q=AA*Q+AB*DQR
      DEP(5)=-DK*P+(TVC*R+DSUM)*Q
      DEQ(5)=DK*Q-DSUM*P
      RETURN
      END
