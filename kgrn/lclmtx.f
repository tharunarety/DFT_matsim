      SUBROUTINE lclmtx(drs,drsd,chi,e,l,ita,it,is,km)
C   ******************************************************************
C   *                                                                *
C   * Calculate the logarithmic derivative of the forward            *
C   * extrapolated free electron solution for complex energy z(lz):  *
C   *                                                                *
C   *   fi^m(r) = fi^m(s^m)*[f^m(r) + g^m(r) * D{fi0^m(s^m)}],       *
C   *                                                                *
C   * where :                                                        *
C   *                                                                *
C   *   f^m = t^m(1)*n + t^m(2)*j and g^m = - t^m(3)*n - t^m(4)*j    *
C   *                                                                *
C   * t^m is the screening matrix for s^m and                        *
C   *                                                                *
C   *                        D{fi(s)} - D{f^m(s)}    f^m(s)          *
C   *   D{fi0^m(s^m)} =  -  ---------------------- * ------          *
C   *                        D{fi(s)} - D{g^m(s)}    g^m(s)          *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE logderivative
      USE partialwaves
      USE potential
      USE radialmesh
      USE slope
      IMPLICIT NONE
      INTEGER :: l, ita, it, is
      REAL(KIND=8) :: w2, sm, drs, drsd
      REAL(KIND=8), DIMENSION(4,0:1) :: km, ks
      REAL(KIND=8) :: e, z, d, dd, xn, xs, chi
      REAL(KIND=8) :: x1, x2, x3, x4, x1d, x2d, x3d, x4d
C
C     Get new logarithmic derivatives for types with s^m>s
C
      w2=sws*sws
C
      z=(e-vmtzr(it,is))*w2
      sm=wsm(ita,it)/sws
      CALL trmtrx(z,sm,km,l,itrans,0)
C
      sm=hsr(ita,it)/sws
      CALL trmtrx(z,sm,ks,l,itrans,0)
C
      km(1:4,1)=w2*km(1:4,1)
      ks(1:4,1)=w2*ks(1:4,1)
C
      x1= km(1,0)*ks(4,0)-km(2,0)*ks(3,0)
      x2=-km(1,0)*ks(2,0)+km(2,0)*ks(1,0)
      x3= km(3,0)*ks(4,0)-km(4,0)*ks(3,0)
      x4=-km(3,0)*ks(2,0)+km(4,0)*ks(1,0)
C
      x1d= km(1,0)*ks(4,1)-km(2,0)*ks(3,1)
     .    +km(1,1)*ks(4,0)-km(2,1)*ks(3,0)
      x2d=-km(1,0)*ks(2,1)+km(2,0)*ks(1,1)
     .    -km(1,1)*ks(2,0)+km(2,1)*ks(1,0)
      x3d= km(3,0)*ks(4,1)-km(4,0)*ks(3,1)
     .    +km(3,1)*ks(4,0)-km(4,1)*ks(3,0)
      x4d=-km(3,0)*ks(2,1)+km(4,0)*ks(1,1)
     .    -km(3,1)*ks(2,0)+km(4,1)*ks(1,0)
C
      xs=drs*x1+x2
      xn=drs*x3+x4
      d=xs/xn
      dd=(drsd*x1+drs*x1d+x2d)/xn-d*(drsd*x3+drs*x3d+x4d)/xn
C
      drs=d
      drsd=dd
      chi=-0.5*sws/wsm(ita,it)*chi*xn
C
      RETURN
      END
