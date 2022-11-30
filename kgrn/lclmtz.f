      SUBROUTINE lclmtz(tm,zm,nzm)
C   ******************************************************************
C   *                                                                *
C   * Calculate the logarithmic derivative of the forward            *
C   * extrapolated free electron solution for complex energies ZM:   *
C   *                                                                *
C   *   fi^m(r) = fi^m(s^m)*[f^m(r) + g^m(r) * D{fi0^m(s^m)}]        *
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
      INTEGER :: nzm, lz, ita, it, l, is, jsr
      REAL(KIND=8) :: w2, sm
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8), DIMENSION(4,0:lmax,0:1) :: km, ks
      COMPLEX(KIND=8), DIMENSION(nzm,4,0:lmax,mnta,nt,ns) :: tm
      COMPLEX(KIND=8) :: z, d, dd, xn, xs
      COMPLEX(KIND=8) :: x1, x2, x3, x4, x1d, x2d, x3d, x4d
C
C     Get new logarithmic derivatives
C
      w2=sws*sws
      DO 20 it=1,nt
      DO 21 ita=1,nta(it)
      jsr=jsrs(ita,it)
      IF(localmt(ita,it).EQ.0) GO TO 21
C
C     Loop for spin and complex energy
C
      DO 22 is=1,ns
      DO 22 lz=1,nzm
C
C    "m" screening, local kappa
C
      z=(zm(lz)-vmtzr(it,is))*w2
      sm=wsm(ita,it)/sws
      CALL trmtrzr(z,sm,km,itrans)
      km(1:4,0:lmax,1)=w2*km(1:4,0:lmax,1)
C
      tm(lz,1:4,0:lmax,ita,it,is)=km(1:4,0:lmax,0)
C
C     's' screening, local kappa
C
      sm=hsr(ita,it)/sws
      CALL trmtrzr(z,sm,ks,itrans)
      ks(1:4,0:lmax,1)=w2*ks(1:4,0:lmax,1)
C
C     Loop for l
C
      DO 23 l=0,lmax
C
      x1= km(1,l,0)*ks(4,l,0)-km(2,l,0)*ks(3,l,0)
      x2=-km(1,l,0)*ks(2,l,0)+km(2,l,0)*ks(1,l,0)
      x3= km(3,l,0)*ks(4,l,0)-km(4,l,0)*ks(3,l,0)
      x4=-km(3,l,0)*ks(2,l,0)+km(4,l,0)*ks(1,l,0)
C
      x1d= km(1,l,0)*ks(4,l,1)-km(2,l,0)*ks(3,l,1)
     .    +km(1,l,1)*ks(4,l,0)-km(2,l,1)*ks(3,l,0)
      x2d=-km(1,l,0)*ks(2,l,1)+km(2,l,0)*ks(1,l,1)
     .    -km(1,l,1)*ks(2,l,0)+km(2,l,1)*ks(1,l,0)
      x3d= km(3,l,0)*ks(4,l,1)-km(4,l,0)*ks(3,l,1)
     .    +km(3,l,1)*ks(4,l,0)-km(4,l,1)*ks(3,l,0)
      x4d=-km(3,l,0)*ks(2,l,1)+km(4,l,0)*ks(1,l,1)
     .    -km(3,l,1)*ks(2,l,0)+km(4,l,1)*ks(1,l,0)
C
      xs=dfi(lz,ita,it,l,is,1)*x1+x2
      xn=dfi(lz,ita,it,l,is,1)*x3+x4
      d=xs/xn
      dd=(dfi(lz,ita,it,l,is,2)*x1+dfi(lz,ita,it,l,is,1)*x1d+x2d)/xn-
     . d*(dfi(lz,ita,it,l,is,2)*x3+dfi(lz,ita,it,l,is,1)*x3d+x4d)/xn
C
      dfi(lz,ita,it,l,is,1)=d
      dfi(lz,ita,it,l,is,2)=dd
      fi0m(lz,ita,it,l,is)=-0.5*sws/wsm(ita,it)*
     .                      cf(jsr,l,ita,it,is,lz)*xn
C
   23 CONTINUE
   22 CONTINUE
   21 CONTINUE
   20 CONTINUE
C
      RETURN
      END
