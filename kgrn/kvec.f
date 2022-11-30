      SUBROUTINE kvec(lk,pkx,pky,pkz,weight)
C   ******************************************************************
C   *                                                                *
C   *    Transformation from symmetry to rectangular coordinates.    *
C   *                                                                *
C   ******************************************************************
      USE bzmesh
      USE bzmesh2
      USE control_data
      USE control_text
      USE csts
      IMPLICIT NONE
      INTEGER :: lk, lper, lpar, k1, k2, k3
      REAL(KIND=8) :: pkx, pky, pkz, weight
C
      IF(mode.EQ.'3D') THEN
         weight=ww(lk)
         k1=kx(lk)
         k2=ky(lk)
         k3=kz(lk)
         IF(ibz.NE.14) THEN
            pkx=dkx*k1+dhx*IABS(k2)
            pky=dky*k2
            pkz=dkz*k3
         ELSE
            pkx=k1*dkx
            pky=k1*tkx(1)+k2*dky
            pkz=k1*tkx(2)+k2*tkx(3)+k3*dkz
         ENDIF
         pkx=pkx*pi
         pky=pky*pi
         pkz=pkz*pi
      ELSE
         lper=(lk-1)/npar+1
         lpar=lk-(lper-1)*npar
         weight=wk(lpar)*dkz2*dkp
         pkx=akx(lpar)
         pky=aky(lpar)
         pkz=akz(lper)
      ENDIF
C
      RETURN
      END
