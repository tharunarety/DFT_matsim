      FUNCTION THFPOT(R,Z,WA)
C   ******************************************************************
C   *                                                                *
C   *    Thomas-Fermi potential.                                     *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    R     : Radius.                                             *
C   *    WA    : Number of electrons -Z-1                            *
C   *    Z     : Atomic number.                                      *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    THFPOT  : Thomas-Fermi potential in Hartree                 *
C   *                                                                *
C   ******************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      S=Z
      WC=SQRT((R*(S+WA)**(1.D0/3.D0))/0.8853D0)
      WD=WC*(0.60112D0*WC+1.81061D0)+1.D0
      WE=WC*(WC*(WC*(WC*(0.04793D0*WC+0.21465D0)
     1  +0.77112D0)+1.39515D0)+1.81061D0)+1.D0
      WC=(Z+WA)*(WD/WE)**2-WA
      THFPOT=-WC/R
      RETURN
      END
