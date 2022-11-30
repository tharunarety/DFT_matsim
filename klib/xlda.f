      SUBROUTINE xlda(rho,ex,exd,exdd,exddd)
C   ******************************************************************
C   *    Calculate the derivatives of the exchange energy            *
C   *    per electron within the local density approximation:        *
C   *                                                                *
C   *    Ex[n] = int n(r)*ex(n) dr.                                  *
C   *                                                                *
C   *    The Ry atomic units are used. Non spin polarized.           *
C   *    The spin polarized values can be obtained by:               *
C   *                                                                *
C   *    Ex[n1,n2] = {Ex[2n1]+Ex[2n2]}/2  or                         *
C   *                                                                *
C   *    ex(n1,n2) = {n1*ex(2n1)+n2*ex(2n2)}/n.                      *
C   *                                                                *
C   *   *On entri:                                                   *
C   *       n(r) = RHO the charge density.                           *
C   *   *On exit:                                                    *
C   *       ex, ex' , ex'' , ex''' within LDA                        *
C   ******************************************************************
      IMPLICIT NONE
      REAL(KIND=8) :: rho, ex, exd, exdd, exddd
      REAL(KIND=8), PARAMETER :: ax = -0.73855876638202d0
C
      IF(rho.LT.1.d-18) THEN
        ex   = 0.d0
        exd  = 0.d0
        exdd = 0.d0
        exddd= 0.d0
        RETURN
      ENDIF
C
      ex   = ax*(rho**(1./3.))
      exd  = 1.d0/3.d0*ax*(rho**(-2./3.))
      exdd =-2.d0/9.d0*ax*(rho**(-5./3.))
      exddd= 10.d0/27.d0*ax*(rho**(-8./3.))
C
      ex   = 2.d0*ex
      exd  = 2.d0*exd
      exdd = 2.d0*exdd
      exddd= 2.d0*exddd
C
      RETURN
      END
