      SUBROUTINE clda(rho,ec,ecd,ecdd,ecddd)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the derivatives of the correlation energy         *
C   *    per electron within the local density approximation.        *
C   *    The Ry atomic units are used.                               *
C   *                                                                *
C   *   *On entri:                                                   *
C   *       n(r) = RHO the charge density.                           *
C   *   *On exit:                                                    *
C   *       ec, ec' , ec'' , ec''' within LDA (Perdew-Wang 1992)     *
C   ******************************************************************
      IMPLICIT NONE
      REAL(KIND=8) :: rho, ec, ecd, ecdd, ecddd
      REAL(KIND=8) :: rs, rsd, rsdd, rsddd, rs12, rs32, rs52, ec0,
     .                ec0d, ec0dd, ec0ddd, q0, q0rs, q1, q2, q3, q4,
     .                q5, q12, q12d, fourpi
      REAL(KIND=8), PARAMETER :: ax = -0.73855876638202d0,
     .              a = 0.0310907d0, a1 = 0.21370d0, b1 = 7.5957d0,
     .              b2 = 3.5876d0, b3 = 1.6382d0, b4 = 0.49294d0
C
      IF(rho.LT.1.d-18) THEN
        ec   = 0.d0
        ecd  = 0.d0
        ecdd = 0.d0
        ecddd= 0.d0
        RETURN
      ENDIF
C
      fourpi=4.d0*DACOS(-1.d0)
      rs   = 1.d0/((fourpi*rho/3.d0)**(1./3.))
      rsd  = -rs/3.d0/rho
      rsdd = 4.d0/9.d0*rs/rho/rho
      rsddd= -28.d0/27.d0*rs/rho/rho/rho
C
      q0   = -2.d0*a*(1.d0+a1*rs)
      q0rs = -2.d0*a*a1
      rs12 = SQRT(rs)
      rs32 = rs12**3
      q1   = 2.d0*a*(b1*rs12+b2*rs+b3*rs32+b4*rs*rs)
      q2   = LOG(1.d0+1.d0/q1)
      ec0  = q0*q2
      q3   = a*(b1/rs12+2.d0*b2+3.d0*b3*rs12+4.d0*b4*rs)
      q12  = q1*(1.d0+q1)
      ec0d = q0rs*q2-q0*q3/q12
      q12d = 1.d0+2.d0*q1
      q4   = a*(-0.5d0*b1/rs32+1.5d0*b3/rs12+4.d0*b4)
      ec0dd= -2.d0*q0rs*q3/q12+q0*q3*q3/q12/q12*q12d-q0*q4/q12
      rs52 = rs12**5
      q5   = a*(0.75d0*b1/rs52-0.75d0*b3/rs32)
      ec0ddd= -3.d0*q0rs*q4/q12+3.d0*q0rs*q3*q3/q12/q12*q12d-
     .       2.d0*q0*q3*q3*q3/q12/q12/q12*(1.d0+3.d0*q1+3.d0*q1*q1)+
     .       3.d0*q0*q3*q4/q12/q12*q12d-q0*q5/q12
C
      ec   = ec0
      ecd  = ec0d*rsd
      ecdd = ec0dd*rsd*rsd+ec0d*rsdd
      ecddd= ec0ddd*rsd*rsd*rsd+3.d0*ec0dd*rsd*rsdd+ec0d*rsddd
C
      ec   = 2.d0*ec
      ecd  = 2.d0*ecd
      ecdd = 2.d0*ecdd
      ecddd= 2.d0*ecddd
C
      RETURN
      END
