      FUNCTION fty(m,mp)
C   ******************************************************************
C   *                                                                *
C   *    Set the matrix used to transform from the spherical to the  *
C   *    real harmonic representation.                               *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      COMPLEX(KIND=8) :: fty, fty0
      REAL(KIND=8), PARAMETER :: r2=0.7071067812d0
      INTEGER :: m, mp, id
C
      id=ABS(m)-ABS(mp)
      IF(id.EQ.0) THEN
        IF(m.LT.0.AND.mp.LT.0) THEN
          fty0=CMPLX(0.d0,1.d0,8)
        ELSEIF(m.LT.0.AND.mp.GT.0) THEN
          fty0=((-1.d0)**m)*CMPLX(0.d0,-1.d0,8)
        ELSEIF(m.GT.0.AND.mp.LT.0) THEN
          fty0=CMPLX(1.d0,0.d0,8)
        ELSEIF(m.GT.0.AND.mp.GT.0) THEN
          fty0=((-1.d0)**m)*CMPLX(1.d0,0.d0,8)
        ENDIF
        IF(IABS(m).GT.0) THEN
          fty=fty0*r2
        ELSE
          fty=CMPLX(1.d0,0.d0,8)
        ENDIF
      ELSE
        fty=(0.d0,0.d0)
      ENDIF
C
      RETURN
      END
