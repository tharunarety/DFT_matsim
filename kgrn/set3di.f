      SUBROUTINE set3di
C   *****************************************************************
C   *                                                               *
C   *    Set the Gauss integration mesh used to obtain the          *
C   *    kinetic force components.                                  *
C   *                                                               *
C   *****************************************************************
      USE control_data ; USE csts ; USE gaussi ; USE message
      IMPLICIT NONE
C
      nth=nfi
      IF(nth.LE.0) THEN
         WRITE(m6,101) 'NTH=',nth
         STOP
      ENDIF
      IF(nfi.LE.0) THEN
         WRITE(m6,101) 'NFI=',nfi
         STOP
      ENDIF
C
      ALLOCATE(thvec(nth),wth(nth))
      ALLOCATE(fivec(nfi),wfi(nfi))
C
      CALL wagaus(0.d0,pi,thvec,wth,nth)
      CALL wagaus(0.d0,twopi,fivec,wfi,nfi)
      WRITE(m6,100) nth,nfi
C
  100 FORMAT(/,' SET3DI:   Angular integration mesh generated with ',
     1       'NTH =',i4,' NFI =',i4)
  101 FORMAT(/,' SET3DI:**',a,i3,' must be non-zero')
      RETURN
      END
