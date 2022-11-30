      SUBROUTINE screen
C   ******************************************************************
C   *                                                                *
C   *    Read the hard sphere radii a(L,R).                          *
C   *    Set up the transformation matrix and energy derivatives.    *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE factorial
      USE control_data
      USE control_text
      USE lattice
      USE screening
      USE message
      IMPLICIT NONE
      CHARACTER(LEN=1), DIMENSION(6) :: txtl=(/"s","p","d","f","g","h"/)
      INTEGER :: iq, l, i, j
      REAL(KIND=8) :: amx, a, ahrd
C
      ALLOCATE(sigma(0:lmax,nq))
      ALLOCATE(tmat(4,0:lmax,nq,0:nder))
      ALLOCATE(bigd(2,0:lmax,nq,0:nder))
      ALLOCATE(cd(2,0:lmax,nq,0:nder))
      ALLOCATE(dfac(0:lmax,nq))
      ALLOCATE(alwats(0:lmaxw,nq,0:nder))
      ALLOCATE(ifib(0:nder,0:nder))
C
      WRITE(m6,100) txtl(1:nl)
      WRITE(m6,*)
C
C     Read the hard sphere radii in units of atomic sphere radii
C
      DO 10 iq=1,nq
      READ(5,'(10x,6f5.2)') sigma(0:lmax,iq)
      WRITE(m6,101) iq,sigma(0:lmax,iq)
   10 CONTINUE
C
C     Test for the hard sphere radii
C
      DO 20 iq=1,nq
      DO 20 l=0,lmax
      ahrd=sigma(l,iq)*wst(iq)*ws 
      IF(ahrd.GT.wi(iq)) THEN
         WRITE(m6,102) iq,l,ahrd,wi(iq)
         STOP
      ENDIF
   20 CONTINUE
C
      ifib(0:nder,0:nder)=0
      ifib(0,0)=1
      DO 21 j=1,nder
      ifib(j,0)=1
      DO 21 i=1,j
   21 ifib(j,i)=ifib(j-1,i-1)+ifib(j-1,i)
      CALL trmtrx
C
      IF(wats.NE.'Y') RETURN
      amx=0.d0
      DO 30 iq=1,nq
      a=sigma(0,iq)*wst(iq)
      IF(amx.LT.a) amx=a
   30 CONTINUE
      DO 31 iq=1,nq
   31 rwats(iq)=rwats(iq)/ws+amx+dwats
      CALL trwats
C
  100 FORMAT(/,' SCREEN:   Hard sphere radii in units of ws:',
     .       //,11x,'site',9x,5(a,7x))
  101 FORMAT(11x,i3,3x,6f8.4)
  102 FORMAT(/,' SCREEN:   Hard sphere radius a(',i2,',',i2,') =',
     .           f9.6,' is bigger than Si =',f9.6)
      RETURN
      END
