      SUBROUTINE gotora(dosreg,dosimg)
C   ******************************************************************
C   *                                                                *
C   * Find the density of states by means of Cauchy relations from   *
C   * the Green's function, determined for the complex energies.     *
C   *                                                                *
C   ******************************************************************
      USE energymesh
      USE control_data
      USE message
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(nzm) :: dosreg, dosimg
      REAL(KIND=8), DIMENSION(nzm) :: dregdx, dimgdx
      REAL(KIND=8), PARAMETER :: dely = 7.d-3, zimmin = 1.d-8
      REAL(KIND=8)  :: delx, signz, dely1, aimz, aimzp
      INTEGER       :: iexitr, lz
C
      delx=REAL(zm(2)-zm(1),8)
      signz=AIMAG(zm(1))/ABS(AIMAG(zm(1)))
      aimz=AIMAG(zm(1))-signz*dely
      dely1=dely
      iexitr=0
   20 CONTINUE
C
C     Get the derivatives
C
      CALL der4(nzm,dosreg,dregdx,delx)
      CALL der4(nzm,dosimg,dimgdx,delx)
C
C     Use the Cauchy relations
C
      DO 21 lz=1,nzm
      dosreg(lz)=dosreg(lz)+signz*dely1*dimgdx(lz)
   21 dosimg(lz)=dosimg(lz)-signz*dely1*dregdx(lz)
      aimzp=aimz
      aimz=aimz-signz*dely
      IF(ABS(aimz).GE.zimmin.AND.aimz*signz.GE.0) THEN
         GOTO 20
      ELSEIF(iexitr.EQ.0) THEN
         aimz=zimmin*signz
         dely1=ABS(aimz-aimzp)
         iexitr=1
         GOTO 20
      ELSE
         GOTO 22
      ENDIF
   22 CONTINUE
C
      RETURN
      END
