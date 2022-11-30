      SUBROUTINE arrayz(all)
C   ******************************************************************
C   *                                                                *
C   * Allocate different arrays.                                     *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE energymesh
      USE slope
      IMPLICIT NONE
      INTEGER :: all
C
      IF(all.EQ.1) THEN
         ALLOCATE(taz(nzm,4,0:lmax,nt,ns))
         ALLOCATE(tsz(nzm,4,0:lmax,mnta,nt,ns,0:1))
         ALLOCATE(tmz(nzm,4,0:lmax,mnta,nt,ns))
         ALLOCATE(tax(nx,4,0:lmax,nt,ns))
         ALLOCATE(tsx(nx,4,0:lmax,mnta,nt,ns,0:1))
         ALLOCATE(tmx(nx,4,0:lmax,mnta,nt,ns))
         taz=zero
         tsz=zero
         tmz=zero
         tax=zero
         tsx=zero
         tmx=zero
      ELSE
         DEALLOCATE(taz,tsz,tmz,tax,tsx,tmx)
      ENDIF
C
      RETURN
      END
