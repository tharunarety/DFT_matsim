      SUBROUTINE bessinit(lmin,lmax)
C   ******************************************************************
C   *                                                                *
C   *    Initialize factorials used in BESSL and BESSLZ              *
C   *                                                                *
C   ******************************************************************
      USE bessl_fact
      IMPLICIT NONE
      INTEGER :: lmax, lmin, lmx, l
C
      lmx=MAX0(lmax,-lmin,2)
      IF(lmin.GT.lmax) THEN
         DEALLOCATE(facbl)
         RETURN
      ENDIF
      ALLOCATE(facbl(-lmx:lmx+1))
      lmaxi=lmx
      lmini=lmin
C
C     Construct table of double factorials facbl(l)=(2l-1)!!
C
      facbl(0)=1.d0
      DO l=1,lmx+1
         facbl(l)=facbl(l-1)*(l+l-1)
      ENDDO
      DO l=-1,lmin,-1
         facbl(l)=facbl(l+1)/(l+l+1)
      ENDDO
C
      RETURN
      END
