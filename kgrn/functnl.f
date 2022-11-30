      SUBROUTINE functnl
C   ******************************************************************
C   *                                                                *
C   *                                                                *
C   ******************************************************************
      USE bzmesh
      USE control_data
      USE control_text
      USE moments
      USE radialmesh
      USE realgaunt
      IMPLICIT NONE
      INTEGER :: maxlmax, lmaxhp
C
      IF(func.NE.'ASA') THEN
         ngaunt=0
         CALL mgaunt(lmax,lmax,lmax2)
         diml=nlm2
      ELSE
         diml=1
      ENDIF
      CALL rotfact
      maxlmax=lmax2
      IF(fcd.EQ.'Y') THEN
         maxlmax=MAX0(lmax2,lmaxh)
         CALL set3di
         lmaxhp=lmaxh+1
         CALL realhinit(lmaxhp)
      ENDIF
      CALL rotm3d(ibz,maxlmax)
      CALL setsym(ibz)
C
C     Set the Madelung matrix
C
      CALL madmtr
C
      RETURN
      END
