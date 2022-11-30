      SUBROUTINE setibz
C   ******************************************************************
C   *                                                                *
C   * Point group:                IBZR(irot,1) = 1.                  *
C   *                                                                *
C   * Improper symmetry elements: IBZR(irot,1) = 0;                  *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE message
      USE symmetry
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: nugam
      REAL(KIND=8), PARAMETER :: err = 1.d-8
      REAL(KIND=8)  :: smm, diffg
      INTEGER, DIMENSION(48) :: ibzr0
      INTEGER  :: irot, jrot, krot, eq, l, m, mp, mpp
C
      ALLOCATE(imsymm(48),nugam(0:lmax,-lmax:lmax,-lmax:lmax))
C
      imsymm(1:48)=0
      imsymm(1)=1           !! Identity
      ibzr0(1:48)=iwrot(1:48)
C
      DO 20 irot=2,nrot
      IF(ibzr0(irot).EQ.0.AND.ibzrot(irot).EQ.1) THEN
         imsymm(irot)=1
         DO 21 jrot=1,nrot
         IF(iwrot(jrot).EQ.1) THEN
            DO 22 l=0,lmax
            DO 22 m=-l,l
            DO 22 mp=-l,l
            smm=0.d0
            DO 23 mpp=-l,l
   23       smm=smm+ugam(jrot,l,m,mpp)*ugam(irot,l,mpp,mp) 
   22       nugam(l,m,mp)=smm
            DO 24 krot=1,nrot
            IF(ibzr0(krot).EQ.0.AND.ibzrot(krot).EQ.1) THEN
               eq=1
               DO 25 l=0,lmax
               DO 25 m=-l,l
               DO 25 mp=-l,l
               diffg=ABS(ugam(krot,l,m,mp)-nugam(l,m,mp))
               IF(diffg.GT.err) eq=0
   25          CONTINUE
               IF(eq.EQ.1) ibzr0(krot)=1
            ENDIF
   24       CONTINUE
         ENDIF
   21    CONTINUE
      ENDIF
   20 CONTINUE
      WRITE(m6,100) SUM(imsymm(1:48)),SUM(iwrot(1:48)),
     .              SUM(ibzrot(1:48))
C
      DEALLOCATE(nugam)
  100 FORMAT(/,' SETIBZ:   dimension of the point group:',/,
     .          11x,'IMSYMM: ',i3,' IBZR: ',i3,' IBZROT: ',i3)
      RETURN
      END
