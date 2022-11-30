      SUBROUTINE alltmp(nrlmw,inf)
C   ******************************************************************
C   *                                                                *
C   *    Allocate temporary arrays.                                  *
C   *                                                                *
C   ******************************************************************
      USE control_text
      USE control_data
      USE temporary
      IMPLICIT NONE
      INTEGER :: nrlmw, inf
C
      IF(inf.EQ.1) THEN
         ALLOCATE(sb(nrlmw,nrlmw),sbd(nrlmw,nrlmw,nder))
         ALLOCATE(work(nrlmw),ipvt(nrlmw))
         ALLOCATE(dett(nrlmw),sc(0:nder))
         IF(wats.EQ.'Y') THEN
            ALLOCATE(sbare(nlm,nlmw,0:nder))
         ELSE
            ALLOCATE(sbare(nlm,nlm,0:nder))
         ENDIF
      ELSEIF(inf.EQ.2) THEN
         DEALLOCATE(sb,sbd)
         DEALLOCATE(work,ipvt)
         DEALLOCATE(dett,sc)
         DEALLOCATE(sbare)
      ENDIF
C
      RETURN
      END
