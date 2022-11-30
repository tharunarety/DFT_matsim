      SUBROUTINE slfarr(all)
C   ******************************************************************
C   *                                                                *
C   * Allocate and deallocate arrays for the self-consistent loops.  *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE energymesh
      USE greenfunc
      USE logderivative
      USE slope
      IMPLICIT NONE
      INTEGER :: all
C
      IF(all.EQ.1) THEN
         ALLOCATE(hgh(nzm,nq,ns,lmax+1:lmaxt))
         ALLOCATE(hghx(nx,nq,ns,lmax+1:lmaxt))
         ALLOCATE(slop(nlmqt,nlmqt),tmp(nlmqt,nlmq))
         slop=zero
         tmp=zero
         ALLOCATE(ga(nzm,nq,nlm,nlm,ns),gax(nx,nq,nlm,nlm,ns))
         ALLOCATE(dtilz(nzm,nq,nlm,nlm,ns),dtilx(nx,nq,nlm,nlm,ns))
         ALLOCATE(bg0(nzm,nq,nlm,ns),bg0x(nx,nq,nlm,ns))
      ELSEIF(all.EQ.0) THEN
         DEALLOCATE(ga,gax,bg0,bg0x,dtilz,dtilx)
         DEALLOCATE(hgh,hghx,slop,tmp)
      ENDIF
C
      RETURN
      END
