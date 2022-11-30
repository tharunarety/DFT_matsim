      SUBROUTINE linarr(all,lmaxact,nzlin)
C   ******************************************************************
C   *                                                                *
C   * Allocate and deallocate arrays for the dyson loops.            *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE energymesh
      USE greenfunc
      USE logderivative
      USE partialwaves
      USE radialmesh
      IMPLICIT NONE
      INTEGER :: all, nzlin, lmaxact
C
      IF(all.EQ.1) THEN
         ALLOCATE(gd(nzlin,nq,nlm,nlm,ns),bgdos(nzlin,nq,nlm))
         gd=zero
         ALLOCATE(cf(dimr,0:lmaxact,mnta,nt,ns,nzlin))
         ALLOCATE(cfr(dimr,0:lmax,mnta,nt,ns,dimecr))
         cf=zero
         cfr=0.d0
         ALLOCATE(dfi(nzlin,mnta,nt,0:lmax,ns,2))
         ALLOCATE(fi0m(nzlin,mnta,nt,0:lmax,ns))
         dfi=zero
         ALLOCATE(gsng(nzlin,mnta,nt,0:lmax,ns))
         ALLOCATE(gi(nzlin,mnta,nq,nlm,nlm,ns))
         ALLOCATE(dos(nzlin,mnta,nt,0:lmax),tnosx(nx,mnta,nt,0:lmax,ns))
         tnosx=zero
         IF(func.NE.'ASA') THEN
            ALLOCATE(cfm(0:lmax,0:lmax,0:lmax2,mnta,nt,ns,nzlin))
            ALLOCATE(cfrm(0:lmax,0:lmax2,mnta,nt,ns,dimecr))
            cfm=zero
            cfrm=0.d0
         ENDIF
      ELSEIF(all.EQ.0) THEN
         DEALLOCATE(gd,bgdos,cf,cfr,dfi,fi0m,gsng,gi,dos,tnosx)
         IF(func.NE.'ASA') DEALLOCATE(cfm,cfrm)
      ENDIF
C
      RETURN
      END
