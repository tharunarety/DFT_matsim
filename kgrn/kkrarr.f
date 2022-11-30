      SUBROUTINE kkrarr(nzm,lmaxact,allc)
C   ******************************************************************
C   *                                                                *
C   * Allocate different arrays for the big itterations.             *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE greenfunc
      USE logderivative
      USE partialwaves
      USE radialmesh
      USE taylor
      IMPLICIT NONE
      INTEGER :: nzm, lmaxact, allc
C
      IF(allc.EQ.1) THEN
         ALLOCATE(cf(dimr,0:lmaxact,mnta,nt,ns,nzm))
         ALLOCATE(cfr(dimr,0:lmax,mnta,nt,ns,dimecr))
         cf=zero
         cfr=0.d0
         ALLOCATE(dfi(nzm,mnta,nt,0:lmax,ns,2))
         ALLOCATE(fi0m(nzm,mnta,nt,0:lmax,ns))
         ALLOCATE(gsng(nzm,mnta,nt,0:lmax,ns),dos(nzm,mnta,nt,0:lmax))
         ALLOCATE(gi(nzm,mnta,nq,nlm,nlm,ns))
         ALLOCATE(tayl(nzm,mder,ns),tayld(nzm,mder,ns),clkw0(nzm))
         tayl=zero
         tayld=zero
         IF(func.NE.'ASA') THEN
            ALLOCATE(cfm(0:lmax,0:lmax,0:lmax2,mnta,nt,ns,nzm))
            ALLOCATE(cfrm(0:lmax,0:lmax2,mnta,nt,ns,dimecr))
            cfm=zero
            cfrm=0.d0
         ENDIF
      ELSE
         DEALLOCATE(cf,cfr,dfi,fi0m,gsng)
         DEALLOCATE(dos,gi)
         DEALLOCATE(tayl,tayld,clkw0)
         IF(func.NE.'ASA') DEALLOCATE(cfm,cfrm)
      ENDIF
C
      RETURN
      END
