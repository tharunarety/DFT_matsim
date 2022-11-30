      SUBROUTINE efxmom(lin,gi,hghx,tnosx,nzlin,efg,dn,ires)
C   ******************************************************************
C   *                                                                *
C   *   Find the contribution to the state density                   *
C   *   by interpolation close to the Fermi level.                   *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE dosmom ; USE energymesh
      USE message ; USE partialwaves ; USE totalenergy
      IMPLICIT NONE
      INTEGER      :: lin, nzlin
      INTEGER      :: ires, lx, nfx, nxx1, nxx2, iq, it, ita, is
      INTEGER      :: lm, l1, l2, l, i
      COMPLEX(KIND=8), DIMENSION(nzlin,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(nx,nq,ns,lmax+1:lmaxt) :: hghx
      COMPLEX(KIND=8), DIMENSION(nx,mnta,nt,0:lmax,ns)  :: tnosx
      COMPLEX(KIND=8), DIMENSION(nx) :: wxx
      REAL(KIND=8), PARAMETER :: tolerr = 0.1d0
      REAL(KIND=8),    DIMENSION(nx) :: nel, deln
      REAL(KIND=8) :: efg, dn, dhx, elnx0, x0, dfx, mdos, smm
      REAL(KIND=8) :: ex, e, weight, nex, err
C
C     Take only 50% from the new Fermi level
C
      dhx=hx
C
      nel(1)=0.d0
      DO 10 lx=2,nx
      smm=0.d0
      DO 11 is=1,ns
      DO 11 l=0,lmax
      DO 11 it=1,nt
      DO 11 ita=1,nta(it)
      nex=AIMAG(tnosx(lx-1,ita,it,l,is))+AIMAG(tnosx(lx,ita,it,l,is))
   11 smm=smm+conc(ita,it)*mmt(it)*nex
   10 nel(lx)=nel(lx-1)+smm
      nel=dhx*nel
C
      elnx0=nel(nx0)-(dn+elt)
      nel=nel-elnx0
      deln=nel-elt
      DO 21 lx=2,nx
      IF(deln(lx-1)*deln(lx).LT.0.d0) THEN
         mdos=SIGN(1.d0,deln(lx)-deln(lx-1))
         nfx=lx-1
         IF(mdos.GT.0.d0) GO TO 30
      ENDIF
   21 CONTINUE
      WRITE(m6,101) efg,hx
      IF(lin.EQ.0) THEN
         WRITE(m6,102)
         DO 22 lx=1,nx
   22    WRITE(m6,103) REAL(zx(lx),8),nel(lx)
      ENDIF
C
C  a) Fermi level is not in the interval from HX(1) to HX(NX)
C     Set up IRES according to the number of electrons on the contour
C     IRES = 1 stop the code (restart with another quess for Ef);
C     IRES = 2 Ef should be somewhat bellow hx(1); 
C     IRES = 3 Ef should be somewhat above hx(nx); 
C             if called from KKREQ : set up new contour and continue iterations 
C             if called from KKREQ0: stop Dyson loops and continue iterations 
C
      ires=1
      IF(nel(1).GT.elt) THEN
         err=(nel(1)-elt)/elt
         IF(err.LT.tolerr) ires=2
         IF(lin.EQ.0) efg=REAL(zx(1),8)
      ELSEIF(nel(nx).LT.elt) THEN
         err=(elt-nel(nx))/elt
         IF(err.LT.tolerr) ires=3
         IF(lin.EQ.0) efg=REAL(zx(nx),8)
      ENDIF
      RETURN
C
C  b) Fermi level is in the interval from HX(1) to HX(NX)
C     Interpolate to Fermi level
C
   30 x0=REAL(zx(nfx),8)
      efg=x0-hx*deln(nfx)/(deln(nfx+1)-deln(nfx))
      IF(dn.LT.0.d0) THEN
         wxx=CMPLX(0.d0,0.d0,8)
         DO 31 lx=nx0+1,nfx
         wxx(lx-1)=wxx(lx-1)+CMPLX(dhx,0.d0,8)
   31    wxx(lx)=wxx(lx)+CMPLX(dhx,0.d0,8)
         dfx=dhx*(efg-x0)/hx
         wxx(nfx)=wxx(nfx)+CMPLX(dfx,0.d0,8)
         wxx(nfx+1)=CMPLX(dfx,0.d0,8)
         nxx1=nx0
         nxx2=nfx+1
      ELSE
         wxx=CMPLX(0.d0,0.d0,8)
         DO 32 lx=nfx+2,nx0
         wxx(lx-1)=wxx(lx-1)+CMPLX(-dhx,0.d0,8)
   32    wxx(lx)=wxx(lx)+CMPLX(-dhx,0.d0,8)
         dfx=dhx*(x0+hx-efg)/hx
         wxx(nfx)=CMPLX(-dfx,0.d0,8)
         wxx(nfx+1)=wxx(nfx+1)+CMPLX(-dfx,0.d0,8)
         nxx1=nfx
         nxx2=nx0
      ENDIF
C
      nfcd=nxx2-nxx1+1
      zfcd(1:nfcd)=zx(nxx1:nxx2)
      wfcd(1:nfcd)=wxx(nxx1:nxx2)
C
      dn=0.d0
      DO 40 is=1,ns
      DO 40 l=0,lmax
      DO 40 it=1,nt
      DO 40 ita=1,nta(it)
C
      tnosx(nxx1:nxx2,ita,it,l,is)=wxx(nxx1:nxx2)*
     .                          tnosx(nxx1:nxx2,ita,it,l,is)
      tnos(ita,it,l,is)=tnos(ita,it,l,is)+
     .              SUM(AIMAG(tnosx(nxx1:nxx2,ita,it,l,is)))
      tnosx(nxx1:nxx2,ita,it,l,is)=zx(nxx1:nxx2)*
     .                          tnosx(nxx1:nxx2,ita,it,l,is)
      emom(ita,it,l,is)=emom(ita,it,l,is)+
     .              SUM(AIMAG(tnosx(nxx1:nxx2,ita,it,l,is)))
      dn=dn+conc(ita,it)*mmt(it)*tnos(ita,it,l,is)
   40 CONTINUE
C
      CALL greenf(lin,nxx1,nxx2,gi,hghx,zx,wxx,nx,nzlin)
C
C     Check precision on Fermi level determination.
C
      IF(ABS(dn-elt).GT.1.d-6) THEN
         WRITE(m6,101) efg,hx
         WRITE(m6,104) elt,dn
         WRITE(m6,105) x0,nfx,nx0,nxx1,nxx2
         WRITE(m6,106)
         DO 50 lx=1,nx
   50    WRITE(m6,107) REAL(zx(lx),8),REAL(wxx(lx),8),nel(lx),deln(lx)
         STOP
      ENDIF
C
      IF(ns.EQ.2) THEN
         tmag=0.d0
         amag=0.d0
         DO 60 it=1,nt
         DO 60 ita=1,nta(it)
         DO 61 l=0,lmax
   61    amag(ita,it)=amag(ita,it)+tnos(ita,it,l,2)-tnos(ita,it,l,1)
         tmag=tmag+conc(ita,it)*mmt(it)*amag(ita,it)
   60    CONTINUE
      ENDIF
C
  101 FORMAT(/,' EFXMOM:** Fermi level not found. EFGS =',f10.6,
     1       ' HX =',f8.4)
  102 FORMAT(/,19x,'ZX          ELN',/)
  103 FORMAT(11x,10f12.4)
  104 FORMAT(/,11x,'ELT=',f10.6,' NOS(Ef)=',f10.6)
  105 FORMAT(/,11x,f10.6,4i3)
  106 FORMAT(/,11x,'    ZX        WXX        ELN        DELN',/)
  107 FORMAT(/,11x,4f10.6)
      RETURN
      END
