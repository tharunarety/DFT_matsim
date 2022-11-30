      SUBROUTINE kkreq(efg)
C   ******************************************************************
C   *                                                                *
C   *    Find the Fermi level by solving the kink equation :         *
C   *                                                                *
C   *      K*v = a*[D{fi0^a(a)} - S^a]*v = 0                         *
C   *                                                                *
C   *    and construct the state density and the new charge density. *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts ; USE density ; USE dosmom ; USE energymesh
      USE greenfunc ; USE logderivative ; USE message
      USE slope
      IMPLICIT NONE
      INTEGER, PARAMETER :: prnt = 0, ntry = 3
      REAL(KIND=8) :: efg, dn
      INTEGER :: is, it, ita, ires, itry, l
C
      CALL kkrarr(nzm,lmaxt,1)
C
C     For incorrect initial guess for the Fermi level
C     the KKR eq. is solved on a new contour
C
      itry=0
   10 CONTINUE
C
C     1. Solve kink equation for the complex contour ZM
C
      CALL setexp(zm,nzm,nz2)
      CALL expans(zm,nzm)
C
C     Solve the Dirac equation and make the logarithmic derivatives
C
      CALL logder(zm,nzm)
C
C     Set up the logarithmic derivative at s^m if requisted
C
      IF(lclmff.NE.0) CALL lclmtz(tmz,zm,nzm)
C
C     Integrate the single-site Green's function on the real axis
C
      CALL singlg(efg,prnt)
C
C     Make the logarithmic derivative of the backward extrapolated
C     free electron solution : fi^a = f^a + g^a * D{fi0^a(a)}/(-a*d)
C
      CALL screen(taz,tsz,zm,nzm)
      CALL logder0(taz,tsz,tmz,zm,nzm)
C
C     Initial guess for the coherent potential function
C
      CALL logdcpa(dtilz,nzm)
C
C     Make the higher basis functions
C
      CALL logderh(zm,nzm,lmaxt)
C
C     Find the k-integrated Green's function 1/(D - S) on the contour
C
      dn=0.d0
      DO 20 is=1,ns
      CALL pathop(is,ga,bg0,gi,hgh,dtilz,nzm)
C
      CALL doscpa(is,bg0,gi,dos,zm,nzm)
C
C     Find the number of electrons on the contour
C
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
C
      DO 22 l=0,lmax
   22 dos(1:nzm,ita,it,l)=wgm(1:nzm)*dos(1:nzm,ita,it,l)
      DO 23 l=0,lmax
      IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
         tnos0(ita,it,l,is)=SUM(AIMAG(dos(1:nz2,ita,it,l)))
      ENDIF
   23 tnos(ita,it,l,is)=SUM(AIMAG(dos(1:nzm,ita,it,l)))
C
      DO 24 l=0,lmax
   24 dos(1:nzm,ita,it,l)=zm(1:nzm)*dos(1:nzm,ita,it,l)
      DO 25 l=0,lmax
      IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
         emom0(ita,it,l,is)=SUM(AIMAG(dos(1:nz2,ita,it,l)))
      ENDIF
   25 emom(ita,it,l,is)=SUM(AIMAG(dos(1:nzm,ita,it,l)))
      dn=dn+conc(ita,it)*mmt(it)*SUM(tnos(ita,it,0:lmax,is))
   21 CONTINUE
   20 CONTINUE
      tnos=spinfc*tnos/pi
      emom=spinfc*emom/pi
      dn=spinfc*dn/pi-elt
      IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
         tnos0=spinfc*tnos0/pi
         emom0=spinfc*emom0/pi
      ENDIF
      IF(prnt.EQ.1) CALL prinfo('KKREQ ')
      IF(crt.EQ.'I') THEN
         IF(msgl.EQ.1) WRITE(msgio,110) dn+elt,elt,efg
         WRITE(m6,110) dn+elt,elt,efg
      ENDIF
C
C     Calculate the imaginar part of the r-diagonal Green's function 
C
      chdn=0.d0
      IF(func.NE.'ASA') grnfm=0.d0
      chdh=0.d0
      CALL greenf(0,1,nzm,gi,hgh,zm,wgm,nzm,nzm)
C
      CALL kkrarr(nzm,lmaxt,0)
C
      CALL kkrarr(nx,lmaxt,1)
      ALLOCATE(tnosx(nx,mnta,nt,0:lmax,ns))
C
C     2. Solve kink equation for the energies ZX close to the Fermi level
C
      CALL setexp(zx,nx,0)
      CALL expans(zx,nx)
C
C     Solve the Dirac equation and make the logarithmic derivatives
C
      CALL logder(zx,nx)
C
C     Set up the logarithmic derivative at s^m if requisted
C
      IF(lclmff.NE.0) CALL lclmtz(tmx,zx,nx)
C
C     Make the logarithmic derivative of the backward extrapolated
C     free electron solution : fi^a = f^a + g^a * D{fi0^a(a)}/(-a*d)
C
      CALL screen(tax,tsx,zx,nx)
      CALL logder0(tax,tsx,tmx,zx,nx)
C
C     Initial guess for the coherent potential function
C
      CALL logdcpa(dtilx,nx)
C
C     Make the higher basis functions
C
      CALL logderh(zx,nx,lmaxt)
C
C     Find the k-integrated Green's function 1/(D - S)
C
      tnosx=zero
      DO 30 is=1,ns
      dos=zero
      CALL pathop(is,gax,bg0x,gi,hghx,dtilx,nx)
C
      CALL doscpa(is,bg0x,gi,dos,zx,nx)
C
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      DO 30 l=0,lmax
   30 tnosx(1:nx,ita,it,l,is)=dos(1:nx,ita,it,l)
      tnosx=spinfc*tnosx/pi
C
C     Find the number of electrons close to EF
C
      ires=0
      CALL efxmom(0,gi,hghx,tnosx,nx,efg,dn,ires)
      IF(ires.EQ.1) THEN
         WRITE(m6,'(/,a)') ' EFXMOM called from KKREQ'
         CLOSE(11,STATUS='DELETE')
         STOP
      ELSEIF(ires.EQ.2.OR.ires.EQ.3) THEN
         itry=itry+1
         DEALLOCATE(zm,wgm,zx)
         CALL zmesh(efg,1)
         DEALLOCATE(tnosx)
         CALL kkrarr(nx,lmaxt,0)
         CALL kkrarr(nzm,lmaxt,1)
         IF(itry.LT.ntry) GO TO 10
         WRITE(m6,100) ntry
         STOP
      ENDIF
C
C     Add the higher tails contribution to the density
C
      IF(lmaxt.GT.lmax) chdn=chdn+chdh
C
      CALL kkrarr(nx,lmaxt,0)
      DEALLOCATE(tnosx)
C
  100 FORMAT(' KKREQ:** ITRY exceeded NTRY =',i3)
  110 FORMAT(/,' KKREQ:  NOS(Ef) =',f10.6,' ELT =',f8.4,' Ef =',f8.4)
      RETURN
      END
