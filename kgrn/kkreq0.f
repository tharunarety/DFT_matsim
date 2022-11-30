      SUBROUTINE kkreq0(efg,nzlin,lin,contl)
C   ******************************************************************
C   *                                                                *
C   *   Find an approximate solution by means of the Dyson equation  *
C   *   for the k-integrated Green's function.                       *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts ; USE density ; USE dosmom ; USE energymesh
      USE greenfunc ; USE logderivative ; USE message
      USE partialwaves ; USE radialmesh ; USE slope ; USE totalenergy
      IMPLICIT NONE
      INTEGER, PARAMETER :: prnt = 0
      REAL(KIND=8) :: efg, dn
      INTEGER :: is, it, ita, ires, nzlin, lin, l, contl
C
C     1. Solve kink equation for the complex contour ZM
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
      CALL logder0(taz,tsz,tmz,zm,nzm)
C
C     Find the k-integrated Green's function by solving the Dyson equation
C
      dn=0.d0
      DO 20 is=1,ns
      dos=zero
      CALL pathop0(is,ga,gi,bg0,bgdos,dtilz,nzm,nzlin)
C
      CALL doscpa0(is,ga,gi,bgdos,dos,zm,nzm,nzlin)
C
C     Find the number of electrons on the contour
C
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
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
      IF(prnt.EQ.1) CALL prinfo('KKREQ0')
C
C     Calculate the imaginar part of the r-diagonal Green's function 
C
      chdn=0.d0
      IF(func.NE.'ASA'.AND.lin.EQ.1) grnfm=0.d0
      CALL greenf(lin,1,nzm,gi,hgh,zm,wgm,nzm,nzlin)
C
C     2. Solve kink equation for the energies ZX close to the Fermi level
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
      CALL logder0(tax,tsx,tmx,zx,nx)
C
C     Find the k-integrated Green's function by solving the Dyson equation
C
      tnosx=zero
C
      DO 31 is=1,ns
      dos=zero
      CALL pathop0(is,gax,gi,bg0x,bgdos,dtilx,nx,nzlin)
C
      CALL doscpa0(is,gax,gi,bgdos,dos,zx,nx,nzlin)
C
      DO 31 it=1,nt
      DO 31 ita=1,nta(it)
      DO 31 l=0,lmax
   31 tnosx(1:nx,ita,it,l,is)=dos(1:nx,ita,it,l)
      tnosx=spinfc*tnosx/pi
C
C     Find the number of electrons close to EF
C
      ires=0
      CALL efxmom(lin,gi,hghx,tnosx,nzlin,efg,dn,ires)
      IF(ires.EQ.1) THEN
         WRITE(m6,'(/,a,i3)') ' EFXMOM called from KKREQ0, LIN=',lin
         CLOSE(11,STATUS='DELETE')
         STOP
      ENDIF
      contl=ires
C
C     Add the higher tails contribution to the density
C
      IF(lmaxt.GT.lmax) chdn=chdn+chdh
C
      RETURN
      END
