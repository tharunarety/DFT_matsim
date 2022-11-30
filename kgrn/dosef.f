      SUBROUTINE dosef(efg)
C   ******************************************************************
C   *                                                                *
C   *    Compute the density of states at the Fermi level by         *
C   *    calculating                                                 *
C   *                                                                *
C   *    NE1= -Int{n(z)[f(z)*log(f(z))+(1-f(z))*log(1-f(z))]}dz      *
C   *                                                                *
C   *    and                                                         *
C   *                                                                *
C   *    N2 = -Int{[f(z)*log(f(z))+(1-f(z))*log(1-f(z))]}dz          *
C   *                                                                *
C   * where f(z) is the Fermi function, n(z) the density of states.  *
C   * The energy integral is performed on a small elliptic contour   *
C   * from (-nres*2*pi*k_B*T,0) to (nres*2*pi*k_B*T,0), and with     *
C   * imaginary axis (pi*k_B*T/2) (i.e. no Matsubara frecvencies     *
C   * within the complex contour).                                   *
C   *                                                                *
C   * DOS(Ef) = NE1/N2                                               *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE csts ; USE dosmom ; USE energymesh
      USE greenfunc  ; USE logderivative ; USE message
      USE slope ; USE totalenergy
      IMPLICIT NONE
      INTEGER, PARAMETER :: prnt = 0
      REAL(KIND=8) :: efg, deflf
      COMPLEX(KIND=8) :: z, fz, ferm
      INTEGER      :: is, it, ita, l, lz
C
      CALL kkrarr(nzm,lmaxt,1)
      ALLOCATE(ga(nzm,nq,nlm,nlm,ns))
      ALLOCATE(dtilz(nzm,nq,nlm,nlm,ns))
      ALLOCATE(bg0(nzm,nq,nlm,ns))
C
C     1. Solve kink equation for the complex contour ZM
C
      CALL setexp(zm,nzm)
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
C     Find the k-integrated Green's function 1/(D - S) on the contour,
C     and construct the g*S, S*g*S matrices
C
      DO 20 is=1,ns
      dos=zero
      CALL dospth(is,ga,bg0,gi,dtilz,nzm)
C
      CALL doscpa(is,bg0,gi,dos,zm,nzm)
C
C     Find the number of electrons on the contour
C
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      DO 21 l=0,lmax
C
      DO 22 lz=1,nzm
      z=zm(lz)-efg*zone
      z=z/tfermi
      fz=zone/(zone+EXP(z))
      ferm=fz*LOG(fz)+(zone-fz)*LOG(zone-fz)
      dos(lz,ita,it,l)=wgm(lz)*dos(lz,ita,it,l)*ferm
   22 CONTINUE
      tdos(ita,it,l,is)=SUM(AIMAG(dos(1:nzm,ita,it,l)))
   21 CONTINUE
   20 CONTINUE
      tdos=-spinfc*tdos/pi
C
      deflf=tfermi*pi*pi/3.d0
      tdos=tdos/deflf
C
      CALL kkrarr(nzm,lmaxt,0)
      DEALLOCATE(ga,dtilz,bg0)
      RETURN
      END
