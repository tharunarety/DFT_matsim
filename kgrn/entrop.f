      SUBROUTINE entrop(efg)
C   ******************************************************************
C   *                                                                *
C   * Compute electronic entropy:                                    *
C   *                                                                *
C   *        -k_B*Int{n(z)[f(z)*log(f(z))+(1-f(z))*log(1-f(z))]}dz   *
C   *                                                                *
C   * where f(z) is the Fermi function, n(z) the density of states.  *
C   * The energy integral is performed on a small elliptic contour   *
C   * from (-nres*2*pi*k_B*T,0) to (nres*2*pi*k_B*T,0), and with     *
C   * imaginary axis (pi*k_B*T/2) (i.e. no Matsubara frecvencies     *
C   * within the complex contour).                                   *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE csts ; USE dosmom ; USE energymesh
      USE greenfunc  ; USE logderivative ; USE message
      USE slope ; USE totalenergy
      IMPLICIT NONE
      INTEGER, PARAMETER :: prnt = 0
      REAL(KIND=8) :: efg, tse, deflf
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
C
      entr(ita,it,l,is)=SUM(AIMAG(dos(1:nzm,ita,it,l)))
   21 CONTINUE
   20 CONTINUE
      entr=-spinfc*entr/pi
C
C     Get T*S
C
      entr=tfermi*entr
C
      deflf=tfermi*pi*pi/(3.d0*LOG(2.d0))
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      tse=SUM(tdos(ita,it,0:lmax,1:ns))*deflf
      ents(ita,it)=ents(ita,it)+SUM(entr(ita,it,0:lmax,1:ns))
      WRITE(m6,101) it,ita,ents(ita,it),tse*tfermi
   30 CONTINUE
C
      CALL kkrarr(nzm,lmaxt,0)
      DEALLOCATE(ga,dtilz,bg0)
  101 FORMAT(/,' ENTROP: IT =',i3,' ITA =',i3,' T*S =',f10.6,
     .         ' T*S(~) =',f10.6)
      RETURN
      END
