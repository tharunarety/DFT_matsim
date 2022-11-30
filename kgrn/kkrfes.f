      SUBROUTINE kkrfes(ef)
C   ******************************************************************
C   *                                                                *
C   * Calculate Fermi surface.                                       *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts ; USE energymesh ; USE greenfunc ; USE logderivative
      USE message ; USE partialwaves ; USE radialmesh ; USE slope
      USE taylor ; USE text ; USE bzmesh ; USE dosmom
      IMPLICIT NONE
      CHARACTER(LEN=80)  :: fermiprn
      REAL(KIND=8) :: ef, dn
      INTEGER, PARAMETER :: prnt = 0
      INTEGER      :: is, lk, it, ita, l, i6, ispace
C
C     Set up k-mesh for the cross-section
C
      CALL fermis
C
C     Open output file
C
      IF(fsts.EQ.0) THEN
         i6=index(for006,' ')-4
         fermiprn(1:i6+3)=for006(1:i6)//'fes' 
         DO ispace=i6+4,80
            fermiprn(ispace:ispace)=' '
         END DO
         OPEN(12,FILE=fermiprn,STATUS='UNKNOWN',FORM='FORMATTED')
         WRITE(m6,100) ibz,nfs
         WRITE(m6,101) fermiprn
         IF(msgl.NE.0) WRITE(msgio,100) ibz,nfs
         IF(msgl.NE.0) WRITE(msgio,101) fermiprn
      ENDIF
C
      CALL kkrarr(nzm,lmax,1)
      ALLOCATE(ga(nzm,nq,nlm,nlm,ns))
      ALLOCATE(dtilz(nzm,nq,nlm,nlm,ns))
      ALLOCATE(bg0(nzm,nq,nlm,ns))
C
C     1. Solve kink equation for the complex energies ZM
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
      CALL singlg(ef,prnt)
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
C     Find the k-resolved Green's function 1/(D - S) on the contour
C
      DO 30 is=1,ns
      CALL fespth(is,bg0,gi,dtilz,zm,wgm,nzm)
C
      IF(fsts.EQ.1) THEN
         dn=0.d0
         CALL doscpa(is,bg0,gi,dos,zm,nzm)
         DO 31 it=1,nt
         DO 32 ita=1,nta(it)
C
         DO 33 l=0,lmax
         dos(1:nzm,ita,it,l)=wgm(1:nzm)*dos(1:nzm,ita,it,l)
         tdos(ita,it,l,is)=SUM(AIMAG(dos(1:nzm,ita,it,l)))
         dn=dn+conc(ita,it)*mmt(it)*
     .      SUM(AIMAG(dos(1:nzm,ita,it,l)))
   33    CONTINUE
   32    CONTINUE
   31    CONTINUE
         dn=spinfc*dn/pi
         WRITE(m6,110) dn
         IF(msgl.NE.0) WRITE(msgio,110) dn
      ENDIF
   30 CONTINUE
      tdos=spinfc*tdos/pi
      IF(prnt.EQ.1) CALL prinfo('KKRFES')
C
      CALL kkrarr(nzm,lmax,0)
      DEALLOCATE(ga,bg0,dtilz)
      IF(fsts.EQ.0) CLOSE(12)
C
  100 FORMAT(/,' KKRFES: Fermi surface calculation for IBZ =',i3,
     .         ' and NKVEC =',i6)
  101 FORMAT(/,9x,'Fermi surface stored in: ',a)
  110 FORMAT(/,9x,'DOS(EF) =',f10.6)
      RETURN
      END
