      SUBROUTINE kkrdos(efg)
C   ******************************************************************
C   *                                                                *
C   * Calculate density of states, Green's function is first         *
C   * calculated at the line, parallel to the real axis. After that  *
C   * Cauchy relations are used in order to obtain DOS for real      *
C   * energy points.                                                 *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts ; USE energymesh ; USE greenfunc ; USE logderivative
      USE message ; USE partialwaves ; USE radialmesh ; USE slope
      USE taylor ; USE text ; USE bzmesh
      IMPLICIT NONE
      CHARACTER(LEN=7) :: spnch
      CHARACTER(LEN=80)  :: dosprn
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dosreg, dosimg, egraf
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dosit, grlold
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: tns, grl
      REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: grlm
      REAL(KIND=8)  :: efg, sz, dost, dostold, nost, stepe
      INTEGER       :: is, iq, it, ita, lz, l, m, l1, l2, lm, i6, ispace
      INTEGER :: ijob
      INTEGER, PARAMETER :: prnt = 0
C
C     Open output file
C
      i6=index(for006,' ')-4
      dosprn(1:i6+3)=for006(1:i6)//'dos' 
      DO ispace=i6+4,80
         dosprn(ispace:ispace)=' '
      END DO
      OPEN(13,FILE=dosprn,STATUS='UNKNOWN',FORM='FORMATTED')
      WRITE(m6,100) ibz,nzm
      ijob=INDEX(dosprn,' ')-1
      WRITE(m6,101) dosprn(1:ijob)
      IF(msgl.NE.0) WRITE(msgio,100) ibz,nzm
      IF(msgl.NE.0) WRITE(msgio,101) dosprn(1:ijob)
C
      spnch='UP+DOWN'
C
      CALL kkrarr(nzm,lmax,1)
      ALLOCATE(dosreg(nzm),dosimg(nzm),egraf(nzm))
      ALLOCATE(grlm(nzm,mnta,nt,0:lmax))
      ALLOCATE(dosit(nt),tns(0:lmax,0:nzm))
      ALLOCATE(grl(0:lmax,nt),grlold(0:lmax))
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
C     Find the k-integrated Green's function 1/(D - S) on the contour
C
      DO 20 is=1,ns
      grlm=zero
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
      dosreg(1:nzm)=REAL(dos(1:nzm,ita,it,l),8)
      dosimg(1:nzm)=AIMAG(dos(1:nzm,ita,it,l))
C
C     Use Cauchy relations, going to the real axis
C
      IF(ABS(AIMAG(zm(1))).GE.1.d-2)
     .    CALL gotora(dosreg,dosimg)
C
      DO 21 lz=1,nzm
      egraf(lz)=REAL(zm(lz),8)-efg
      sz=AIMAG(zm(lz))/ABS(AIMAG(zm(lz)))
      grlm(lz,ita,it,l)=grlm(lz,ita,it,l)-
     .                  spinfc*sz/pi*dosimg(lz)
   21 CONTINUE
C
C     Output
C
      IF(ns.EQ.2.AND.is.EQ.1) spnch='DOWN   '
      IF(ns.EQ.2.AND.is.EQ.2) spnch='UP     '
C
C     Total DOS and DOS(IT)
C
      stepe=egraf(2)-egraf(1)
      WRITE(13,109) spnch,(it,it=1,nt)
      WRITE(13,*)
      nost=0.d0
      dostold=0.d0
      DO 30 lz=1,nzm
      dost=0.d0
      dosit=0.d0
      DO 31 it=1,nt
      DO 31 ita=1,nta(it)
      DO 31 l=0,lmax
      dosit(it)=dosit(it)+conc(ita,it)*grlm(lz,ita,it,l)
   31 dost=dost+conc(ita,it)*mmt(it)*grlm(lz,ita,it,l)
      nost=nost+(dostold+dost)/2.d0*stepe
      dostold=dost
      WRITE(13,110) egraf(lz),dost/nq,nost/nq,dosit(1:nt)
   30 CONTINUE
C
C     Loop for IT and ITA
C
      DO 40 it=1,nt
      DO 40 ita=1,nta(it)
C
C     DOS(ITA,IT) and DOS(ITA,IT,IL)
C
      WRITE(13,120) it,ttxt(ita,it),spnch,txtl(1:nl)
      WRITE(13,*)
      grlold(0:lmax)=0.d0
      tns(0:lmax,0)=0.d0
      DO 41 lz=1,nzm
      dosit(it)=0.d0
      DO 42 l=0,lmax
      grl(l,it)=0.d0
      grl(l,it)=grl(l,it)+conc(ita,it)*grlm(lz,ita,it,l)
   43 dosit(it)=dosit(it)+conc(ita,it)*grlm(lz,ita,it,l)
      tns(l,lz)=tns(l,lz-1)+
     .          (grlold(l)+grl(l,it))/2.*stepe*mmt(it)
   42 grlold(l)=grl(l,it)
   41 WRITE(13,111) egraf(lz),dosit(it),grl(0:lmax,it)
C
C     NOS(ITA,IT),NOS(ITA,IT,IL)
C
      WRITE(13,121)
      DO 50 lz=1,nzm
   50 WRITE(13,111) egraf(lz),SUM(tns(0:lmax,lz)),tns(0:lmax,lz)
      WRITE(13,*) '   '
C
   40 CONTINUE
C
   20 CONTINUE
C
      CALL kkrarr(nzm,lmax,0)
      DEALLOCATE(dosreg,dosimg,egraf,grlm)
      DEALLOCATE(dosit,tns,grl,grlold)
      DEALLOCATE(ga,dtilz,bg0)
C
  100 FORMAT(/,' KKRDOS: DOS calculation for IBZ =',i3,
     .         ' and NZD =',i6)
  101 FORMAT(/,9x,'DOS stored in: ',a)
  109 FORMAT(/,5x,'Total DOS and NOS and partial (IT) DOS',a,
     .       //,5x,'E',5x,'Total',4(4x,i2))
  110 FORMAT(1x,f8.4,10(1x,f8.3))
  111 FORMAT(1x,f8.4,10(1x,f8.3))
  120 FORMAT(/,5x,'Sublattice',i3,' Atom ',a,
     .       ' spin ',a,//,5x,'E',5x,'Total',4(7x,a))
  121 FORMAT(/,'  TNOS ',/)
C
      RETURN
      END
