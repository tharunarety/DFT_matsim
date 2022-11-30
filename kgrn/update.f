      SUBROUTINE update(ef)
C   ******************************************************************
C   *                                                                *
C   *    Save potential parameters, renormalised charge densities,   *
C   *    total charge densities, potentials, and total energies      *
C   *    on FOR002                                                   *
C   *                                                                *
C   ******************************************************************
      USE atomicdens   ; USE botomtop     ; USE bzmesh
      USE control_data ; USE control_text ; USE density
      USE energymesh   ; USE message      ; USE moments
      USE pota         ; USE potential    ; USE potparam
      USE radialmesh   ; USE softcore     ; USE text
      USE totalenergy
      IMPLICIT NONE
      REAL(KIND=8) :: ef
      INTEGER, PARAMETER :: mo=30
      INTEGER :: iq, it, ita, is, ip, l, ixch, jri, jrn
C
C     Save potential parameters
C
      IF(iter.GT.1) THEN
         OPEN(2,FILE=for002,STATUS='UNKNOWN',FORM='UNFORMATTED')
         OPEN(9,FILE=for009,STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDIF
      WRITE(2) 'KGRN',nq,nl,conv,func,txch,afm
      WRITE(2) sws,wst(1:nq),wsi(1:nq),wsc(1:nq),dx
      WRITE(2) eb,ef,ns,nt,itq(1:nq),nta(1:nt)
      WRITE(2) ebt,etp,0.d0,dexch
      ixch=ixc-1
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
   20 WRITE(2) ws(ita,it),hsr(ita,it),r1(ita,it),jwss(ita,it),
     .         jris(ita,it),ttxt(ita,it)
C
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      txtp(ita,it)='Potential parameters from KGRN for'//ttxt(ita,it)
      txtp(ita,it)(59:69)=dato//'  '
      WRITE(2) txtp(ita,it)
      jri=jris(ita,it)
      WRITE(2) 'CHRD:KGRN ',dato
      WRITE(2) eln(ita,it),qtr(ita,it),qs(ita,it),qcpa(ita,it),
     .         qsca(it),nz(ita,it),ion(ita,it),ixch,jri
      WRITE(2) nqns(ita,it,1:mo),nks(ita,it,1:mo),nels(ita,it,1:mo),
     .         ncorbs(ita,it,1:mo),symbols(ita,it),configs(ita,it),
     .         dens(ita,it,1:ns,1:mo),dq1s(ita,it,1:ns,1:mo),
     .         izs(ita,it),norbs(ita,it),ions(ita,it),eonec(ita,it)
      WRITE(2) eone(ita,it),emadl(ita,it),vint(ita,it),enuc(ita,it),
     .         ecor(ita,it),2.d0*eval(ita,it)+ecor(ita,it),exct(ita,it),
     .         excc(ita,it),ekin(ita,it),etot(ita,it),okae(ita,it)
      WRITE(2) cor(ita,it,1:jri)
      jrn=jrsm(ita,it)
      DO 21 is=1,ns
      WRITE(2) pots(ita,it,is),potw(ita,it,is),vmtz(is),vmtzr(it,is)
      WRITE(2) chde(ita,it,is,1:jri)
      WRITE(2) chdo(ita,it,is,1:jrn)
   21 WRITE(2) v(1:jrn,ita,it,is)
C
      DO 22 l=0,lmax
      DO 22 it=1,nt
      DO 22 ita=1,nta(it)
      DO 22 is=1,ns
      DO 22 ip=1,pan
   22 WRITE(2) eny(l,ita,it,is,ip),cc(l,ita,it,is,ip),
     .         bot(l,ita,it,is,ip),top(l,ita,it,is,ip)
C
      DO 23 iq=1,nq
   23 WRITE(2) qlmo(iq,1:diml)
C
C     Save the complex contour
C
      WRITE(9) zmsh
      WRITE(9) eps,nz1,nz2,nz3,nres,nzm,nx,hx
      WRITE(9) zm(1:nzm),wgm(1:nzm)
      WRITE(9) ibz,nkx,nky,nkz,ibz2,nkx2,nky2,nkz2
C
      CLOSE(2)
      CLOSE(9)
C
      WRITE(m6,'(/,2a,4x,a)')
     .      ' Potential parameters, charge densities,',
     .      ' and potential stored on:',DATO
      IF(msgl.EQ.1) WRITE(msgio,'(/,a)')
     .         ' UPDATE: Potential etc. stored on FOR002'
      RETURN
      END
