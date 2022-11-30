      SUBROUTINE oldin(efgs,ixch,in)
C   ******************************************************************
C   *                                                                *
C   *    Read potential parameters, renormalised charge densities,   *
C   *    total charge densities, potentials, and total energies      *
C   *    generated earlier by BGRN.                                  *
C   *                                                                *
C   ******************************************************************
      USE atomicdens   ; USE botomtop     ; USE bzmesh
      USE control_data ; USE control_text ; USE density
      USE energymesh   ; USE lattice      ; USE message
      USE moments      ; USE pota         ; USE potential
      USE potparam     ; USE radialmesh   ; USE softcore
      USE text         ; USE totalenergy 
      IMPLICIT NONE
      CHARACTER(LEN=3) :: funco
      CHARACTER(LEN=4) :: tprg, tdum
      INTEGER, PARAMETER :: mo = 30
      INTEGER :: in, ixch, nqd, nld, iq, it, ita, is, ip, jri, jrn, l
      INTEGER :: jq, jt, ntait, ntajt
      REAL(KIND=8) :: efgs, efg1, dum
C
      REWIND in
      READ(in) tprg,nqd,nld,conv,funco
      IF(tprg.EQ.'ATOM') THEN
         tpot='N'
      ELSE
         tpot='Y'
         IF(funco.NE.func) THEN
            WRITE(m6,104) funco,func
            STOP
         ENDIF
      ENDIF
      IF(nqd.NE.nq) THEN
         WRITE(m6,102) nqd,nq
         STOP
      ENDIF
      IF(nld.NE.nl) THEN
         WRITE(m6,103) nld,nl
         STOP
      ENDIF
      WRITE(m6,101) tprg
      READ(in) sws,wst(1:nq),wsi(1:nq),wsc(1:nq),dx
      READ(in) eb,efg1,ns,nt,itq(1:nq),nta(1:nt)
      IF(in.EQ.2) efgs=efg1
      WRITE(m6,100) nt,ns,efgs
      READ(in) ebt,etp,dum,dexch
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
   20 READ(in) ws(ita,it),hsr(ita,it),r1(ita,it),jwss(ita,it),
     .         jris(ita,it),ttxt(ita,it)
C
      CALL setrms
C 
      CALL arrays(1)
C
      DO 22 it=1,nt
      DO 22 ita=1,nta(it)
      READ(in) txtp(ita,it)
      READ(in) tdum
      READ(in) eln(ita,it),qtro(ita,it),qs(ita,it),qcpa(ita,it),
     .         qsca(it),nz(ita,it),ion(ita,it),ixch,jri
      IF(tpot.EQ.'Y') THEN
         READ(in) nqns(ita,it,1:mo),nks(ita,it,1:mo),nels(ita,it,1:mo),
     .            ncorbs(ita,it,1:mo),symbols(ita,it),configs(ita,it),
     .            dens(ita,it,1:ns,1:mo),dq1s(ita,it,1:ns,1:mo),
     .            izs(ita,it),norbs(ita,it),ions(ita,it),eonec(ita,it)
         READ(in) eone(ita,it),emadl(ita,it),vint(ita,it),enuc(ita,it),
     .            ecor(ita,it),dum,exct(ita,it),excc(ita,it),
     .            ekin(ita,it),etot(ita,it),okae(ita,it)
      ELSE
         READ(in) dum
      ENDIF
      READ(in) core(ita,it,1:jri)
      IF(tpot.EQ.'N') THEN
         jrn=jri
      ELSE
         jrn=jrsm(ita,it)
      ENDIF
      DO 22 is=1,ns
      IF(tpot.EQ.'Y') THEN
         READ(in) pots(ita,it,is),potw(ita,it,is),vmtz(is),
     .            vmtzr(it,is)
         READ(in) chde(ita,it,is,1:jri)
      ENDIF
      READ(in) chdo(ita,it,is,1:jrn)
   22 READ(in) v(1:jrn,ita,it,is)
C
      IF(tpot.EQ.'Y') THEN
         DO 23 l=0,lmax
         DO 23 it=1,nt
         DO 23 ita=1,nta(it)
         DO 23 is=1,ns
         DO 23 ip=1,pan
         READ(in) eny(l,ita,it,is,ip),cc(l,ita,it,is,ip),
     .            bot(l,ita,it,is,ip),top(l,ita,it,is,ip)
   23    CONTINUE
C
         DO 24 iq=1,nq
   24    READ(in) qlmo(iq,1:diml)
      ENDIF
C
      IF(msgl.EQ.1) WRITE(msgio,'(/,2a)')
     .         ' OLDIN: Potential etc. generated by ',tprg
C
C     Determine the number MMT(IT) of type t atoms.
C
      DO 25 it=1,nt
      mmt(it)=0
      DO 25 iq=1,nq
      IF(itq(iq).EQ.it) mmt(it)=mmt(it)+1
   25 CONTINUE
C
      DO 26 iq=1,nq
      it=itq(iq)
      IF(mmt(it).GT.1) THEN
         ntait=nta(it)
         IF(ntait.GT.1) THEN
            DO 27 jq=1,nq
            jt=itq(jq)
            IF(jq.NE.iq.AND.jt.EQ.it) THEN
               ntajt=nta(jt)
               IF(ntajt.NE.ntait) THEN
                  WRITE(m6,110) iq,it,ntait,jq,jt,ntajt 
                  STOP
               ENDIF
               DO 28 ita=1,ntait
               IF(nz(ita,it).NE.nz(ita,jt)) THEN
                  WRITE(m6,111) iq,it,ita,nz(ita,it),
     .                          jq,jt,ita,nz(ita,it)
                  STOP
               ENDIF
   28          CONTINUE
            ENDIF
   27       CONTINUE
         ENDIF
      ENDIF
   26 CONTINUE
C
      IF(tpot.EQ.'Y') THEN
         IF(strt.NE.'N') THEN
            READ(9) zmsh
            READ(9) eps,nz1,nz2,nz3,nres,nzm,nx,hx
            READ(9) dum
            READ(9) ibz,nkx,nky,nkz,ibz2,nkx2,nky2,nkz2
         ENDIF
      ENDIF
C
  100 FORMAT(/,' OLDIN:    NT     =',i3,' NS   =',i3,//,11x,
     .       'EFGS   =',f10.6)
  101 FORMAT(/,11x,'Self-consistent bulk program: ',a,/)
  102 FORMAT(/,' OLDIN:**  NQD =',i3,' Must be equal to NQ =',i3)
  103 FORMAT(/,' OLDIN:**  NLD =',i3,' Must be equal to NL =',i3)
  104 FORMAT(/,' OLDIN:**  funco =',a,' func =',a)
  110 FORMAT(/,' OLDIN:**  IQ =',i3,' IT =',i3,' NTA =',i3,//,
     .                11x,'JQ =',i3,' JT =',i3,' NTA =',i3,//
     .                11x,' inconsistent setup !!')
  111 FORMAT(/,' OLDIN:**  IQ =',i3,' IT =',i3,' ITA =',i3,' Z =',i3,//,
     .                11x,'JQ =',i3,' JT =',i3,' ITA =',i3,' Z =',i3,//,
     . 11x,' the order of sorts on the sublattices should be the same ')
      RETURN
      END