      SUBROUTINE atomp(efgs)
C   ******************************************************************
C   *                                                                *
C   *    Prepare and perform atomic calculations.                    *
C   *                                                                *
C   ******************************************************************
      USE atomb
      USE atomicdens
      USE control_data
      USE density
      USE lattice
      USE potential
      USE softcore
      USE radialmesh
      USE text
      USE totalenergy
      IMPLICIT NONE
      INTEGER, PARAMETER :: mo=30
      INTEGER :: iq, it, ita, jwsam
      REAL(KIND=8) :: efgs, xwsa, dxa, dr1a, etotc
C
      ALLOCATE(hqtr(mnta,nt),kitq(nq),knta(nt))
      ALLOCATE(hws(mnta,nt),hhsr(mnta,nt),hwst(nq),hwsi(nq),hwsc(nq))
      ALLOCATE(etotcore(mnta,nt))
C
C     Transfer default values to atomic common block /ATOB/
C
      knt=nt
      knta(1:nt)=nta(1:nt)
      knq=nq
      knl=nl
      kns=ns
      hsws=sws
      hefgs=efgs
      hdexch=dexch
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      hws(ita,it)=ws(ita,it)
      hhsr(ita,it)=hsr(ita,it)
      hwst(iq)=wst(iq)
      hwsi(iq)=wsi(iq)
      hwsc(iq)=wsc(iq)
      hqtr(ita,it)=qtro(ita,it)
   20 kitq(iq)=itq(iq)
C
      CALL atomd(txch,dxa,dr1a)
C
C     Generate radial mesh
C
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      IF(nz(ita,it).GT.0) THEN
         r1(ita,it)=dr1a/nz(ita,it)
      ELSE
         r1(ita,it)=dr1a/20.d0
      ENDIF
      jwsam=LOG(ws(ita,it)/r1(ita,it))/dxa
      jwsam=2*(jwsam/2)+2
      xwsa=jwsam*dxa
      r1(ita,it)=ws(ita,it)*EXP(-xwsa)
      jwss(ita,it)=jwsam+1
      jris(ita,it)=jwss(ita,it)+2
   21 CONTINUE
C
      DO 22 it=1,nt
      DO 22 ita=1,nta(it)
      etotc=0.d0
      CALL atomi(ws(ita,it),r1(ita,it),jwss(ita,it),nz(ita,it),it,ita)
      CALL atomc(ws(ita,it),it,ita,ns,etotc)
      etotcore(ita,it)=etotc
      CALL atoms(it,ita,nt,mnta,jwss,jris,r1)
   22 CONTINUE
C
      RETURN
      END
