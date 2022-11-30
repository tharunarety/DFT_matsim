      SUBROUTINE chdren(ef,lin)
C   ******************************************************************
C   *                                                                *
C   *   Renormalize the exact spherical density (chdo) for the ASA   *
C   *   Madelung potential needed for the next iteration.            *
C   *   The renormalized density (chde) is used for the ASA total    *
C   *   energy as well.                                              *
C   *                                                                *
C   *   For each sublattice there is a net charge which is treated   *
C   *   by the usual Madelung term (qlmo). The deviation between     *
C   *   the actual net charge (qsit) and qlmo is due to the SCA      *
C   *   and it is treated by the MADC term.                          *
C   *   For disordered alloys the difference between the actual      *
C   *   net charge for the alloy component (qs) and the net charge   *
C   *   on the sublattice (qsit) is treated by the SIM correction.   *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE density ; USE message
      USE moments    ; USE radialmesh   ; USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(mnta,nt) :: zp, add
      REAL(KIND=8) :: facns, zwc, addz, r, qtot, qadd, qtrit, qsit, ef
      INTEGER :: lin, iq, it, ita, ntait, is, jrn, jws, jsr, jri, ir
C
C     Enforce normalization for the whole unit cell (not per atom !!!)
C
      facns=spinfc/2.d0
      qtot=0.d0
      addz=0.d0
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jws=jwss(ita,it)
      jrn=jrsm(ita,it)
      jsr=jsrs(ita,it)
      IF(ns.EQ.1) THEN
         fi(1:jrn)=chdo(ita,it,1,1:jrn)*ri(1:jrn,ita,it)
      ELSE
         fi(1:jrn)=(chdo(ita,it,1,1:jrn)+chdo(ita,it,2,1:jrn))*
     .             ri(1:jrn,ita,it)
      ENDIF
      CALL simpn(fi,dx,jws,zwc)
      zp(ita,it)=zwc
      qs(ita,it)=zwc-nz(ita,it)
C
      qtot=qtot+conc(ita,it)*mmt(it)*zn0(ita,it)
      addz=addz+conc(ita,it)*mmt(it)*zp(ita,it)
   20 CONTINUE
      addz=(znt-addz)/qtot*facns
C
C     Set up the renormalized density
C
      IF(ABS(addz).GT.1.d-16) THEN
         add(1:mnta,1:nt)=0.d0
         qadd=0.d0
C
C        Renormalize only those sites which have Z <> 0 or
C        the renormalized density is positive
C
         DO 21 it=1,nt
         DO 21 ita=1,nta(it)
         IF(vac(ita,it).NE.1) THEN
            IF(nz(ita,it).GT.0.OR.
     .         (qs(ita,it)+addz*zn0(ita,it)).GT.0.d0) THEN
               add(ita,it)=addz
               qadd=qadd+conc(ita,it)*mmt(it)*zn0(ita,it)
            ENDIF
         ENDIF
   21    CONTINUE
         add(1:mnta,1:nt)=add(1:mnta,1:nt)*qtot/qadd
C
         DO 22 is=1,ns
         DO 22 it=1,nt
         DO 22 ita=1,nta(it)
         jrn=jrsm(ita,it)
   22    chde(ita,it,is,1:jrn)=chdo(ita,it,is,1:jrn)+
     .                         add(ita,it)*chdr(ita,it,1:jrn)
      ENDIF
C
C     Calculate the new monopole moments
C
      qtot=0.d0
      DO 23 it=1,nt
      DO 23 ita=1,nta(it)
      jws=jwss(ita,it)
      IF(ns.EQ.1) THEN
         fi(1:jws)=chde(ita,it,1,1:jws)*ri(1:jws,ita,it)
      ELSE
         fi(1:jws)=(chde(ita,it,1,1:jws)+
     .              chde(ita,it,2,1:jws))*ri(1:jws,ita,it)
      ENDIF
      CALL simpn(fi,dx,jws,zwc)
      qtr(ita,it)=zwc-nz(ita,it)
      qtot=qtot+conc(ita,it)*mmt(it)*qtr(ita,it)
   23 CONTINUE
C
C     Check for the charge neutrality
C
      IF(ABS(qtot).GT.1.d-6) THEN
         WRITE(m6,100) 'SUM(QTR)=',qtot,' QTR =',qtr
         IF(msgl.NE.0) WRITE(msgio,100) 'SUM(QTR)=',qtot,
     .                                  ' QTR =',qtr
         STOP
      ENDIF
C
C     qs(ita,it)  net charge within the ws(ita,it) sphere;
C     qtr(ita,it) net charge within the ws(ita,it) sphere after the
C                 density was renormalized;
C     Ordered alloy:
C       qsit     total charge (relative to Z) within the ws(it) sphere
C       qsca(it) extra (missing) charge within the ws(it) sphere which
C                needs to be redistributed on the coordination spheres;
C                it is due to the SCA to the Poisson's eq.
C     Disordered alloy:
C       qtrit average charge transfer on the sublattice it (used in Madelung)
C       qcpa(ita,it) extra (missing) charge within the ws(ita,it) sphere
C                    relative to qsit and it is treated by the SIM;
C                    this term is due to the CPA.
C
      DO 24 it=1,nt
      ntait=nta(it)
      qsit=SUM(conc(1:ntait,it)*qs(1:ntait,it))   !! Qsca
      qtrit=SUM(conc(1:ntait,it)*qtr(1:ntait,it)) !! Qlm
      qsca(it)=qsit-qtrit                         !! Qsca - Qlm
      qcpa(1:ntait,it)=qs(1:ntait,it)-qsit        !! Qsim
   24 CONTINUE
C
C     Average zero order moments
C
      DO 30 iq=1,nq
      it=itq(iq)
      ntait=nta(it)
      qtrit=SUM(conc(1:ntait,it)*qtr(1:ntait,it))
   30 qlmo(iq,1)=qtrit
C
  100 FORMAT(/,' CHDREN:  ',a,f10.6,a,10f10.6)
C
      RETURN
      END
