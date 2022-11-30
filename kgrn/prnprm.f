      SUBROUTINE prnprm(key,ef,etotal)
C   ******************************************************************
C   *                                                                *
C   *    Prints potential parameters etc.                            *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE botomtop   ; USE control_data
      USE control_text ; USE density ; USE dosmom ; USE greenfunc
      USE message    ; USE potential    ; USE potparam
      USE radialmesh ; USE slope        ; USE text
      USE totalenergy
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(mnta,nt,0:lmax) :: stoi
      REAL(KIND=8), DIMENSION(mnta,nt)        :: chdw, ewcorr
      REAL(KIND=8) :: tolmag = 1.d-3
      REAL(KIND=8) :: etotal, ef, cmdliq, okacorr, ewacorr
      REAL(KIND=8) :: hopfield, fl, flp1, ttdos
      INTEGER      :: key, iq, it, ita, l, lp1, ntait, is, ip, jws
C
C     Write potential parameters
C
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      DO 20 ip=1,pan
      WRITE(m6,101) ttxt(ita,it),hsr(ita,it),sws
      WRITE(m6,102) txch,dato
      WRITE(m6,103) ip,((thead(l,is),l=0,lmax),is=1,ns)
      WRITE(m6,104) ((eny(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,105) ((dny(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,106) ((d2(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,107) ((omm(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,108) ((signfi(l,ita,it,is,ip)
     .              *sfim(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,109) ((fmofp(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,110) ((amy(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,111) ((tl(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,112) ((top(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,113) ((cc(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,114) ((bot(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,115) ((vl(l,ita,it,is,ip),l=0,lmax),is=1,ns)
      WRITE(m6,141) ((logdl(l,ita,it,is),l=0,lmax),is=1,ns)
      WRITE(m6,116) ebt,etp
      WRITE(m6,117) potm(ita,it)
      WRITE(m6,118) (pots(ita,it,is),is=1,ns)
   20 WRITE(m6,119) (potw(ita,it,is),is=1,ns)
C
      IF(key.EQ.1) RETURN
C
C     Calculate spin moments
C
      IF(ns.EQ.2) THEN
         DO 21 it=1,nt
         DO 21 ita=1,nta(it)
         DO 21 l=0,lmax
         IF(ABS(amag(ita,it)).GT.tolmag) THEN
            stoi(ita,it,l)=(cc(l,ita,it,1,pan)-cc(l,ita,it,2,pan))/
     .                     amag(ita,it)
         ELSE
            stoi(ita,it,l)=0.d0
         ENDIF
   21    CONTINUE
      ENDIF
C
C     Print moments
C
      WRITE(m6,120) ef
      WRITE(m6,121) vmtz(1:ns)
      DO 22 it=1,nt
      DO 22 ita=1,nta(it)
      WRITE(m6,122) ttxt(ita,it)
      DO 23 ip=1,pan
      WRITE(m6,103) ip,((thead(l,is),l=0,lmax),is=1,ns)
      WRITE(m6,*)
      IF(ip.EQ.pan) THEN
         WRITE(m6,123) ((tnos(ita,it,l,is),l=0,lmax),is=1,ns),
     .                 SUM(tnos(ita,it,0:lmax,1:ns))
         WRITE(m6,130) ((tdos(ita,it,l,is),l=0,lmax),is=1,ns),
     .                 SUM(tdos(ita,it,0:lmax,1:ns))
      ELSE
         WRITE(m6,123) ((tnos0(ita,it,l,is),l=0,lmax),is=1,ns),
     .                 SUM(tnos0(ita,it,0:lmax,1:ns))
      ENDIF
      IF(ns.EQ.2) THEN
         WRITE(m6,124) amag(ita,it),tmag
         WRITE(m6,125) (stoi(ita,it,l),l=0,lmax)
      ENDIF
   23 CONTINUE
C
      DO is=1,ns
      hopfield=0.d0
      ttdos=SUM(tdos(ita,it,0:lmax,is))
      IF(ABS(ttdos).GT.1.d-6) THEN
         DO l=0,lmax-1
         lp1=l+1
         fl=tdos(ita,it,l,is)/ttdos
         flp1=tdos(ita,it,lp1,is)/ttdos
         hopfield=hopfield+lp1*epmatrix(l,ita,it,is)**2.*
     .         fl*flp1/(2.d0*l+1.d0)/(2.d0*l+3.d0)
         ENDDO
         hopfield=2.d0*ttdos*hopfield
         WRITE(m6,140) hopfield,hopfield*48.58725d0     ! Ry/Bohr^2-->eV/A^2
      ENDIF
      ENDDO
   22 CONTINUE
C
      DO 24 it=1,nt
      DO 24 ita=1,nta(it)
      jws=jwss(ita,it)
      chdw(ita,it)=chdo(ita,it,1,jws)
      IF(ns.EQ.2) chdw(ita,it)=chdw(ita,it)+chdo(ita,it,2,jws)
   24 CONTINUE
C
      ewcorr=0.d0
      DO 25 iq=1,nq
      it=itq(iq)
      ntait=nta(it)
      cmdliq=1.8d0+SUM(vmdl(iq,1:nq))/2.d0
      DO 25 ita=1,ntait
   25 ewcorr(ita,it)=ewcorr(ita,it)+
     .      cmdliq*ws(ita,it)/9.d0*chdw(ita,it)*chdw(ita,it)
      DO 26 it=1,nt
   26 ewcorr(1:ntait,it)=ewcorr(1:ntait,it)/mmt(it)
C
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
   30 WRITE(m6,126) eone(ita,it),vint(ita,it),ekin(ita,it),ecor(ita,it),
     .            enuc(ita,it),emadl(ita,it),eval(ita,it),exct(ita,it),
     .            excc(ita,it),etot(ita,it),okae(ita,it),
     .            etot(ita,it)+okae(ita,it),ents(ita,it),
     .            etot(ita,it)-ents(ita,it)
      IF(softc.NE.'Y') THEN
         WRITE(m6,127)
      ELSE
         WRITE(m6,128)
      ENDIF
C
      ewacorr=0.d0
      okacorr=0.d0
      DO 31 it=1,nt
      IF(it.EQ.1) WRITE(m6,129)
      DO 31 ita=1,nta(it)
      WRITE(m6,'(/,11x,a,f14.6)') 'Kinetic energy   ',ekin(ita,it)
      WRITE(m6,'(11x,a,f14.6)') 'Madelung energy  ',emadl(ita,it)
      WRITE(m6,'(11x,a,f14.6)') 'El-ion energy    ',enuc(ita,it)
      WRITE(m6,'(11x,a,f14.6)') 'El-el energy     ',
     .      ecor(ita,it)+eval(ita,it)
      WRITE(m6,'(11x,a,f14.6)') 'Exc energy       ',
     .      exct(ita,it)-excc(ita,it)
      WRITE(m6,'(30x,a)') '------------'
      WRITE(m6,'(11x,a,f14.6)') 'Total            ',etot(ita,it)
      WRITE(m6,'(11x,a,f14.6)') 'Total+OKAE       ',
     .      etot(ita,it)+okae(ita,it)
      WRITE(m6,'(11x,a,f14.6,a,f10.6)') 'Total+Ewald      ',
     .      etot(ita,it)+ewcorr(ita,it),' 4pi SS n(S)=',chdw(ita,it)
      okacorr=okacorr+conc(ita,it)*mmt(it)*okae(ita,it)
      ewacorr=ewacorr+conc(ita,it)*mmt(it)*ewcorr(ita,it)
   31 CONTINUE
      IF(nt.GT.1.OR.nta(nt).GT.1) THEN
         WRITE(m6,'(/,11x,a,2f14.6)') 'Total energy     ',etotal
         WRITE(m6,'(/,11x,a,2f14.6)') 'Total energy+OKA ',
     .                    etotal+okacorr/nq
         WRITE(m6,'(/,11x,a,f14.6)') 'Total+Ewald      ',
     .                    etotal+ewacorr/nq
      ENDIF
      WRITE(m6,'(30x,a)') '============'
C
      IF(nt.EQ.1.AND.conv.EQ.'Y') THEN
      DO ita=1,nta(nt)
      WRITE(m6,'(/,a,f10.6,f14.6,5f10.6)') 'result ',ws(ita,nt),
     .          etot(ita,nt),vmtz(1:ns)-ef,ef,vmtz(1:ns)
      ENDDO
      ENDIF
C
  101 FORMAT(/,' Atom:',a4,12x,'S =',f12.6,8x,'SWS =',f12.6)
  102 FORMAT(' Exch:',a3,30x,'Date : ',a)
  103 FORMAT(/,' Panel	    =',i3,3x,a,9(6x,a))
  104 FORMAT(/,' E-nu       =',8f12.6)
  105 FORMAT(' D-nu       =',8f12.6)
  141 FORMAT(' D(Ef)      =',8f12.6)
  106 FORMAT(' Dnud       =',8f12.6)
  107 FORMAT(  ' Omega-     =',8f12.6)
  108 FORMAT(' S*Phisq-   =',8f12.6)
  109 FORMAT(' Phi-/Phi+  =',8f12.6)
  110 FORMAT(/,' Mu         =',8f12.6)
  111 FORMAT(' Tau        =',8f12.6)
  112 FORMAT(' A          =',8f12.6)
  113 FORMAT(' C          =',8f12.6)
  114 FORMAT(' B          =',8f12.6)
  115 FORMAT(' V          =',8f12.6)
  116 FORMAT(/,' Bot, Top   =',8f12.6)
  117 FORMAT(' V-Madelung =',f12.6)
  118 FORMAT(' V(S) Up,Dwn=',f12.6,36x,f12.6)
  119 FORMAT(' V(W) Up,Dwn=',f12.6,36x,f12.6)
  120 FORMAT(/,' EF         =',f12.6)
  121 FORMAT(/,' VMTZ Up,Dwn=',f12.6,36x,f12.6)
  122 FORMAT(/,' Atom:',a4)
  123 FORMAT(' Nos(Ef)    =',8f12.6)
  124 FORMAT(' Magn. mom. =',8f12.6)
  125 FORMAT(' Stoner-I   =',8f12.6)
  126 FORMAT(/,' PRNPRM:   EONE =',f14.6,' VINT =',f14.6,' EKIN =',
     .        f14.6,/,11x,'ECOR =',f14.6,' ENUC =',f14.6,' EMADL=',
     .        f14.6,/,11x,'EVAL =',f14.6,' EXCT =',f14.6,' EXCC =',
     .        f14.6,/,11x,'ETOT =',f14.6,' OKAE =',f14.6,' E[n] =',
     .        f14.6,/,11x,' T*S =',f14.6,' E-TS =',f14.6)
  127 FORMAT(/,11x,'Valence energies in the frozen core approximation',
     .       ' (Rydbergs)')
  128 FORMAT(/,11x,'Total energies in the soft core approximation',
     .       ' (Rydbergs)')
  129 FORMAT(11x,'calculated from output charge density:')
  130 FORMAT(' Dos(Ef)    =',8f12.6)
  140 FORMAT(' Hopfield   =',f12.6,' (Ry/Bohr^2), ',f12.6,' (eV/AA^2)')
C
      RETURN
      END
