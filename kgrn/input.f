      SUBROUTINE input(efgs,lmaxh0,nlin)
C   ******************************************************************
C   *                                                                *
C   *    Read control data, charge densities and potentials for a    *
C   *    new run. A cold start performs atomic calculations to       *
C   *    initiate the self-consistent band calculations. In case     *
C   *    of poor convergence the calculations may be continued from  *
C   *    the data stored in the potential file.                      *
C   *                                                                *
C   ******************************************************************
      USE atomb        ; USE atomicdens   ; USE bzmesh
      USE control_data ; USE control_text ; USE density
      USE diracparam
      USE energymesh   ; USE gaussi       ; USE greenfunc
      USE lattice      ; USE message      ; USE pota
      USE potential    ; USE potparam     ; USE radialmesh
      USE slope        ; USE softcore     ; USE text
      USE volume
      IMPLICIT NONE
      CHARACTER(LEN=1) :: kmsh, fixsh
      CHARACTER(LEN=4) :: txth
      CHARACTER(LEN=5) :: label
      INTEGER, PARAMETER :: mo=30
      INTEGER :: lmaxh0, nlin, in, iq, it, ita, jta, ntait
      INTEGER :: nzit, nz1h, oomtz
      REAL(KIND=8), PARAMETER :: boltzk = 6.336d-6, conctol = 1.d-4
      REAL(KIND=8) :: efgs, tol1
      REAL(KIND=8) :: sri, srm, qtrh, splith, dsws, wsh, concit
C
      READ(5,'(a)') label
      IF(label.NE.'Band:') THEN
         WRITE(m6,'(a)') ' INPUT: Header must be Band:'
         STOP
      ENDIF
      READ(5,'(6(7x,i3))') niter,nlin,nprn,ncpa,nt,mnta
      READ(5,'(8x,a,5(9x,a))') mode,frc,dosc,ops,afm,crt
      READ(5,'(5(7x,i3),9x,a)') lmaxh,lmaxt,nfi,oomtz,shf,softc
      IF(lmaxh.GT.lmaxh0) lmaxh=lmaxh0
      IF(lmaxt.GT.lmaxh0) lmaxt=lmaxh0
      nlmh=(lmaxh+1)*(lmaxh+1)
      nlmqh=nlmh*nq
      IF(afm.EQ.'F') THEN
         ns=2
      ELSEIF(afm.EQ.'M'.OR.afm.EQ.'m') THEN
         ns=2
      ELSEIF(afm.EQ.'A') THEN
         ns=2
      ELSEIF(afm.EQ.'P') THEN
         ns=1
      ELSE
         WRITE(m6,100) afm
         STOP
      ENDIF
C
C     Allocate arrays of nq,nt,ns,...
C
      ALLOCATE(itq(nq),nta(nt),conc(mnta,nt),vmtz(ns),vmtzr(nt,ns))
      ALLOCATE(localmt(mnta,nt))
      ALLOCATE(ttxt(mnta,nt),txtp(mnta,nt))
      ALLOCATE(mmt(nt),qtro(mnta,nt),split(mnta,nt),fixst(mnta,nt))
      ALLOCATE(hsr(mnta,nt),nz(mnta,nt),wsm(mnta,nt))
      ALLOCATE(jris(mnta,nt),jwss(mnta,nt),r1(mnta,nt),ws(mnta,nt))
      ALLOCATE(asc(nlmq))
      ALLOCATE(eln(mnta,nt),ion(mnta,nt),qtr(mnta,nt),qti(mnta,nt))
      ALLOCATE(nqns(mnta,nt,mo),nks(mnta,nt,mo),nels(mnta,nt,mo))
      ALLOCATE(ncorbs(mnta,nt,mo),symbols(mnta,nt))
      ALLOCATE(configs(mnta,nt),dens(mnta,nt,ns,mo),dq1s(mnta,nt,ns,mo))
      ALLOCATE(izs(mnta,nt),norbs(mnta,nt),ions(mnta,nt),eonec(mnta,nt))
      ALLOCATE(vac(mnta,nt))
      ALLOCATE(hdexch(mnta,nt),dexch(mnta,nt),dexcho(mnta,nt))
      ALLOCATE(fixsjta(mnta,nt),splitoo(mnta,nt),splito(mnta,nt))
      vmtz=0.d0
      vmtzr=0.d0
      fixsjta=0
      hdexch=0.d0
      dexch=0.d0
      dexcho=0.d0
      vac=0
      conc=0.d0
      qtro=0.d0
      split=0.d0
      hsr=0.d0
      nz=0
      jris=0
      jwss=0
      r1=0.d0
      ws=0.d0
      eln=0.d0
      ion=0
      qtr=0.d0
      qti=0.d0
      nqns=0
      nks=0
      nels=0
      ncorbs=0
      dens=0.d0
      dq1s=0.d0
      izs=0
      norbs=0
      ions=0
      eonec=0.d0
      localmt=0
C
C     Read integers specifying k and z mesh
C
      READ(5,'(9x,a,4(7x,i3),9x,a)') kmsh,ibz,nkx,nky,nkz,fllbz
      IF(ibz.NE.lat) THEN
         WRITE(m6,'(/,a,i3,a,i3,a)') ' INPUT:** IBZ =',ibz,
     .        ' should be the same as LAT =',lat,' from KSTR'
         STOP
      ENDIF
      READ(5,'(9x,a,5(7x,i3))') kmsh,ibz2,nkx2,nky2,nkz2
      READ(5,'(9x,a,4(7x,i3),6x,i4)') zmsh,nz1,nz2,nz3,nres,nzd
      IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
         pan=2
      ELSE
         pan=1
      ENDIF
      IF(zmsh.EQ.'c'.OR.zmsh.EQ.'e'.OR.zmsh.EQ.'m') THEN
         nz1h=nz1/2
         IF(nz1.NE.2*nz1h) THEN
            WRITE(m6,'(/,a,i3,2a)') ' INPUT:** nz1 =',nz1,
     .               ' must be even for ZMSH =',zmsh
            STOP
         ENDIF
      ENDIF
C
      ALLOCATE(eny(0:lmax,mnta,nt,ns,pan),cc(0:lmax,mnta,nt,ns,pan))
      ALLOCATE(top(0:lmax,mnta,nt,ns,pan),bot(0:lmax,mnta,nt,ns,pan))
C
      READ(5,'(4(8x,f7.3))') depth,imagz,eps,elim
      READ(5,'(4(8x,f7.3))') amix,efmix,vmtz0,mmom
      IF(afm.EQ.'F') THEN
         dexch=0.5d0*mmom
         IF(ABS(mmom).GT.1.d-06) WRITE(m6,101) 0.5d0*mmom
      ELSEIF(afm.EQ.'M') THEN
         dexch=0.5d0*mmom
         WRITE(m6,102) mmom
      ENDIF
      READ(5,'(4(8x,f7.4))') tole,tolef,tolcpa,tfermi
      READ(5,'(10x,f10.6,7x,i3,8x,f7.2,8x,f7.4)') sws,nsws,dsws,alphmd
      IF(nsws.GT.1) THEN
         WRITE(m6,'(a)')  
     .    ' INPUT:*** NSWS > 1 is not implemented'
         STOP
      ENDIF
C
      lclmff=0
      IF(oomtz.EQ.0) THEN
         fixvmtz=0
         WRITE(m6,120)
      ELSEIF(oomtz.EQ.1) THEN
         fixvmtz=1
         WRITE(m6,121) vmtz0
      ELSEIF(oomtz.EQ.2) THEN
         fixvmtz=2
         WRITE(m6,122)
      ELSEIF(oomtz.EQ.3) THEN
         fixvmtz=3
         WRITE(m6,123) vmtz0
      ELSEIF(oomtz.EQ.4) THEN
         fixvmtz=2
         lclmff=1
         WRITE(m6,124)
      ELSE
         WRITE(m6,'(a)')
     .    ' INPUT:*** FIXG > 4 is not implemented'
         STOP
      ENDIF
C
      CALL setsws(dsws)
C
      tfermi=tfermi*boltzk
      IF(lmaxt.GT.lmax) THEN
         WRITE(m6,109) lmaxt
      ELSE
         lmaxt=lmax
      ENDIF
      nlmt=(lmaxt+1)*(lmaxt+1)
      nlmqt=nlmt*nq
      IF(fcd.EQ.'Y') WRITE(m6,110) lmaxh
      softz=0
      IF(softc.EQ.'Y') THEN
         softz=2
         WRITE(m6,111)
      ELSEIF(softc.EQ.'Z') THEN
         softz=1
         WRITE(m6,117)
      ENDIF
      IF(softz.EQ.1) softc='N'
      IF(frc.EQ.'Y') WRITE(m6,112) shf
      IF(dosc.EQ.'D') THEN
         WRITE(m6,113)
      ELSEIF(dosc.EQ.'S') THEN
         WRITE(m6,114)
      ELSEIF(dosc.EQ.'N') THEN
         CONTINUE
      ELSE
         WRITE(m6,116) dosc
         STOP
      ENDIF
      IF(crt.NE.'M'.AND.crt.NE.'m'.AND.crt.NE.'I') THEN
         WRITE(m6,115) crt
         STOP
      ENDIF
      READ(5,'(a)') label
      IF(label.NE.'Setup') THEN
         WRITE(m6,'(A)') ' INPUT: Header must be Setup'
         STOP
      ENDIF
C
C     Set and test parameters for the z-mesh
C
      READ(5,'(2(8x,f7.3),2(7x,i3),9x,a)') efgs,hx,nx,nz0,stmp
      IF(zmsh.EQ.'c'.OR.zmsh.EQ.'e'.OR.zmsh.EQ.'m') THEN
         nz1h=nz0/2
         IF(nz0.NE.2*nz1h) THEN
            WRITE(m6,'(/,a,i3,2a)') ' INPUT:** nz0 =',nz0,
     .               ' must be even for ZMSH =',zmsh
            STOP
         ENDIF
      ENDIF
      vmtz(1:ns)=efgs-vmtz0
      WRITE(m6,103) tole,tolef,tolcpa,depth,imagz,eps,amix,efmix,
     .              alphmd,
     .              nq,nt,nl,ns,mnta,kmsh,ibz,nkx,nky,nkz,zmsh,
     .              nz1,nz2,nz3,nres,niter,nlin,ncpa,afm,crt
      IF(2*(nx/2).EQ.nx) THEN
         WRITE(m6,104) nx
         STOP
      ENDIF
      READ(5,'()')
C
      empt='N'
C
C     Read setup
C
      nta=0
   20 CONTINUE
      READ(5,'(a)') label
      BACKSPACE(5)
      IF(label.NE.'Atom:') THEN
         READ(5,'(a,2x,3i3,i4,4f7.3,2f5.2,2x,a)') txth,iq,it,
     .       ita,nzit,concit,srm,sri,wsh,qtrh,splith,fixsh
         IF(it.GT.nt) THEN
            WRITE(m6,105) 'IT',it,'NT',nt
            STOP
         ENDIF
         IF(ita.GT.nta(it)) nta(it)=ita
         IF(nta(it).GT.mnta) THEN
            WRITE(m6,105) 'NTA(IT)',nta(it),'MNTA',mnta
            STOP
         ENDIF
         nz(ita,it)=nzit
         IF(nzit.EQ.0) empt='Y'
         itq(iq)=it
         IF(txth.EQ.'Va  ') vac(ita,it)=1
         ttxt(ita,it)=txth
         concit=MAX(0.d0,concit)
         conc(ita,it)=concit
C
         srm=MAX(1.d0,srm)
         IF(lclmff.EQ.0) srm=1.d0
         wsm(ita,it)=srm
         hsr(ita,it)=sri
         ws(ita,it)=wsh
         qtro(ita,it)=qtrh
         split(ita,it)=splith
         fixst(ita,it)=fixsh
         GO TO 20
      ENDIF
C
C     Test for partial DLM calculation
C
      IF(afm.EQ.'F') fixst='N'
      splito=split
      splitoo=split
      IF(afm.EQ.'m') THEN
         DO it=1,nt
         IF(nta(it).GT.1) THEN
            DO ita=1,nta(it)
            DO jta=1,nta(it)
            IF(nz(ita,it).EQ.nz(jta,it)) THEN
               IF(split(ita,it)*split(jta,it).LT.0.d0) THEN
                 IF(fixst(ita,it).NE.'Y') THEN
                   IF(conc(ita,it).LE.conc(jta,it)) THEN
                   IF(fixsjta(jta,it).NE.ita) fixsjta(ita,it)=jta
                   ENDIF
                 ENDIF
               ENDIF
            ENDIF
            ENDDO
            ENDDO
         ENDIF
         ENDDO
      ENDIF
C
      IF(empt.EQ.'Y') WRITE(m6,125)
      srde=0
      IF(srde.EQ.1) WRITE(m6,126)
C
      CALL sites
C
      IF(strt.EQ.'A') THEN
C
C        Determine the number MMT(IT) of type atoms.
C
         DO 21 it=1,nt
         mmt(it)=0
         DO 21 iq=1,nq
         IF(itq(iq).EQ.0) THEN
            WRITE(m6,106) iq,itq(iq)
            STOP
         ELSEIF(itq(iq).EQ.it) THEN
            mmt(it)=mmt(it)+1
         ENDIF
   21    CONTINUE
         DO 22 it=1,nt
         IF(mmt(it).EQ.0) THEN
            WRITE(m6,107) it,mmt(it)
            WRITE(m6,108) itq(1:nq)
            STOP
         ENDIF
   22    CONTINUE
      ENDIF
C
C     Renormalize the concentration
C
      cpa='N'
      DO 23 it=1,nt
      ntait=nta(it)
      concit=SUM(conc(1:ntait,it))
      IF(ABS(concit-1.d0).GT.conctol) THEN
         WRITE(m6,127) it
         WRITE(m6,128) conc(1:ntait,it)
      ENDIF
      conc(1:ntait,it)=conc(1:ntait,it)/concit
      DO 24 ita=1,ntait
      concit=conc(ita,it)
      IF(concit.GT.conctol.AND.concit.LT.(1.d0-conctol)) cpa='Y'
   24 CONTINUE
   23 CONTINUE
C
C     Compute meps, the relative machine precision
C
      meps=1.0d0
   30 meps=meps/2.0d0
      tol1=1.0d0+meps
      IF(tol1.GT.1.0d0) GO TO 30
C
      thead=RESHAPE((/'s -1/2','p -1/2','d -1/2','f -1/2',
     .      'g -1/2','h -1/2','s  1/2','p  1/2','d  1/2',
     .      'f  1/2','g  1/2','h  1/2'/),(/ 6,2 /))
C
  100 FORMAT(/,' INPUT:**  AFM =',a,' not implemented.',
     .       ' Must be P for paramagnetic, F for ferromagnetic,',
     .       ' A for antiferromagnetic, M for fix total spin,',
     .       ' m for fixed individual spin')
  101 FORMAT(/,' INPUT:**  Warning: You have introduced an external',
     .       ' exchange field',/,11x,'into the calculations, DEXC =',
     .       f10.6)
  102 FORMAT(/,' INPUT:    Fixed spin calculation, MMOM = ',f10.6)
  103 FORMAT(/,11x,
     .       'TOLE   =',e10.4,' TOLEF =',e10.4,' TOLCPA=',e10.4,/,
     .   11x,'DEPTH  =',f10.7,' IMAGZ =',f10.7,' EPS   =',f10.7,/,
     .   11x,'AMIX   =',f10.7,' EFMIX =',f10.7,' ALPCPA=',f10.7,//,
     .   11x,'NQ     =',i3,' NT    =',i3,' NL    =',i3,' NS    =',i3,
     .       ' MNTA  =',i3,/,
     .   11x,'KMSH   =  ',a,' IBZ   =',i3, ' NKX   =',i3,' NKY   =',i3,
     .       ' NKZ   =',i3,/,
     .   11x,'ZMSH   =  ',a,' NZ1   =',i3,' NZ2   =',i3,' NZ3   =',i3,
     .       ' NRES  =',i3,/,
     .   11x,'NITER  =',i3,' NLIN  =',i3,' NCPA  =',i3,' AFM   =  ',a,
     .       ' CRT   =  ',a)
  104 FORMAT(/,' INPUT:**  NX =',i3,'. Must be odd')
  105 FORMAT(/,' INPUT:**  ',a,' =',i3,' must be smaller than ',
     .         a,' =',i3)
  106 FORMAT(/,' INPUT:**  ITQ(',i2,') = ',i2,'. Must be non-zero')
  107 FORMAT(/,' INPUT:**  IT =',i3,' MMT =',i3)
  108 FORMAT(11x,'ITQ(IQ) =',10i3)
  109 FORMAT(/,11x,'Higher tails included in the charge density',
     .             ' Lmaxt =',i3)
  110 FORMAT(/,11x,'Full charge density is calculated for      ',
     .             ' Lmaxh =',i3)
  111 FORMAT(/,11x,'Soft core approximation')
  112 FORMAT(/,11x,'Force calculation included, shift = ',i3)
  113 FORMAT(/,11x,'Density of states (DOS) is calculatated')
  114 FORMAT(/,11x,'Fermi surface is calculated')
  115 FORMAT(/,' INPUT:**  CRT =',a,' not implemented.',
     .       ' Must be M or m for metals and I for insulators')
  116 FORMAT(/,' INPUT:**  DOS =',a,' not implemented.',
     .       ' Must be N, or D for DOS and S for Fermi surface')
  117 FORMAT(/,11x,'All electron calculation with fixed core')
  120 FORMAT(/,' INPUT:    The f and g equations are solved',
     .         ' self-consistently')
  121 FORMAT(/,' INPUT:    The f equation is solved self-consistently',//,
     .  11x,'The g is fixed above its non-overlapping value by ',f10.6)
  122 FORMAT(/,' INPUT:    Non-overlapping values are used for f and g')
  123 FORMAT(/,' INPUT:    Non-overlapping values is used for f',//,
     .  11x,'The g is fixed above its non-overlapping value by ',f10.6)
  124 FORMAT(/,' INPUT:    Non-overlapping values are used for f and g',
     .  //,11x,'Local muffin-tin zero is turned on')
  125 FORMAT(/,11x,'There are empty spheres in the structure')
  126 FORMAT(/,11x,'The scalar relativistic D_dot is used')
  127 FORMAT(/,11x,'WARNING: sublattice IT =',i3,
     .             ' has non-integer occupation:',/)
  128 FORMAT(11x,8f8.3)
C
      RETURN
      END
