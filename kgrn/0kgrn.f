C   ******************************************************************
C   *                                                                *
C   *                         *** KGRN ***                           *
C   *                                                                *
C   *         Self-consistent EMTO(CPA) Green's function             *
C   *                      bulk calculations                         *
C   *                                                                *
C   *                           L. Vitos                             *
C   *                  Applied Materials Physics,                    *
C   *       Department of Materials Science and Engineering          *
C   *                Royal Institute of Technology                   *
C   *                 SE-100 44 Stockholm, Sweden                    *
C   *                                                                *
C   *                        H. L. Skriver                           *
C   *                     Department of Physics                      *
C   *                Technical University of Denmark                 *
C   *                        DK-2800 Lyngby                          *
C   *                                                                *
C   ******************************************************************
      PROGRAM kgrn
      USE control_data ; USE energymesh ; USE message ; USE text
      USE control_text ; USE radialmesh ; USE volume
      IMPLICIT NONE
      CHARACTER(LEN=1) :: zmshold
      REAL(KIND=8) :: ef, efgs, hxd, hx0, vmix
      INTEGER :: nlin, nz10, isws, newsph, lmaxh0, lmaxmax
      INTEGER :: ijob
C
      version='LV 24 Oct 2013'
C
      CALL fhndlr
      CALL setcst
      ijob=INDEX(job,' ')-1
      IF(msgl.EQ.1) WRITE(msgio,'(/,5a)') ' KGRN: Job = ',
     .       job(1:ijob),',   Functional = ',func
      CALL trnsfm(1,lmaxh0)
      CALL input(efgs,lmaxh0,nlin)
      CALL atoma
      CALL setzmp
C
      nz10=nz1
      hx0=hx
C
      ALLOCATE(zfcd(nx),wfcd(nx))
C
C     Loop for volumes
C
      DO 20 isws=1,nsws
      CALL functnl
      CALL nxtsws(isws,efgs,hx)
      IF(conv.EQ.'N') THEN
         nz1=nz0
         hxd=hx0
         newsph=0
         vmix=0.1d0
      ELSE
         nz1=nz10
         hxd=hx
         newsph=1
         vmix=1.d0
      ENDIF
      CALL setdrc(dx)
      CALL setpsn(dx)
      CALL setxcp(ixc,txch,exchf,ns)
      lmaxmax=MAX0(lmaxt,lmaxh)+1
      CALL bessinit(-1,lmaxmax)
      CALL start(efgs)
      CALL wavefc(0,efgs)
      CALL ebtop(1)
      CALL setkms
      CALL zmesh(efgs,1)
      efl(isws)=efgs
      ef=efgs
C
C     Generate the k-dependent slope matrix and store it if stmp = Y
C     Test for the Taylor expansion
C
      CALL blochs
C
C     Self-consistent loops
C
      CALL slfcns(isws,nlin,nz10,newsph,ef,vmix,hxd,lmaxh0)
C
C     Save potential for restart
C
      CALL update(ef)
C
C     Calculate density of state
C
      IF(dosc.EQ.'D') THEN
         DEALLOCATE(zm,wgm,zx)
         zmshold=zmsh
         zmsh='D'
         CALL zmesh(ef,1)
         CALL arrayz(1)
         CALL kkrdos(ef)
         CALL arrayz(0)
      ELSEIF(dosc.EQ.'S') THEN
         DEALLOCATE(zm,wgm,zx)
         zmshold=zmsh
         zmsh='S'
         CALL zmesh(ef,1)
         CALL arrayz(1)
         CALL kkrfes(ef)
         CALL arrayz(0)
      ENDIF
C
C     Calculate the full non-spherical charge density and forces
C
      CALL arrays(0)
      IF(fcd.NE.'Y') GO TO 20
      WRITE(m6,'(/,a,i3)') 
     .  ' Calculation of the full charge density, lmaxh =',lmaxh
      IF(msgl.NE.0) WRITE(msgio,'(/,a,i3)') 
     .  ' Calculation of the full charge density, lmaxh =',lmaxh
      IF(dosc.NE.'N') THEN
         DEALLOCATE(zm,wgm,zx)
         zmsh=zmshold
         CALL zmesh(ef,1)
      ENDIF
C
C     Complete the zmesh close to the Fermi level
C
      ALLOCATE(zm0(nzm),wgm0(nzm))
      zm0=zm
      wgm0=wgm
      DEALLOCATE(zm,wgm)
      nzm=nzm+nfcd
      ALLOCATE(zm(nzm),wgm(nzm))
      zm(1:nzm-nfcd)=zm0
      wgm(1:nzm-nfcd)=wgm0
      zm(nzm-nfcd+1:nzm)=zfcd(1:nfcd)
      wgm(nzm-nfcd+1:nzm)=wfcd(1:nfcd)
      DEALLOCATE(zfcd,wfcd)
C
      CALL arrayz(1)
      CALL kkrfcd(ef,0)
      CALL eldens
C
C     Force calculation
C
      IF(frc.EQ.'Y') CALL gtforce
C
C     Calculate SCA density
C
      WRITE(m6,'(/,a,i3)') 
     .  ' Calculation of the SCA full charge density, lmax =',lmax
      IF(msgl.NE.0) WRITE(msgio,'(/,a,i3)') 
     .  ' Calculation of the SCA full charge density, lmax =',lmax
C
      CALL kkrfcd(ef,1)
      CALL eldens
C
C     Save full-charge density and kinetic force components
C
   30 CONTINUE
      CALL fcdupd(ef)
C
   20 CONTINUE
      CALL htimer(dato,clock)
      IF(msgl.NE.0) 
     .   WRITE(msgio,'(/,4a)') ' KGRN calculation finished at ',
     .              clock,' on ',dato
      WRITE(m6,'(/,4a)') ' KGRN calculation finished at ',
     .              clock,' on ',dato
C
      STOP
      END
