      SUBROUTINE slfcns(isws,nlin,nz10,newsph,ef,vmix,hxd,lmaxh0)
C   ******************************************************************
C   *                                                                *
C   *  Self-consistent loops                                         *
C   *                                                                *
C   ******************************************************************
      USE bzmesh ; USE control_data ; USE control_text ; USE energymesh 
      USE message      ; USE potential ; USE text
      use radialmesh
      IMPLICIT NONE
      CHARACTER(LEN=1) :: convrg, zmshold
      REAL(KIND=8) :: ef, vmix, hxd, mixcore = 1.d0, tef = 0.01d0
      REAL(KIND=8) :: ei, efo, efi, efn, erren, erref, etotal, delef
      INTEGER :: irpt = 0, prnt = 0, on = 1, off = 0
      INTEGER :: isws, nlin, nz10, newsph, contl
      INTEGER :: ncore, lin, mixkey, nzlin, lmaxh0
C
      ei=0.d0
      efi=0.d0
      mixkey=ns
      ncore=1
      convrg='N'
C
      DO 20 iter=1,niter
      CALL arrayz(on)
      CALL slfarr(on)
C
C     Construct Green's function on the contour and integrate over k.
C
      efn=ef
      CALL kkreq(ef)
      erref=ABS(efi-ef)
      efi=ef
C
C     Get total energy within SCA
C
      CALL totale(etotal)
      erren=ABS(ei-etotal)
      ei=etotal
      IF(iter.EQ.1) CALL prnprm(off,ef,etotal)
C
C     Test for convergence
C
      IF(crt.EQ.'M'.OR.crt.EQ.'m') THEN
         IF(erren.LT.tole.AND.erref.LT.tolef) convrg='Y'
      ELSEIF(crt.EQ.'I') THEN
         IF(erren.LT.tole) convrg='Y'
      ENDIF
      IF(convrg.EQ.'Y') THEN
         IF(irpt.EQ.0) THEN
            irpt=1
            convrg='N'
         ENDIF
      ELSE
         irpt=0
      ENDIF
C
      IF(convrg.EQ.'Y') THEN
         ef=efn
C
C        Set up the proper potential spheres and continue iterations
C
         IF(ops.EQ.'Y'.AND.newsph.EQ.0) THEN
            newsph=1
            CALL potsph
            GO TO 24
         ENDIF
         CALL totale(etotal)
         CALL wavefc(on,ef)
         CALL ebtop(on)
         CALL slfarr(off)
C
         IF((dosc.NE.'N'.OR.fcd.EQ.'Y').AND.stmp.EQ.'Y') 
     .       CALL trnsfm(2,lmaxh0)
C
C        DOS
C
         DEALLOCATE(zm,wgm,zx)
         zmshold=zmsh
         zmsh='S'
         fsts=1
         CALL zmesh(ef,off)
         CALL kkrfes(ef)
         zmsh=zmshold
         fsts=0
C
C        Alternative DOS (used for Hopfield)
C
         DEALLOCATE(zm,wgm,zx)
         zmshold=zmsh
         zmsh='T'
         CALL zmesh(ef,on)
         CALL dosef(ef)
         zmsh=zmshold
C
C        Entropy
C
         IF(zmsh.EQ.'F'.OR.zmsh.EQ.'f') THEN
            DEALLOCATE(zm,wgm,zx)
            zmshold=zmsh
            zmsh='T'
            CALL zmesh(ef,on)
            CALL entrop(ef)
            zmsh=zmshold
         ENDIF
C
         DEALLOCATE(zm,wgm,zx)
         CALL zmesh(ef,on)
         conv='Y'
         CALL prnprm(off,ef,etotal)
         CALL htimer(dato,clock)
         WRITE(m6,100) iter,clock,dato
         IF(msgl.NE.0) WRITE(msgio,100) iter,clock,dato
         CALL arrayz(off)
         GO TO 22
      ENDIF
   24 CONTINUE
C
C     Iterate on a fixed contour using the k-integrated Green's function
C
      IF(nlin.GT.0) THEN
         nzlin=MAX0(nzm,nx)
         CALL linarr(on,lmax,nzlin)
         efo=ef
         DO 21 lin=1,nlin
         IF(softc.EQ.'Y') THEN
            CALL softp(mixcore,ncore)
         ENDIF
         CALL newchd(prnt)
         CALL mltpm(lin)
         CALL mixchd(amix,mixkey,on)
         CALL chdren(ef,lin)
         CALL newpot(prnt,lin)
         IF(fixvmtz.LT.2) CALL optpot(fixvmtz,off,ef,vmix)
         CALL kkreq0(ef,nzlin,on,contl)
         delef=ABS(efo-ef)*1000.d0
         IF(delef.LT.tolef.AND.ncore.EQ.1) GO TO 23
         efo=ef
         IF(contl.NE.0) GO TO 25
   21    CONTINUE
   25    WRITE(m6,110) ' Linear loop not converged'
   23    CONTINUE
         CALL linarr(off,lmax,nzlin)
      ELSE
         lin=1
         contl=0
         IF(softc.EQ.'Y') THEN
            CALL softp(mixcore,ncore)
         ENDIF
         CALL newchd(prnt)
         CALL mltpm(lin)
         CALL mixchd(amix,mixkey,lin)
         CALL chdren(ef,lin)
         CALL newpot(prnt,lin)
      ENDIF
C
C     Prepare parameters for the next iteration
C
      nz1=MIN0(nz1+2,nz10)
      IF(contl.EQ.0) THEN
         IF(crt.EQ.'m') THEN
            hxd=hxd/2.d0
            hx=MAX(hxd,0.0001d0)
         ELSEIF(crt.EQ.'M') THEN
            IF(iter.LE.10) THEN
               hxd=hxd-0.001d0
            ELSE
               hxd=hxd/2.d0
            ENDIF
            hx=MAX(hxd,0.001d0)
         ELSEIF(crt.EQ.'I') THEN
            hx=hxd
            IF(iter.GT.3.AND.erref.LT.tef) efmix=0.d0
         ENDIF
         ef=(1.d0-efmix)*efi+efmix*ef
      ENDIF
      DEALLOCATE(zm,wgm,zx)
      CALL zmesh(ef,off)
C
      vmix=MIN(2.d0*vmix,1.d0)
      CALL optpot(fixvmtz,on,ef,vmix)
C
      CALL wavefc(on,ef)
      CALL wmsg(ef,etotal,erren,erref,lin,isws)
C
      CALL arrayz(off)
      CALL slfarr(off)
C
      CALL update(ef)
   20 CONTINUE
C
C     Print non-converged results for information
C
      CALL totale(etotal)
      CALL wavefc(on,ef)
      CALL ebtop(on)
      CALL prnprm(off,ef,etotal)
      conv='N'
      WRITE(m6,110) ' Not converged'
      IF(msgl.NE.0) WRITE(msgio,110) ' Not converged'
C
      IF((dosc.NE.'N'.OR.fcd.EQ.'Y').AND.stmp.EQ.'Y') 
     .    CALL trnsfm(2,lmaxh0)
C
   22 CONTINUE
C
      CLOSE(11,STATUS='DELETE')
C
  100 FORMAT(/,' Converged in',i3,' iterations at ',a,2x,a)
  110 FORMAT(/,a)
      RETURN
      END
