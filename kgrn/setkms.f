      SUBROUTINE setkms
C   ******************************************************************
C   *                                                                *
C   *    Construct the k-mesh to be used in the calculation of       *
C   *    the Green's function.                                       *
C   *                                                                *
C   ******************************************************************
      USE bzmesh ; USE bzmesh2 ; USE control_data
      USE control_text ; USE message
      IMPLICIT NONE
      INTEGER, PARAMETER :: mkm = 100000, mpar = 100, mper = 101
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fbx, fby, fbz, fbw
      REAL(KIND=8) :: wtot, sumw, pkx, pky, pkz, weight
      INTEGER :: lk, lkt, fblk, fbnkvec, sumnk
C
      IF(mode.EQ.'3D') THEN
         sumnk=nkx+nky+nkz
         IF(sumnk.EQ.1) GO TO 22
C
         ALLOCATE(kx(mkm),ky(mkm),kz(mkm),ww(mkm))
         ALLOCATE(fbx(48),fby(48),fbz(48),fbw(48))
C
         CALL kmesh(ibz,nkx,nky,nkz,nkvec,dkx,dky,dkz,dhx,
     .        boa,coa,alf,bet,gam,m6,mkm,tkx,tky,tkz,kx,ky,kz,ww)
         WRITE(m6,106) nkvec,nkx,nky,nkz,dkx,dky,dkz
         IF(msgl.NE.0) WRITE(msgio,106) nkvec,nkx,nky,nkz,dkx,dky,dkz
C
C        Normalize k-point weights
C
         wtot=SUM(ww(1:nkvec))
         ww(1:nkvec)=ww(1:nkvec)/wtot
C
         IF(fllbz.NE.'Y') CALL setibz
         lkt=0
         DO 20 lk=1,nkvec
         CALL kvec(lk,pkx,pky,pkz,weight)
         CALL fullbz(pkx,pky,pkz,weight,fbx,fby,fbz,fbw,fbnkvec)
         DO 20 fblk=1,fbnkvec
   20    lkt=lkt+1
C
         ALLOCATE(fkx(lkt),fky(lkt),fkz(lkt),fkw(lkt))
C
         lkt=0
         sumw=0.d0
         DO 21 lk=1,nkvec
         CALL kvec(lk,pkx,pky,pkz,weight)
         CALL fullbz(pkx,pky,pkz,weight,fbx,fby,fbz,fbw,fbnkvec)
         DO 21 fblk=1,fbnkvec
         lkt=lkt+1
         fkx(lkt)=fbx(fblk)
         fky(lkt)=fby(fblk)
         fkz(lkt)=fbz(fblk)
         fkw(lkt)=fbw(fblk)
   21    sumw=sumw+fkw(lkt)
         fbnkvec=lkt
         IF(ABS(sumw-1.d0).GT.1.d-8) THEN
            WRITE(m6,102) sumw
            STOP
         ENDIF
C
C        Select the non-equivalent k-vectors
C
         CALL nekvec(fbnkvec)
C
         IF(nkvec.LT.fbnkvec.AND.fllbz.NE.'Y') WRITE(m6,100) fbnkvec
         nkvec=fbnkvec
         IF(fllbz.EQ.'Y') WRITE(m6,101) nkvec
         IF(msgl.NE.0) WRITE(msgio,103) nkvec
C
         DEALLOCATE(fbx,fby,fbz,fbw)
         RETURN
C
C        Gamma point calculation for molecules
C
  22     CONTINUE
         fllbz='Y'
         lkt=1
         ALLOCATE(fkx(lkt),fky(lkt),fkz(lkt),fkw(lkt))
         nkvec=1
         fkx(1)=0.d0
         fky(1)=0.d0
         fkz(1)=0.d0
         fkw(1)=1.d0
         IF(fllbz.EQ.'Y') WRITE(m6,101) nkvec
         WRITE(m6,104)
      ELSEIF(mode.EQ.'2D') THEN
         ALLOCATE(akx(mpar),aky(mpar),akz(mper),wk(mpar))
         CALL primkv(nq)
         CALL kmesh2(ibz2,nkx2,nky2,nkz2,nkvec,mpar,mper)
C
         ALLOCATE(fkx(nkvec),fky(nkvec),fkz(nkvec),fkw(nkvec))
C
         DO 30 lk=1,nkvec
         CALL kvec(lk,pkx,pky,pkz,weight)
         fkx(lk)=pkx
         fky(lk)=pky
         fkz(lk)=pkz
   30    fkw(lk)=weight
         DEALLOCATE(akx,aky,akz,wk)
      ELSE
         WRITE(m6,105) mode
         STOP
      ENDIF
C
  100 FORMAT(/,' SETKMS:   IBZ increased, NKVEC = ',i4)
  101 FORMAT(/,' SETKMS:   Integration performed in the full-BZ,',
     .                ' NKVEC = ',i4)
  102 FORMAT(/,' SETKMS:**  SUMW = ',f12.8,' it should be 1.')
  103 FORMAT(/,' SETKMS:   NKVEC = ',i4)
  104 FORMAT(/,11x,'Gamma point calculation')
  105 FORMAT(/,' SETKMS:**  MODE = ',a,'. Must be 2D or 3D.')
  106 FORMAT(/,' KMESH:    NKVEC =',i6,' NKX =',i3,' NKY =',i3,
     .       ' NKZ =',i3,/,
     .       11x,'DKX =',f10.6,' DKY =',f10.6,' DKZ =',f10.6)
      RETURN
      END
