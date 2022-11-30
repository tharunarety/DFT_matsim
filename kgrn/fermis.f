      SUBROUTINE fermis
C   ******************************************************************
C   *                                                                *
C   * Set up the k-mesh for cross section in BZ.                     *
C   * Units of pi/a are used.                                        *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE control_text ; USE bzmesh ; USE message
      IMPLICIT NONE
      REAL(KIND=8)  :: tw, pkx, pky, pkz, weight
      INTEGER       :: lx, ly, lz, lk
      fs=1
      IF(ibz.EQ.4) fs=2
C
      IF(fsts.EQ.1) THEN
         nfs=nkvec
         GO TO 90
      ELSE
         IF(stmp.NE.'Y') THEN
            DEALLOCATE(fsx,fsy,fsz,fsw)
         ELSE
            DEALLOCATE(kx,ky,kz,ww)
         ENDIF
      ENDIF
C
C     Here we should write data for cross-sections for each BZ
C     There are nfs points and they are stored in a similar way
C     as in kmesh: 
C     kx, ky, kz points (integers) and dkx, dky, dkz, dhx steps
C    
      IF(ibz.LE.3) THEN
C
C        fcc (010)+(110)
C
         nky=400
         nkx=2*nky-1
         nkz=nky
         nfs=nkz*nkx
         ALLOCATE(kx(nfs),ky(nfs),kz(nfs),ww(nfs))
C
         tw=2.d0
         IF(ibz.EQ.1) tw=1.d0
         dkx=2.d0*tw/(nkx-1)
         dky=tw/(nky-1)
         dkz=tw/(nkz-1)
C
         lk=0
         DO 10 lz=1,nkz
         DO 11 lx=-(nky-1),nky-1
         lk=lk+1
         IF(lx.LE.0) THEN
C
C           (010)
C
            kx(lk)=lx
            ky(lk)=0
            kz(lk)=lz-1
         ELSE
C
C           (110)
C
            kx(lk)=lx
            ky(lk)=lx
            kz(lk)=lz-1
         ENDIF
         ww(lk)=1.d0
   11    CONTINUE
   10    CONTINUE
      ELSEIF(ibz.EQ.4) THEN
         nkx=100
         nky=nkx
         nkz=nkx
         IF(fs.EQ.1) THEN
            nfs=nkz*(2*nkx-1)
         ELSEIF(fs.EQ.2) THEN
            nfs=(2*nkx-1)*(2*nky-1)
         ENDIF
         ALLOCATE(kx(nfs),ky(nfs),kz(nfs),ww(nfs))
C
         dkx=4.d0/3.d0/(nkx-1)
         IF(fs.EQ.1) THEN
            dky=4.d0/3.d0/SQRT(3.d0)/(nky-1)
         ELSEIF(fs.EQ.2) THEN
            dky=2.d0/SQRT(3.d0)/(nky-1)
         ENDIF
         dkz=1.d0/coa/(nkz-1)
         dhx=0.d0
C
         lk=0
         IF(fs.EQ.1) THEN
C
C           hcp GKHA+GMLA
C
            DO 20 lz=1,nkz
            DO 21 lx=-(nkx-1),nkx-1
            lk=lk+1
            IF(lx.LE.0) THEN
C
C              GKA
C
               kx(lk)=lx
               ky(lk)=0
               kz(lk)=lz-1
            ELSE
C
C              GMA
C
               kx(lk)=lx
               ky(lk)=lx
               kz(lk)=lz-1
            ENDIF
            ww(lk)=1.d0
   21       CONTINUE
   20       CONTINUE
         ELSEIF(fs.EQ.2) THEN
C
C           hcp GKM
C
            DO 22 ly=-(nky-1),nky-1
            DO 23 lx=-(nkx-1),nkx-1
            lk=lk+1
            kx(lk)=lx
            ky(lk)=ly
            kz(lk)=0
            ww(lk)=1.d0
   23       CONTINUE
   22       CONTINUE
         ENDIF
      ELSE
         WRITE(m6,100)
         STOP
      ENDIF
C
      GO TO 90
C
   89 CONTINUE
      WRITE(m6,110) nky
      STOP
C
   90 CONTINUE
      ALLOCATE(fsx(nfs),fsy(nfs),fsz(nfs),fsw(nfs))
      DO 91 lk=1,nfs
      IF(nkvec.GT.1) THEN
         CALL kvec(lk,pkx,pky,pkz,weight)
      ELSE
C
C        Gamma point
C
         pkx=0.d0
         pky=0.d0
         pkz=0.d0
         weight=1.d0
      ENDIF
      fsx(lk)=pkx
      fsy(lk)=pky
      fsz(lk)=pkz
      fsw(lk)=weight
   91 CONTINUE
      IF(fsts.EQ.1.AND.nkvec.GT.1) DEALLOCATE(kx,ky,kz,ww)
C
  100 FORMAT(/,' FERMIS: no cross section implemented for IBZ =',i3)
  110 FORMAT(/,' FERMIS: too small NKY =',i3)
      RETURN
      END
