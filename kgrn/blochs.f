      SUBROUTINE blochs
C   ******************************************************************
C   *                                                                *
C   *    Calculate the Bloch transformed slope matrix                *
C   *    and store it in a temporary file.                           *
C   *                                                                *
C   *    Test for the accuracy of the Taylor expansion within        *
C   *    the actual energy window (Gamma point).                     *
C   *                                                                *
C   ******************************************************************
      USE bzmesh       ; USE control_data  ; USE control_text
      USE energymesh   ; USE message       ; USE potential
      USE radialmesh   ; USE slope         ; USE taylor
      IMPLICIT NONE
      INTEGER :: lk, case, jd, lz, i, j, clkw0lz, dt
      INTEGER :: error, k, gkp
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ss, ssp
      REAL(KIND=8) :: ss1, ss2
      REAL(KIND=8), PARAMETER :: tol =1.d-10
      REAL(KIND=8) :: avmtz, zs2, enrg, er1, er2, err
      REAL(KIND=8) :: pkx, pky, pkz, pka
C
      IF(stmp.EQ.'Y') THEN
         OPEN(11,FILE=for011,STATUS='UNKNOWN',FORM='UNFORMATTED',
     .                           IOSTAT=error)
         IF(error.NE.0) THEN
            WRITE(m6,*) 'File 11 could not be opened'
            STOP
         ENDIF
      ENDIF
C
      ALLOCATE(sa(nlmqt,nlmq,0:nder),STAT=error)
      IF(error.NE.0) THEN
         WRITE(m6,*) 'Space for Sa could not be allocated'
         STOP
      ENDIF
      ALLOCATE(sa0(0:nder),ss(0:nder))
      IF(expan.EQ.'D') THEN
         ALLOCATE(sap(nlmqt,nlmq,0:nderp),STAT=error)
         IF(error.NE.0) THEN
            WRITE(m6,*) 'Space for Sap could not be allocated'
            STOP
         ENDIF
         ALLOCATE(sap0(0:nderp),ssp(0:nderp))
      ENDIF
C
      dt=1
      gkp=0
      DO 20 lk=1,nkvec
      pkx=fkx(lk)
      pky=fky(lk)
      pkz=fkz(lk)
      pka=pkx*pkx+pky*pky+pkz*pkz
      IF(pka.LT.1.d-8.OR.stmp.EQ.'Y') THEN
         CALL smtrx(nlmt,pkx,pky,pkz)
      ENDIF
C
      IF(stmp.EQ.'Y') THEN
         DO 21 i=1,nlmqt
         DO 21 j=1,nlmq
         WRITE(11) sa(i,j,0:nder)
         IF(expan.EQ.'D') WRITE(11) sap(i,j,0:nderp)
   21    CONTINUE
      ENDIF
C
C     Test for the Taylor expansion for the diagonal (dt) terms
C     in Gamma (k=0) point
C
      IF(pka.LT.1.d-8) THEN
         gkp=gkp+1
C
C        Save diagonal (dt) terms for latter use
C
         sa0(0:nder)=sa(dt,dt,0:nder)
         IF(expan.EQ.'D') sap0(0:nderp)=sap(dt,dt,0:nderp)
C
         avmtz=SUM(vmtz(1:ns))/ns
         case=1
         DO 22 lz=1,nzm
         enrg=REAL(zm(lz),8)
         zs2=(enrg-avmtz)*sws*sws 
         IF(expan.EQ.'S') THEN
            zs2=zs2-kw20
            ss(0)=REAL(sa0(0),8)
            DO 23 jd=1,nder
   23       ss(jd)=REAL(sa0(jd),8)*zs2**jd/facd(jd)
            ss1=SUM(ss(0:nder))
C
            IF(ABS(ss1).GT.tol) THEN
               err=ABS(ss(nder)/ss1)
            ELSE
               err=ABS(ss(nder))
            ENDIF
C
            clkw0lz=1
         ELSEIF(expan.EQ.'D') THEN
            zs2=zs2-kw20
            ss(0)=REAL(sa0(0),8)
            DO 24 jd=1,nder
   24       ss(jd)=REAL(sa0(jd),8)*zs2**jd/facd(jd)
            ss1=SUM(ss(0:nder))
            IF(ABS(ss1).GT.tol) THEN
               er1=ABS(ss(nder)/ss1)
            ELSE
               er1=ABS(ss(nder))
            ENDIF
C
            zs2=zs2-kw20p
            ssp(0)=REAL(sap0(0),8)
            DO 25 jd=1,nderp
   25       ssp(jd)=REAL(sap0(jd),8)*zs2**jd/facd(jd)
            ss2=SUM(ssp(0:nderp))
            IF(ABS(ss2).GT.tol) THEN
               er2=ABS(ssp(nderp)/ss2)
            ELSE
               er2=ABS(ssp(nderp))
            ENDIF
C
C           Choose the more accurate one
C
            IF(er1.LT.er2) THEN
               err=er1
               clkw0lz=1
            ELSE
               err=er2
               clkw0lz=2
            ENDIF 
         ENDIF
C
         IF(err.GT.1.d-4) THEN
            IF(case.EQ.1) WRITE(m6,100)
            case=case+1
            IF(clkw0lz.EQ.1) THEN
               WRITE(m6,101) lz,clkw0lz,enrg,ss1,ss(1:nder)
            ELSE
               WRITE(m6,101) lz,clkw0lz,enrg,ss2,ssp(1:nderp)
            ENDIF
         ENDIF
C
   22    CONTINUE
      ENDIF
C
   20 CONTINUE
      DEALLOCATE(sa,ss)
      IF(expan.EQ.'D') DEALLOCATE(sap,ssp)
C
      IF(stmp.EQ.'Y') THEN
         DEALLOCATE(salpl,STAT=error)
         IF(error.NE.0) THEN
            WRITE(m6,*) 'Space for Salpl could not be deallocated'
            STOP
         ENDIF
         IF(expan.EQ.'D') THEN
            DEALLOCATE(salplp,STAT=error)
            IF(error.NE.0) THEN
               WRITE(m6,*) 'Space for Salplp could not be deallocated'
               STOP
            ENDIF
         ENDIF
         REWIND(11)
      ENDIF
      IF(gkp.NE.1) THEN
         WRITE(m6,102) gkp
         STOP
      ENDIF
C
  100 FORMAT(/,' WARNING:  Taylor expansion might not be converged:',
     .       //,' LZ',' exp',3x,'E',5x,'S00',5x,'S00_Dots(1,2,...)',/)
  101 FORMAT(i3,1x,i2,1x,f6.2,30f8.4)
  102 FORMAT(/,' WARNING:  No Gamma point in k-mesh, gkp= ',i3)
      RETURN
      END
