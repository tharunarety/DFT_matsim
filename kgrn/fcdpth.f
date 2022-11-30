      SUBROUTINE fcdpth(is,gah,bg0,gi,dtil,nzm)
C   ******************************************************************
C   *                                                                *
C   * 1.Calculate the path operator g^a(z) defined as:               *
C   *                                                                *
C   *   K^a(z) * g^a(z) = a*[S^a - D^a] * g^a(z) = 1                 *
C   *                                                                *
C   * 2.Calculate the density of state G(z) defined as:              *
C   *                                                                *
C   *   G(z) = g^a(z) * K^dot^a(z)                                   *
C   *                                                                *
C   *   and integrate over the 3D BZ.                                *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE bzmesh ; USE control_data ; USE control_text
      USE kinkmatrix   ; USE logderivative ; USE message
      USE partialwaves ; USE potential     ; USE slope
      USE taylor
      IMPLICIT NONE
      INTEGER :: is, i, lk, jd, lz, nzm, iq, iq0, it, ita, ntait
      INTEGER :: icon, icont, convt, doh, minit
      INTEGER :: l, m, lm, lmp, j, ndera, twol
      COMPLEX(KIND=8), DIMENSION(nzm,nq,ns,nlmf,nlmf) :: gah
      COMPLEX(KIND=8), DIMENSION(nzm,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: bg
      COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ga
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,ns) :: bg0
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: dtil
      COMPLEX(KIND=8), DIMENSION(nlm,nlm) :: ma, mb, mc
      COMPLEX(KIND=8) :: diag, ty, tyd, prod
      REAL(KIND=8)    :: pkx, pky, pkz, weight, e, ah
      REAL(KIND=8)    :: rconv, iconv
      INTEGER, DIMENSION(:), ALLOCATABLE :: conv
C
      ALLOCATE(sa(nlmqf,nlmq,0:nder))
      ALLOCATE(sal(nlmq,nlmq,0:nder))
      IF(expan.EQ.'D') THEN
         ALLOCATE(sap(nlmqf,nlmq,0:nderp))
         ALLOCATE(sapl(nlmq,nlmq,0:nderp))
      ENDIF
      ALLOCATE(kinkmd(nlmq,nlmq))
      ALLOCATE(bg(nzm,nq,nlm,nlm),ga(nzm,nq,nlm,nlm,ns))
      ALLOCATE(conv(nzm))
C
C     Solve CPA equation
C
      doh=0
      conv(1:nzm)=0
      icpa=0
   10 icpa=icpa+1
C
      DO 11 lz=1,nzm
      IF(conv(lz).EQ.0) THEN
         bg(lz,1:nq,1:nlm,1:nlm)=zero
         ga(lz,1:nq,1:nlm,1:nlm,is)=zero
      ENDIF
   11 CONTINUE
C
      DO 20 lk=1,nkvec
      pkx   =fkx(lk)
      pky   =fky(lk)
      pkz   =fkz(lk)
      weight=fkw(lk)
C
C     Get the Bloch transformed slope matrices
C
      CALL smtrx(nlmf,pkx,pky,pkz)
      CALL mslope(nlmf)
C
      DO 30 lz=1,nzm
      IF(conv(lz).GT.0) GO TO 30
C
C     Taylor expansion for the slope matrix S^a and S^a_dot
C
      IF(clkw0(lz).EQ.1) THEN
         kinkm(1:nlmq,1:nlmq)=sal(1:nlmq,1:nlmq,0)
         ndera=nder
      ELSE
         kinkm(1:nlmq,1:nlmq)=sapl(1:nlmq,1:nlmq,0)
         ndera=nderp
      ENDIF
      kinkmd=zero
      DO 31 jd=1,ndera
      ty=tayl(lz,jd,is)
      tyd=tayld(lz,jd,is)
      IF(clkw0(lz).EQ.1) THEN
         kinkm(1:nlmq,1:nlmq)=kinkm(1:nlmq,1:nlmq)+
     .                        ty*sal(1:nlmq,1:nlmq,jd)
         kinkmd(1:nlmq,1:nlmq)=kinkmd(1:nlmq,1:nlmq)+
     .                         tyd*sal(1:nlmq,1:nlmq,jd)
      ELSE
         kinkm(1:nlmq,1:nlmq)=kinkm(1:nlmq,1:nlmq)+
     .                        ty*sapl(1:nlmq,1:nlmq,jd)
         kinkmd(1:nlmq,1:nlmq)=kinkmd(1:nlmq,1:nlmq)+
     .                         tyd*sapl(1:nlmq,1:nlmq,jd)
      ENDIF
   31 CONTINUE
C
C     kinkm = a*S^a - a*D^a
C
      DO 32 i=1,nlmq
      ah=asc(i)
      kinkm (i,1:nlmq)=ah*kinkm (i,1:nlmq) 
   32 kinkmd(i,1:nlmq)=ah*kinkmd(i,1:nlmq) 
C
C     kinkm = S^a - D^a
C
      DO 33 iq=1,nq
      iq0=(iq-1)*nlm
      DO 33 lm=1,nlm
      i=lm+iq0
      DO 33 lmp=1,nlm
      j=lmp+iq0
   33 kinkm(i,j)=kinkm(i,j)-dtil(lz,iq,lm,lmp,is)
C
      CALL mtrxinv(nlmq,kinkm,unit,gak,info)
      IF(info.NE.0) STOP 'FCDPTH-1: INFO <> 0'
C
      gak=weight*gak
      IF(doh.EQ.1) THEN
         IF(clkw0(lz).EQ.1) THEN
            CALL hghfcd(is,gah,sa,nder,lz,nzm)
         ELSE
            CALL hghfcd(is,gah,sap,nderp,lz,nzm)
         ENDIF
      ENDIF
C
C     Find the path operator and the density of states
C
      DO 34 iq=1,nq
      iq0=(iq-1)*nlm
      DO 34 lm=1,nlm
      i=lm+iq0
      DO 34 lmp=1,nlm
      j=lmp+iq0
      ga(lz,iq,lm,lmp,is)=ga(lz,iq,lm,lmp,is)+gak(i,j)
C
      diag=SUM(gak(i,1:nlmq)*kinkmd(1:nlmq,j))
      bg(lz,iq,lm,lmp)=bg(lz,iq,lm,lmp)+diag
C
   34 CONTINUE
   30 CONTINUE
   20 CONTINUE
C
C     Rotate the path operator and the density of states
C     obtained in the IBZ for the FBZ
C
      IF(fllbz.NE.'Y') THEN
         DO 40 lz=1,nzm
         IF(conv(lz).GT.0) GO TO 40
         CALL intgf(ga(lz,1:nq,1:nlm,1:nlm,is),
     .             bg(lz,1:nq,1:nlm,1:nlm))
   40    CONTINUE
      ENDIF
C
C     Eq. 1b Set up the Green's functions for the alloy components
C
      CALL dyson(is,ga,gi,dtil,nzm,nzm,0)
C
C     Eq. 1c
C
      convt=0
      DO 50 lz=1,nzm
      IF(conv(lz).GT.0) GO TO 50
      icont=0
      DO 51 iq=1,nq
      it=itq(iq)
      ntait=nta(it)
      icon=0
C
C     Find new coherent potential function for sublattice iq
C
      IF(ntait.GT.1) THEN
         ma=zero
         DO 52 ita=1,ntait
C
         ma(1:nlm,1:nlm)=ma(1:nlm,1:nlm)+
     .                   conc(ita,it)*gi(lz,ita,iq,1:nlm,1:nlm,is)
   52    CONTINUE
C
C        [sum(sort) c(sort) gi(z)(sort)]^-1
C
         CALL mtrxinv(nlm,ma,unil,mb,info)
         IF(info.NE.0) STOP 'FCDPTH-2: INFO <> 0'
         ma(1:nlm,1:nlm)=ga(lz,iq,1:nlm,1:nlm,is)
C
C        g0(z)^-1
C
         CALL mtrxinv(nlm,ma,unil,mc,info)
         IF(info.NE.0) STOP 'FCDPTH-3: INFO <> 0'
C
C        New coherent potential function (100 % mixing)
C
         ma(1:nlm,1:nlm)=mb(1:nlm,1:nlm)-mc(1:nlm,1:nlm)
         dtil(lz,iq,1:nlm,1:nlm,is)=dtil(lz,iq,1:nlm,1:nlm,is)-
     .                              ma(1:nlm,1:nlm)
C
C        Test for the convergency
C
         DO 53 lm=1,nlm
         DO 53 lmp=1,nlm
         prod=ma(lm,lmp)
         rconv=ABS(REAL(prod,8))
         iconv=ABS(AIMAG(prod))
C
         IF(lm.EQ.lmp) THEN
            rconv=rconv/ABS(REAL(dtil(lz,iq,lm,lmp,is),8))
            iconv=iconv/ABS(AIMAG(dtil(lz,iq,lm,lmp,is)))
         ENDIF
         IF(rconv.LE.tolcpa.AND.iconv.LE.tolcpa) icon=icon+1
   53    CONTINUE
c         IF(msgl.NE.0) WRITE(msgio,100) icpa,lz,icon,nlm*nlm
      ELSE
         icon=nlm*nlm
      ENDIF
      icont=icont+nlm*nlm-icon
   51 CONTINUE
C
      IF(icont.EQ.0) conv(lz)=icpa
      convt=convt+icont
   50 CONTINUE
C
      IF(doh.EQ.0) THEN
         IF(convt.GT.0) THEN
            IF(icpa.LT.ncpa) GO TO 10
            minit=ncpa
            DO 54 lz=1,nzm
            IF(conv(lz).NE.0.AND.conv(lz).LT.minit) minit=conv(lz)
   54       CONTINUE
            IF(mnta.GT.1) THEN
            IF(msgl.NE.0) WRITE(msgio,110) is,minit,ncpa
            WRITE(m6,110) is,minit,ncpa
            ENDIF
         ELSE
            IF(mnta.GT.1) THEN
            IF(msgl.NE.0) WRITE(msgio,120) is,MINVAL(conv),MAXVAL(conv)
            WRITE(m6,120) is,MINVAL(conv),MAXVAL(conv)
            ENDIF
         ENDIF
C
         doh=1
         conv(1:nzm)=0
         GO TO 10
      ENDIF
C
C     Save the matrices
C
      DO 60 iq=1,nq
      DO 60 lm=1,nlm
   60 bg0(1:nzm,iq,lm,is)=bg(1:nzm,iq,lm,lm)
C
      DEALLOCATE(sa,sal)
      IF(expan.EQ.'D') DEALLOCATE(sap,sapl)
      DEALLOCATE(bg,kinkmd,conv,ga)
C
  100 FORMAT(/,' PATHOP: ICPA =',i3,' LZ =',i3,' ICON =',i3,
     .         ' ITOT =',i3)
  110 FORMAT(/,' PATHOP: CPA equation for IS =',i2,' not converged',
     .   /,9x,'number of iterations performed ',i2,' to ',i2)
  120 FORMAT(/,' FCDPTH: CPA equation for IS =',i3,
     .         ' converged in ',i2,'-',i2,' iterations')
      RETURN
      END
