      SUBROUTINE fespth(is,bg0,gi,dtil,zm,wgm,nzm)
C   ******************************************************************
C   *                                                                *
C   *    1. Calculate the path operator g(z) and coherent potential  *
C   *       function D(z) by solving the CPA equations:              *
C   *                                                                *
C   *       a)  K(z) * g(z) = [a*S(z) - D(z)] * g(z) = 1             *
C   *                                                                *
C   *       b)  gi(z) = g(z) + g(z) * [Di(z) - D(z)] *gi(z)          *
C   *                                                                *
C   *       c)  g(z)  = sum(sort) c(sort) * gi(z)(sort).             *
C   *                                                                *
C   * Notation:                                                      *
C   *                                                                *
C   *    ga    ----->  g(z)                                          *
C   *    dtil  ----->  D(z)                                          *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE csts
      USE bzmesh ; USE control_data
      USE control_text ; USE kinkmatrix ; USE logderivative
      USE message ; USE partialwaves
      USE potential ; USE slope
      USE taylor
      IMPLICIT NONE
      INTEGER :: is, i, j, lk, jd, lz, nzm, iq, iq0, lm, lmp, ndera
      INTEGER :: it, ita, ntait, icon, icont, convt, minit, l, m
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: ga
      COMPLEX(KIND=8), DIMENSION(nzm,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,ns) :: bg0
      COMPLEX(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: bg
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: dtil
      COMPLEX(KIND=8), DIMENSION(nlm,nlm) :: ma, mb, mc
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm, wgm
      COMPLEX(KIND=8) :: diag, ty, tyd, prod, dosd, dosk
      REAL(KIND=8)    :: pkx, pky, pkz, ah, weight, sz, dosr
      REAL(KIND=8)    :: pku, pkv, pkvo
      REAL(KIND=8)    :: rconv, iconv, twol, epol, npol
      INTEGER, DIMENSION(:), ALLOCATABLE :: conv
C
C     Allocate temporary blocks
C
      ALLOCATE(sa(nlmq,nlmq,0:nder))
      IF(expan.EQ.'D') THEN
         ALLOCATE(sap(nlmq,nlmq,0:nderp))
      ENDIF
      ALLOCATE(conv(nzm))
C
C     Solve CPA equation
C
      conv(1:nzm)=0
      icpa=0
   10 icpa=icpa+1
C
C     Eq. 1a
C
      DO 11 lz=1,nzm
      IF(conv(lz).EQ.0) ga(lz,1:nq,1:nlm,1:nlm,is)=zero
   11 CONTINUE
C
      DO 20 lk=1,nkvec
C
C     Get the Bloch transformed slope matrices
C
      weight=fkw(lk)
      pkx   =fkx(lk)
      pky   =fky(lk)
      pkz   =fkz(lk)
      CALL smtrx(nlm,pkx,pky,pkz)
C
      DO 30 lz=1,nzm
      IF(conv(lz).GT.0) GO TO 30
C
C     Taylor expansion for the slope matrix S^a and S^a_dot
C
      IF(clkw0(lz).EQ.1) THEN
         kinkm(1:nlmq,1:nlmq)=sa(1:nlmq,1:nlmq,0)
         ndera=nder
      ELSE
         kinkm(1:nlmq,1:nlmq)=sap(1:nlmq,1:nlmq,0)
         ndera=nderp
      ENDIF
      DO 31 jd=1,ndera
      ty=tayl(lz,jd,is)
      IF(clkw0(lz).EQ.1) THEN
         kinkm(1:nlmq,1:nlmq)=kinkm(1:nlmq,1:nlmq)+
     .                        ty*sa(1:nlmq,1:nlmq,jd)
      ELSE
         kinkm(1:nlmq,1:nlmq)=kinkm(1:nlmq,1:nlmq)+
     .                        ty*sap(1:nlmq,1:nlmq,jd)
      ENDIF
   31 CONTINUE
C
C     kinkm = a*S^a - a*D^a
C
      DO 32 i=1,nlmq
      ah=asc(i)
      kinkm(i,1:nlmq)=ah*kinkm(i,1:nlmq) 
   32 CONTINUE
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
      IF(info.NE.0) STOP 'PATHOP-1: INFO <> 0'
C
      gak=weight*gak
C
C     Find the path operator
C
      DO 34 iq=1,nq
      iq0=(iq-1)*nlm
      DO 34 lm=1,nlm
      i=lm+iq0
      DO 34 lmp=1,nlm
      j=lmp+iq0
      ga(lz,iq,lm,lmp,is)=ga(lz,iq,lm,lmp,is)+gak(i,j)
   34 CONTINUE
C
   30 CONTINUE
   20 CONTINUE
C
C     Rotate the path operator
C     obtained in the IBZ for the FBZ
C
      IF(fllbz.NE.'Y') THEN
         DO 40 lz=1,nzm
         IF(conv(lz).GT.0) GO TO 40
         CALL rotgf(ga(lz,1:nq,1:nlm,1:nlm,is))
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
         IF(info.NE.0) STOP 'PATHOP-2: INFO <> 0'
         ma(1:nlm,1:nlm)=ga(lz,iq,1:nlm,1:nlm,is)
C
C        g0(z)^-1
C
         CALL mtrxinv(nlm,ma,unil,mc,info)
         IF(info.NE.0) STOP 'PATHOP-3: INFO <> 0'
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
C
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
      IF(convt.GT.0) THEN
         IF(icpa.LT.ncpa) GO TO 10
         minit=ncpa
         DO 54 lz=1,nzm
         IF(conv(lz).NE.0.AND.conv(lz).LT.minit) minit=conv(lz)
   54    CONTINUE
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
C     Eq. 1b Set up the Green's functions for the alloy components
C
      IF(icpa.GT.1) CALL dyson(is,ga,gi,dtil,nzm,nzm,0)
C
C     Start calculating the Fermi surface
C
      ALLOCATE(kinkmd(nlmq,nlmq),bg(nzm,nq,nlm,nlm))
      bg=zero
      DO 60 lz=1,nzm
      sz=wgm(lz)*spinfc/pi
C
      DO 61 lk=1,nfs
C
C     Get the Bloch transformed slope matrices
C
      weight=fsw(lk)
      pkx   =fsx(lk)
      pky   =fsy(lk)
      pkz   =fsz(lk)
      CALL smtrx(nlm,pkx,pky,pkz)
C
C     Taylor expansion for the slope matrix S^a and S^a_dot
C
      IF(clkw0(lz).EQ.1) THEN
         kinkm(1:nlmq,1:nlmq)=sa(1:nlmq,1:nlmq,0)
         ndera=nder
      ELSE
         kinkm(1:nlmq,1:nlmq)=sap(1:nlmq,1:nlmq,0)
         ndera=nderp
      ENDIF
      kinkmd=zero
      DO 62 jd=1,ndera
      ty=tayl(lz,jd,is)
      tyd=tayld(lz,jd,is)
      IF(clkw0(lz).EQ.1) THEN
         kinkm(1:nlmq,1:nlmq)=kinkm(1:nlmq,1:nlmq)+
     .                        ty*sa(1:nlmq,1:nlmq,jd)
         kinkmd(1:nlmq,1:nlmq)=kinkmd(1:nlmq,1:nlmq)+
     .                         tyd*sa(1:nlmq,1:nlmq,jd)
      ELSE
         kinkm(1:nlmq,1:nlmq)=kinkm(1:nlmq,1:nlmq)+
     .                        ty*sap(1:nlmq,1:nlmq,jd)
         kinkmd(1:nlmq,1:nlmq)=kinkmd(1:nlmq,1:nlmq)+
     .                         tyd*sap(1:nlmq,1:nlmq,jd)
      ENDIF
   62 CONTINUE
C
C     kinkm = a*S^a - a*D^a
C
      DO 63 i=1,nlmq
      ah=asc(i)
      kinkm(i,1:nlmq)=ah*kinkm(i,1:nlmq)
      kinkmd(i,1:nlmq)=ah*kinkmd(i,1:nlmq)
   63 CONTINUE
C
C     kinkm = S^a - D^a
C
      DO 64 iq=1,nq
      iq0=(iq-1)*nlm
      DO 64 lm=1,nlm
      i=lm+iq0
      DO 64 lmp=1,nlm
      j=lmp+iq0
      kinkm(i,j)=kinkm(i,j)-dtil(lz,iq,lm,lmp,is)
   64 CONTINUE
C
      CALL mtrxinv(nlmq,kinkm,unit,gak,info)
      IF(info.NE.0) STOP 'PATHOP-1: INFO <> 0'
C
      gak=weight*gak
C
C     Find the path operator
C
      IF(fsts.EQ.1) THEN
         DO 65 iq=1,nq
         iq0=(iq-1)*nlm
         DO 65 lm=1,nlm
         i=lm+iq0
         DO 65 lmp=1,nlm
         j=lmp+iq0
         diag=SUM(gak(i,1:nlmq)*kinkmd(1:nlmq,j))
         bg(lz,iq,lm,lmp)=bg(lz,iq,lm,lmp)+diag
   65    CONTINUE
      ENDIF
C
      dosk=zero
      DO 66 iq=1,nq
      iq0=(iq-1)*nlm
      DO 67 l=0,lmax
      twol=2*l+1
      DO 68 m=-l,l
      lm=l*l+l+m+1
      i=lm+iq0
      it=itq(iq)
      dosd=SUM(gak(i,1:nlmq)*kinkmd(1:nlmq,i))
      DO 69 ita=1,nta(it)
      diag=dosd-gi(lz,ita,iq,lm,lm,is)*dfi(lz,ita,it,l,is,2)
C
C     Extract the poles of D_dot
C
      diag=diag-twol*gsng(lz,ita,it,l,is)
      DO 70 j=1,necr(l,ita,it,is)
      epol=ecr(l,ita,it,is,j)
      npol=nocr(l,ita,it,is,j)
   70 diag=diag+twol*npol/(zm(lz)-epol)
C 
      dosk=dosk+conc(ita,it)*diag
   69 CONTINUE
   68 CONTINUE
   67 CONTINUE
   66 CONTINUE
C
      IF(fsts.EQ.0) THEN
         dosr=sz*AIMAG(dosk)
         IF(fs.EQ.1) THEN
            pkv=pkz
            IF(pkx.LE.0) pku=pkx
            IF(pkx.GT.0) pku=SQRT(pkx*pkx+pky*pky)
         ELSEIF(fs.EQ.2) THEN
            pku=pkx
            pkv=pky
         ENDIF
         IF(lk.EQ.1) THEN
            WRITE(12,140) pku,pkv,dosr
         ELSE
            IF(ABS(pkv-pkvo).LT.1.d-5) THEN
               WRITE(12,140) pku,pkv,dosr
            ELSE
               WRITE(12,141) pku,pkv,dosr
            ENDIF
         ENDIF
         pkvo=pkv
      ENDIF
   61 CONTINUE
      IF(fsts.EQ.1) THEN
         CALL rotgf(bg(lz,1:nq,1:nlm,1:nlm))
         DO 71 iq=1,nq
         DO 72 lm=1,nlm
         bg0(1:nzm,iq,lm,is)=bg(1:nzm,iq,lm,lm)
   72    CONTINUE
   71    CONTINUE
      ENDIF
   60 CONTINUE
C
      DEALLOCATE(sa)
      IF(expan.EQ.'D') THEN
         DEALLOCATE(sap)
      ENDIF
      DEALLOCATE(kinkmd,bg,conv)
C
  110 FORMAT(/,' FESPTH: CPA equation for IS =',i2,' not converged',
     .   /,9x,'number of iterations performed ',i2,' to ',i2)
  120 FORMAT(/,' FESPTH: CPA equation for IS =',i3,
     .         ' converged in ',i2,'-',i2,' iterations')
  140 FORMAT(2f10.6,f16.6)
  141 FORMAT(/,2f10.6,f16.6)
      RETURN
      END
