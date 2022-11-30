      SUBROUTINE trnsfm(case,lmaxh0)
C   ******************************************************************
C   *                                                                *
C   *    Read lattice translation vectors, the Madelung matrices,    *
C   *    and the corresponding slope matrices generated in KSTR.     *
C   *                                                                *
C   *    If expan = D it reads a second set of slope matrix for      *
C   *    a big negative kappa values.                                *
C   *                                                                *
C   *    If expan = M it reads a second set of slope matrix for      *
C   *    a big negative kappa values and constructs one slope        *
C   *    with large extended derivatives.                            *
C   *                                                                *
C   ******************************************************************
      USE bzmesh     ; USE control_data ; USE control_text
      USE lattice    ; USE message      ; USE potential
      USE radialmesh ; USE slope        ; USE taylor
      USE text
      IMPLICIT NONE
      CHARACTER(LEN=1)  :: modes, high
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: am
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: bm
      INTEGER, DIMENSION(:), ALLOCATABLE :: kpvt
      REAL(KIND=8) :: dkp, tdi
      INTEGER :: case, iq, iv, iv1, iv2, lmp, lm, npa, jd
      INTEGER :: nlmh0, nlh0, lmaxh0, nlmact, error, nder0, mm
C
      CHARACTER(LEN=1)  :: modesp, highp
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vmdlp, sigmap
      REAL(KIND=8) :: cmdlp, dum
      INTEGER :: nqp, nlp, nvp, nlh0p, itransp, idum, info, i, j
C
      IF(case.EQ.1) THEN
         oexpan=expan
      ELSEIF(case.EQ.2) THEN
         expan=oexpan
         OPEN(1,FILE=for001,STATUS='OLD',FORM='UNFORMATTED')
      ENDIF
C
      IF(expan.EQ.'S') THEN
         WRITE(m6,100)
      ELSEIF(expan.EQ.'D') THEN
         WRITE(m6,101)
      ELSEIF(expan.EQ.'M') THEN
         WRITE(m6,102)
      ELSE
         WRITE(m6,103)
         STOP
      ENDIF
C
      READ(1) txts,modes,high
      WRITE(m6,110) txts
      IF(mode.EQ.'2D'.AND.modes.NE.'S') THEN
         WRITE(m6,111) mode
         STOP
      ENDIF
      IF(mode.EQ.'3D'.AND.modes.NE.'B') THEN
         WRITE(m6,112) mode
         STOP
      ENDIF
      IF(high.NE.'Y') THEN
         WRITE(m6,113) high
         STOP
      ENDIF
      READ(1) nq,nl,nlm,nlmq,nv,nlh0,nlmh0,lat
      lmaxh0=nlh0-1
      nlmact=nlmh0
      lmax=nl-1
      lmax2=MAX0(2*lmax,1)
      nl2=lmax2+1
      nlm2=nl2*nl2
C
      IF(case.EQ.1) THEN
         ALLOCATE(nviq(0:nq),iqbas(nv),jqbas(nv))
         ALLOCATE(qx(nq),qy(nq),qz(nq),tx(nv),ty(nv),tz(nv))
         ALLOCATE(wst(nq),wsi(nq),wsc(nq))
         ALLOCATE(vmdl(nq,nq))
      ENDIF
C
      READ(1) nviq(0:nq),iqbas(1:nv),jqbas(1:nv)
      READ(1) npa
      IF(case.EQ.1) THEN
         READ(1) wsa,wst(1:nq),wsi(1:nq),wsc(1:nq)
      ELSE
         READ(1) dum
      ENDIF
      READ(1) bsx,bsy,bsz,qx(1:nq),qy(1:nq),qz(1:nq)
      READ(1) tkx,tky,tkz,boa,coa,alf,bet,gam
      READ(1) tx(1:nv),ty(1:nv),tz(1:nv)
      READ(1) kw20,nder,itrans
      nder0=nder
      IF(expan.EQ.'M') nder0=nder+1+nder
      IF(itrans.NE.3) THEN
         WRITE(m6,114) itrans
         STOP
      ENDIF
      WRITE(m6,115) itrans,kw20,nder
C
      IF(case.EQ.1) THEN
         ALLOCATE(sigma(0:lmax,nq))
         ALLOCATE(tmat(4,0:lmax,nq,0:nder0),facd(0:nder0))
      ENDIF
C
C     Factorial used in the Taylor expansion
C
      facd(0)=1.d0
      DO 20 jd=1,nder0
   20 facd(jd)=jd*facd(jd-1)
C
      READ(1) sigma(0:lmax,1:nq)
      DO 21 iq=1,nq
   21 sigma(0:lmax,iq)=sigma(0:lmax,iq)*wst(iq)
      READ(1) tmat(1:4,0:lmax,1:nq,0:nder)
C
      READ(1) cmdl,vmdl(1:nq,1:nq)
C
      ALLOCATE(salpl(nlmact,nlm,nv,0:nder0),STAT=error)
      IF(error.NE.0) THEN
         WRITE(m6,*) 'Space for Salpl could not be allocated'
         STOP
      ENDIF
      DO 30 iq=1,nq
      DO 30 jd=0,nder
      iv1=nviq(iq-1)+1
      iv2=nviq(iq)
      DO 30 iv=iv1,iv2
      READ(1) ((salpl(lmp,lm,iv,jd),lmp=1,nlmact),lm=1,nlm)
   30 CONTINUE
      CLOSE(1)
C
      mder=nder
      IF(expan.EQ.'S') RETURN
C
C     Read second slope matrix and check for the structure and setup
C
      OPEN(1,FILE=for001p,STATUS='OLD',FORM='UNFORMATTED')
C
      READ(1) txts,modesp,highp
      WRITE(m6,110) txts
      IF(modesp.NE.modes) THEN
         WRITE(m6,120) 'MODESp',modesp,'MODES',modes
         STOP
      ENDIF
      IF(highp.NE.high) THEN
         WRITE(m6,120) 'HIGHp',highp,'HIGH',high
         STOP
      ENDIF
      READ(1) nqp,nlp,nlm,nlmq,nvp,nlh0p,nlmh0
      IF(nqp.NE.nq) THEN
         WRITE(m6,121) 'NQp',nqp,'NQ',nq
         STOP
      ENDIF
      IF(nlp.NE.nl) THEN
         WRITE(m6,121) 'NLp',nlp,'NL',nl
         STOP
      ENDIF
      IF(nvp.NE.nv) THEN
         WRITE(m6,121) 'NVp',nvp,'NV',nv
         STOP
      ENDIF
      IF(nlh0p.NE.nlh0) THEN
         WRITE(m6,121) 'NLH0p',nlh0p,'NLH0',nlh0
         STOP
      ENDIF
C
      READ(1) idum
      READ(1) idum
      READ(1) dum
      READ(1) dum
      READ(1) dum
      READ(1) dum
      READ(1) kw20p,nderp,itransp
      IF(itransp.NE.itrans) THEN
         WRITE(m6,111) 'ITRANSp',itransp,'ITRANS',itrans
         STOP
      ENDIF
      WRITE(m6,115) itrans,kw20p,nderp
      mder=MAX0(nder,nderp)
C
      ALLOCATE(vmdlp(nq,nq))
      ALLOCATE(sigmap(0:lmax,nq))
      READ(1) sigmap(0:lmax,1:nq)
      DO 40 iq=1,nq
   40 sigmap(0:lmax,iq)=sigmap(0:lmax,iq)*wst(iq)
      READ(1) dum
      READ(1) cmdlp,vmdlp(1:nq,1:nq)
      IF(ABS(cmdlp-cmdl).GT.1.d-10) THEN
         WRITE(m6,122) 'CMDLp',cmdlp,'CMDL',cmdl
         STOP
      ENDIF
      IF(ABS(SUM(vmdlp)-SUM(vmdl)).GT.1.d-10) THEN
         WRITE(m6,123) 'VMDLp','VMDL'
         STOP
      ENDIF
      IF(ABS(SUM(sigmap)-SUM(sigma)).GT.1.d-10) THEN
         WRITE(m6,123) 'SIGMAp','SIGMA'
         STOP
      ENDIF
      DEALLOCATE(vmdlp,sigmap)
C
      ALLOCATE(salplp(nlmact,nlm,nv,0:nderp),STAT=error)
      IF(error.NE.0) THEN
         WRITE(m6,*) 'Space for Salplp could not be allocated'
         STOP
      ENDIF
      DO 50 iq=1,nq
      DO 50 jd=0,nderp
      iv1=nviq(iq-1)+1
      iv2=nviq(iq)
      DO 50 iv=iv1,iv2
      READ(1) ((salplp(lmp,lm,iv,jd),lmp=1,nlmact),lm=1,nlm)
   50 CONTINUE
      CLOSE(1)
C
      IF(expan.EQ.'D') RETURN
C
C     Set up the modified TwoC expansion for salpl
C
      IF(nder0.LT.(nder+nderp+1)) THEN
         WRITE(m6,104) nder,nderp
         STOP
      ENDIF
C
      mm=nderp+1
      ALLOCATE(am(mm,mm),bm(mm),kpvt(mm))
C
      dkp=kw20p-kw20
      DO i=0,nderp
      DO j=0,nderp
      am(i+1,j+1)=dkp**(nder+1+j-i)/facd(nder+1+j-i)
      ENDDO
      ENDDO
      CALL dgefa(am,mm,mm,kpvt,info)
      IF(info.NE.0) STOP 'info<>0'
C
      DO 60 iq=1,nq
      iv1=nviq(iq-1)+1
      iv2=nviq(iq)
      DO 60 iv=iv1,iv2
      DO 60 lmp=1,nlmact
      DO 60 lm=1,nlm
C
      DO i=0,nderp
      tdi=0.d0
      DO j=i,nder
      tdi=tdi+salpl(lmp,lm,iv,j)*dkp**(j-i)/facd(j-i)
      ENDDO
      bm(i+1)=salplp(lmp,lm,iv,i)-tdi
      ENDDO
      CALL dgesl(am,mm,mm,kpvt,bm,0)
C
      DO j=0,nderp
      salpl(lmp,lm,iv,nder+1+j)=bm(j+1)
      ENDDO
C
   60 CONTINUE
      nder=nder+1+nderp
      mder=nder
      WRITE(m6,105) nder
      expan='S'
      DEALLOCATE(salplp,am,bm,kpvt)
C
  100 FORMAT(/,' TRNSFM:   Single Taylor expansion for the',
     .         ' slope matrix')
  101 FORMAT(/,' TRNSFM:   Double Taylor expansion for the',
     .         ' slope matrix')
  102 FORMAT(/,' TRNSFM:   Modified TwoC expansion for the',
     .         ' slope matrix')
  103 FORMAT(/,' TRNSFM: **EXPAN must be S, D or M')
  104 FORMAT(/,' TRNSFM: **For EXPAN = M, the nderp > nder',2i4)
  105 FORMAT(/,11x,'Modified TwoC expansion with NDER =',i3)
  110 FORMAT(/,11x,a)
  111 FORMAT(/,' TRNSFM:** MODE =',a,' but BSTR-MODE is not S.')
  112 FORMAT(/,' TRNSFM:** MODE =',a,' but BSTR-MODE is not B.')
  113 FORMAT(/,' TRNSFM:** Incorrect transfer matrices,',
     .         ' HIGH = ',a)
  114 FORMAT(/,' TRNSFM:** ITRANS =',i3,' is not implemented')
  115 FORMAT(/,11x,'ITRANS =',i3,' (kappa0*w)^2 =',f8.3,' NDER =',i3)
  120 FORMAT(/,' TRNSFM:** ',a,'=',a,' is different of ',a,'=',a)
  121 FORMAT(/,' TRNSFM:** ',a,'=',i3,' is different of ',a,'=',i3)
  122 FORMAT(/,' TRNSFM:** ',a,'=',f10.6,
     .         ' is different of ','=',f10.6,a)
  123 FORMAT(/,' TRNSFM:** ',a,' is different of ',a)
C
      RETURN
      END
