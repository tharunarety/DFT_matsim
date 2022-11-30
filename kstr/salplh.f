      SUBROUTINE salplh(iq,ir1,ir2,nrq,nrlmw,nrlm)
C   ******************************************************************
C   *                                                                *
C   *  Evaluate the high-low subblock of the slope matrix Sa and     *
C   *  energy derivatives using the blowing-up technique:            *
C   *                                                                *
C   *    Sa(R'H',RL) = S0(R'H',RL)*t1 - S0(R'H',i)*t3*Sa(i,RL)       *
C   *                                                                *
C   *  where i = {R''L''} and {Watson L''} .                         *
C   *  Note: Higers are not calculated on the Watson sphere.         *
C   *        Only R' = Q' are calculated.                            *
C   *                                                                *
C   ******************************************************************
      USE factorial
      USE control_data
      USE control_text
      USE lattice
      USE lmtolm
      USE screening
      USE message
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: fac, max, min
      INTEGER :: iq, ir1, ir2, nrq, nrlmw, irl, irp, jqp, 
     .           lmp, lp, irlm, jrlm, id, jd, nrlm
C
      ALLOCATE(sbd(nlmh,nrlmw,0:nder),sc(0:nder))
      ALLOCATE(sah(nlmh,nrlmw,0:nder))
      ALLOCATE(siq(nlm,nrlm,0:nder))
C
C     Save the first row of Sa(L',L)
C
      siq(1:nlm,1:nrlm,0:nder)=sa(1:nlm,1:nrlm,0:nder)
C
      IF(msgl.EQ.1) WRITE(msgio,'(/,a)')
     .          ' Initialize unscreened high-low subblock'
      CALL mbareh(iq,ir1,ir2,nrlmw)
C
C     t1-t3*Sa
C
      irl=-nlm
      DO 20 irp=ir1,ir2
      jqp=jqbas(irp)
      irl=irl+nlm
      do 20 lmp=1,nlm
      lp=llx(lmp)
      irlm=irl+lmp
      DO 21 jrlm=1,nrlmw
C
      sc(0:nder)=0.d0
      DO 22 jd=0,nder
      DO 22 id=0,jd
   22 sc(jd)=sc(jd)+ifib(jd,id)*tmat(3,lp,jqp,jd-id)*
     .                          sa(irlm,jrlm,id)
   21 sa(irlm,jrlm,0:nder)=-sc(0:nder)
   20 sa(irlm,irlm,0:nder)=tmat(1,lp,jqp,0:nder)+
     .                     sa(irlm,irlm,0:nder)
      IF(wats.EQ.'Y') THEN
        irl=nrq*nlm
        do 23 lmp=1,nlmw
        lp=llx(lmp)
        irlm=irl+lmp
        DO 24 jrlm=1,nrlmw
C
        sc(0:nder)=0.d0
        DO 25 jd=0,nder
        DO 25 id=0,jd
   25   sc(jd)=sc(jd)-ifib(jd,id)*alwats(lp,iq,jd-id)*
     .                          sa(irlm,jrlm,id)
   24   sa(irlm,jrlm,0:nder)=-sc(0:nder)
   23   sa(irlm,irlm,0)=1.d0+sa(irlm,irlm,0)
      ENDIF
C
C     1/t4*S0
C
      DO 26 lmp=1,nlm
      lp=llx(lmp)
      DO 26 jrlm=1,nrlmw
      sc(0:nder)=0.d0
      DO 27 jd=0,nder
      DO 27 id=0,jd
   27 sc(jd)=sc(jd)+ifib(jd,id)*cd(2,lp,iq,jd-id)*
     .       sbd(lmp,jrlm,id)
   26 sbd(lmp,jrlm,0:nder)=sc(0:nder)
C
C     t2/t4
C
      IF(msgl.EQ.1) WRITE(msgio,'(a,i5)')
     .          ' Construct high-low subblock'
      sah(1:nlmh,1:nrlmw,0:nder)=0.d0
      DO 28 lmp=1,nlm
      lp=llx(lmp)
   28 sah(lmp,lmp,0:nder)=bigd(2,lp,iq,0:nder)
C
C     Calculate 1/t4*S0 * (t1-t3*Sa)
C
      DO 29 jd=0,nder
      DO 29 id=0,jd
      fac=ifib(jd,id)
C
C     alpha*A(m,k)*B(k,n) + beta*C(m,n) --> C(m,n)
C
C     dgemm('N','N',m,n,k,alpha, A,dim_A, B,dim_B, beta, C,dim_C)
C
      CALL dgemm('N','N',nlmh,nrlmw,nrlmw,fac,sbd(1,1,jd-id),
     .           nlmh,sa(1,1,id),nrlmw,1.d0,sah(1,1,jd),nlmh)
   29 CONTINUE
C
C     Symmetrize structure constant
C
      CALL symsth(ir1,ir2)
C
C     Print structure constant
C
      CALL prnsth(iq,ir1,ir2)
C
      IF(store.EQ.'Y') CALL storeh(nrq)
C
      siq(1:nlm,1:nrlm,0:nder)=siq(1:nlm,1:nrlm,0:nder)-
     .                         sah(1:nlm,1:nrlm,0:nder)
      max=DABS(MAXVAL(siq))
      min=DABS(MINVAL(siq))
      IF(max.GT.errp.OR.min.GT.errp) THEN
        WRITE(m6,100) max,MAXLOC(siq),min,MINLOC(siq)
        IF(msgl.NE.0) WRITE(msgio,100) 
     .                max,MAXLOC(siq),min,MINLOC(siq)
      ENDIF
      DEALLOCATE(sbd,sc)
      DEALLOCATE(sah,siq)
C
  100 FORMAT(/,' SALPLH:   error in the lower part of the ',
     .         'slope matrix after the blowing up',//,
     .         11x,'max =',f10.6,' with location ',3i5,/,
     .         11x,'min =',f10.6,' with location ',3i5)
      RETURN
      END
