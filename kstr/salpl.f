      SUBROUTINE salpl
C   ******************************************************************
C   *                                                                *
C   *  Evaluate the slope matrix Sa and energy derivatives for the   *
C   *  screening parameters {t1,t2,t3,t4} calculated in "screen":    *
C   *                                                                *
C   *    (t3*t4 + t3*S0*t3) * (t1/t3 - Sa) = (t1*t4-t2*t3)           *
C   *                                                                *
C   *  After inversion only the Sa(R',R) block is retained, then     *
C   *  Sa(R',W), Sa(W ,R) and Sa(W ,W) blocks are dropped.           *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE factorial
      USE control_data
      USE control_text
      USE lattice
      USE lmtolm
      USE screening
      USE message
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: fac
      INTEGER :: iq, ir1, ir2, nrq, nrlm, nrlmw, jd, id, i, lmp, 
     .           info, irl, irp, jqp, irlm, lp
C
C     Loop for sites
C
      DO 20 iq=1,nq
      IF(msgl.EQ.1) WRITE(msgio,'(/,a)')
     .          ' Initialize unscreened low-low subblock'
      ir1=nriq(iq-1)+1
      ir2=nriq(iq)
      nrq=ir2-ir1+1
      nrlm=nlm*nrq
C
C     Increase the dimension for the Watson sphere
C
      nrlmw=nrlm+nlmw
C
      ALLOCATE(sa(nrlmw,nrlmw,0:nder))
      CALL alltmp(nrlmw,1)
C
      sa(1:nrlmw,1:nrlmw,0)=0.d0
      sb(1:nrlmw,1:nrlmw)=0.d0
      sbd(1:nrlmw,1:nrlmw,1:nder)=0.d0
C
C     Set up B = t3*t4 + t3*S0*t3 and B^(n)
C
      CALL mbares(iq,ir1,ir2,nrlm)
C
C     Invert  B
C
      IF(msgl.EQ.1) WRITE(msgio,'(a,i5)')
     .          ' Construct low-low subblock'
      CALL dsifa(sb,nrlmw,nrlmw,ipvt,info)
      IF(info.NE.0) STOP 'info <> 0'
C
C     A = - B^-1 * [t1*t4 - t2*t3]
C
      jd=0
      DO 21 i=1,nrlmw
      sa(i,i,jd)=-dett(i)
   21 CALL dsisl(sb,nrlmw,nrlmw,ipvt,sa(1,i,jd))
C
C     -B^(n)*A for n=1,nder
C
      IF(msgl.EQ.1.AND.nder.ge.1) 
     .             WRITE(msgio,'(a)') ' First order derivative'
      DO 22 id=1,nder
      fac=-1.d0
      CALL dsymm('L','U',nrlmw,nrlmw,fac,sbd(1,1,id),
     .           nrlmw,sa(1,1,jd),nrlmw,0.d0,sa(1,1,id),nrlmw)
      IF(id.EQ.1) THEN
C
C        A^(1) = B^-1*(-B^(1)*A) first order derivative
C
         DO 23 i=1,nrlmw
   23    CALL dsisl(sb,nrlmw,nrlmw,ipvt,sa(1,i,id))   
      ENDIF
   22 CONTINUE
C
C     sum(i=0,n-1) -ifib(n,i)*B^(n-i)*A^(i) for n=2,nder
C
      IF(msgl.EQ.1.AND.nder.gt.1) 
     .             WRITE(msgio,'(a)') ' Higher order derivatives'
      DO 24 jd=2,nder
      DO 25 id=1,jd-1
      fac=-ifib(jd,id)
   25 CALL dsymm('L','U',nrlmw,nrlmw,fac,sbd(1,1,jd-id),
     .           nrlmw,sa(1,1,id),nrlmw,1.d0,sa(1,1,jd),nrlmw)
      DO 26 i=1,nrlmw
   26 CALL dsisl(sb,nrlmw,nrlmw,ipvt,sa(1,i,jd))
   24 CONTINUE
C
C     Sa = A + t1/t3
C
      irl=-nlm
      DO 27 irp=ir1,ir2
      jqp=jqbas(irp)
      irl=irl+nlm
      DO 27 lmp=1,nlm
      lp=llx(lmp)
      irlm=irl+lmp
   27 sa(irlm,irlm,0:nder)=bigd(1,lp,jqp,0:nder)+
     .                   sa(irlm,irlm,0:nder)
C
C     Symmetrize structure constant
C
      CALL symstr(ir1,ir2)
C
      CALL alltmp(nrlmw,2)
C
C     Print structure constant
C
      CALL prnstr(iq,ir1,ir2)
C
      IF(store.EQ.'Y') CALL stores(nrq)
C
C     Calculate the high-l structure constants if requested
C
      IF(high.EQ.'Y') CALL salplh(iq,ir1,ir2,nrq,nrlmw,nrlm)
      DEALLOCATE(sa)
   20 CONTINUE
C
      RETURN
      END
