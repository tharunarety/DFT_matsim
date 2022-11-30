      SUBROUTINE trmtrx
C   ******************************************************************
C   *   Makes the transformation matrix tmat = {t1,t2,t3,t4}.        *
C   *   Input: itrans: characterizes the transformation:             * 
C   *      0:old transformation with empirical alpha values:         *
C   *             |N^a(kappa)> identical to |N^0(kappa)>             *
C   *             |J^a(kappa)> =|J^0(kappa)>-alpha|K^0(kappa)>       *
C   *             |J^a(kappa)> = 0 at hard core radius a=sigma*wst   *
C   *             sigma = wst*[2(2l+1)*alpha]^(1/(2l+1))             * 
C   *             the Wronskian W{N^a,J^a} equals W{N^0,J^0}=w/2     *
C   *      1:old transformation:                                     *
C   *             |N^a(kappa)> identical to |N^0(kappa)>             *
C   *             |J^a(kappa)> =|J^0(kappa)>-alpha|K^0(kappa)>       *
C   *             |J^a(kappa)> = 0 at hard core radius a=sigma*wst   *
C   *                          sigma is input parameter              * 
C   *             the Wronskian W{N^a,J^a} equals W{N^0,J^0}=w/2     *
C   *      2:head |N^a(kappa)> has value 1 and derivative 0          *
C   *             |J^a(kappa)> has value 0 and derivative 1/a*a      *
C   *                          at hard core radius a=sigma*wst       *
C   *                          sigma is input parameter              * 
C   *      3:head |N^a(kappa)> has value 1 and derivative 0          *
C   *             |J^a(kappa)> has value 0 and derivative -1/a       *
C   *                          at hard core radius a=sigma*wst       *
C   *                          sigma is input parameter              * 
C   *   Output:                                                      *
C   *          tmat  :transformation matrix for head and tail        *
C   *                 and  (kappa*ws)^2-derivative                   *
C   *          bigd  :t1/t3 and t2/t4 and energy derivatives         *
C   *          cd    :t3*t4 and 1/t4  and energy derivatives         *
C   *          dfac  :t1*t4 - t2*t3   detrminant of tmat             *
C   *                                                                *
C   *   Remarks: The transformation matrix is:                       *
C   *                                                                *
C   *            |N^a>=tmat(1)*|N>+tmat(2)*|J>                       *
C   *            |J^a>=tmat(3)*|N>+tmat(4)*|J>                       *
C   *                                                                *
C   *        and ||N^a>=|N^a>-|J^a>S^a                               *
C   ******************************************************************
      USE basis
      USE factorial
      USE control_data
      USE lattice
      USE screening
      USE message
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: tcd, cm1, dm1
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: ac, bd, jntmp
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: al, nm1, dn, dj
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: t
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: fi, gi, j, n
      REAL(KIND=8) :: r2, rfct, rfct2, hcr, fct1, fct2
      INTEGER      :: i, iq, k, l, jder, nder0
C
      ALLOCATE(tcd(0:nder),cm1(0:nder),dm1(0:nder),t(4,0:nder))
      ALLOCATE(ac(0:nder),bd(0:nder),jntmp(0:nder),al(0:nder))
      ALLOCATE(nm1(0:nder),dn(0:nder),dj(0:nder),j(0:nder),n(0:nder))
C
      nder0=MAX0(nder,1)
      ALLOCATE(fi(-nder0:nder0+lmax),gi(-nder0:nder0+lmax))
C
      IF(itrans.EQ.0) THEN
         DO 10 iq=1,nq
         IF(lmax.EQ.0) THEN
            sigma(0,iq)=2.d0*0.2143D0
         ELSEIF(lmax.EQ.1) THEN
            sigma(0,iq)=2.d0*0.2872D0
            sigma(1,iq)=(2.d0*(2.d0*1.+1)*0.02582D0)**(1./3.)
         ELSEIF(lmax.EQ.2) THEN
            sigma(0,iq)=2.d0*0.3485D0
            sigma(1,iq)=(2.d0*(2.d0*1.+1)*0.05303D0)**(1./3.)
            sigma(2,iq)=(2.d0*(2.d0*2.+1)*0.010714D0)**(1./5.)
         ELSEIF(lmax.EQ.3) THEN
            sigma(0,iq)=2.d0*0.385057D0
            sigma(1,iq)=(2.d0*(2.d0*1.+1)*0.073209D0)**(1./3.)
            sigma(2,iq)=(2.d0*(2.d0*2.+1)*0.0224806D0)**(1./5.)
            sigma(3,iq)=(2.d0*(2.d0*3.+1)*0.0060691D0)**(1./7.)
         ELSE
            WRITE(m6,*) ' TRMTRX: itrans=0 not implemented for',
     .                  ' lmax=',lmax
            STOP
         ENDIF
         sigma(0:lmax,iq)=sigma(0:lmax,iq)/wst(iq)
   10    CONTINUE
      ELSEIF(itrans.LT.0.OR.itrans.GT.3) THEN
         WRITE(m6,*) ' TRMTRX: bad itrans  ',itrans
         itrans=2
         WRITE(m6,103) itrans
      ENDIF
C
      DO 20 iq=1,nq
      DO 20 l=0,lmax
      hcr=sigma(l,iq)*wst(iq)
      r2=hcr*hcr
      CALL bessh(kap2*r2,-nder0,l+nder0,fi(-nder0),gi(-nder0))
C
      rfct=hcr**(-l-1)
      rfct2=hcr**l
C
C     make n_l and its energy derivatives
C
      n(0)   = rfct*gi(l)
      DO 21 jder=1,nder
   21 n(jder)= n(jder-1)*(r2/2.)*(gi(l-jder)/gi(l-jder+1))/
     .         (2*(l-jder)+1)
C
C     make n_(l+1) and its energy derivs
C
      dn(0) = rfct/hcr*gi(l+1)
      DO 22 jder=1,nder
   22 dn(jder)= dn(jder-1)*(r2/2.)*(gi(l-jder+1)/gi(l-jder+2))/
     .          (2*(l-jder)+3)
C
C     make energy derives of radial deriv of n_l
C     nb.. not logarithmic deriv
C
      DO 23 jder=0,nder
   23 dn(jder)=(l*n(jder)/hcr -(2*l+1)*dn(jder))/ws
C
C     make inverse of n_l, and its energy derivs
C
      nm1(0)=1./n(0)
      jntmp(1:nder)=0.d0
      nm1(1:nder)=0.d0
      DO 24 jder=1,nder
      DO 25 i=1,jder
   25 jntmp(jder)=jntmp(jder)+n(i)*nm1(jder-i)*ifib(jder-1,i-1)
      DO 26 i=1,jder
   26 nm1(jder)=nm1(jder)-jntmp(i)*nm1(jder-i)*ifib(jder-1,i-1)
   24 CONTINUE
C
C     make j_l and its energy derivs
C
      j(0)   = rfct2*fi(l)
      DO 27 jder=1,nder
   27 j(jder)= -j(jder-1)*(r2/2.)*(fi(l+jder)/fi(l+jder-1))/
     .         (2*(l+jder)-1)
C
C     make j_(l-1) and its energy derivs
C
      dj(0) = rfct2/hcr*fi(l-1)
      DO 28 jder=1,nder
   28 dj(jder)= -dj(jder-1)*(r2/2.)*(fi(l+jder-1)/fi(l+jder-2))/
     .          (2*(l+jder)-3)
C
C     make enrgy derivs of radial derivs of j_l
C     nb.. not logarithmic deriv
C
      DO 29 jder=0,nder
   29 dj(jder)=(-(l+1)*j(jder)/hcr +(2*l-1)*dj(jder))/ws
C
C     make alpha and its derivs
C
      al(0)=j(0)*nm1(0)
      al(1:nder)=0.d0
      DO 30 jder=1,nder
      DO 31 i=0,jder
   31 al(jder)=al(jder)+j(i)*nm1(jder-i)*ifib(jder,i)
   30 CONTINUE
C
C     itrans dependant part
C
      fct1=2.*hcr*hcr*ws
      IF(itrans.EQ.0.OR.itrans.EQ.1) THEN
         dfac(l,iq) = 1.d0
      ELSEIF(itrans.EQ.2) THEN
         fct2 = 2./ws
         dfac(l,iq) = 2./ws
      ELSEIF(itrans.EQ.3) THEN
         fct2 = -2.*hcr
         dfac(l,iq) = -2.*hcr
      ENDIF
C
      IF(itrans.EQ.0.OR.itrans.EQ.1) THEN
         t(1,0) = 1.d0
         t(1,1:nder) = 0.d0
         t(2,0:nder) = 0.d0
         t(3,0:nder) = -al(0:nder)
         t(4,0) = 1.d0
         t(4,1:nder) = 0.d0
         tcd(0:nder) = -al(0:nder)
      ELSE
         t(1,0:nder) = fct1*dj(0:nder)
         t(2,0:nder) = -fct1*dn(0:nder)
         t(3,0:nder) = -fct2*j(0:nder)
         t(4,0:nder) = fct2*n(0:nder)
         tcd(0)=t(3,0)*t(4,0)
         tcd(1:nder)=0.d0
         DO 32 jder=1,nder
         DO 33 i=0,jder
   33    tcd(jder)=tcd(jder)+ifib(jder,i)*t(3,jder-i)*t(4,i)
   32    CONTINUE
      ENDIF
C
C     end itrans dependant part
C
C     make inverse of c and its energy derivs
C
      cm1(0)=1./t(3,0)
      jntmp(1:nder)=0.d0
      cm1(1:nder)=0.d0
      DO 34 jder=1,nder
      DO 35 i=1,jder
   35 jntmp(jder)=jntmp(jder)+t(3,i)*cm1(jder-i)*
     .            ifib(jder-1,i-1)
      DO 36 i=1,jder
   36 cm1(jder)=cm1(jder)-jntmp(i)*cm1(jder-i)*ifib(jder-1,i-1)
   34 CONTINUE
C
C     make a/c  and energy derivs thereof
C
      jntmp(0) = t(1,0)*cm1(0)
      ac(0)=jntmp(0)
      jntmp(1:nder)=0.d0
      DO 37 jder=1,nder
      DO 38 i=0,jder
   38 jntmp(jder)=jntmp(jder)+ifib(jder,i)*t(1,jder-i)*cm1(i)
   37 ac(jder)=jntmp(jder)
C
C     make inverse of d and its energy derivs
C
      dm1(0)=1./t(4,0)
      jntmp(1:nder)=0.d0
      dm1(1:nder)=0.0d0
      DO 39 jder=1,nder
      DO 40 i=1,jder
   40 jntmp(jder)=jntmp(jder)+t(4,i)*dm1(jder-i)*
     .            ifib(jder-1,i-1)
      DO 41 i=1,jder
   41 dm1(jder)=dm1(jder)-jntmp(i)*dm1(jder-i)*ifib(jder-1,i-1)
   39 CONTINUE
C
C     make b/d  and energy derivs thereof
C
      jntmp(0) = t(2,0)*dm1(0)
      jntmp(1:nder)=0.d0
      bd(0)=jntmp(0)
      DO 42 jder=1,nder
      DO 43 i=0,jder
   43 jntmp(jder)=jntmp(jder)+ifib(jder,i)*t(2,jder-i)*
     .            dm1(i)
   42 bd(jder)=jntmp(jder)
C
C     copy to arrays
C
      tmat(1:4,l,iq,0:nder)=t(1:4,0:nder)
      bigd(1,l,iq,0:nder)=ac(0:nder)
      bigd(2,l,iq,0:nder)=bd(0:nder)
      cd(1,l,iq,0:nder)=tcd(0:nder)
      cd(2,l,iq,0:nder)=dm1(0:nder)
   20 CONTINUE
      IF(nprn.EQ.1) THEN
         WRITE(m6,100)
         WRITE(m6,101)((iq,l,(tmat(k,l,iq,0),k=1,4),
     .                        l=0,lmax),iq=1,nq)
         WRITE(m6,102)
         DO 50 jder=1,nder
   50    WRITE(m6,101)((iq,l,(tmat(k,l,iq,jder),k=1,4),
     .                         l=0,lmax),iq=1,nq)
      ENDIF
C
      DEALLOCATE(tcd,cm1,dm1,t)
      DEALLOCATE(ac,bd,jntmp,al)
      DEALLOCATE(nm1,dn,dj,j,n)
      DEALLOCATE(fi,gi)
C
100   FORMAT(/' TRMTRX:   IQ  L',6x,'A',12x,'B',12x,'C',12x,'D',/)
101   FORMAT(1000(10x,2i3,4f13.8/))
102   FORMAT(' TRMTRX:   IQ  L',6x,'ADOT',9x,'BDOT',9x,'CDOT',9x,
     .       'DDOT',/)
103   FORMAT( ' TRMTRX:   itrans set to',i2)
      RETURN
      END
