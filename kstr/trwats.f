      SUBROUTINE trwats
C   ******************************************************************
C   *   Makes the screening constants of the Watson orbitals:        *
C   *                                                                *
C   *     N^al_W(rwats)=N^0_W(rwats)-alwats*J^0_W(rwats)= 0          *
C   *                                                                *
C   *   Output:                                                      *
C   *          alwats:screening constants for Watson-sphere          *
C   *                 and  (kappa*ws)^2-derivative                   *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE factorial
      USE control_data
      USE lattice
      USE screening
      USE message
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: fi, gi, j, n
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: jm1, jntmp, al
      REAL(KIND=8) :: rwats1,r2,rfac,rfac2
      INTEGER :: i, jder, iq, l
C
      ALLOCATE(j(0:nder),n(0:nder),jm1(0:nder))
      ALLOCATE(jntmp(0:nder),al(0:nder))
      ALLOCATE(fi(-nder:nder+lmaxw),gi(-nder:nder+lmaxw))
C
      DO 20 iq=1,nq
      rwats1=rwats(iq)
      r2=rwats1*rwats1
      rfac=1.d0
      rfac2=1.d0/rwats1
C
      CALL bessh(kap2*r2,-nder,lmaxw+nder,fi(-nder),gi(-nder))
C
      DO  20 l=0,lmaxw
      rfac=rfac/rwats1
      rfac2=rfac2*rwats1
C
C     make n_l and its energy derivatives
C
      n(0)   = rfac*gi(l)
      DO 21 jder=1,nder
   21 n(jder)= n(jder-1)*(r2/2.)*(gi(l-jder)/gi(l-jder+1))/
     .         (2*(l-jder)+1)        
C
C     make j_l and its energy derivs
C
      j(0)   = rfac2*fi(l)
      DO 22 jder=1,nder
   22 j(jder)= -j(jder-1)*(r2/2.)*(fi(l+jder)/fi(l+jder-1))/
     .         (2*(l+jder)-1)
C
C     make inverse of j_l, and its energy derivs
C
      jm1(0)=1./j(0)
      jntmp(1:nder)=0.d0
      jm1(1:nder)=0.d0
      DO 23 jder=1,nder
      DO 24 i=1,jder
   24 jntmp(jder)=jntmp(jder)+j(i)*jm1(jder-i)*ifib(jder-1,i-1)
      DO 23 i=1,jder
   23 jm1(jder)=jm1(jder)-jntmp(i)*jm1(jder-i)*ifib(jder-1,i-1)
C
C     make alpha and its derivs
C
      al(0)=n(0)*jm1(0)
      al(1:nder)=0.d0
      DO 25 jder=1,nder
      DO 25 i=0,jder
   25 al(jder)=al(jder)+n(i)*jm1(jder-i)*ifib(jder,i)

      alwats(l,iq,0:nder)=al(0:nder)
   20 CONTINUE
      IF(nprn.EQ.1) THEN
         WRITE(m6,100)
         DO 30 iq=1,nq
         DO 30 l=0,lmaxw
         WRITE(m6,101) iq,l,rwats(iq),alwats(l,iq,0:nder)
   30    CONTINUE
      ENDIF

100   FORMAT(/,' TRWATS:   IQ L  RWATS/WS',4x,
     .         ' Screening constants of the Watson orbitals',/)
101   FORMAT(11x,2i2,f7.3,f8.4,f11.4,2e7.1e1,3e8.1e2)
      RETURN
      END
