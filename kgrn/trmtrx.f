      SUBROUTINE trmtrx(kap2,sigma,t,l,itrans,npr)
C   ******************************************************************
C   *   Makes the transformation matrix t = {t1,t2,t3,t4}.           *
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
C   *             the Wronskian W{N^a,J^a} equals W{N^0,J^0}=w/2     *
C   *      2:head |N^a(kappa)> has value 1 and derivative 0          *
C   *             |J^a(kappa)> has value 0 and derivative 1/a*a      *
C   *                          at hard core radius a=sigma*wst       *
C   *      3:head |N^a(kappa)> has value 1 and derivative 0          *
C   *             |J^a(kappa)> has value 0 and derivative -1/a       *
C   *                          at hard core radius a=sigma*wst       *
C   *   Output:                                                      *
C   *          tmat  :transformation matrix for head and tail        *
C   *          dfac  :t1*t4 - t2*t3   detrminant of tmat             *
C   *                                                                *
C   *   Remarks: The transformation matrix is:                       *
C   *                                                                *
C   *            |N^a>=tmat(1)*|N>+tmat(2)*|J>                       *
C   *            |J^a>=tmat(3)*|N>+tmat(4)*|J>                       *
C   *                                                                *
C   *        and ||N^a>=|N^a>-|J^a>S^a                               *
C   ******************************************************************
      USE control_data
      USE lattice
      USE message
      USE radialmesh
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(-1:lmax+1) :: fi, gi
      REAL(KIND=8), DIMENSION(4,0:1)     :: t
      REAL(KIND=8)  :: j, n, dj, dn, jd, nd, djd, dnd
      REAL(KIND=8)  :: kappa, det, dett, kap2, sigma
      REAL(KIND=8)  :: r2, rfct, rfct2, hcr, fct1, fct2
      REAL(KIND=8), PARAMETER ::  tol = 1.d-7
      INTEGER       :: l, npr, itrans
C
      IF(npr.EQ.1) WRITE(m6,100)
      hcr=sigma
      r2=hcr*hcr
      kappa=kap2*r2
      CALL bessl(kappa,-1,l+1,fi(-1),gi(-1))
      rfct=hcr**(-l-1)
      rfct2=hcr**l
C
      n  = rfct*gi(l)
      nd =n*(r2/2.)*(gi(l-1)/gi(l))/(2*(l-1)+1)
      dn= rfct/hcr*gi(l+1)
      dnd=dn*(r2/2.)*(gi(l)/gi(l+1))/(2*(l-1)+3)
      dn=(l*n/hcr -(2*l+1)*dn)/sws
      dnd=(l*nd/hcr -(2*l+1)*dnd)/sws
      j = rfct2*fi(l)
      jd=-j*(r2/2.)*(fi(l+1)/fi(l))/(2*(l+1)-1)
      dj= rfct2/hcr*fi(l-1)
      djd=-dj*(r2/2.)*(fi(l)/fi(l-1))/(2*(l+1)-3)
      dj=(-(l+1)*j/hcr +(2*l-1)*dj)/sws
      djd=(-(l+1)*jd/hcr +(2*l-1)*djd)/sws
C
      IF(itrans.EQ.2) THEN
         fct2 = 2.d0/sws
         det  = 2.d0/sws
      ELSEIF(itrans.EQ.3) THEN
         fct2 = -2.d0*hcr
         det  = -2.d0*hcr
      ENDIF
C
      IF(itrans.EQ.0.OR.itrans.EQ.1) THEN
         det  = 1.d0
         t(1,0) = 1.d0
         t(1,1) = 0.d0
         t(2,0:1) = 0.d0
         t(3,0) = -j/n
         t(3,1) = -j/n*(jd/j-nd/n)
         t(4,0) = 1.d0
         t(4,1) = 0.d0
      ELSE
         fct1 =2.*hcr*hcr*sws
         t(1,0) = fct1*dj
         t(1,1) = fct1*djd
         t(2,0) = -fct1*dn
         t(2,1) = -fct1*dnd
         t(3,0) = -fct2*j
         t(3,1) = -fct2*jd
         t(4,0) = fct2*n
         t(4,1) = fct2*nd
      ENDIF
      dett=t(1,0)*t(4,0)-t(2,0)*t(3,0)
      IF(ABS(dett-det).GT.tol) WRITE(m6,'(4f10.6)') dett,det
      dett=t(1,0)*t(4,1)+t(1,1)*t(4,0)-t(2,0)*t(3,1)-t(2,1)*t(3,0)
      IF(ABS(dett).GT.tol) WRITE(m6,'(4f10.6)') dett
C
      IF(npr.EQ.1) WRITE(m6,101) l,t
   20 CONTINUE
C
100   FORMAT(/' TRMTRX:   IT  L',8x,'A',13x,'B',13x,'C',13x,'D',/)
101   FORMAT(10x,i3,8f7.4,/,16x,8f7.4,/)
103   FORMAT( ' TRMTRX:   itrans set to',i2)
      RETURN
      END
