      SUBROUTINE trmtrz(kap2,sigma,tmat,itrans,dimnta,npr)
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
      INTEGER          :: dimnta, it, ita, l, lmin, npr, itrans, ntait
      COMPLEX(KIND=8), DIMENSION(-1:lmax+1) :: fi, gi
      COMPLEX(KIND=8), DIMENSION(4,0:1)     :: t
      COMPLEX(KIND=8), DIMENSION(4,0:lmax,dimnta,nt,0:1)  :: tmat
      COMPLEX(KIND=8)  :: j, n, dj, dn, jd, nd, djd, dnd
      COMPLEX(KIND=8)  :: kappa, det, dett, kap2
      REAL(KIND=8),    DIMENSION(0:lmax,mnta,nt)   :: sigma
      REAL(KIND=8)     :: r2, rfct, rfct2, hcr, fct1, fct2
      REAL(KIND=8),    PARAMETER ::  tol = 1.d-7
C
      lmin=-1
C
      IF(npr.EQ.1) WRITE(m6,100)
      ntait=1
      DO 20 it=1,nt
      IF(dimnta.GT.1) ntait=nta(it)
      DO 20 ita=1,ntait
      DO 20 l=0,lmax
      hcr=sigma(l,ita,it)
      r2=hcr*hcr
      kappa=kap2*r2
      CALL besslz(kappa,lmin,l+1,fi(lmin),gi(lmin))
      rfct=hcr**(-l-1)
      rfct2=hcr**l
C
      n  = rfct*gi(l)                                !! Eq. (4)
      nd =n*(r2/2.)*(gi(l-1)/gi(l))/(2*(l-1)+1)      !! Eq. (8) and (4)
      dn= rfct/hcr*gi(l+1)                           !! Eq. (4) for (l+1)
      dnd=dn*(r2/2.)*(gi(l)/gi(l+1))/(2*(l-1)+3)     !! Eq. (8) for (l+1)
      dn=(l*n/hcr -(2*l+1)*dn)/sws                   !! Eq. (12)
      dnd=(l*nd/hcr -(2*l+1)*dnd)/sws                !! Eq. (12) for .
      j = rfct2*fi(l)                                !! Eq. (3)
      jd=-j*(r2/2.)*(fi(l+1)/fi(l))/(2*(l+1)-1)      !! Eq. (7) and (3)
      dj= rfct2/hcr*fi(l-1)                          !! Eq. (3) for (l-1)
      djd=-dj*(r2/2.)*(fi(l)/fi(l-1))/(2*(l+1)-3)    !! Eq. (7) for (l-1)
      dj=(-(l+1)*j/hcr +(2*l-1)*dj)/sws              !! Eq. (11)
      djd=(-(l+1)*jd/hcr +(2*l-1)*djd)/sws           !! Eq. (11) for .
C
      IF(itrans.EQ.2) THEN
         fct2 = 2.d0/sws
         det  = CMPLX(2.d0/sws,0.d0,8)
      ELSEIF(itrans.EQ.3) THEN
         fct2 = -2.d0*hcr
         det  = CMPLX(-2.d0*hcr,0.d0,8)
      ENDIF
C
      IF(itrans.EQ.0.OR.itrans.EQ.1) THEN
         det  = CMPLX(1.d0,0.d0,8)
         t(1,0) = CMPLX(1.d0,0.d0,8)
         t(1,1) = CMPLX(0.d0,0.d0,8)
         t(2,0:1) = CMPLX(0.d0,0.d0,8)
         t(3,0) = -j/n
         t(3,1) = -j/n*(jd/j-nd/n)
         t(4,0) = CMPLX(1.d0,0.d0,8)
         t(4,1) = CMPLX(0.d0,0.d0,8)
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
      IF(ABS(dett-det).GT.tol)
     .                WRITE(m6,'(a,4f10.6)') 'TRMTRZ: ',dett,det
      dett=t(1,0)*t(4,1)+t(1,1)*t(4,0)-t(2,0)*t(3,1)-t(2,1)*t(3,0)
      IF(ABS(dett).GT.tol) WRITE(m6,'(a,4f10.6)') 'TRMTRZ: ',dett
C
      IF(npr.EQ.1) WRITE(m6,101) it,l,t
      tmat(1:4,l,ita,it,0:1)=t(1:4,0:1)
   20 CONTINUE
C
100   FORMAT(/' TRMTRZ:   IT  L',8x,'A',13x,'B',13x,'C',13x,'D',/)
101   FORMAT(10x,2i3,8f7.4,/,16x,8f7.4,/)
103   FORMAT( ' TRMTRZ:   itrans set to',i2)
      RETURN
      END
