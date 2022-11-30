      SUBROUTINE bessl(y,lmin,lmax,fi,gi)
C   ******************************************************************
C   *                                                                *
C   *   Calculate spherical Bessel functions j, Neumann functions n. *
C   *                                                                *
C   *  *On entry:                                                    *
C   *                                                                *
C   *   y     : = e*r**2 = (kappa*r)**2 = x**2                       *
C   *   lmin  : minimum l, may be negative                           *
C   *   lmax  : maximum l                                            *
C   *                                                                *
C   *  *On exit:                                                     *
C   *                                                                *
C   *   fi    : fi(l)=j_l(x)/x^l*(2l-1)!!/2     ; l=lmin,.,lmax,     *
C   *   gi    : gi(l)=-n_l(x)*x^(l+1)/(2l-1)!!  ; l=lmin,.,lmax,     *
C   *                                                                *
C   *   The Bessel functions are calculated for lmax and lmax-1 by   *
C   *   an expansion in powers of x^2=y:                             *
C   *                                                                *
C   *                 (-x^2/2)^k                                     *
C   *   fi =  Sum_k  --------------  = dum(lmx+1-l)                  *
C   *                k! (2l+1+2k)!!                                  *
C   *                                                                *
C   *   The remaining Bessel functions are obtained by recursion:    *
C   *                                                                *
C   *   j_{l-2}(x)=(2l-1)j_{l-1}(x)/x-j{l}(x); l=lmx-2,..,-lmx-1     *
C   *                                                                *
C   *   and the Neumann function are then given by                   *
C   *                                                                *
C   *   n_l(x)=j_{-l-1}*(-1)^{l+1}                                   *
C   *                                                                *
C   * Warning:  For x > 40 this algorithm is numerically unstable !! *
C   * Warning:  for lmax < -1 program must be checked                *
C   *                                                                *
C   * Routine based on  LMASA-52. Cleaned by HLS 12/1-98.            *
C   *                                                                *
C   ******************************************************************
      USE bessl_fact
      IMPLICIT NONE
      INTEGER :: lmin, lmax, nf, tlp1, ll1, ll2
      INTEGER :: i, isn, j1,j 2, k, l, lmx, lmxp1, lmxp2
      REAL(KIND=8) :: y, dt, dt2, my, t
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dum
      REAL(KIND=8), DIMENSION(lmin:*)         :: fi, gi
C
      lmx=MAX0(lmax,-lmin,2)
      if(lmx.gt.lmaxi.or.lmin.lt.lmini) write(*,*)
     .  ' BESSL:!!!',lmax,lmx,lmaxi,lmin,lmini
C
C     Treat kappa = 0 as special case
C
      IF(DABS(y).LT.tol) THEN
        DO l=lmin,lmax
           fi(l)=1.d0/(4*l+2)
           gi(l)=1.d0
        ENDDO
        RETURN
      ENDIF
C
      ALLOCATE(dum(lmx*2+2))
C
C     Construct fi_(lmax) = dum(1)=j_(lmax)(x)/x^{lmax} from series 
C     expansion
C
      my=-y
      tlp1=lmx+lmx+1
      dt=1.d0
      t=1.d0
      i=0
      DO WHILE (DABS(dt).GE.tol)
         i=i+2
         dt2=i+tlp1
         dt=dt*my/(i*dt2)
         t=t+dt
      ENDDO
      dum(1)=t/facbl(lmx+1)
C
C     Construct fi_(lmax-1) = dum(2)=j_(lmax-1)(x)/x^{lmax-1} from series 
C     expansion
C
      tlp1=tlp1-2
      dt=1.d0
      t=1.d0
      i=0
      DO WHILE (DABS(dt).GE.tol)
         i=i+2
         dt2=i+tlp1
         dt=dt*my/(i*dt2)
         t=t+dt
      ENDDO
      dum(2)=t/facbl(lmx)
C
C     Recursion for dum(k)=j_{lmx+1-k}(x)/x^{lmx+1-k}=fi(lmx+1-k)
C
      ll1=lmx+lmx+1
      ll2=ll1+1
      nf=ll1
      DO k=3,ll2
         nf=nf-2
         dum(k)=nf*dum(k-1)-y*dum(k-2)
      ENDDO
C
C     Get fi from dum: fi_l = dum(lmax+1-l) 
C     and gi from fi : gi_l = fi_(-l-1)*(-1)**(l+1) 
C  
      lmxp1=lmx+1
      lmxp2=lmx+2
      isn=(-1)**lmin
      DO k=lmin,lmax
         j1=lmxp1-k
         j2=lmxp2+k
         fi(k)=dum(j1)
         gi(k)=dum(j2)*isn
         isn=-isn
      ENDDO
C
C     Normalize Bessel and Neumann functions
C           
      DO l=lmin,lmax
         fi(l)=fi(l)*facbl(l)*0.5d0
         gi(l)=gi(l)/facbl(l)
      ENDDO
C
      DEALLOCATE(dum)
      RETURN
      END
