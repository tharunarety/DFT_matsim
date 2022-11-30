      SUBROUTINE gtneum(lmax2,r,ux,uy,uz)
C   ******************************************************************
C   *                                                                *
C   *  Generate the Neumann function and energy derivatives.         *
C   *                                                                *
C   *  Input: lmax2     maximum l ;                                  *
C   *         r         radial distance in units of ws ;             *
C   *         u         unit vector .                                *
C   *  Output:hank(lm,0)  Nl(k^2,r)*Ylm*(2l-1)!!                     *
C   *         hank(lm,j)  the jth order derivative.                  *
C   *                     for Nl(k^2,r) see s0lplk.                  *
C   ******************************************************************
      USE bessel
      USE control_data
      USE factorial
      USE realharmonic
      IMPLICIT NONE
      REAL(KIND=8) :: r, ux, uy, uz, rfac, r2
      INTEGER :: lmax2, lmxnd, nlm2, l, m, lm, j
C
      lmxnd=lmax2+nder
      nlm2=(lmax2+1)*(lmax2+1)
      ALLOCATE(hank(nlm2,0:nder))
      ALLOCATE(fi(-nder:lmxnd),gi(-nder:lmxnd),gdot(0:nder))
      ALLOCATE(ylm(nlm2))
C
      CALL bessh(kap2*r*r,-nder,lmxnd,fi(-nder),gi(-nder))
      CALL realhr(ylm,lmax2,ux,uy,uz)
C
      rfac=1.d0
      r2=r*r/2.d0
      DO 20 l=0,lmax2
      gdot(0)=gi(l)
      DO 21 j=1,nder
   21 gdot(j)=gdot(j-1)*r2*(gi(l-j)/gi(l-j+1))/(2.d0*(l-j)+1.d0)
      rfac=rfac/r
      DO 20 m=-l,l
      lm=l*l+l+m+1
      hank(lm,0:nder)=gdot(0:nder)*ylm(lm)*rfac*fac2(l)
   20 CONTINUE
C
      DEALLOCATE(fi,gi,gdot)
      DEALLOCATE(ylm)
C
      RETURN
      END
