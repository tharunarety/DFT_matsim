      SUBROUTINE realhinit(lmax)
C   ******************************************************************
C   *                                                                *
C   *    Initialize normalization of real harmonics.                 *
C   *                                                                *
C   ******************************************************************
      USE realhr_norm
      IMPLICIT NONE
      INTEGER :: lmax, l, m, nmax
      REAL(KIND=8), PARAMETER :: pi=3.1415926535897932d0
C
      lmaxi=lmax
      nmax=2*lmax
      ALLOCATE(ylnorm(0:lmax,0:lmax))
      ALLOCATE(al(0:nmax),alinv(lmax),tlp1(0:lmax))
C
      DO 20 l=0,lmax
   20 ylnorm(0,l)=DSQRT((2*l+1)/(4.d0*pi))
C
      DO 21 m=1,lmax
      DO 21 l=m,lmax
   21 ylnorm(m,l)=ylnorm(m-1,l)/dsqrt((l-m+1.d0)*(l+m))
C
      al(0)=0.d0
      DO 22 l=1,nmax
   22 al(l)=l
      DO 23 l=1,lmax
   23 alinv(l)=1.d0/al(l)
C
      DO 24 l=0,lmax
   24 tlp1(l)=2*l+1
C
      RETURN
      END
