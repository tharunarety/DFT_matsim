      SUBROUTINE poisson(rchd,nz,epot,potb,s,it,ita,jsr,lp)
C   ******************************************************************
C   *                                                                *
C   *    Solves the l-dependent Poisson equation. Generalized from   *
C   *    Loucks 'AUGMENTED PLANE WAVE METHOD', Benjamin, New York    *
C   *    (1967)  P.98                                                *
C   *                                                                *
C   *    INPUT  : lp   = l+1.                                        *
C   *             it   = Type of sublattice.                         *
C   *             ita  = Type of atom.                               *
C   *             nz   = Atomic number.                              *
C   *             potb = Potential at the sphere radius S.           *
C   *             rchd = 4.*Pi*Radius**2*Charge density.             *
C   *    OUTPUT : epot = Electrostatic potential from nucleus and    *
C   *                    electron charge density subjected to the    *
C   *                    dipole boundary condition.                  *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE poissonparam
      USE radialmesh
      IMPLICIT NONE
      INTEGER :: nz, jsr, jri, it, ita, lp, itop, j, jv, ir
      REAL(KIND=8), DIMENSION(*) :: rchd, epot
      REAL(KIND=8) :: s, twoz, d, dv, r, potb
C
      e(1)=0.d0
      f(1)=f1(lp)
      jri=jsr+2
      twoz=2.d0*nz
C
C     Solve Poisson
C
      itop=jri-1
      DO 21 j=2,itop
      d=c*SQRT(ri(j,ita,it))*(edl*rchd(j+1)+10.d0*rchd(j)+rchd(j-1)/edl)
      f(j)=c2(lp)-1.d0/f(j-1)
   21 e(j)=(d/a(lp)+e(j-1))/f(j)
C
      epot(jsr)=potb*SQRT(s)
C
      DO 22 j=1,jsr-1
      jv=jsr-j
   22 epot(jv)=e(jv)+epot(jv+1)/f(jv)
      DO 23 j=jsr,jri-1
   23 epot(j+1)=(epot(j)-e(j))*f(j)
C
C     Impose proper boundary condition
C
      dv=twoz/s
      DO 24 ir=1,jri
      r=ri(ir,ita,it)
   24 epot(ir)=epot(ir)/SQRT(r)-twoz/r+dv
C
      RETURN
      END
