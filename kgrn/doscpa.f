      SUBROUTINE doscpa(is,bg0,gi,dos,zm,nzm)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the second part of the density of state G(z):     *
C   *                                                                *
C   *    G(z) = g(z) * S^dot(z) - sum(sort) c(sort) gi(z) * Di^dot   *
C   *                                                                *
C   * Notation:                                                      *
C   *                                                                *
C   *   bg0  ---->  g(z) * S^dot(z)                                  *
C   *   gi   ---->  gi(z) Green's function for alloy components.     *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE logderivative
      USE partialwaves
      IMPLICIT NONE
      INTEGER :: is, nzm, i, j, lz, iq
      INTEGER :: it, ita, l, m, lm, twol
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,ns) :: bg0
      COMPLEX(KIND=8), DIMENSION(nzm,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(nzm,mnta,nt,0:lmax) :: dos
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8) :: diag, dfd
      REAL(KIND=8)    :: epol, npol
C
C     Find density of state from bg0(z) and gi(z)
C
      dos=zero
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 l=0,lmax
      DO 20 m=-l,l
      lm=l*l+l+m+1
      DO 20 lz=1,nzm
      DO 20 ita=1,nta(it)
      dfd=dfi(lz,ita,it,l,is,2)
      diag=bg0(lz,iq,lm,is)-gi(lz,ita,iq,lm,lm,is)*dfd
C
      dos(lz,ita,it,l)=dos(lz,ita,it,l)+diag/mmt(it)
   20 CONTINUE
C
C     Extract the poles of D_dot
C
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      DO 30 l=0,lmax
      twol=2*l+1
      DO 30 lz=1,nzm
   30 dos(lz,ita,it,l)=dos(lz,ita,it,l)-twol*gsng(lz,ita,it,l,is)
C
      DO 31 it=1,nt
      DO 31 ita=1,nta(it)
      DO 31 l=0,lmax
      twol=2*l+1
      DO 31 j=1,necr(l,ita,it,is)
      epol=ecr(l,ita,it,is,j)
      npol=nocr(l,ita,it,is,j)
      DO 31 lz=1,nzm
      diag=npol/(zm(lz)-epol)
   31 dos(lz,ita,it,l)=dos(lz,ita,it,l)+twol*diag
C
      RETURN
      END
