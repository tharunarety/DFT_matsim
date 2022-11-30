      SUBROUTINE doscpa0(is,g0,gi,bgdos,dos,zm,nzm,nzlin)
C   ******************************************************************
C   *                                                                *
C   *   Calculate the density of state G(z) defined as:              *
C   *                                                                *
C   *   G(z) = g^a(z) * K^dot^a(z)  (not matrix product !!!)         *
C   *                                                                *
C   *   using the following equation:                                *
C   *                                                                *
C   *      G(z) = g^a(z) * [G0(z)/g^a0(z) + a*(D^a0 - D^a)^dot].     *
C   *                                                                *
C   *                      (not matrix product !!!)                  *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text ; USE logderivative
      USE message ; USE partialwaves
      IMPLICIT NONE
      INTEGER :: is, nzlin, i, j, lz, nzm, iq
      INTEGER :: it, ita, l, m, lm, twol
      COMPLEX(KIND=8), DIMENSION(nzlin,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: g0
      COMPLEX(KIND=8), DIMENSION(nzlin,nq,nlm) :: bgdos
      COMPLEX(KIND=8), DIMENSION(nzlin,mnta,nt,0:lmax) :: dos
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8) :: diag, dfd
      REAL(KIND=8)    :: epol, npol
C
C     Find density of state
C
      dos=zero
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 l=0,lmax
      DO 20 ita=1,nta(it)
      DO 20 lz=1,nzm
      dfd=dfi(lz,ita,it,l,is,2)
      DO 20 m=-l,l
      lm=l*l+l+m+1
      diag=bgdos(lz,iq,lm)-gi(lz,ita,iq,lm,lm,is)*dfd
C
   20 dos(lz,ita,it,l)=dos(lz,ita,it,l)+diag/mmt(it)
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
