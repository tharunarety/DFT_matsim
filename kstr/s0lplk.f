      SUBROUTINE s0lplk(nlmp)
C   ******************************************************************
C   * Calculate one-center expansion coefficents to J of N           *
C   * and their energy derivative. The expansion theorem:            *
C   *                                                                *
C   *       N_{RL}(k^2,r) = - J_{R'L'}(k^2,r) * S_{R'L',RL}          *
C   *                                                                *
C   * where : N_RL(r)     = N_L(r-R)                                 *
C   *         N_L(k^2,r)  = - (kw)^(l+1)/(2l-1)!! n_l(kr) Y_L(r)     *
C   *         J_RL(r)     = J_L(r-R)                                 *
C   *         J_L(k^2,r)  =   (kw)^(-l) *(2l-1)!!/2 j_l(kr) Y_L(r)   *
C   *                                                                *
C   * n_l(x) and j_l(x) are the Neumann (Hankel for k^2<0) and Bessel*
C   * functions of fractional order and                              *
C   *                                                                *
C   * S_R'L'RL = -8 pi sum_L'' I_{LL'L""} * [-(kw)^2]^{(l'+l-l'')/2} *
C   *                                                                *
C   *     (-1)^l*(2l''-1)!!/(2l'-1)!!/(2l-1)!! * N_L''(k^2,R-R').    *
C   *                                                                *
C   * I_{LL'L""} <> 0 for even l+l'+l'' => S_R'L'RL are real !!      *
C   *                                                                *
C   ******************************************************************
      USE bessel
      USE factorial
      USE control_data
      USE lmtolm
      USE realgaunt
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: rgnt
      INTEGER :: nlmp, i, lmp, lp, lm, l, lmpp, ip, j, k
C
      sbare(1:nlmp,1:nlm,0:nder)=0.d0
C
C     Upper triangular part
C
      DO 20 i=1,ngaunt
      lmp=lmpg(i)
      if(lmp.gt.nlmp) GO TO 20
      lp=llx(lmp)
      lm=lmg(i)
      l=llx(lm)
      lmpp=lmppg(i)
      rgnt=gnt(i)
      ip=kpow(i)
C
      sc(0:nder)=0.d0
      DO 21 j=0,nder
      DO 21 k=0,j
   21 sc(j)=sc(j)+ifib(j,k)*edot(ip,k)*hank(lmpp,j-k)*rgnt
      sbare(lmp,lm,0:nder)=sbare(lmp,lm,0:nder)+
     .                     fack(lp,l)*sc(0:nder)
   20 CONTINUE
C
C     Lower triangular part
C
      DO 22 lmp=1,nlm
      lp=llx(lmp)
      DO 22 lm=lmp+1,nlm
      l=llx(lm)
   22 sbare(lm,lmp,0:nder)=sbare(lmp,lm,0:nder)*sig(l)*sig(lp)
C
      DEALLOCATE(hank)
      RETURN
      END
