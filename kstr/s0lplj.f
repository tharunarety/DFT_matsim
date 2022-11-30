      SUBROUTINE s0lplj(r)
C   ******************************************************************
C   * Calculate one-center expansion coefficents to J of J           *
C   * and their energy derivative. The expansion theorem:            *
C   *                                                                *
C   *       J_{RL}(k^2,r) = - J_{R'L'}(k^2,r) * S_{R'L',RL}          *
C   *                                                                *
C   * where : J_RL(r)     = J_L(r-R)                                 *
C   *         J_L(k^2,r)  =   (kw)^(-l) *(2l-1)!!/2 j_l(kr) Y_L(r)   *
C   *                                                                *
C   * j_l(x) are the Bessel functions of fractional order and        *
C   *                                                                *
C   * S_R'L'RL = -8 pi sum_L'' I_{LL'L""} * [-(kw)^2]^{(l'-l+l'')/2} *
C   *                                                                *
C   *    (-1)^l''*(2l-1)!!/(2l'-1)!!/(2l''-1)!! * J_L''(k^2,R-R').   *
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
      REAL(KIND=8) :: r,rgnt
      INTEGER :: i, lmp, lp, lm, l, lmpp, ip, j, k
C
      sbare(1:nlm,1:nlmw,0:nder)=0.d0
C
      IF(r.LT.err) THEN
        DO 10 lmp=1,nlm
   10   sbare(lmp,lmp,0)=-1.d0
        RETURN
      ENDIF
      DO 20 i=1,ngauntw
      lmp=lmpgw(i)
      lp=llx(lmp)
      lm=lmgw(i)
      l=llx(lm)
      lmpp=lmppgw(i)
      rgnt=gntw(i)
      ip=jpow(i)
C
      sc(0:nder)=0.d0
      DO 21 j=0,nder
      DO 21 k=0,j
   21 sc(j)=sc(j)+ifib(j,k)*edot(ip,k)*bess(lmpp,j-k)*rgnt
      sbare(lmp,lm,0:nder)=sbare(lmp,lm,0:nder)+
     .                     facj(lp,l)*sc(0:nder)
   20 CONTINUE
C
      DEALLOCATE(bess)
      RETURN
      END
