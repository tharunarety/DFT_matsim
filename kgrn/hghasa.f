      SUBROUTINE hghasa(is,hgh,ss,ndera,lz,nzm,weight)
C   ******************************************************************
C   *                                                                *
C   * Set up the following matrices:                                 *
C   *                                                                *
C   *            g_L L ;                                             *
C   *            g_H H = SUM_L'L(S_H L'*g_L'L*S_LH);                 *
C   *                                                                *
C   *   where g is the path operator and S the slope matrix.         *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE kinkmatrix ; USE slope ; USE taylor
      IMPLICIT NONE
      INTEGER         :: nzm, is, lz, l, lm, lmp, m, ndera
      INTEGER         :: iq, il, ilh, jq, jl, jlh, kq, kl, klh
      COMPLEX(KIND=8), DIMENSION(nzm,nq,ns,lmax+1:lmaxt) :: hgh
      COMPLEX(KIND=8), DIMENSION(nlmqt,nlmq,0:ndera) :: ss
      COMPLEX(KIND=8) :: smm
      REAL(KIND=8) :: weight
C
      IF(lmaxt.LE.lmax) RETURN
C
C     Taylor expansion for the slope matrix S^a
C
      DO 10 iq=1,nq
      il=(iq-1)*nlm
      ilh=(iq-1)*nlmt
      DO 10 lm=1,nlm
      DO 10 jq=1,nq
      jlh=(jq-1)*nlmt
      DO 10 lmp=nlm+1,nlmt
      slop(jlh+lmp,ilh+lm)=ss(jlh+lmp,il+lm,0)
      slop(jlh+lmp,ilh+lm)=slop(jlh+lmp,ilh+lm)+
     .   SUM(tayl(lz,1:ndera,is)*ss(jlh+lmp,il+lm,1:ndera))
   10 CONTINUE
      DO 11 iq=1,nq
      ilh=(iq-1)*nlmt
      DO 11 lm=nlm+1,nlmt
      DO 11 jq=1,nq
      jl=(jq-1)*nlm
      jlh=(jq-1)*nlmt
      DO 11 lmp=1,nlm
      slop(jlh+lmp,ilh+lm)=CONJG(ss(ilh+lm,jl+lmp,0))
      slop(jlh+lmp,ilh+lm)=slop(jlh+lmp,ilh+lm)+
     .   SUM(tayl(lz,1:ndera,is)*CONJG(ss(ilh+lm,jl+lmp,1:ndera)))
   11 CONTINUE
C
      DO 20 iq=1,nq
      il =(iq-1)*nlm
      DO 20 lm=1,nlm
      DO 20 jq=1,nq
      jlh=(jq-1)*nlmt
      DO 20 lmp=nlm+1,nlmt
      smm=zero
      DO 21 kq=1,nq
      kl =(kq-1)*nlm
      klh=(kq-1)*nlmt
   21 smm=smm+SUM(slop(jlh+lmp,klh+1:klh+nlm)*gak(kl+1:kl+nlm,il+lm))
   20 tmp(jlh+lmp,il+lm)=smm
C
      DO 30 iq=1,nq
      il =(iq-1)*nlm
      ilh=(iq-1)*nlmt
      DO 30 l=lmax+1,lmaxt
      smm=zero
      DO 31 m=-l,l
      lm=l*l+l+m+1
      DO 31 kq=1,nq
      kl =(kq-1)*nlm
      klh=(kq-1)*nlmt
   31 smm=smm+SUM(tmp(ilh+lm,kl+1:kl+nlm)*slop(klh+1:klh+nlm,ilh+lm))
   30 hgh(lz,iq,is,l)=hgh(lz,iq,is,l)+weight*smm
C
      RETURN
      END
