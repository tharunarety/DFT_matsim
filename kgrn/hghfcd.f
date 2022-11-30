      SUBROUTINE hghfcd(is,gah,ss,ndera,lz,nzm)
C   ******************************************************************
C   *                                                                *
C   * Set up the following matrices:                                 *
C   *                                                                *
C   *            g_L L ;                                             *
C   *            g_L H = SUM_L'(g_LL'*S_L'H);                        *
C   *            g_H L = SUM_L'(S_HL'*g_L'L);                        *
C   *            g_H'H = SUM_L'L(S_H'L'*g_L'L*S_LH);                 *
C   *                                                                *
C   *   where g is the path operator and S the slope matrix.         *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE kinkmatrix ; USE message
      USE slope        ; USE taylor
      IMPLICIT NONE
      INTEGER         :: nzm, is, lz, jd, lm, lmp, ndera
      INTEGER         :: iq, il, ilh, jq, jl, jlh, kq, kl, klh
      COMPLEX(KIND=8), DIMENSION(nzm,nq,ns,nlmf,nlmf) :: gah
      COMPLEX(KIND=8), DIMENSION(nlmqf,nlmq,0:ndera) :: ss
      COMPLEX(KIND=8) :: smm
C
C     Taylor expansion for the slope matrix S^a
C
      DO 10 iq=1,nq
      il=(iq-1)*nlm
      ilh=(iq-1)*nlmf
      DO 10 lm=1,nlm
      DO 10 jq=1,nq
      jlh=(jq-1)*nlmf
      DO 10 lmp=1,nlmf
      slop(jlh+lmp,ilh+lm)=ss(jlh+lmp,il+lm,0)
      slop(jlh+lmp,ilh+lm)=slop(jlh+lmp,ilh+lm)+
     .   SUM(tayl(lz,1:ndera,is)*ss(jlh+lmp,il+lm,1:ndera))
   10 CONTINUE
      DO 11 iq=1,nq
      ilh=(iq-1)*nlmf
      DO 11 lm=nlm+1,nlmf
      DO 11 jq=1,nq
      jl=(jq-1)*nlm
      jlh=(jq-1)*nlmf
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
      jlh=(jq-1)*nlmf
      DO 20 lmp=nlm+1,nlmf
      smm=zero
      DO 21 kq=1,nq
      kl =(kq-1)*nlm
      klh=(kq-1)*nlmf
   21 smm=smm+SUM(slop(jlh+lmp,klh+1:klh+nlm)*gak(kl+1:kl+nlm,il+lm))
   20 tmp(jlh+lmp,il+lm)=smm
C
      DO 30 iq=1,nq
      il =(iq-1)*nlm
      ilh=(iq-1)*nlmf
C
      DO 31 lm=1,nlm
      DO 31 lmp=1,nlm
   31 gah(lz,iq,is,lmp,lm)=gah(lz,iq,is,lmp,lm)+
     .                     gak(il+lmp,il+lm)
C
      DO 32 lm=1,nlm
      DO 32 lmp=nlm+1,nlmf
   32 gah(lz,iq,is,lmp,lm)=gah(lz,iq,is,lmp,lm)+tmp(ilh+lmp,il+lm)
C
      DO 33 lm=nlm+1,nlmf
      DO 33 lmp=1,nlm
      smm=zero
      DO 34 kq=1,nq
      kl =(kq-1)*nlm
      klh=(kq-1)*nlmf
   34 smm=smm+SUM(gak(il+lmp,kl+1:kl+nlm)*slop(klh+1:klh+nlm,ilh+lm))
   33 gah(lz,iq,is,lmp,lm)=gah(lz,iq,is,lmp,lm)+smm
C
      DO 35 lm=nlm+1,nlmf
      DO 35 lmp=nlm+1,nlmf
      smm=zero
      DO 36 kq=1,nq
      kl =(kq-1)*nlm
      klh=(kq-1)*nlmf
   36 smm=smm+SUM(tmp(ilh+lmp,kl+1:kl+nlm)*slop(klh+1:klh+nlm,ilh+lm))
   35 gah(lz,iq,is,lmp,lm)=gah(lz,iq,is,lmp,lm)+smm
C
   30 CONTINUE
C
      RETURN
      END
