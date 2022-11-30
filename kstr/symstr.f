      SUBROUTINE symstr(ir1,ir2)
C   ******************************************************************
C   *                                                                *
C   *  Symmetrize structure constants S(R'L'RL) = S(RLR'L') .        *
C   *                                                                *
C   ******************************************************************
      USE factorial
      USE control_data
      USE lmtolm
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: siglpl, temp
      INTEGER :: ir1, ir2, jrl, jr, l1, l2, lp, l, i, j, jd
C
      return
      IF(itrans.EQ.3) RETURN
      jrl=-nlm
      DO 20 jr=ir1,ir2
      jrl=jrl+nlm
      l2=jrl+nlm
      DO 20 lp=1,nlm
      l1=jrl+lp+1
      i=lp
      DO 20 j=l1,l2
      l=j-jrl
      siglpl=sig(llx(l))*sig(llx(lp))
      DO 20 jd=0,nder
      temp=(sa(i,j,jd)+siglpl*sa(l,jrl+lp,jd))/2.d0
      sa(i,j,jd)=temp
   20 sa(l,jrl+lp,jd)=siglpl*temp
C
      RETURN
      END
