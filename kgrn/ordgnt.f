      SUBROUTINE ordgnt(nlmf)
C   ******************************************************************
C   *                                                                *
C   * Order the real Gaunt numbers after (lp,l) and lpp.             *
C   *                                                                *
C   ******************************************************************
      USE realgaunt
      IMPLICIT NONE
      INTEGER :: nlmf, ign, lmp, lm, dimg
C
      ALLOCATE(ngnt(nlmf,nlmf))
C
      ngnt=0
      DO 20 ign=1,ngaunt
      lmp=lmpg(ign)
      lm=lmg(ign)
      ngnt(lmp,lm)=ngnt(lmp,lm)+1
   20 CONTINUE
C
      dimg=MAXVAL(ngnt)
      ALLOCATE(ignt(nlmf,nlmf,dimg))
C
      ngnt=0
      DO 21 ign=1,ngaunt
      lmp=lmpg(ign)
      lm=lmg(ign)
      ngnt(lmp,lm)=ngnt(lmp,lm)+1
      ignt(lmp,lm,ngnt(lmp,lm))=ign
   21 CONTINUE
C
      RETURN
      END
