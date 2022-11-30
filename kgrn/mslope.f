      SUBROUTINE mslope(nlmact)
C   ******************************************************************
C   *                                                                *
C   *    Fill up the lower part of the slope matrix.                 *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE control_text
      USE slope        ; USE taylor
      IMPLICIT NONE
      INTEGER :: nlmact, iq, il, jq, jl, jlh, lm, lmp
C
      DO 20 iq=1,nq
      il =(iq-1)*nlm
      DO 20 jq=1,nq
      jl =(jq-1)*nlm
      jlh=(jq-1)*nlmact
      DO 20 lm=1,nlm
      DO 20 lmp=1,nlm
      sal(jl+lmp,il+lm,0:nder)=sa(jlh+lmp,il+lm,0:nder)
   20 CONTINUE
C
      IF(expan.EQ.'S') RETURN
      DO 30 iq=1,nq
      il =(iq-1)*nlm
      DO 30 jq=1,nq
      jl =(jq-1)*nlm
      jlh=(jq-1)*nlmact
      DO 30 lm=1,nlm
      DO 30 lmp=1,nlm
      sapl(jl+lmp,il+lm,0:nderp)=sap(jlh+lmp,il+lm,0:nderp)
   30 CONTINUE
C
      RETURN
      END
