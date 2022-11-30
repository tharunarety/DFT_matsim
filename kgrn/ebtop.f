      SUBROUTINE ebtop(iprn)
C   ******************************************************************
C   *                                                                *
C   *    Establish the position and width of the energy bands based  *
C   *    on the potential parameters.                                *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE botomtop
      USE control_data
      USE message
      USE potparam
      IMPLICIT NONE
      REAL(KIND=8) :: ebot, etop
      INTEGER      :: iprn, it, ita, is, ip, l
C
      ebot=100.
      etop=-100.
      DO 21 is=1,ns
      DO 21 ip=1,pan
      DO 21 l=0,lmax
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      IF(conc(ita,it).GT.1.d-8) THEN
         IF(bot(l,ita,it,is,ip).LT.ebot) ebot=bot(l,ita,it,is,ip)
         IF(top(l,ita,it,is,ip).GT.etop) etop=top(l,ita,it,is,ip)
      ENDIF
   21 CONTINUE
      ebt=ebot
      etp=etop
      IF(iprn.EQ.1) THEN
         WRITE(m6,101)
         WRITE(m6,102) ebt,etp
      ENDIF
C
  101 FORMAT(/,' EBTOP:      Bot       Top',/)
  102 FORMAT(11x,2f10.6)
      RETURN
      END
