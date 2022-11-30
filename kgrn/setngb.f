      SUBROUTINE setngb
C   ******************************************************************
C   *                                                                *
C   *    Set up the Madelung potential correction defined as:        *
C   *                                                                *
C   *    MADC = alpha*2/S (Ry units)                                 *
C   *                                                                *
C   *    where:                                                      *
C   *           w      is the atomic sphere radius                   *
C   *                                                                *
C   *           alpha = NN(fcc)/w .                                  *
C   *                                                                *
C   *           NN(fcc) nearest neighbour distance in fcc lattice.   *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE moments ; USE message ; USE radialmesh
      IMPLICIT NONE
      INTEGER, PARAMETER :: npr = 1
      REAL(KIND=8), PARAMETER :: alpha = 0.5526694571d0
      INTEGER :: iq, it, ita
C
      ALLOCATE(madc(mnta,nq))
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
   20 madc(ita,iq)=2.d0*alpha/ws(ita,it)
      IF(npr.EQ.1) THEN
         WRITE(m6,100)
         WRITE(m6,101) (madc(1:nta(itq(iq)),iq),iq=1,nq)
      ENDIF
C
  100 FORMAT(/,' SETNGB:   Madelung correction for each site and sort')
  101 FORMAT(/,11x,8f10.6)
      RETURN
      END
