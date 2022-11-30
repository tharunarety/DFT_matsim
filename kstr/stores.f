      SUBROUTINE stores(nrq)
C   ******************************************************************
C   *                                                                *
C   *  Store structure constants S(R'L'RL) for R' from the unit cell.*
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE message
      USE temporary
      IMPLICIT NONE
      INTEGER :: jd, jrl, ir, nrq, l1, l2, lp, l
C
      DO 20 jd=0,nder
      jrl=-nlm
      DO 20 ir=1,nrq
      jrl=jrl+nlm
      l1=jrl+1
      l2=jrl+nlm
   20 WRITE(1) ((sa(lp,l,jd),lp=1,nlm),l=l1,l2)
      IF(msgl.EQ.1) WRITE(msgio,'(/,a,i5)')
     .         ' Store low-low slope matrices'
C
      WRITE(m6,100)
  100 FORMAT(/,' STORES:   Slope matrices low-low stored on',
     .       ' FOR001')
      RETURN
      END
