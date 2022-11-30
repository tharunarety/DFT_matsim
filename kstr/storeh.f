      SUBROUTINE storeh(nrq)
C   ******************************************************************
C   *                                                                *
C   *  Store structure constants S(R'H'RL) for R' from the unit cell.*
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
   20 WRITE(2) ((sah(lp,l,jd),lp=1,nlmh),l=l1,l2)
      IF(msgl.EQ.1) WRITE(msgio,'(/,a,i5)')
     .         ' Store high-low slope matrices'
C
      WRITE(m6,100)
  100 FORMAT(/,' STOREH:   Slope matrices high-low stored on',
     .       ' FOR002')
      RETURN
      END
