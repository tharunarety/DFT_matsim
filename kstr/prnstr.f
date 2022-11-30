      SUBROUTINE prnstr(iq,ir1,ir2)
C   ******************************************************************
C   *                                                                *
C   *    Print structure constants S(R'L'RL) for R' in the           *
C   *    unit cell.                                                  *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE lattice
      USE message
      USE temporary
      IMPLICIT NONE
      INTEGER :: iq, ir1, ir2, jrl, jr, l1, l2, lp, l, jd
C
      IF(nprn.NE.1) RETURN
C
      DO 20 jd=0,nder
      jrl=-nlm
      DO 20 jr=ir1,ir2
      jrl=jrl+nlm
      WRITE(m6,101) jd,jr,iq,rx(jr),ry(jr),rz(jr)
      l1=jrl+1
      l2=jrl+nlm
      DO 20 lp=1,nlm
   20 WRITE(m6,102) (sa(lp,l,jd),l=l1,l2)
C
  101 FORMAT(/,20x,'Slope matrix (Dot =',i2,') for',/,22x,
     1      'vector number ',i4,', IQ =',i4,/,29x,'RX0,RY0,RZ0',/,20x,
     2       3f10.6)
  102 FORMAT(9f8.4)
      RETURN
      END
