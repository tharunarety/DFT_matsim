      SUBROUTINE PRIMKR
C   ******************************************************************
C   *                                                                *
C   *    Generate primitive vectors of reciprocal space, real space  *
C   *    cell volume, and atomic Wigner Seitz radius.                *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    BKX                                                         *
C   *    BKY   : Primitive translations in reciprocal space in       *
C   *    BKZ     units of one over the lattice spacing.              *
C   *    VOL   : Cell volume in units of the lattice spacing cubed.  *
C   *    WS    : Atomic Wigner Seitz radius in units of the lattice  *
C   *            spacing.                                            *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE csts
      USE control_data
      USE message
      IMPLICIT NONE
      INTEGER :: I, I1, I2
  101 FORMAT(' ',10X,'(',F10.5,',',F10.5,',',F10.5,'  )')
  102 FORMAT(//,11X,'Primitive vectors of reciprocal space',
     1       /,11X,'in units of 2*pi/a',/)
  103 FORMAT(//,' PRIMKR:',3X,'Wigner Seitz radius:',F10.6,
     1       ' VOL   =',F10.6)
C
C     Find primitive vectors (BKX,BKY,BKZ) of reciprocal space
C
      DO 22 I=1,3
      I1=1+MOD(I,3)
      I2=1+MOD(I1,3)
      CALL CROSS(BSX(I1),BSY(I1),BSZ(I1),BSX(I2),BSY(I2),BSZ(I2),
     1           BKX(I),BKY(I),BKZ(I))
   22 CONTINUE
      VOL=ABS(BSX(1)*BKX(1)+BSY(1)*BKY(1)+BSZ(1)*BKZ(1))
C
C     Atomic Wigner-Seitz radius in units of a
C
      SWS=(3.D0*VOL/4.D0/PI/NQ3)**(1.D0/3.D0)
      WRITE(M6,103) SWS,VOL
C
      WRITE(M6,102)
      DO 23 I=1,3
      BKX(I)=BKX(I)/VOL
      BKY(I)=BKY(I)/VOL
      BKZ(I)=BKZ(I)/VOL
   23 WRITE(M6,101) BKX(I),BKY(I),BKZ(I)
      DO 24 I=1,3
      BKX(I)=BKX(I)*TWOPI
      BKY(I)=BKY(I)*TWOPI
   24 BKZ(I)=BKZ(I)*TWOPI
      RETURN
      END
