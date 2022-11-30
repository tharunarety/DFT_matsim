      SUBROUTINE primkr
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
      USE control_data
      USE csts
      USE message
      IMPLICIT NONE
      INTEGER :: i, i1, i2
C
C     Find primitive vectors (BKX,BKY,BKZ) of reciprocal space
C
      DO 22 i=1,3
      i1=1+MOD(i,3)
      i2=1+MOD(i1,3)
      CALL cross(bsx(i1),bsy(i1),bsz(i1),bsx(i2),bsy(i2),bsz(i2),
     .           bkx(i),bky(i),bkz(i))
   22 CONTINUE
      vol=ABS(bsx(1)*bkx(1)+bsy(1)*bky(1)+bsz(1)*bkz(1))
C
C     Atomic Wigner-Seitz radius in units of a
C
      ws=(3.d0*vol/4.d0/pi/NQ3)**(1.d0/3.d0)
      WRITE(m6,103) ws,vol
C
      WRITE(m6,102)
      DO 23 i=1,3
      bkx(i)=bkx(i)/vol
      bky(i)=bky(i)/vol
      bkz(i)=bkz(i)/vol
   23 WRITE(m6,101) bkx(i),bky(i),bkz(i)
      DO 24 i=1,3
      bkx(i)=bkx(i)*twopi
      bky(i)=bky(i)*twopi
   24 bkz(i)=bkz(i)*twopi
C
  101 FORMAT(' ',10x,'(',f10.5,',',f10.5,',',f10.5,'  )')
  102 FORMAT(//,11x,'Primitive vectors of reciprocal space',
     .       /,11x,'in units of 2*pi/a',/)
  103 FORMAT(//,' PRIMKR:',3x,'Wigner Seitz radius:',f10.6,
     .       ' VOL   =',f10.6)
      RETURN
      END
