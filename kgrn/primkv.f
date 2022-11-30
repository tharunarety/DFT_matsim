      SUBROUTINE PRIMKV(NQ)
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
      USE lattice
      USE message
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ONE=1)
      COMMON/KSPC/BKX(3),BKY(3),BKZ(3)
  101 FORMAT(//,' PRIMKR:',3X,'Wigner Seitz radius:',F10.6)
  102 FORMAT(/,11X,'Primitive translations of real space',/,
     1       11X,'in units of the lattice spacing a:',/)
  103 FORMAT(/,11X,'Basis vectors in the primitive cell',/,11X,
     1       'in units of the lattice spacing a:',/)
  104 FORMAT(/,11X,'Primitive 3D reciprocal translations',/,
     1       11X,'in units of 2*pi/a:',/)
  105 FORMAT(10X,'(',F10.5,',',F10.5,',',F10.5,'  )')
C
      PI=ACOS(-ONE)
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
      WS=(3.*VOL/4./PI/NQ)**(1./3.)
C
      WRITE(M6,101) WS
C
C     Print primitive real and reciprocal space vectors from <BSTR>.
C
      WRITE(M6,102)
      DO 120 I=1,3
  120 WRITE(M6,105) BSX(I),BSY(I),BSZ(I)
      WRITE(M6,103)
      DO 121 IQ=1,NQ
  121 WRITE(M6,105) QX(IQ),QY(IQ),QZ(IQ)
C
      WRITE(M6,104)
      DO 23 I=1,3
      BKX(I)=BKX(I)/VOL
      BKY(I)=BKY(I)/VOL
      BKZ(I)=BKZ(I)/VOL
   23 WRITE(M6,105) BKX(I),BKY(I),BKZ(I)
      TWOPI=2.*PI
      DO 24 I=1,3
      BKX(I)=BKX(I)*TWOPI
      BKY(I)=BKY(I)*TWOPI
   24 BKZ(I)=BKZ(I)*TWOPI
C
      RETURN
      END
