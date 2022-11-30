      SUBROUTINE MADL3
C   ******************************************************************
C   *                                                                *
C   *    Generate the 3D Madelung matrices.                          *
C   *                                                                *
C   ******************************************************************
      USE madelung_lattice
      USE message
      IMPLICIT NONE
  101 FORMAT(' Calculate and store 3D Madelung matrices')
  102 FORMAT(/,' ***       Calculate and store 3D Madelung matrices')
  103 FORMAT(/,' ***       End of Madelung calculation')
C
      READ(5,'(3(10X,F10.6))') ALAMDA,AMAX,BMAX
      RMAX=AMAX/ALAMDA
      GMAX=2.D0*ALAMDA*BMAX
      IF(MSGL.EQ.1) WRITE(MSGIO,101)
      WRITE(M6,102)
      CALL LATT3M
      CALL MADL3D
      WRITE(M6,103)
      RETURN
      END
