      SUBROUTINE setdrc(dx)
C   ******************************************************************
C   *                                                                *
C   *     Set constants for DIRAC.                                   *
C   *                                                                *
C   ******************************************************************
      USE diracparam
      IMPLICIT NONE
      REAL(KIND=8) :: dx
C
      TEST=1.D09
      CLIGHT=274.07446D0
      CSQ=CLIGHT*CLIGHT
      EXPDXH=EXP(DX/2.D0)
      DXD8=DX/8.0D0
      A1=DX*3.3D0
      A2=-DX*4.2D0
      A3=DX*7.8D0
      A4=DX*14.D0/45.D0
      A5=DX*64.D0/45.D0
      A6=DX*24.D0/45.D0
      FAC3=7.D0/8.D0*EXPDXH**2
      FAC4=.5D0*EXPDXH**2
      RETURN
      END
