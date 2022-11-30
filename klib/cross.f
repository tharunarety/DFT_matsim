      SUBROUTINE cross(ax,ay,az,bx,by,bz,cx,cy,cz)
C   ******************************************************************
C   *                                                                *
C   *   Cross product (CX,CY,CZ)=(AX,AY,AZ)x(BX,BY,BZ)               *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      REAL(KIND=8) :: ax, ay, az, bx, by, bz, cx, cy, cz
C
      cx=ay*bz-by*az
      cy=bx*az-ax*bz
      cz=ax*by-bx*ay
C
      RETURN
      END

