      SUBROUTINE volumk(x,y,z,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,v)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the volume and centre of gravity of a tetrahedron *
C   *    spanned by vectors:                                         *
C   *                                                                *
C   *        (X1,Y1,Z1);(X2,Y2,Z2);(X3,Y3,Z3);(X4,Y4,Z4).            *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      REAL(KIND=8) :: x, y, z, v, p1, p2, p3
      REAL(KIND=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
C
      x=0.25d0*(x1+x2+x3+x4)
      y=0.25d0*(y1+y2+y3+y4)
      z=0.25d0*(z1+z2+z3+z4)
C
      CALL cross(x2-x1,y2-y1,z2-z1,x3-x1,y3-y1,z3-z1,p1,p2,p3)
      v=DABS((x4-x1)*p1+(y4-y1)*p2+(z4-z1)*p3)/6.d0
C
      RETURN
      END
