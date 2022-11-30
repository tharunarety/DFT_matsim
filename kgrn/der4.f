      SUBROUTINE der4(n,x,y,del)
C   ******************************************************************
C   *                                                                *
C   * Find the derivatives (4-order method)                          *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      INTEGER      :: n, n2, i
      REAL(KIND=8), DIMENSION(n) :: x, y
      REAL(KIND=8) :: del
C
      n2=n-2
      DO 20 i=3,n2
   20 y(i)=((x(i-2)-x(i+2))-8.d0*(x(i-1)-x(i+1)))/(12.d0*del)
      y(2)=(x(5)-6.d0*x(4)+18.d0*x(3)-1.d1*x(2)-3.d0*x(1))/(12.d0*del)
      y(n-1)=(3.d0*x(n)+1.d1*x(n-1)-18.d0*x(n-2)+6.d0*x(n-3)-x(n-4))/
     *(12.d0*del)
      y(1)=(-3.d0*x(1)+4.d0*x(2)-x(3))/(2.d0*del)
      y(n)=(3.d0*x(n)-4.d0*x(n-1)+x(n-2))/(2.d0*del)
C
      RETURN
      END
