C     Last change:  HLS  20 Jan 98    8:33 am
      SUBROUTINE htimer(dato,clock)
C   ******************************************************************
C   *                                                                *
C   *    Returns dd-mm-yy in dato*9 and hr:min in clock*5            *
C   *                                                                *
C   ******************************************************************
      CHARACTER*9                :: dato
      CHARACTER*5                :: clock
      CHARACTER*3, DIMENSION(12) :: month
      CHARACTER*10               :: dd, tt
      INTEGER                    :: dt(8)
      DATA month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     1           'Oct','Nov','Dec'/
C
      CALL date_and_time(DATE=dd, TIME=tt, VALUES=dt)
C
      dato=dd(7:8)//'-'//month(dt(2))//'-'//dd(3:4)
      clock=tt(1:2)//':'//tt(3:4)
      RETURN
      END
