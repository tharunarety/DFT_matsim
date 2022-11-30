      SUBROUTINE atoma
C   ******************************************************************
C   *                                                                *
C   *   Open scratch file for atomic calculation.                    *
C   *                                                                *
C   ******************************************************************
      USE control_data
      IMPLICIT NONE
      CHARACTER(LEN=120) :: string
      INTEGER :: n, i, it, ita
C
      OPEN(8,STATUS='SCRATCH')
      n=5
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
   20 n=n+6
      DO 21 i=1,n
      READ(5,'(a)') string
      WRITE(8,'(a)') string
   21 CONTINUE
      RETURN
      END
