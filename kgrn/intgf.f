      SUBROUTINE intgf(ga,bg)
C   ******************************************************************
C   *                                                                *
C   *   Calculate the path operator and density of states for the    *
C   *   full Brillouin zone.                                         *
C   *                                                                *
C   ******************************************************************
      USE control_data
      IMPLICIT NONE
      COMPLEX(KIND=8), DIMENSION(nq,nlm,nlm) :: ga
      COMPLEX(KIND=8), DIMENSION(nq,nlm,nlm) :: bg
C
C     Rotate and sum up the matrix elements for the star of k
C
      CALL rotgf(ga(1:nq,1:nlm,1:nlm))
      CALL rotgf(bg(1:nq,1:nlm,1:nlm))
C
      RETURN
      END
