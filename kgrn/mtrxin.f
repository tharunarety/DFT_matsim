      SUBROUTINE mtrxinv(n,a,b,c,info)
C
      IMPLICIT NONE
      INTEGER :: n, lda, ldb, nrhs, info
      COMPLEX(KIND=8), DIMENSION(n,n) :: a, b, c
      INTEGER, DIMENSION(n) :: ipiv
C
      lda=n
      ldb=n
      nrhs=n
C
      c=b
      CALL zgesv( N, NRHS, A, LDA, IPIV, C, LDB, INFO )
      IF(info.NE.0) STOP 'MTRXINV: INFO <> 0'
C
      RETURN
      END
