      SUBROUTINE gensim(f,d,n,fi)
C     ***************************************************************
C     *                                                             *
C     *    Integrates F from X(1) to X(n) where X is Louck's mesh.  *
C     *                                                             *
C     ***************************************************************
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(*) :: f
      REAL(KIND=8) :: d, fi, ff, corr
      INTEGER :: n, n1, j, n2
C
      IF(n.GE.3) THEN
         IF(MOD(n,2).EQ.0) THEN
           n1=n-1
           ff=-f(n-2)+8.d0*f(n-1)+5.d0*f(n)
           corr=ff*d/12.d0
         ELSE
           n1=n
           corr=0.d0
         ENDIF
         n2=n1-1
         ff=0.d0
         DO 20 j=2,n2,2
   20    ff=ff+(f(j-1)+4.d0*f(j)+f(j+1))
         fi=d*ff/3.d0+corr
      ELSEIF(n.EQ.2) THEN
         ff=5.d0*f(1)+8.d0*f(2)-f(3)
         fi=d*ff/12.d0
      ELSEIF(n.LE.1) THEN
         fi=0.d0
      ENDIF
C
      RETURN
      END
