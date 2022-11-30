      SUBROUTINE softp(mix,ncore)
C   ******************************************************************
C   *                                                                *
C   *    Prepare and perform atomic calculations for the soft core.  *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE control_data
      USE density
      USE pota
      USE softcore
      USE radialmesh
      USE temporary
      IMPLICIT NONE
      INTEGER, PARAMETER :: mw = 999
      REAL(KIND=8), DIMENSION(mw) :: corn
      REAL(KIND=8) :: mix, mix1
      INTEGER :: ncore, icore, it, ita, is, ir, jws, jri
C
      mix1=1.d0-mix
      ncore=1
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jws=jwss(ita,it)
      jri=jws+2
C
C     Set up the new core density
C
      eonec(ita,it)=0.d0
      IF(nz(ita,it).EQ.0) GO TO 20
      DO 21 is=1,ns
      IF(is.EQ.1) corn(1:mw)=0.d0
C
      CALL softi(ws(ita,it),r1(ita,it),jws,nz(ita,it),it,ita,is)
      CALL softc(ws(ita,it),corn,icore,it,ita,is,ns)
      ncore=icore*ncore
   21 CONTINUE
      cor(ita,it,1:jri)=mix1*cor(ita,it,1:jri)+mix*corn(1:jri)
C
   20 CONTINUE
C
      RETURN
      END
