      SUBROUTINE setsws(dsws)
C   ******************************************************************
C   *                                                                *
C   *    Set the average Wigner-Seitz radius or radii if the total   *
C   *    energy minimum is to be determined.                         *
C   *                                                                *
C   ******************************************************************
      USE message
      USE radialmesh
      USE volume
      IMPLICIT NONE
      INTEGER :: nswsh, isw0, i, isws
      REAL(KIND=8) :: dsws, sws0
C
      nswsh=nsws/2
      ALLOCATE(swsl(2*nswsh+1))
C
      IF(nsws.EQ.1) THEN
         swsl(1)=sws
      ELSE
         isw0=sws/dsws+0.05d0
         sws0=isw0*dsws
         isws=0
         DO 20 i=-nswsh,nswsh
         isws=isws+1
         swsl(isws)=sws0-i*dsws
   20    CONTINUE
         nsws=isws
         sws=swsl(1)
         WRITE(m6,101) swsl(1:nsws)
      ENDIF
      ALLOCATE(efl(nsws),etol(nsws))
C
  101 FORMAT(/,' SETSWS:   SWS',6f10.6)
      RETURN
      END
