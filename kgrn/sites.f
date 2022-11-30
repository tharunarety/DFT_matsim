      SUBROUTINE sites
C   ******************************************************************
C   *                                                                *
C   *    Check if the setup in the data file is consistent with      *
C   *    the setup in the structure progrem, KSTR.                   *
C   *    For different local enviroments IT should be different.     *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE message
      USE radialmesh
      IMPLICIT NONE
      INTEGER :: iq, it, jq, jt
      REAL(KIND=8) :: tol = 1.d-6
C
      DO 20 iq=1,nq
      it=itq(iq) 
      DO 20 jq=1,nq
      jt=itq(jq) 
      IF(it.EQ.jt) THEN
         IF(ABS(wst(iq)-wst(jq)).GT.tol) THEN
            WRITE(m6,100) iq,it,wst(iq),jq,jt,wst(jq)
            STOP
         ENDIF
         IF(ABS(wsi(iq)-wsi(jq)).GT.tol) THEN
            WRITE(m6,101) iq,it,wsi(iq),jq,jt,wsi(jq)
            STOP
         ENDIF
         IF(ABS(wsc(iq)-wsc(jq)).GT.tol) THEN
            WRITE(m6,102) iq,it,wsc(iq),jq,jt,wsc(jq)
            STOP
         ENDIF
      ENDIF
   20 CONTINUE
C
  100 FORMAT(/,' SITES:*** IQ =',i3,' IT=',i3,' WST=',f10.6,
     .                   ' JQ =',i3,' JT=',i3,' WST=',f10.6)
  101 FORMAT(/,' SITES:*** IQ =',i3,' IT=',i3,' WSI=',f10.6,
     .                   ' JQ =',i3,' JT=',i3,' WSI=',f10.6)
  102 FORMAT(/,' SITES:*** IQ =',i3,' IT=',i3,' WSC=',f10.6,
     .                   ' JQ =',i3,' JT=',i3,' WSC=',f10.6)
      RETURN
      END
