      SUBROUTINE wmsg(ef,etotal,erren,erref,lin,isws)
C   ******************************************************************
C   *                                                                *
C   *    Write messages to screen and output file                    *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE control_text ; USE density
      USE dosmom ; USE message ; USE moments
      IMPLICIT NONE
      REAL(KIND=8) :: ef, etotal, erren, erref, fcpa
      INTEGER      :: lin, isws, it
      erren=MIN(1.d0,erren)
      erref=MIN(1.d0,erref)
C
      fcpa=SUM(1.d0*nta(1:nt))/nt
      IF(msgl.EQ.1) THEN
         WRITE(msgio,110) iter,etotal,erren
         WRITE(msgio,120) lin,ef,erref
         IF(afm.NE.'P') WRITE(msgio,130) tmag
         IF(nt.GT.1) THEN
            WRITE(msgio,140) (qtr(1:nta(it),it),it=1,nt)
         ENDIF
         WRITE(msgio,150) qsca(1:nt)
         IF(fcpa.GT.1.d0) THEN
            WRITE(msgio,160) (qcpa(1:nta(it),it),it=1,nt)
         ENDIF
      ENDIF
      WRITE(m6,110) iter,etotal,erren
      IF(afm.NE.'P') WRITE(m6,130) tmag
      WRITE(m6,120) lin,ef,erref
      IF(nt.GT.1) THEN
         WRITE(m6,140) (qtr(1:nta(it),it),it=1,nt)
      ENDIF
      WRITE(m6,150) qsca(1:nt)
      IF(fcpa.GT.1.d0) THEN
         WRITE(m6,160) (qcpa(1:nta(it),it),it=1,nt)
      ENDIF
C
  110 FORMAT(/,' KGRN:  Iteration no. ',i3,' Etot =',f16.6,
     .         ' erren =',f11.8)
  120 FORMAT(/,8x,'Dyson loops   ',i3,' EF   =',f16.6,
     .         ' erref =',f11.8)
  130 FORMAT(/,8x,'Magn. mom. =',f7.3)
  140 FORMAT(/,' KGRN:  QTR  =',6f10.6)
  150 FORMAT(/,' KGRN:  QSCA =',6f10.6)
  160 FORMAT(/,' KGRN:  QCPA =',6f10.6)
      RETURN
      END
