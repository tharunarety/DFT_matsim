      SUBROUTINE madmtr
C   ******************************************************************
C   *                                                                *
C   *    Read the Madelung matrix.                                   *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE message
      USE moments
      USE potential
      USE text
      IMPLICIT NONE
      REAL(KIND=8) :: dvv
      INTEGER      :: nqmad, nlmad, iq, jq, lm, lmp, nlm0
C
C     Read the Madelung matrix for Q'l'm';Qlm
C
      READ(4) txt
      READ(4) nqmad,nlmad
      nlmmad=nlmad*nlmad
      WRITE(m6,100) txt,nlmmad
C
      ALLOCATE(vmad(nq,nq,nlmmad,nlmmad))
C
      IF(nl2.GT.nlmad) THEN
         WRITE(m6,'(''NL2 > NLMAD'')')
         STOP
      ENDIF
      IF(nqmad.NE.nq) THEN
         WRITE(m6,'(''NQ in Madelung <> NQ in KGRN'')')
         STOP
      ENDIF
C
      DO 20 iq=1,nq
      DO 20 jq=1,nq
      READ(4) ((vmad(iq,jq,lm,lmp),lm=1,nlmmad),lmp=1,nlmmad)
   20 CONTINUE
C
C     Check for the Madelung matrix
C
      DO 30 iq=1,nq
      DO 30 jq=1,nq
      dvv=ABS(vmad(iq,jq,1,1)-vmdl(iq,jq))
      IF(dvv.GT.1.d-6) THEN
         WRITE(m6,110) iq,jq,vmad(iq,jq,1,1),vmdl(iq,jq)
         STOP
      ENDIF
   30 CONTINUE
C
      IF(nprn.EQ.2) THEN
         IF(func.EQ.'ASA') THEN
            WRITE(m6,120)
            DO 40 iq=1,nq
   41       WRITE(m6,122) vmad(iq,1:nq,1,1)
   40       CONTINUE
         ELSE
            nlm0=MIN0(diml,9)
            DO 42 iq=1,nq
            DO 42 jq=iq,nq
            WRITE(m6,121) iq,jq
            DO 43 lmp=1,nlm0
   43       WRITE(m6,122) vmad(iq,jq,lmp,1:lmp)
   42       CONTINUE
         ENDIF
      ENDIF
      IF(nt.EQ.1.OR.afm.EQ.'A') vmad(1:nq,1:nq,1,1)=0.d0
C
  100 FORMAT(/,' MADMTR:   ',a,//,11x,'NLMAD  = ',i3)
  110 FORMAT(/,' MADMTR:** Inconsistency in Madelung matrix:',//,
     1       11x,'IQ=',i2,' JQ=',i2,' VMAD=',f10.6,' VMDL=',f10.6)
  120 FORMAT(/,11x,'Madelung matrix',/)
  121 FORMAT(/,11x,'Madelung matrix for IQ = ',i2,' JQ = ',i2,/)
  122 FORMAT(10x,16f5.2)
C
      RETURN
      END
