      SUBROUTINE prinfo(pro)
C   ******************************************************************
C   *                                                                *
C   * Print site-,l- and spin projected NOS and EONE.                *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE dosmom ; USE greenfunc ; USE message ; USE bzmesh
      IMPLICIT NONE
      CHARACTER(LEN=6) :: pro
      REAL(KIND=8), DIMENSION(0:lmax) :: tnost0, emomt0
      REAL(KIND=8), DIMENSION(0:lmax) :: tnost, emomt
      REAL(KIND=8) :: stnos
      INTEGER      :: is, it, ita, l
C
      IF(fsts.EQ.1) GO TO 30
C
C     NOS and EMOM printout
C
      stnos=0.d0
      DO 20 is=1,ns
C
      tnost=0.d0
      emomt=0.d0
      tnost0=0.d0
      emomt0=0.d0
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
C
C     Print to MSGIO
C
      IF(msgl.EQ.1) THEN
         IF(zmsh.EQ.'m'.OR.zmsh.EQ.'M'.OR.zmsh.EQ.'f') THEN
            WRITE(msgio,102) pro,tnos0(ita,it,0:lmax,is),
     .                       SUM(tnos0(ita,it,0:lmax,is))
            WRITE(msgio,112) emom0(ita,it,0:lmax,is),
     .                       SUM(emom0(ita,it,0:lmax,is))
         ENDIF
         WRITE(msgio,100) pro,tnos(ita,it,0:lmax,is),
     .                    SUM(tnos(ita,it,0:lmax,is))
         WRITE(msgio,110) emom(ita,it,0:lmax,is),
     .                    SUM(emom(ita,it,0:lmax,is))
      ENDIF
C
C     Print to M6
C
      IF(zmsh.EQ.'m'.OR.zmsh.EQ.'M'.OR.zmsh.EQ.'f') THEN
         WRITE(m6,102) pro,tnos0(ita,it,0:lmax,is),
     .                 SUM(tnos0(ita,it,0:lmax,is))
         WRITE(m6,112) emom0(ita,it,0:lmax,is),
     .                 SUM(emom0(ita,it,0:lmax,is))
         DO 22 l=0,lmax
         tnost0(l)=tnost0(l)+conc(ita,it)*mmt(it)*tnos0(ita,it,l,is)
         emomt0(l)=emomt0(l)+conc(ita,it)*mmt(it)*emom0(ita,it,l,is)
   22    CONTINUE
      ENDIF
      WRITE(m6,100) pro,tnos(ita,it,0:lmax,is),
     .              SUM(tnos(ita,it,0:lmax,is))
      WRITE(m6,110) emom(ita,it,0:lmax,is),
     .              SUM(emom(ita,it,0:lmax,is))
      DO 23 l=0,lmax
      tnost(l)=tnost(l)+conc(ita,it)*mmt(it)*tnos(ita,it,l,is)
      emomt(l)=emomt(l)+conc(ita,it)*mmt(it)*emom(ita,it,l,is)
   23 CONTINUE
   21 CONTINUE
C
C     Print total TNOS and EMOM
C
      IF(msgl.EQ.1) THEN
         IF(zmsh.EQ.'m'.OR.zmsh.EQ.'M'.OR.zmsh.EQ.'f') THEN
            WRITE(msgio,103) tnost0(0:lmax),SUM(tnost0(0:lmax))
            WRITE(msgio,113) emomt0(0:lmax),SUM(emomt0(0:lmax))
         ENDIF
         WRITE(msgio,101) tnost(0:lmax),SUM(tnost(0:lmax))
         WRITE(msgio,111) emomt(0:lmax),SUM(emomt(0:lmax))
      ENDIF
      IF(zmsh.EQ.'m'.OR.zmsh.EQ.'M'.OR.zmsh.EQ.'f') THEN
         WRITE(m6,103) tnost0(0:lmax),SUM(tnost0(0:lmax))
         WRITE(m6,113) emomt0(0:lmax),SUM(emomt0(0:lmax))
      ENDIF
      WRITE(m6,101) tnost(0:lmax),SUM(tnost(0:lmax))
      WRITE(m6,111) emomt(0:lmax),SUM(emomt(0:lmax))
      stnos=stnos+SUM(tnost(0:lmax))
   20 CONTINUE
C
      IF(msgl.EQ.1) WRITE(msgio,120) stnos,elt
      WRITE(m6,120) stnos,elt
C
      RETURN
C
C     DOS printout
C
   30 CONTINUE
      stnos=0.d0
      DO 40 is=1,ns
C
      tnost=0.d0
      DO 41 it=1,nt
      DO 41 ita=1,nta(it)
C
C     Print to MSGIO
C
      IF(msgl.EQ.1) THEN
         WRITE(msgio,130) pro,tdos(ita,it,0:lmax,is),
     .                    SUM(tdos(ita,it,0:lmax,is))
      ENDIF
C
C     Print to M6
C
      WRITE(m6,130) pro,tdos(ita,it,0:lmax,is),
     .              SUM(tdos(ita,it,0:lmax,is))
      DO 42 l=0,lmax
      tnost(l)=tnost(l)+conc(ita,it)*mmt(it)*tdos(ita,it,l,is)
   42 CONTINUE
   41 CONTINUE
C
C     Print total DOS
C
      IF(msgl.EQ.1) THEN
         WRITE(msgio,131) tnost(0:lmax),SUM(tnost(0:lmax))
      ENDIF
      WRITE(m6,131) tnost(0:lmax),SUM(tnost(0:lmax))
      stnos=stnos+SUM(tnost(0:lmax))
   40 CONTINUE
C
      IF(msgl.EQ.1) WRITE(msgio,140) stnos
      WRITE(m6,140) stnos
C
      RETURN
C
  100 FORMAT(/,' ',a,': Tnos (s,p,...,tot) : ',7f8.4)
  101 FORMAT(/,'         TNOS (s,p,...,tot) : ',7f8.4)
  102 FORMAT(/,' ',a,': Tnos (lower panel) : ',7f8.4)
  103 FORMAT(/,'         TNOS (lower panel) : ',7f8.4)
  110 FORMAT(/,'         Eone (s,p,...,tot) : ',7f8.4)
  111 FORMAT(/,'         EONE (s,p,...,tot) : ',7f8.4)
  112 FORMAT(/,'         Eone (lower panel) : ',7f8.4)
  113 FORMAT(/,'         EONE (lower panel) : ',7f8.4)
  120 FORMAT(/,'         NOS(Ef) =',f10.6,' ELT =',f8.4)
  130 FORMAT(/,' ',a,': dos (s,p,...,tot) : ',7f8.4)
  131 FORMAT(/,'         DOS (s,p,...,tot) : ',7f8.4)
  140 FORMAT(/,'         DOS(Ef) =',f10.6)
      END
