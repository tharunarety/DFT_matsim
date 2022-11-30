      SUBROUTINE mltpm(lin)
C   ******************************************************************
C   *                                                                *
C   * Calculate multipole moments within the ASA spheres.            *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts       ; USE density      ; USE dosmom
      USE message    ; USE moments      ; USE radialmesh
      USE realgaunt  ; USE symmetry     ; USE temporary
      IMPLICIT NONE
      CHARACTER(LEN=5) :: tp(16)
      DATA TP/'  s  ','  y  ','  z  ','  x  ','  xy ','  yz ','3zz-1',
     .        '  xz ','xx-yy','3xx y',' xyz ','5zz y','5zz z','5zz x',
     .        'xx  z','3yy x'/
      REAL(KIND=8) :: totnos, smm, wq, qmom, gnti, facns, amagit
      INTEGER, PARAMETER :: fbzmom = 0
      INTEGER      :: lin, iq, it, is
      INTEGER      :: l, m, lm, mp, lmp, lpp, lmpp, ign
      INTEGER      :: irot, iqrot, nlm0
C
      IF(func.EQ.'ASA') RETURN
C
      qlmn=0.d0
      DO 20 is=1,ns
      DO 21 ign=1,ngaunt
      lpp=lppg(ign) 
      gnti=gnt(ign)/(2.d0*lpp+1.d0)
      lmp=lmpg(ign)
      lm=lmg(ign)
      lmpp=lmppg(ign)
   21 qlmn(1:nq,lmpp)=qlmn(1:nq,lmpp)+gnti*grnfm(1:nq,is,lmp,lm,lpp)
   20 CONTINUE
      facns=spinfc*sqfpi/pi
      qlmn=facns*qlmn
C
      IF(fllbz.NE.'Y'.AND.fbzmom.EQ.1) THEN
C
C        Calculate the total multipole moments from Q(IBZ) as
C
C               Q(FBZ) = U(T)*QIBZ 
C
C        where U(T) is the rotation matrix generated in ROTM3D.
C
         qlm0=0.d0
C
C        Loop for site
C
         DO 30 iq=1,nq
         wq=wqst(iq)
C
C        Loop for lm
C
         DO 31 l=0,lmax2
         DO 31 m=-l,l
         lm=l*l+l+m+1
C
C        Sum for symmetry operators
C
         DO 32 irot=1,nrot
         IF(ibzr(irot,iq).NE.0) THEN
C
C           If symmetry element IROT transform IQ to IQ' rearrange
C           correspondent partial moments
C
            iqrot=iprmt(irot,iq)
            smm=0.d0
C
C           Sum for l'm'
C
            DO 33 mp=-l,l
            lmp=l*l+l+mp+1
   33       smm=smm+ugam(irot,l,m,mp)*qlmn(iq,lmp)*wq
            qlm0(iqrot,lm)=qlm0(iqrot,lm)+smm
         ENDIF
   32    CONTINUE
   31    CONTINUE
   30    CONTINUE
C
         qlmn=qlm0
      ENDIF
C
C     Print non-zero multipole moments if requested
C
      IF(lin.EQ.1) THEN
         totnos=SUM(qlmn(1:nq,1))
         nlm0=MIN0(nlm2,16)
         qmom=SUM(ABS(qlmn(1:nq,1:nlm0)))
         IF(qmom.GT.1.d-6) THEN
            IF(ns.EQ.1) THEN
               WRITE(m6,101) tp(1:nlm0)
               WRITE(m6,'(A)') ' '
               DO 50 iq=1,nq
   50          WRITE(m6,102) iq,qlmn(iq,1:nlm0)
               WRITE(m6,103) totnos
            ELSE
               WRITE(m6,101) tp(1),'Spin ',tp(2:nlm0)
               WRITE(m6,'(A)') ' '
               DO 51 iq=1,nq
               it=itq(iq)
               amagit=SUM(conc(1:nta(it),it)*amag(1:nta(it),it))
   51          WRITE(m6,102) iq,qlmn(iq,1),amagit,qlmn(iq,2:nlm0)
               WRITE(m6,103) totnos
            ENDIF
         ENDIF
         qmom=SUM(ABS(qlmn(1:nq,nlm0+1:nlm2)))
         IF(qmom.GT.1.d-6) THEN
            DO 52 iq=1,nq
            WRITE(m6,104) iq
            DO 52 lm=nlm0+1,nlm2
            IF(ABS(qlmn(iq,lm)).GT.1.d-6) WRITE(m6,105) lm,qlmn(iq,lm)
   52       CONTINUE
         ENDIF
      ENDIF
C
  101 FORMAT(/,' MLTPM:    Multipole moments',//,'  IQ ',10(1x,a,1x))
  102 FORMAT(1x,i3,1x,18f7.4)
  103 FORMAT(' Tot',f8.4)
  104 FORMAT(/,' MLTPM:    Non vanishing higher moments for IQ = ',i2,/)
  105 FORMAT(11x,i3,f10.6)
      RETURN
      END
