      SUBROUTINE eldens
C   ******************************************************************
C   *                                                                *
C   * Calculate the full charge density n(r) = sum_L n_L(r) * Y_L(r).*
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text
      USE csts       ; USE density      ; USE dosmom
      USE message    ; USE moments      ; USE pota
      USE radialmesh ; USE softcore     ; USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: qlmnt
      REAL(KIND=8) :: facns, totnos, qmom, sqfl, swsl, r, fint, smm
      INTEGER :: lm, nlm0, l, m, lp
      INTEGER :: iq, it, ita, is, jri, jrc, jws, ir
      CHARACTER*5 :: tp(16)
      DATA TP/'  s  ','  y  ','  z  ','  x  ','  xy ','  yz ','3zz-1',
     .        '  xz ','xx-yy','3xx y',' xyz ','5zz y','5zz z','5zz x',
     .        'xx  z','3yy x'/
C
      IF(msgl.NE.0) WRITE(msgio,100)
      ALLOCATE(qlmn(nq,nlmf),qlmnt(mnta,nq,nlmf))
C
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
   20 qti(ita,it)=SUM(tnos(ita,it,0:lmax,1:ns))
C
C     Check for the charge neutrality
C
      smm=0.d0
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      qti(ita,it)=qti(ita,it)-eln(ita,it)
      smm=smm+conc(ita,it)*mmt(it)*qti(ita,it)
   21 CONTINUE
      IF(ABS(smm).GT.1.d-6) THEN
         WRITE(m6,107) smm,qti
         IF(msgl.NE.0) WRITE(msgio,107) smm,qti
      ENDIF
C
      facns=spinfc/pi
      chdl=facns*chdl
C
      IF(fllbz.NE.'Y') CALL rotchd
C
C     Add the soft core correction to the charge-density
C
      facns=spinfc/2.d0
      IF(softc.EQ.'Y') THEN
         DO 22 is=1,ns
         DO 22 iq=1,nq
         it=itq(iq)
         DO 22 ita=1,nta(it)
         jws=jwss(ita,it)
   22    chdl(ita,iq,is,1,1:jws)=chdl(ita,iq,is,1,1:jws)+
     .                       facns*cor(ita,it,1:jws)/sqfpi
      ENDIF
C
C     Set up the average multipole moments
C
      qlmn=0.d0
      qlmnt=0.d0
      DO 30 iq=1,nq
      it=itq(iq)
      DO 31 ita=1,nta(it)
C
      jws=jwss(ita,it)
      DO 32 l=0,lmaxf
      swsl=sws**l
      sqfl=sqfpi/(2.d0*l+1.d0)
      DO 32 m=-l,l
      lp=l+1
      lm=l*l+l+m+1
C
      DO 33 ir=1,jws
      r=ri(ir,ita,it)
   33 fi(ir)=SUM(chdl(ita,iq,1:ns,lm,ir))*(r**lp)
      CALL simpn(fi,dx,jws,fint)
      qlmnt(ita,iq,lm)=qlmnt(ita,iq,lm)+fint*sqfl/swsl
   32 qlmn(iq,lm)=qlmn(iq,lm)+conc(ita,it)*fint*sqfl/swsl
   31 CONTINUE
   30 CONTINUE
C
      totnos=SUM(qlmn(1:nq,1))
C
C     Print non-zero multipole moments
C
      nlm0=MIN0(nlmf,16)
      qmom=SUM(ABS(qlmn(1:nq,1:nlm0)))
      IF(qmom.GT.1.d-6) THEN
         IF(ns.EQ.1) THEN
            WRITE(m6,101) tp(1:nlm0)
            WRITE(m6,'(A)') ' '
            DO 40 iq=1,nq
            it=itq(iq)
            DO 40 ita=1,nta(it)
   40       WRITE(m6,102) iq,qlmnt(ita,iq,1:nlm0)
            WRITE(m6,104) totnos
         ELSE
            WRITE(m6,101) tp(1),'Spin ',tp(2:nlm0)
            WRITE(m6,'(A)') ' '
            DO 41 iq=1,nq
            it=itq(iq)
            DO 41 ita=1,nta(it)
   41       WRITE(m6,103) iq,qlmnt(ita,iq,1),amag(ita,it),
     .                    qlmnt(ita,iq,2:nlm0)
            WRITE(m6,104) totnos
         ENDIF
      ENDIF
      qmom=SUM(ABS(qlmn(1:nq,nlm0+1:nlmf)))
      IF(qmom.GT.1.d-6) THEN
         DO 42 iq=1,nq
         it=itq(iq)
         DO 42 ita=1,nta(it)
         WRITE(m6,105) iq,ita
         DO 42 lm=nlm0+1,nlmf
         IF(ABS(qlmnt(ita,iq,lm)).GT.1.d-6) 
     .      WRITE(m6,106) lm,qlmnt(ita,iq,lm)
   42    CONTINUE
      ENDIF
C
C     Zero order moments
C
      qlmn(1:nq,1)=0.d0
      DO 50 iq=1,nq
      it=itq(iq)
      DO 50 ita=1,nta(it)
   50 qlmn(iq,1)=qlmn(iq,1)+conc(ita,it)*qtr(ita,it)
C
      DEALLOCATE(qlmnt)
C
  100 FORMAT(/,' ELDENS: Multipole moments')
  101 FORMAT(/,' ELDENS:   Multipole moments',//,'  IQ ',10(1x,a,1x))
  102 FORMAT(1x,i3,1x,16f7.4)
  103 FORMAT(1x,i3,1x,17f7.4)
  104 FORMAT(' Tot',f8.4)
  105 FORMAT(/,' ELDENS:   Non vanishing higher moments for IQ = ',i2,
     .         ' ITA =',i2,/)
  106 FORMAT(11x,i3,f10.6)
  107 FORMAT(/,' ELDENS: delQ =',f10.6,' QTI =',10f10.6)
  110 FORMAT(/,' ELDENS:  IQ = ',i3,' IS = ',i3,
     .      //,15x,'r',15x,'n(r)',/)
  111 FORMAT(11x,f10.6,3f16.6)
      RETURN
      END
