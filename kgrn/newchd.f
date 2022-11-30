      SUBROUTINE newchd(prnt)
C   ******************************************************************
C   *                                                                *
C   * Calculate the spherical part of the new charge density.        *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text ; USE csts
      USE density    ; USE dosmom       ; USE greenfunc
      USE message    ; USE pota         ; USE radialmesh
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: facns, r, smm, addz, zwc
      INTEGER :: prnt, iq, it, ita, is, jrn, jws, ir
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
         WRITE(m6,100) 'QTI(dos)=',qti
         IF(msgl.NE.0) WRITE(msgio,100) 'QTI(dos)=',qti
         STOP
      ENDIF
C
C     Devide the density with number of atoms per type (MMT)
C     and with pi (from the contour integral)
C     Multiply by 2 if NS = 1
C     
      DO 22 is=1,ns
      DO 22 it=1,nt
      facns=spinfc/mmt(it)/pi
      DO 22 ita=1,nta(it)
      jrn=jrsm(ita,it)
   22 chdn(ita,it,is,1:jrn)=chdn(ita,it,is,1:jrn)*facns
C
C     Add core charge-density
C
      facns=spinfc/2.d0
      DO 23 is=1,ns
      DO 23 it=1,nt
      DO 23 ita=1,nta(it)
      jws=jwss(ita,it)
   23 chdn(ita,it,is,1:jws)=chdn(ita,it,is,1:jws)+
     .                      facns*cor(ita,it,1:jws)
C
      IF(prnt.EQ.1) THEN
         DO 30 it=1,nt
         DO 30 ita=1,nta(it)
         jrn=jrsm(ita,it)
         DO 30 is=1,ns
         WRITE(m6,110) it,ita,is
         fi(1:jrn)=chdn(ita,it,is,1:jrn)/fourpi
         CALL diffn(fi,rhop,rhopp,jrn,dx)
         rhopp(1:jrn)=(rhopp(1:jrn)-rhop(1:jrn))/ri2(1:jrn,ita,it)
         rhop(1:jrn)=rhop(1:jrn)/ri(1:jrn,ita,it)
         DO 30 ir=1,jrn
         r=ri(ir,ita,it)
   30    WRITE(m6,120) r,fi(ir),fi(ir)-facns*cor(ita,it,ir)/fourpi,
     .                 facns*cor(ita,it,ir)/fourpi,
     .                 facns*cor(ita,it,ir)/fourpi/fi(ir)
      ENDIF
C
  100 FORMAT(/,' NEWCHD:  ',a,10f10.6)
  110 FORMAT(/,' NEWCHD:  IT = ',i3,' ITA =',i3,' IS = ',i3,
     .      //,15x,'r',15x,'n(r)',/)
  120 FORMAT(11x,f10.6,8f16.6)
C
      RETURN
      END
