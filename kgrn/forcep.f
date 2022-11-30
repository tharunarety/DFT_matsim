      SUBROUTINE forcep
C   ******************************************************************
C   *                                                                *
C   *    Print the force components.                                 *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE force
      USE message    ; USE radialmesh   ; USE text
      IMPLICIT NONE
      CHARACTER(LEN=4) :: a
      REAL(KIND=8) :: ftotx, ftoty, ftotz, ftotxc, ftotyc, ftotzc
      REAL(KIND=8) :: tol = 1.d-5, si
      INTEGER :: iq, it, ita, jsi, loop
C
      loop=1
   10 CONTINUE
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      jsi=jwsi(ita,it)-shf
      si=ri(jsi,ita,it)
      a=ttxt(ita,it)
      IF(loop.EQ.1) THEN
         IF(iq.EQ.1) WRITE(m6,100)
         IF(nz(ita,it).EQ.0) GO TO 20
      ELSE
         IF(iq.EQ.1) WRITE(m6,101)
         IF(nz(ita,it).NE.0) GO TO 20
      ENDIF
C
      IF(ABS(fix(iq)+fiy(iq)+fiz(iq)).GT.tol) THEN
         WRITE(m6,110) iq,fix(iq),fiy(iq),fiz(iq)
      ENDIF
C
      ftotx=fex(iq)+fix(iq)+fxx(iq)+fsx(iq)+fnx(iq)
      ftotxc=ftotx+fqx(iq)+fcx(iq)
      IF(DABS(ftotx).GT.tol) WRITE(m6,120) a,si,' x',fex(iq),
     .   fxx(iq),fsx(iq)+fnx(iq),ftotx,fqx(iq),fcx(iq),ftotxc
      ftoty=fey(iq)+fiy(iq)+fxy(iq)+fsy(iq)+fny(iq)
      ftotyc=ftoty+fqy(iq)+fcy(iq)
      IF(DABS(ftoty).GT.tol) WRITE(m6,120) a,si,' y',fey(iq),
     .   fxy(iq),fsy(iq)+fny(iq),ftoty,fqy(iq),fcy(iq),ftotyc
      ftotz=fez(iq)+fiz(iq)+fxz(iq)+fsz(iq)+fnz(iq)
      ftotzc=ftotz+fqz(iq)+fcz(iq)
      IF(DABS(ftotz).GT.tol) WRITE(m6,120) a,si,' z',fez(iq),
     .   fxz(iq),fsz(iq)+fnz(iq),ftotz,fqz(iq),fcz(iq),ftotzc
   20 CONTINUE
      IF(loop.EQ.1) THEN
         loop=2
         GO TO 10
      ENDIF
C
  100 FORMAT(/,' FORCEP: forces in Ry/a.u. :',//,
     . 3x,'Atom',4x,'Sf',11x,'ele',7x,'Fxc',7x,'kin',7x,'tot',
     . 7x,'srf',7x,'vol',7x,'tot'/)
  101 FORMAT(/,' FORCEP: forces in Ry/a.u. :',//,
     . 3x,'Empty',3x,'Sf',11x,'ele',7x,'Fxc',7x,'kin',7x,'tot',
     . 7x,'srf',7x,'vol',7x,'tot'/)
  110 FORMAT(/,' non zero self interaction term IQ=',i3,3f10.6)
  120 FORMAT(3x,a,f8.4,2x,a,1x,8f10.6)
      RETURN
      END
