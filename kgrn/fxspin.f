      SUBROUTINE fxspin(ita,it,lin)
C   ******************************************************************
C   *                                                                *
C   * Get the exchange splitting that generates the input moment.    *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE dosmom ; USE potential ; USE message
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: a0 = 0.01d0, a1 = 0.5d0, a2 = 0.2d0
      REAL(KIND=8) :: ek0, ek1, ek2, ammom
      INTEGER :: ita,it,lin
C
      ek0=mmom-tmag
      ek1=mmom-tmago
      ek2=mmom-tmagoo
      IF(ABS(mmom).GT.2.d0) THEN
         ammom=ABS(mmom)
         ek0=ek0/ammom
         ek1=ek1/ammom
         ek2=ek2/ammom
      ENDIF
      ek0=ek0/nq
      ek1=ek1/nq
      ek2=ek2/nq
      dexch(ita,it)=dexcho(ita,it)+
     .              a0*(ek0-ek1+a1*ek0+a2*(ek0-2.d0*ek1+ek2))
C
      tmagoo=tmago
      tmago=tmag
      dexcho(ita,it)=dexch(ita,it)
      IF(lin.EQ.1) THEN
         IF(msgl.NE.0) WRITE(msgio,100) mmom,tmag,dexch(ita,it)
         WRITE(m6,100) mmom,tmag,dexch(ita,it)
      ENDIF
  100 FORMAT(/,' FXSPIN:   MMOM =',f7.3,' TMAG =',f7.3,' DEXC =',f10.6)
C
      RETURN
      END
