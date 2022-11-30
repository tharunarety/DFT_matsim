      SUBROUTINE expans(zm,nzm)
C   ******************************************************************
C   *                                                                *
C   *   Generate the Taylor expansion's coefficients.                *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE control_text ; USE potential
      USE radialmesh   ; USE taylor
      IMPLICIT NONE
      INTEGER         :: nzm, is, jd, lz
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8) :: zs2
C
      DO 20 is=1,ns
      DO 20 lz=1,nzm
      zs2=(zm(lz)-vmtz(is))*sws*sws
      IF(clkw0(lz).EQ.1) THEN
         zs2=zs2-kw20
      ELSE
         zs2=zs2-kw20p
      ENDIF
C
      DO 21 jd=1,mder
      tayl(lz,jd,is)=zs2**jd/facd(jd)
   21 tayld(lz,jd,is)=zs2**(jd-1)/facd(jd-1)
   20 CONTINUE
      tayld=tayld*sws*sws
C
      RETURN
      END
