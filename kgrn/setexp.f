      SUBROUTINE setexp(zm,nzm,nzm2)
C   ******************************************************************
C   *                                                                *
C   *    Set up the closest expansion center for each energy point.  *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE control_text ; USE potential
      USE radialmesh   ; USE slope        ; USE taylor
      IMPLICIT NONE
      INTEGER :: nzm, nzm2, jd, lz
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8), DIMENSION(0:mder) :: ss
      COMPLEX(KIND=8) :: ss1, ss2
      REAL(KIND=8), PARAMETER :: tol =1.d-10
      REAL(KIND=8) :: avmtz, zs2, zs20, s1, s2, enrg
C
      IF(expan.EQ.'S') THEN
         clkw0(1:nzm)=1
         RETURN
      ENDIF
C
      avmtz=SUM(vmtz(1:ns))/ns
      DO 20 lz=1,nzm
C
      IF(zmsh.NE.'M'.AND.zmsh.NE.'m') THEN
         enrg=REAL(zm(lz),8)
         zs20=(enrg-avmtz)*sws*sws 
C
C        First center
C
         zs2=zs20-kw20
         ss(0)=sa0(0)
         DO 21 jd=1,nder
   21    ss(jd)=sa0(jd)*zs2**jd/facd(jd)
         ss1=SUM(ss(0:nder))
         s1=ABS(ss1)
         IF(s1.GT.tol) s1=ABS(ss(nder)/s1)
C
C        Second center
C
         zs2=zs20-kw20p
         ss(0)=sap0(0)
         DO 25 jd=1,nderp
   25    ss(jd)=sap0(jd)*zs2**jd/facd(jd)
         ss2=SUM(ss(0:nderp))
         s2=ABS(ss2)
         IF(s2.GT.tol) s2=ABS(ss(nderp)/s2)
C
C        Choose the more accurate one
C
         IF(s1.LT.s2) THEN
            clkw0(lz)=1
         ELSE
            clkw0(lz)=2
         ENDIF
      ELSE
         IF(lz.LE.nzm2) THEN
            clkw0(lz)=2
         ELSE
            clkw0(lz)=1
         ENDIF
      ENDIF
C
   20 CONTINUE
C
      RETURN
      END
