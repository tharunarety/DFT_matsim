      SUBROUTINE GAUNTY(LAM,MY,LAMP,MYP,LAMPP,MYPP,GNTA)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the Gaunt coefficient :                           *
C   *                                                                *
C   *    I(l,m;l',m';l'',m'')=int[Y(l,m)*Y(l',m')*Y(l'',m'')]        *
C   *                                                                *
C   *    where Y(l,m) are the real harmonics.                        *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      COMPLEX(KIND=8) :: SUM, TY, TYP, TYPP, FTY
      REAL(KIND=8), PARAMETER :: SQFPI = 3.544907702D0
      REAL(KIND=8) :: GNTA, FACL, FACLM, GAUNTC, GNTI
      INTEGER :: LAM, MY, LAMP, MYP, LAMPP, MYPP, MYA, MYS, MYPA, MYPS,
     1           MYPPA, MPPABS, M, MP, MPP
  104 FORMAT(/,' GAUNTY:** Gaunt coefficient should be real.',/,11X,
     1       'GNT =',2E15.6)
C
      FACL=DSQRT(1.D0+2*LAMPP)
      SUM=(0.D0,0.D0)
      MYA=ABS(MY)
      MYS=2*MYA
      IF(MYA.EQ.0) MYS=1
      MYPA=ABS(MYP)
      MYPS=2*MYPA
      IF(MYPA.EQ.0) MYPS=1
      MYPPA=ABS(MYPP)
      DO 21 M=-MYA,MYA,MYS
      TY=FTY(MY,M)
      FACLM=FACL*(-1)**M
      DO 22 MP=-MYPA,MYPA,MYPS
      TYP=FTY(MYP,MP)
      MPP=MP+M
      MPPABS=ABS(MPP)
      IF(MPPABS.EQ.MYPPA) THEN
         TYPP=FTY(MYPP,MPP)
         SUM=SUM+FACLM*GAUNTC(LAMPP,LAMP,MP,LAM,-M)
     1                      *CONJG(TY)*CONJG(TYP)*TYPP
      ENDIF
   22 CONTINUE
   21 CONTINUE
      GNTI=AIMAG(SUM)
      IF(GNTI.GT.1.D-06) THEN
         WRITE(*,104) SUM
         STOP
      ENDIF
      GNTA=SUM/SQFPI
      RETURN
      END
C
