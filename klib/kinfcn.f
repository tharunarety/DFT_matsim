C     Last change:  HLS  29 Jun 98   11:19 am
      SUBROUTINE KINFCN(FUNC,IL,RHO,RHOD,RHODD,TS)
C   ******************************************************************
C   *                                                                *
C   *    Set up the kinetic energy density: Ts[n]=int t[n]dr (Ry).   *
C   *                                                                *
C   *                                               L.V. 27.03.1997. *
C   ******************************************************************
      USE csts
      IMPLICIT NONE
      CHARACTER(LEN=3) :: FUNC
      REAL(KIND=8) :: RHO, RHOD, RHODD, TS, TFC, C0, TS0, FK, S, S2,
     1                S4, S6, TS2, C4, A1, A2, TS4, FMR, X, FMY, PX,
     2                UK, UM, FPBE, ALPH, BETA, B, BS, P0, P1, FBC,
     3                C, FTK, FPW, A, D, FDK, GAM, FTK1, FTK2, FTK3,
     4                FTK4, FB, FC, FD, CD, TSN
      INTEGER :: IL
C
      IF(RHO.LT.1.D-16) THEN
         TS=0.D0
         RETURN
      ENDIF
      TFC=(3.D0*PI*PI)**(2.D0/3.D0)
      C0=3.D0/5.D0*TFC
      TS0=C0*(RHO**(5.D0/3.D0))
      FK=(3.D0*PI*PI*RHO)**(1.D0/3.D0)
      S=DABS(RHOD)/2.D0/FK/RHO
      S2=S*S
      S4=S2*S2
      S6=S4*S2
      TS2=5.D0/27.D0*TS0*S*S
C
C     Thomas-Fermi 
C     from : Dreizler-Gross, Density Functional Theory;Springer-Verlag (1990)
C
      IF(FUNC.EQ.'TS0') THEN
         TS=TS0
         GO TO 30
      ENDIF
C
C     Second order gradient expansion
C     from : Dreizler-Gross, Density Functional Theory;Springer-Verlag (1990)
C
      IF(FUNC.EQ.'TS2') THEN
         TS=TS0+TS2
         GO TO 30
      ENDIF
C
C     Fourth order gradient expansion
C     from : Dreizler-Gross, Density Functional Theory;Springer-Verlag (1990)
C
      IF(FUNC.EQ.'TS4') THEN
         C4=1.D0/270.D0/((3.D0*PI*PI)**(2.D0/3.D0))
         A1=RHODD/RHO
         A2=RHOD/RHO
         TS4=C4*(RHO**(1.D0/3.D0))*(A1*A1
     1      -9.D0/8.D0*A1*A2*A2+1.D0/3.D0*A2*A2*A2*A2)
         TS=TS0+TS2+TS4
         GO TO 30
      ENDIF
C
C     Murphy's expansion
C     from : J.P.Perdew; Physics Letters A 165, 79 (1992)
C
      IF(FUNC.EQ.'TMR') THEN
         FMR=1.D0+5.D0/27.D0*S2+8.D0/243.D0*S4
         TS=TS0*FMR
         GO TO 30
      ENDIF
C
C     Locally truncated functional
C     from : E.W.Pearson and R.G.Gordon; J.Chem.Phys. 82, 881 (1985)
C
      IF(FUNC.EQ.'TLT') THEN
         IF(TS2.LT.TS0) THEN
            TS=TS0+TS2
         ELSE
            TS=TS0
         ENDIF
         GO TO 30
      ENDIF
C
C     Pearson smooth cutoff with chi=1
C     from : D.J.Lacks and R.G.Gordon; J.Chem.Phys. 100, 4446 (1993)
C
      IF(FUNC.EQ.'PRS') THEN
         TS   =TS0+TS2/(1.D0+S6)
         GO TO 30
      ENDIF
C
C     Weizsacker functional
C     from : Dreizler-Gross, Density Functional Theory;Springer-Verlag (1990)
C
      IF(FUNC.EQ.'TWZ') THEN
         TS=TS0+9.D0*TS2
         GO TO 30
      ENDIF
C
C     Local functional : t[n] = t^0[n] + t^2[n]   if  t^2[n]/t^0[n] < 1
C                        t[n] = 9 * t^2[n]        if  t^2[n]/t^0[n] > 1
C
      IF(FUNC.EQ.'LOC') THEN
         X=TS2/TS0
         TS=TS0+TS2
         IF(X.GT.1.D0) TS=9.D0*TS2
         GO TO 30
      ENDIF
C
C     Parametrized local functional
C
      IF(FUNC.EQ.'PWZ') THEN
         X=TS2/TS0
         FMY=(1.D0+0.95D0*X+9.D0*0.396D0*X*X*X)/
     1                (1.D0-0.05D0*X+0.396D0*X*X)
         TS=TS0*FMY
         GO TO 30
      ENDIF
C
C     Pade's functional
C     from : A.E.DePristo and J.D.Kress; Phys.Rev. B, 35, 438 (1987)
C
      IF(FUNC.EQ.'PDE') THEN
         X=TS2/TS0
         CALL PADEFN(X,PX)
         TS=TS0*PX
         GO TO 30
      ENDIF
C
C     Perdew-96 functional
C     from : J.P.Perdew, K.Burke and M.Ernzerhof; Phys.Rev.Lett. (1997)
C
      IF(FUNC.EQ.'PBE') THEN
         UK=0.804d0
         UM=0.2195149727645171D0
         FPBE=1.D0+UK-UK/(1.D0+UM/UK*S2)
         TS=TS0*FPBE
         GO TO 30
      ENDIF
C
C     Becke functional
C     from : H.Lee,C.Lee and R.G.Parr; Phys.Rev. A 44,768 (1991)
C     and    A.D.Becke; Phys.Rev. A, 38, 3098 (1988)
C
      IF(FUNC.EQ.'B88') THEN
         ALPH=0.0045D0
         BETA=0.0253D0
         B=(48.D0*PI*PI)**(1.D0/3.D0)
         BS=B*S
         P0=DSQRT(1.D0+BS*BS)
         P1=DLOG(BS+P0)
         FBC=1.D0+ALPH*BS*BS/(1.D0+BETA*BS*P1)
         TS=TS0*FBC
         GO TO 30
      ENDIF
C
C     Thakkar functional
C     from : A.J.Thakkar; Phys.Rev. A, 46, 6920 (1992)
C
      IF(FUNC.EQ.'THK') THEN
         B=(48.D0*PI*PI)**(1.D0/3.D0)
         C=32.D0**(1.D0/3.D0)
         BS=B*S
         P0=DSQRT(1.D0+BS*BS)
         P1=DLOG(BS+P0)
         FTK=1.D0+0.0055D0*BS*BS/(1.D0+0.0253D0*BS*P1)-
     1         0.072D0*BS/(1.D0+C*BS)
         TS=TS0*FTK
         GO TO 30
      ENDIF
C
C     Perdew-86 functional
C     from : J.P.Perdew and W.Yue; Phys.Rev. B 33, 8800 (1986)
C
      IF(FUNC.EQ.'P86') THEN
         FPW=(1.D0+1.296D0*S2+14.D0*S4+0.2D0*S6)**(1.D0/15.D0)
         TS=TS0*FPW
         GO TO 30
      ENDIF
C
C     Perdew-91 functional
C     from : J.P.Perdew et.al.; Phys.Rev. B 46, 6671 (1992)
C
      IF(FUNC.EQ.'P91') THEN
         A=0.19645D0
         B=7.7956D0
         P0=1.D0/DSQRT(1.D0+B*B*S2)
         P1=DLOG(B*S+1.D0/P0)
         FPW=(1.D0+A*S*P1+(0.2743-0.1508*DEXP(-100.D0*S2))*S2)/
     1         (1.D0+A*S*P1+0.004D0*S4)
         TS=TS0*FPW
         GO TO 30
      ENDIF
C
C     Lembarki-Chermette functional
C     from : J.P.Perdew et.al.; Phys.Rev. B 46, 6671 (1992)
C     and    A.Lembarki and H.Chermette; Phys.Rev. A, 50, 5328 (1994)
C
      IF(FUNC.EQ.'TLC') THEN
         A=0.093907D0
         B=76.320D0
         P0=1.D0/DSQRT(1.D0+B*B*S2)
         P1=DLOG(B*S+1.D0/P0)
         FPW=(1.D0+A*S*P1+(0.26608D0-0.0809615D0*
     1         DEXP(-100.D0*S2))*S2)/(1.D0+A*S*P1+0.000057767D0*S4)
         TS=TS0*FPW
         GO TO 30
      ENDIF
C
C     DePristo and Kress
C     from : A.E.DePristo and J.D.Kress; J.Chem.Phys. 96, 1425 (1987)
C
      IF(FUNC.EQ.'DPK') THEN
         B=(48.D0*PI*PI)**(1.D0/3.D0)
         C=324.D0*((18.D0*PI)**(1.D0/3.D0))*PI
         D=7.D0/C
         FDK=1.D0+D*B*B*S2*(1.D0+0.861504D0*B*S)/
     1         (1.D0+0.044286D0*B*B*S2)
         TS=TS0*FDK
         GO TO 30
      ENDIF
C
C     Ou-Yang and Levy
C     from : H.Ou-Yang and M.Levy; Int.J.Quantum Chem. 40, 379 (1991)
C
      IF(FUNC.EQ.'OL1') THEN
         B=(48.D0*PI*PI)**(1.D0/3.D0)
         TS=TS0+TS2+0.00187D0*TS0*B*S
         GO TO 30
      ENDIF
C
C     Ou-Yang and Levy
C     from : H.Ou-Yang and M.Levy; Int.J.Quantum Chem. 40, 379 (1991)
C
      IF(FUNC.EQ.'OL2') THEN
         B=(48.D0*PI*PI)**(1.D0/3.D0)
         C=32.D0**(1.D0/3.D0)
         TS=TS0+TS2+0.0245D0*TS0*B*S/(1.D0+C*B*S)
         GO TO 30
      ENDIF
C
C     Plumer and Stott functional
C     from : RM.L.Plumer and M.J.Stott; J.Phys.C 18, 4143 (1985)
C
      IF(FUNC.EQ.'TPS') THEN
         BETA=1.0D0
         X=BETA/6.D0*S
         GAM=1.D0-8.D0/9.D0/(1.D0+5.D0/3.D0*X*X)
         TS=TS0+9.D0*GAM*TS2
         GO TO 30
      ENDIF
C
C     modified Thakkar functional
C
      IF(FUNC.EQ.'MTH') THEN
         B=(48.D0*PI*PI)**(1.D0/3.D0)
         C=32.D0**(1.D0/3.D0)
         BS=B*S
         P0=DSQRT(1.D0+BS*BS)
         P1=DLOG(BS+P0)
         FTK1=1.D0
         FTK2=BS*BS/(1.D0+0.0253D0*BS*P1)
         FTK3=BS/(1.D0+C*BS)
         FTK4=BS
         FB=.0093431669D0
         FC=-.2474215994D0
         FD=.0062554713D0
         TS=TS0*(FTK1+FB*FTK2+FC*FTK3+FD*FTK4)
         GO TO 30
      ENDIF
C
C     modified DePristo and Kress
C
      IF(FUNC.EQ.'MDK') THEN
         B=(48.D0*PI*PI)**(1.D0/3.D0)
         C=324.D0*((18.D0*PI)**(1.D0/3.D0))*PI
         D=7.D0/C
         FDK=1.D0+D*B*B*S2*(1.D0+0.868396D0*B*S)/
     1         (1.D0+0.024800D0*B*B*S2)
         TS=TS0*FDK
         GO TO 30
      ENDIF
C
      WRITE(*,*) 'KINFCN:  FUNC=',FUNC,' is not implemented, if you',
     1           ' know please do it'
      STOP
C
   30 CONTINUE
      IF(IL.EQ.0) RETURN
C
C     The grad^2 term in case of extended systems
C     from : Dreizler-Gross, Density Functional Theory;Springer-Verlag (1990)
C
      CD =-1.D0/6.D0
      TSN=CD*RHODD
      TS=TS+TSN
C
      RETURN
      END
