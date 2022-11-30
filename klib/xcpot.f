      SUBROUTINE XCPOT(IXC,RHO1,RHO2,RHO,RHOP,RHOPP,RR,V1,V2,EXC)
C   ******************************************************************
C   *                                                                *
C   *    Calculates exchange-correlation potential according to the  *
C   *    value of IXC. The constants used in the various expressions *
C   *    are set in SETXCP and in the data statement below.          *
C   *                                                                *
C   *    HLS: 31-Oct-96                                              *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(TOLD=1.D-20,TOLDD=1.D-20)
      COMMON/SXCP/XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,ALPM,BLPM,
     1            FOURPI,AW,BW,CW,ACA,BCA,CCA,DCA,FCA,OCA,PCA,QCA,RCA,
     2            SCA,TCA,UCA,NSS
      DIMENSION RHOS(2),RHOP(2),RHOPP(2)
      DATA AP,XP0,BP,CP,QP,CP1,CP2,CP3/0.0621814,-0.10498,3.72744,
     1     12.9352,6.1519908,1.2117833,1.1435257,-0.031167608/
      DATA AF,XF0,BF,CF,QF,CF1,CF2,CF3/0.0310907,-0.32500,7.060428,
     2     18.0578,4.7309269,2.9847935,2.7100059,-0.1446006/
  100 FORMAT(/,' XCPOT:**  RR less than 1. D-10 for IXC =',I3)
C
      IF(RHO1.LT.TOLD.OR.RHO2.LT.TOLDD) THEN
         V1=0.D0
         V2=0.D0
         EXC=0.D0
         RETURN
      ENDIF
      RS1=(FOURPI*RHO/3.)**OTH
      RS=1./RS1
      GO TO (1,2,1,3,4,5,6,7,8,9),IXC
C
C     Barth-Hedin  J. PHYS. C5,1629(1972)
C
    1 RSF=RS/XCRF
      RSF2=RSF*RSF
      RSF3=RSF2*RSF
      RSP=RS/XCRP
      RSP2=RSP*RSP
      RSP3=RSP2*RSP
      FCF=(1.D0+RSF3)*LOG(1.D0+1.D0/RSF)+.5D0*RSF-RSF2-OTH
      FCP=(1.D0+RSP3)*LOG(1.D0+1.D0/RSP)+.5D0*RSP-RSP2-OTH
      EPSCP=-XCCP*FCP
      EPSCF=-XCCF*FCF
      EPSXP=-.91633059D0/RS
      CNY=5.1297628D0*(EPSCF-EPSCP)
      X=RHO1/RHO
      FX=(X**FTH+(1.D0-X)**FTH-AA)/BB
      EXC=EPSXP+EPSCP+FX*(CNY+FTH*EPSXP)/5.1297628D0
      ARS=-1.22177412D0/RS+CNY
      BRS=-XCCP*LOG(1.D0+XCRP/RS)-CNY
      TRX1=(2.D0*X)**OTH
      V1=ARS*TRX1+BRS
      TRX2=(2.D0*RHO2/RHO)**OTH
      V2=ARS*TRX2+BRS
      RETURN
C
C     Slater X-Alpha
C
    2 EXC=-0.75D0*XALPHA*(0.5D0*RHO)**OTH
      V1=-XALPHA*(RHO1)**OTH
      V2=-XALPHA*(RHO2)**OTH
      RETURN
C
C     Vosko-Wilk-Nusair  Can. J. Phys. 58,1200(1980)
C
    3 X=SQRT(RS)
      XPX=X*X+BP*X+CP
      XFX=X*X+BF*X+CF
      S=(RHO2-RHO1)/RHO
      SP=1.D0+S
      SM=1.D0-S
      S4=S**4-1.D0
      FS=(SP**FTH+SM**FTH-2.D0)/AA
      BETA=1.D0/(2.74208D0+3.182D0*X+0.09873D0*X*X+0.18268D0*X**3)
      DFS=FTH*(SP**OTH-SM**OTH)/AA
      DBETA=-(0.27402D0*X+0.09873D0+1.591D0/X)*BETA**2
      ATNP=ATAN(QP/(2.D0*X+BP))
      ATNF=ATAN(QF/(2.D0*X+BF))
      ECP=AP*(LOG(X*X/XPX)+CP1*ATNP-CP3*(LOG((X-XP0)**2/XPX)+CP2*ATNP
     1))
      ECF=AF*(LOG(X*X/XFX)+CF1*ATNF-CF3*(LOG((X-XF0)**2/XFX)+CF2*ATNF
     1))
      EC=ECP+FS*(ECF-ECP)*(1.D0+S4*BETA)
      TP1=(X*X+BP*X)/XPX
      TF1=(X*X+BF*X)/XFX
      UCP=ECP-AP/3.D0*(1.D0-TP1-CP3*(X/(X-XP0)-TP1-XP0*X/XPX))
      UCF=ECF-AF/3.D0*(1.D0-TF1-CF3*(X/(X-XF0)-TF1-XF0*X/XFX))
      UC0=UCP+(UCF-UCP)*FS
      UC20=UC0+(ECF-ECP)*SM*DFS
      UC10=UC0-(ECF-ECP)*SP*DFS
      DUC=(UCF-UCP)*BETA*S4*FS+(ECF-ECP)*(-RS/3.D0)*DBETA*S4*FS
      DUC2=DUC+(ECF-ECP)*BETA*SM*(4.D0*S**3*FS+S4*DFS)
      DUC1=DUC-(ECF-ECP)*BETA*SP*(4.D0*S**3*FS+S4*DFS)
      UC1=UC10+DUC1
      UC2=UC20+DUC2
      EPX=-0.91633059D0/RS*(1.D0+FTH*FS/5.1297628D0)
      AMYX2=-1.22177412D0/RS*SP**OTH
      AMYX1=-1.22177412D0/RS*SM**OTH
      EXC=EC+EPX
      V1=UC1+AMYX1
      V2=UC2+AMYX2
      RETURN
C
C     Perdew, Burke, and Ernzerhof, submiited to PRL, May96
C
C     LDA version
C
    4 IF(RR.LT.1.D-10) THEN
         WRITE(6,100) IXC
         STOP
      ENDIF
      LGGA=0
      RHOS(1)=RHO1
      RHOS(2)=RHO2
      CALL PBEGGA(LGGA,RHOS,RHOP,RHOPP,RR,EXC,V1,V2)
      RETURN
C
C     Wigner expression
C
    5 RS78=1.D0/(RS+7.8D0)
      EXC=-0.916D0*RS1-0.88D0*RS78
c
      V1=CW*RS78*RS78-AW*RS1-BW*RS78
      V2=V1
      RETURN
C
C     Ceperley-Alder. Parametrization by Perdew and Zunger.
C
    6 IF(RS.GE.1.D0) THEN
         SQRTRS=SQRT(RS)
         DENOM1=1.D0/(1.D0+ACA*SQRTRS+BCA*RS)
         EX=-0.9164D0*RS1
         EC=-0.2846D0*DENOM1
         EXC=EX+EC
         V1=FCA*EX+EC*(1.D0+CCA*SQRTRS+DCA*RS)*DENOM1
         V2=V1
      ELSE
         RSLOG=LOG(RS)
         RSLN=RS*RSLOG
         EX=-0.9164D0*RS1
         EC=-OCA+PCA*RSLOG-QCA*RS+RCA*RSLN
         EXC=EX+EC
         V1=FCA*EX-SCA+PCA*RSLOG-TCA*RS+UCA*RSLN
         V2=V1
      ENDIF
      RETURN
C
C     Perdew, Burke, and Ernzerhof, PBE 1996
C
C     GGA version
C
    7 IF(RR.LT.1.D-10) THEN
         WRITE(6,100) IXC
         STOP
      ENDIF
      LGGA=1
      RHOS(1)=RHO1
      RHOS(2)=RHO2
      CALL PBEGGA(LGGA,RHOS,RHOP,RHOPP,RR,EXC,V1,V2)
      RETURN
C
C     Local Airy Gas for exchange
C     Ceperley-Alder (by Perdew-Wang) for correlation
C
    8 IF(RR.LT.1.D-10) THEN
         WRITE(6,100) IXC
         STOP
      ENDIF
      RHOS(1)=RHO1
      RHOS(2)=RHO2
      CALL lagexc(RHOS,RHOP,RHOPP,RR,EXC,V1,V2)
      RETURN
C
C     Perdew, et al. PBEsol 2007
C
    9 IF(RR.LT.1.D-10) THEN
         WRITE(6,100) IXC
         STOP
      ENDIF
      LGGA=1
      RHOS(1)=RHO1
      RHOS(2)=RHO2
      CALL PBEsol(LGGA,RHOS,RHOP,RHOPP,RR,EXC,V1,V2)
      RETURN
C
      END
      SUBROUTINE PBEGGA(LGGA,N,ND,NDD,R,EXC,MUXC1,MUXC2)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the exchange-correlation energy density by        *
C   *    Perdew-Burke-Ernzerhof generalized gradient approximation.  *
C   *                                                                *
C   *   *On entry:                                                   *
C   *       n(r)   = the charge density per spin (N) ;               *
C   *       n'(r)  = dn/dr (ND) ;                                    *
C   *       n''(r) = d^2n/dr^2 (NDD) .                               *
C   *       LGGA = 0 LDA (local density approximation) ;             *
C   *            = 1 GGA (generalized gradient approximation) .      *
C   *   *On exit:                                                    *
C   *       exc    = the exchange-correlation energy density,        *
C   *       muxc   = the exchange-correlation potential.             *
C   *                                                                *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 N(2),ND(2),NDD(2),R,EXC,MUXC1,MUXC2
      REAL*8 KF,EXI,EX,NI,NDI,NDDI,NABLA,NABLA2,S,T,U,MUXI
      REAL*8 RS,EC,ZET,FK,SK,G,H,UU,VV,WW,VCUP,VCDN,DVCUP,DVCDN
      INTEGER I,LGGA,LPOT
      COMMON/SXCP/XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,ALPM,BLPM,
     1            FOURPI,AW,BW,CW,ACA,BCA,CCA,DCA,FCA,OCA,PCA,QCA,RCA,
     2            SCA,TCA,UCA,NSS
      DIMENSION VX(2)
C
      PI=0.25D0*FOURPI
      LPOT=1
C
C     Exchange energy and potential
C
      EX=0.D0
      DO 10 I=1,2
      NI=N(I)+N(I)
      EXI=0.D0
      MUXI=0.D0
      NDI=ND(I)+ND(I)
      NDDI=NDD(I)+NDD(I)
      KF=(3.D0*PI*PI*NI)**OTH
      NABLA=DABS(NDI)
      S=0.5D0*NABLA/KF/NI
      NABLA2=2.D0/R*NDI+NDDI
      T=NABLA2/4.D0/KF/KF/NI
      U=NABLA*NDDI/8.D0/KF/KF/KF/NI/NI
      CALL EXCHPBE(NI,S,U,T,LGGA,LPOT,EXI,MUXI)
      VX(I)=MUXI
   10 EX=EX+N(I)*EXI
C
C     Correlation energy and potential
C
      NI=N(1)+N(2)
      EC=0.D0
      VCUP=0.D0
      DVCUP=0.D0
      H=0.D0
      VCDN=0.D0
      DVCDN=0.D0
      NDI=ND(1)+ND(2)
      NDDI=NDD(1)+NDD(2)
      ZET=(N(1)-N(2))/NI
      G=((1.D0+ZET)**(2.D0/3.D0)+(1.D0-ZET)**(2.D0/3.D0))/2.D0
      NABLA=DABS(NDI)
      NABLA2=2.D0/R*NDI+NDDI
      FK=(3.D0*PI*PI*NI)**OTH
      SK=DSQRT(4.D0*FK/PI)
      T=NABLA/2.D0/SK/NI/G
      UU=NABLA*NDDI/((2.D0*SK*G)**3.D0)/NI/NI
      VV=NABLA2/((2.D0*SK*G)**2.D0)/NI
      WW=(NDI*ND(1)-NDI*ND(2)-ZET*NDI*NDI)/((2.D0*SK*G)**2.D0)
     1  /NI/NI
      RS=(3.D0/FOURPI/NI)**OTH
C
      CALL CORPBE(RS,ZET,T,UU,VV,WW,LGGA,LPOT,EC,VCUP,VCDN,
     1                                          H,DVCUP,DVCDN)
C
C     Convert to Rydberg
C
      MUXC1=2.D0*(VX(1)+VCUP+DVCUP)
      MUXC2=2.D0*(VX(2)+VCDN+DVCDN)
      EX=2.D0*EX/(N(1)+N(2))
      EC=2.D0*(EC+H)
C
      EXC=EX+EC
C
      RETURN
      END
      SUBROUTINE lagexc(N,ND,NDD,R,EXC,MUXC1,MUXC2)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the exchange-correlation energy density by        *
C   *    Local Airy Gas approximation for the exchange energy and    *
C   *    Ceperlay-Alder (Perdew-Wang) local correlation functional.  *
C   *                                                                *
C   *   *On entry:                                                   *
C   *       n(r)   = the charge density per spin (N) ;               *
C   *       n'(r)  = dn/dr (ND) ;                                    *
C   *       n''(r) = d^2n/dr^2 (NDD) .                               *
C   *   *On exit:                                                    *
C   *       exc    = the exchange-correlation energy density,        *
C   *       muxc   = the exchange-correlation potential.             *
C   *                                                                *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 N(2),ND(2),NDD(2),R,EXC,MUXC1,MUXC2
      REAL*8 KF,EXI,EX,NI,NDI,NDDI,NABLA,NABLA2,S,T,U,MUXI
      REAL*8 RS,EC,ZET,FK,SK,G,H,UU,VV,WW,VCUP,VCDN,DVCUP,DVCDN
      INTEGER I,POT
      COMMON/SXCP/XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,ALPM,BLPM,
     1            FOURPI,AW,BW,CW,ACA,BCA,CCA,DCA,FCA,OCA,PCA,QCA,RCA,
     2            SCA,TCA,UCA,NSS
      DIMENSION VX(2)
C
      PI=0.25D0*FOURPI
      POT=1
      LGGA=0
C
C     Exchange energy and potential
C
      EX=0.D0
      DO 10 I=1,2
      NI=N(I)+N(I)
      EXI=0.D0
      MUXI=0.D0
      NDI=ND(I)+ND(I)
      NDDI=NDD(I)+NDD(I)
      KF=(3.D0*PI*PI*NI)**OTH
      NABLA=DABS(NDI)
      S=0.5D0*NABLA/KF/NI
      NABLA2=2.D0/R*NDI+NDDI
      T=NABLA2/4.D0/KF/KF/NI
      U=NABLA*NDDI/8.D0/KF/KF/KF/NI/NI
      CALL exchlag(NI,S,U,T,POT,EXI,MUXI)
      VX(I)=MUXI
   10 EX=EX+N(I)*EXI
C
C     Correlation energy and potential
C
      NI=N(1)+N(2)
      EC=0.D0
      VCUP=0.D0
      DVCUP=0.D0
      H=0.D0
      VCDN=0.D0
      DVCDN=0.D0
      NDI=ND(1)+ND(2)
      NDDI=NDD(1)+NDD(2)
      ZET=(N(1)-N(2))/NI
      G=((1.D0+ZET)**(2.D0/3.D0)+(1.D0-ZET)**(2.D0/3.D0))/2.D0
      NABLA=DABS(NDI)
      NABLA2=2.D0/R*NDI+NDDI
      FK=(3.D0*PI*PI*NI)**OTH
      SK=DSQRT(4.D0*FK/PI)
      T=NABLA/2.D0/SK/NI/G
      UU=NABLA*NDDI/((2.D0*SK*G)**3.D0)/NI/NI
      VV=NABLA2/((2.D0*SK*G)**2.D0)/NI
      WW=(NDI*ND(1)-NDI*ND(2)-ZET*NDI*NDI)/((2.D0*SK*G)**2.D0)
     1  /NI/NI
      RS=(3.D0/FOURPI/NI)**OTH
C
      CALL CORPBE(RS,ZET,T,UU,VV,WW,LGGA,POT,EC,VCUP,VCDN,
     1                                       H,DVCUP,DVCDN)
C
C     Convert to Rydberg
C
      MUXC1=2.D0*(VX(1)+VCUP+DVCUP)
      MUXC2=2.D0*(VX(2)+VCDN+DVCDN)
      EX=2.D0*EX/(N(1)+N(2))
      EC=2.D0*(EC+H)
C
      EXC=EX+EC
C
      RETURN
      END
      SUBROUTINE PBEsol(LGGA,N,ND,NDD,R,EXC,MUXC1,MUXC2)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the exchange-correlation energy density by        *
C   *    PBEsol generalized gradient approximation.                  *
C   *                                                                *
C   *   *On entry:                                                   *
C   *       n(r)   = the charge density per spin (N) ;               *
C   *       n'(r)  = dn/dr (ND) ;                                    *
C   *       n''(r) = d^2n/dr^2 (NDD) .                               *
C   *       LGGA = 0 LDA (local density approximation) ;             *
C   *            = 1 GGA (generalized gradient approximation) .      *
C   *   *On exit:                                                    *
C   *       exc    = the exchange-correlation energy density,        *
C   *       muxc   = the exchange-correlation potential.             *
C   *                                                                *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 N(2),ND(2),NDD(2),R,EXC,MUXC1,MUXC2
      REAL*8 KF,EXI,EX,NI,NDI,NDDI,NABLA,NABLA2,S,T,U,MUXI
      REAL*8 RS,EC,ZET,FK,SK,G,H,UU,VV,WW,VCUP,VCDN,DVCUP,DVCDN
      INTEGER I,LGGA,LPOT
      COMMON/SXCP/XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,ALPM,BLPM,
     1            FOURPI,AW,BW,CW,ACA,BCA,CCA,DCA,FCA,OCA,PCA,QCA,RCA,
     2            SCA,TCA,UCA,NSS
      DIMENSION VX(2)
C
      PI=0.25D0*FOURPI
      LPOT=1
C
C     Exchange energy and potential
C
      EX=0.D0
      DO 10 I=1,2
      NI=N(I)+N(I)
      EXI=0.D0
      MUXI=0.D0
      NDI=ND(I)+ND(I)
      NDDI=NDD(I)+NDD(I)
      KF=(3.D0*PI*PI*NI)**OTH
      NABLA=DABS(NDI)
      S=0.5D0*NABLA/KF/NI
      NABLA2=2.D0/R*NDI+NDDI
      T=NABLA2/4.D0/KF/KF/NI
      U=NABLA*NDDI/8.D0/KF/KF/KF/NI/NI
      CALL EXCHsol(NI,S,U,T,LGGA,LPOT,EXI,MUXI)
      VX(I)=MUXI
   10 EX=EX+N(I)*EXI
C
C     Correlation energy and potential
C
      NI=N(1)+N(2)
      EC=0.D0
      VCUP=0.D0
      DVCUP=0.D0
      H=0.D0
      VCDN=0.D0
      DVCDN=0.D0
      NDI=ND(1)+ND(2)
      NDDI=NDD(1)+NDD(2)
      ZET=(N(1)-N(2))/NI
      G=((1.D0+ZET)**(2.D0/3.D0)+(1.D0-ZET)**(2.D0/3.D0))/2.D0
      NABLA=DABS(NDI)
      NABLA2=2.D0/R*NDI+NDDI
      FK=(3.D0*PI*PI*NI)**OTH
      SK=DSQRT(4.D0*FK/PI)
      T=NABLA/2.D0/SK/NI/G
      UU=NABLA*NDDI/((2.D0*SK*G)**3.D0)/NI/NI
      VV=NABLA2/((2.D0*SK*G)**2.D0)/NI
      WW=(NDI*ND(1)-NDI*ND(2)-ZET*NDI*NDI)/((2.D0*SK*G)**2.D0)
     1  /NI/NI
      RS=(3.D0/FOURPI/NI)**OTH
C
      CALL CORsol(RS,ZET,T,UU,VV,WW,LGGA,LPOT,EC,VCUP,VCDN,
     1                                          H,DVCUP,DVCDN)
C
C     Convert to Rydberg
C
      MUXC1=2.D0*(VX(1)+VCUP+DVCUP)
      MUXC2=2.D0*(VX(2)+VCDN+DVCDN)
      EX=2.D0*EX/(N(1)+N(2))
      EC=2.D0*(EC+H)
C
      EXC=EX+EC
C
      RETURN
      END
