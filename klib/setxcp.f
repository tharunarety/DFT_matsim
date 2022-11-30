      SUBROUTINE SETXCP(IXC,TXCH,EXCHF,NS)
C   ******************************************************************
C   *                                                                *
C   *    Initialize constants and text for XCPOT.                    *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    EXCHF : Slater exchange factor equal to 1 for full exchange *
C   *    IXC   :=1        Barth-Hedin                                *
C   *            2        Slater X-Alpha                             *
C   *            3        Barth-Hedin-Janak                          *
C   *            4        Vosko-Wilk-Nusair                          *
C   *            5        Perdew-Wang 92 LDA                         *
C   *            6        Wigner exchange                            *
C   *            7        Perdew-Zunger                              *
C   *            8        PBE, Pe3rdew et al, 1996                   *
C   *            9        Local Airy Gas Approximation               *
C   *           10        PPBEsol, Perdew et al, 2007                *
C   *    NS    : Number of spins                                     *
C   *                                                                *
C   *    HLS:  31-Oct-96                                             *
C   ******************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      PARAMETER(ONE=1)
      CHARACTER TXCH*3
      COMMON/SXCP/XCCP,XCCF,XCRP,XCRF,XALPHA,OTH,FTH,AA,BB,ALPM,BLPM,
     1            FOURPI,AW,BW,CW,ACA,BCA,CCA,DCA,FCA,OCA,PCA,QCA,RCA,
     2            SCA,TCA,UCA,NSS
  101 FORMAT(/,' SETXCP:** Spin polarization not implemented for IXC',
     1       ' =',I3,' Xcpot =',A)
  102 FORMAT(/,' SETXCP:',3X,'Slater exchange alpha =',F10.6)
  103 FORMAT(/,' SETXCP:** IXC =',I3,' not implemented')
C
      PI=ACOS(-ONE)
      FOURPI=4.D0*PI
      OTH=1.D0/3.D0
      NSS=NS
C
      IF(IXC.EQ.1) THEN
C
C        Barth-Hedin J. Phys. C5,1629(1972)
C
         TXCH='B-H'
         XCCP=0.0504D0
         XCCF=0.0254D0
         XCRP=30.D0
         XCRF=75.D0
         FTH=4.D0/3.D0
         AA=0.5D0**OTH
         BB=1.D0-AA
      ELSEIF(IXC.EQ.2) THEN
C
C        Slater X-Alpha
C
         TXCH='X-A'
         XALPHA=6.D0*EXCHF*(3.D0/FOURPI)**OTH
         WRITE(6,102) EXCHF
      ELSEIF(IXC.EQ.3) THEN
C
C        Barth-Hedin-Janak Phys. Rev. B12,1257(1975)
C
         TXCH='BHJ'
         XCCP=0.045D0
         XCCF=0.0225D0
         XCRP=21.D0
         XCRF=53.D0
         FTH=4.D0/3.D0
         AA=0.5D0**OTH
         BB=1.D0-AA
      ELSEIF(IXC.EQ.4) THEN
C
C        Vosko-Wilk-Nusair Can. J. Phys. 58,1200(1980)
C
         TXCH='VWN'
         FTH=4.D0/3.D0
         AA=2.D0**FTH-2.D0
      ELSEIF(IXC.EQ.5) THEN
C
C        Ceperley-Alder. Parametrization by Perdew and Wang 1992
C
         TXCH='P-W'
      ELSEIF(IXC.EQ.6) THEN
C
C        Wigner exchange
C
         TXCH='WXC'
         IF(NS.EQ.2) THEN
            WRITE(6,101) IXC,TXCH
            STOP
         ENDIF
         AW=0.916D0*4.D0/3.D0
         BW=0.88D0*4.D0/3.D0
         CW=0.88D0*7.8D0/3.D0
      ELSEIF(IXC.EQ.7) THEN
C
C        Ceperley-Alder. Parametrization by Perdew and Zunger
C
         TXCH='P-Z'
         IF(NS.EQ.2) THEN
            WRITE(6,101) IXC,TXCH
            STOP
         ENDIF
         ACA=1.0529D0
         BCA=0.3334D0
         CCA=7.D0*ACA/6.D0
         DCA=4.D0*BCA/3.D0
         FCA=4.D0/3.D0
         OCA=0.096D0
         PCA=0.0622D0
         QCA=0.0232D0
         RCA=0.004D0
         SCA=OCA+PCA/3.D0
         TCA=(2.D0*QCA+RCA)/3.D0
         UCA=2.D0*RCA/3.D0
      ELSEIF(IXC.EQ.8) THEN
         TXCH='PBE'
      ELSEIF(IXC.EQ.9) THEN
         TXCH='LAG'
      ELSEIF(IXC.EQ.10) THEN
         TXCH='P07'
      ELSE
         WRITE(6,103) IXC
         STOP
      ENDIF
      RETURN
      END
