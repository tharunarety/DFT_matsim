      SUBROUTINE cdirac(e,k,it,ita,il,is,ispo,beta)
C   ******************************************************************
C   *                                                                *
C   *    Solves the Dirac equation. The constants in /SDRC/ are set  *
C   *    in SETDRC                                                   *
C   *                                                                *
C   *    ISPO                                                        *
C   *     0      Without spin-orbit interaction                      *
C   *    ELSE    With spin-orbit interaction                         *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE diracparam
      USE message
      USE partialwaves
      USE potential
      USE radialmesh
      IMPLICIT NONE
      COMPLEX(KIND=8), DIMENSION(dimr) :: ps, qs
      COMPLEX(KIND=8)  :: dq3, dp3, dp4, qsn, psn, dq4, dq2, vnp1,
     .                 dq1, dp1, vmnp1, dp2, vmh, vh, qp, pp,
     .                 psnp1, qrel, prel, qsnp1, qsnm4, psnm2,
     .                 qsnm1, psnm1, qsnm2, psnm4, qsnm3, psnm3,
     .                 r, vmn, vn, e, dlta, sb2, sa2, c3, qc, 
     .                 pc, sb1, c1, c2, sa1, sb0, alocsq
      REAL(KIND=8) :: twoz, beta
      INTEGER      :: k, it, ita, il, is, ispo
      INTEGER      :: n, n1, nit ,jws, kp1
C
      IF(ISPO.EQ.0.AND.K.GT.0) THEN
         WRITE(m6,'(A,I3)') ' ***DIRAC: No spin-orbit but KAPPA =',K
         STOP
      ENDIF
      JWS=JSRS(ita,IT)
      R=RI(1,ita,IT)
C
C     V = RADIUS**2*POTENTIAL
C
      TWOZ=-V(1,ita,IT,IS)/R
      VN=V(1,ita,IT,IS)
      VMN=VP(1,ita,IT,IS)
      KP1=K+1
      C1=VN/R**2-E
      C2=(1.D0,0.D0)-C1/CSQ
      P(1)=(1.D-5,0.D0)
C
C     Empty sphere case, i.e. TWOZ .LT. 2.
C
      IF(TWOZ.LT.2.) THEN
         ALOCSQ=(0.D0,0.D0)
         Q(1)=C1/C2/(2.D0*(-K-1)+3.D0)*R*P(1)
         BETA=IABS(K)
      ELSEIF(ISPO.EQ.0) THEN
C
C        Without spin-orbit
C
         ALOCSQ=-KP1/CSQ
         BETA=SQRT(KP1*K+1.D0-(TWOZ/CLIGHT)**2)
         SB0=(BETA+K)*CSQ/TWOZ
         SA1=(BETA-1.D0-(TWOZ/CLIGHT)**2)/(2.D0*BETA+1.D0)
         SB1=-SB0+CSQ/TWOZ*(BETA+KP1)*SA1
         SA2=(-2.D0*KP1*SA1+3.D0*KP1+TWOZ/CSQ*(BETA-K+2.D0)*SB1)/
     1        4.D0/(BETA+1.D0)
         SB2=-SB1+CSQ/TWOZ*(BETA+K+2.D0)*SA2
         DLTA=R*CSQ/TWOZ
         Q(1)=(SB0+DLTA*(SB1+DLTA*SB2))/(1.D0+DLTA*(SA1+DLTA*SA2))*
     1         P(1)
      ELSE
C
C           With spin-orbit
C
            ALOCSQ=(0.D0,0.D0)
            BETA=SQRT(K*K-(TWOZ/CLIGHT)**2)
            SB0=(BETA+K)*CSQ/TWOZ
            SA1=(BETA+K-(TWOZ/CLIGHT)**2)/(2.D0*BETA+1.D0)
            SB1=-SB0+CSQ/TWOZ*(BETA+K+1.D0)*SA1
            SA2=(BETA-K+2.D0)/4.D0/(BETA+1.D0)*TWOZ/CSQ*SB1
            SB2=-SB1+CSQ/TWOZ*(BETA+K+2.D0)*SA2
            DLTA=R*CSQ/TWOZ
            Q(1)=(SB0+DLTA*(SB1+DLTA*SB2))/(1.D0+DLTA*(SA1+DLTA*SA2))
     1          *P(1)
      ENDIF
C
      C3=(VMN-2.D0*VN)/C2/C2*ALOCSQ
      PS(1)=R*C2*Q(1)-K*P(1)
      QS(1)=K*Q(1)+(R*C1-C3/R**3)*P(1)
C
      DO 21 N=1,5
      PC=P(N)
      QC=Q(N)
      DP1=DX*(R*C2*QC-K*PC)
      DQ1=DX*(K*QC+(R*C1-C3/R**3)*PC)
      PC=PC+.5D0*DP1
      QC=QC+.5D0*DQ1
      R=R*EXPDXH
      VNP1=V(N+1,ita,IT,IS)
      VMNP1=VP(N+1,ita,IT,IS)
      VH=(VN+VNP1)*.5D0+(VMN-VMNP1)*DXD8
      VMH=1.5D0*(VNP1-VN)/DX-(VMN+VMNP1)/4.D0
      C1=VH/R/R-E
      C2=1.D0-C1/CSQ
      C3=(VMH-2.D0*VH)/C2/C2*ALOCSQ
      DP2=DX*(R*C2*QC-K*PC)
      DQ2=DX*(K*QC+(R*C1-C3/R**3)*PC)
      PC=PC+.5D0*(DP2-DP1)
      QC=QC+.5D0*(DQ2-DQ1)
      DP3=DX*(R*C2*QC-K*PC)
      DQ3=DX*(K*QC+(R*C1-C3/R**3)*PC)
      PC=PC+DP3-.5D0*DP2
      QC=QC+DQ3-.5D0*DQ2
      N1=N+1
      R=RI(N1,ita,IT)
      C1=VNP1/R/R-E
      C2=1.D0-C1/CSQ
      C3=(VMNP1-2.D0*VNP1)/C2/C2*ALOCSQ
      DP4=DX*(R*C2*QC-K*PC)
      DQ4=DX*(K*QC+(R*C1-C3/R**3)*PC)
      P(N1)=P(N)+(DP1+2.D0*(DP2+DP3)+DP4)/6.D0
      Q(N1)=Q(N)+(DQ1+2.D0*(DQ2+DQ3)+DQ4)/6.D0
      PS(N1)=R*C2*Q(N1)-K*P(N1)
      QS(N1)=K*Q(N1)+(R*C1-C3/R**3)*P(N1)
      VN=VNP1
      VMN=VMNP1
   21 CONTINUE
C
      PSN  =PS(6)
      QSN  =QS(6)
      PSNM1=PS(5)
      QSNM1=QS(5)
      PSNM2=PS(4)
      QSNM2=QS(4)
      PSNM3=PS(3)
      QSNM3=QS(3)
      PSNM4=PS(2)
      QSNM4=QS(2)
      DO 22 N=7,JWS
      R=RI(N,ita,IT)
      C1=V(N,ita,IT,IS)/R/R-E
      C2=1.D0-C1/CSQ
      C3=(VP(N,ita,IT,IS)-2.D0*V(N,ita,IT,IS))/C2/C2*ALOCSQ
      PP=P(N-6)+A1*(PSN+PSNM4)+A2*(PSNM1+PSNM3)+A3*PSNM2
      QP=Q(N-6)+A1*(QSN+QSNM4)+A2*(QSNM1+QSNM3)+A3*QSNM2
      NIT=0
   25 PSNP1=R*C2*QP-K*PP
      QSNP1=K*QP+(R*C1-C3/R**3)*PP
      PC=P(N-4)+A4*(PSNP1+PSNM3)+A5*(PSN+PSNM2)+A6*PSNM1
      QC=Q(N-4)+A4*(QSNP1+QSNM3)+A5*(QSN+QSNM2)+A6*QSNM1
      IF(ABS(TEST*(PC-PP)).GT.ABS(PC)) GO TO 24
      IF(ABS(TEST*(QC-QP)).LE.ABS(QC)) GO TO 23
   24 IF(NIT.EQ.5) GO TO 26
      NIT=NIT+1
      PP=PC
      QP=QC
      GO TO 25
   26 PREL=(PC-PP)/PC
      QREL=(QC-QP)/QC
      WRITE(m6,100) IT,ITA,IL,IS,E,N,PREL,QREL
   23 P(N)=PC
      Q(N)=QC
      PS(N)=PSNP1
      QS(N)=QSNP1
      PSNM4=PSNM3
      PSNM3=PSNM2
      PSNM2=PSNM1
      PSNM1=PSN
      PSN  =PSNP1
      QSNM4=QSNM3
      QSNM3=QSNM2
      QSNM2=QSNM1
      QSNM1=QSN
      QSN  =QSNP1
   22 CONTINUE
C
  100 FORMAT(' Hard test in <CDIRAC> : IT=',i2,' ITA=',i2,
     .       ' IL=',i2,' IS=',i2,' E=',2f10.6,' N=',i3,
     .       ' PREL=',2E12.5,' QREL=',2E12.5,'****')
      RETURN
      END
