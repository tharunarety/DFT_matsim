      SUBROUTINE MADL3D
C   ******************************************************************
C   *                                                                *
C   *    Calculates the Madelung potential matrix which is used to   *
C   *    set up the Madelung potential in the self-consistent bulk   *
C   *    programs.                                                   *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE csts
      USE madelung_lattice
      USE madelung_matrix
      USE message
      IMPLICIT NONE
      REAL(KIND=8) :: ALAMSQ, FACTOR, SUM0G, DGSQ, BETASQ, AMDL, SUM0R,
     1                DRI, ALPHA, ERFCA, QPX, QPY, QPZ, QPPX, QPPY,
     2                QPPZ, GDOTQ, EXPB, X, Y, Z, DERFC
      INTEGER :: I, IQ, JQ, NQM, JQP
  101 FORMAT(/,' MADL3D:',3X,'Madelung constant with compensating',
     1       ' background:',//,11X,'CMDL = ',F12.8)
  102 FORMAT(' ',10X,6F12.8)
  103 FORMAT(/,11X,'Madelung constant (ionic) :',//,11X,'AMDL = ',
     1       F12.8,//,11X,'Potentials:',//,11X,'IQ',5X,'Q(IQ)',7X,
     2       'VQ(IQ)',/)
  104 FORMAT(' ',10X,I2,2F12.8,4X,2F12.8)
  105 FORMAT(/,11X,'Madelung monopole potential(s):',//)
C
      ALLOCATE(VMADL(NQ,NQ))
      ALAMSQ=ALAMDA*ALAMDA
      FACTOR=TWOPI/VOL/ALAMSQ/2.
C
C     QPP .EQ. 0
C
      SUM0G=0.D0
      DO 20 I=2,NUMVG
      DGSQ=DG(I)**2
      BETASQ=DGSQ/4./ALAMSQ
      SUM0G=SUM0G+EXP(-BETASQ)/BETASQ
   20 CONTINUE
C
      AMDL=FACTOR*(SUM0G-1.0D0)
C
      SUM0R=0.D0
      DO 21 I=2,NR0
      DRI=DR(I)
      ALPHA=ALAMDA*DRI
      ERFCA=DERFC(ALPHA)
      SUM0R=SUM0R+ERFCA/DRI
   21 CONTINUE
C
      AMDL=AMDL+SUM0R-2.D0*ALAMDA/SQRTPI
      DO 22 IQ=1,NQ
   22 VMADL(IQ,IQ)=AMDL*WS*2.D0
      IF(NQ.EQ.1)  GO TO 26
C
C     QPP .NE. 0
C
      NQM=NQ-1
      DO 23 JQ=1,NQM
      JQP=JQ+1
      QPX=QX(JQ)
      QPY=QY(JQ)
      QPZ=QZ(JQ)
      DO 23 IQ=JQP,NQ
      QPPX=QX(IQ)-QPX
      QPPY=QY(IQ)-QPY
      QPPZ=QZ(IQ)-QPZ
      SUM0G=0.D0
      DO 24 I=2,NUMVG
      DGSQ=DG(I)**2
      GDOTQ=AKX(I)*QPPX+AKY(I)*QPPY+AKZ(I)*QPPZ
      BETASQ=DGSQ/4.D0/ALAMSQ
      EXPB=EXP(-BETASQ)
      SUM0G=SUM0G+EXPB/BETASQ*COS(GDOTQ)
   24 CONTINUE
C
      AMDL=FACTOR*(SUM0G-1.0D0)
C
      SUM0R=0.D0
      DO 25 I=1,NUMVR
      X=ASX(I)-QPPX
      Y=ASY(I)-QPPY
      Z=ASZ(I)-QPPZ
      DRI=SQRT(X*X+Y*Y+Z*Z)
      IF(DRI.GT.RMAX) GO TO 25
      ALPHA=ALAMDA*DRI
      ERFCA=DERFC(ALPHA)
      SUM0R=SUM0R+ERFCA/DRI
   25 CONTINUE
      AMDL=AMDL+SUM0R
      VMADL(IQ,JQ)=2.D0*WS*AMDL
   23 CONTINUE
C
C     Madelung potentials
C
   26 DO 28 JQ=1,NQ
      DO 28 IQ=1,JQ
   28 VMADL(IQ,JQ)=VMADL(JQ,IQ)
C
C     Madelung constant for ASA correction
C
      AMDL=0.D0
      DO 27 JQ=1,NQ
      DO 27 IQ=1,NQ
   27 AMDL=AMDL+VMADL(IQ,JQ)
      CMDL=-0.5D0*AMDL/NQ
      WRITE(M6,101) CMDL
      WRITE(M6,105)
      DO 29 IQ=1,NQ
   29 WRITE(M6,102) (VMADL(IQ,JQ),JQ=1,NQ)
C
      deallocate(asx,asy,asz,dr)
      deallocate(akx,aky,akz,dg)
C
      RETURN
      END
