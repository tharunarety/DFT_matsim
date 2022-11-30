      SUBROUTINE softi(WSATM,R1IT,JWSAIT,NZIT,IT,ita,IS)
C   ******************************************************************
C   *                                                                *
C   *   Read atomic data                                             *
C   *                                                                *
C   ******************************************************************
      USE message
      USE potential
      USE softcore
      USE text
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30,ML=4)
      CHARACTER SYMBOL*4,CONFIG*24,TORBT*4,TEXTA*30,TITLE*4,
     1          TXCH1*20
      COMMON/ADAM/DEP(5),DEQ(5),BE,VC,TVC,DK,DM
      COMMON/ATMM/DEN(MO),DQ1(MO),DFL(MO),DGC(MW,MO),DPC(MW,MO),RMAX,
     2            DR1,FOURPI,RWAT,CHDT(MW),CHDC(MW),DEXV,DEXE,VMIX,
     3            TESTE,TESTY,TESTV,TETS,NODEA(ML),NCORB(MO),NQN(MO),
     4            NQL(MO),NK(MO),NMAX(MO),NEL(MO),NORB,IWAT,IEX,JWSA,
     5            JWOA,JRIA,NITER,IZ,ION,ICEL,NPRNA
      COMMON/DIRA/DV(MW),DR(MW),DP(MW),DQ(MW),DX,Z,NSTOP,NES,TEST,NP
      COMMON/TXTA/TEXTA,TITLE(MO),TXCH1
      DIMENSION TORBT(9),DC(MW)
      DATA TORBT/'s1/2','p3/2','p1/2','d5/2','d3/2','f7/2','f5/2',
     1           'g9/2','g7/2'/
C
      symbol=symbols(ita,it)
      iz=izs(ita,it)
      norb=norbs(ita,it)
      ion=ions(ita,it)
      config=configs(ita,it)
      TEXTA(1:30)=SYMBOL//'  '//CONFIG
      z=iz
      iwat=0
      nqn(1:norb)=nqns(ita,it,1:norb)
      nk(1:norb)=nks(ita,it,1:norb)
      nel(1:norb)=nels(ita,it,1:norb)
      ncorb(1:norb)=ncorbs(ita,it,1:norb)
C
C     Find number of conduction electrons
C
      ICEL=0
      DO 31 J=1,NORB
      IF(NCORB(J).EQ.1) THEN
         ICEL=ICEL+NEL(J)
      ENDIF
   31 CONTINUE
C
C     Generate radial mesh
C
      R1=R1IT
      JWSA=JWSAIT
      NP=1.D0+LOG(RMAX/R1)/DX
      NP=2*(NP/2)+1
      JRIA=JWSA+2
      DR(1)=R1
      DO 20 I=2,NP
  20  DR(I)=DR(1)*EXP((I-1)*DX)
C
      DVAL=Z*Z/(VC*VC)
      DO 22 I=1,NORB
      IF(NCORB(I).EQ.0) THEN
         DEN(I)=dens(ita,it,is,i)
         NQL(I)=IABS(NK(I))
         IF(NK(I).LT.0) NQL(I)=NQL(I)-1
         DFL(I)=NK(I)*NK(I)
         DFL(I)=SQRT(DFL(I)-DVAL)
         L=2*IABS(NK(I))
         J=NQL(I)+IABS(NK(I))
         TITLE(I)=TORBT(J)
      ENDIF
   22 CONTINUE
C
C     Generate potential
C
      DO 23 ir=1,jwsa
      r=dr(ir)
      dv(ir)=v(ir,ita,it,is)/r/r/2.d0
   23 CONTINUE
      dv(jwsa+1:np)=vmtzr(it,is)/2.d0
      dq1(1:norb)=dq1s(ita,it,is,1:norb)
C
      RETURN
      END
