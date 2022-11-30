      SUBROUTINE ATOME(it,ita,etotv)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the valence contribution to the total energy      *
C   *    within the local density approximation.                     *
C   *                                                                *
C   ******************************************************************
      USE softcore
      USE message
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30,ML=4)
      COMMON/ATMM/DEN(MO),DQ1(MO),DFL(MO),DGC(MW,MO),DPC(MW,MO),RMAX,
     2            DR1,FOURPI,RWAT,CHDT(MW),CHDC(MW),DEXV,DEXE,VMIX,
     3            TESTE,TESTY,TESTV,TETS,NODEA(ML),NCORB(MO),NQN(MO),
     4            NQL(MO),NK(MO),NMAX(MO),NEL(MO),NORB,IWAT,IEX,JWSA,
     5            JWOA,JRIA,NITER,IZ,ION,ICEL,NPRNA
      COMMON/DIRA/DV(MW),DR(MW),DP(MW),DQ(MW),DX,Z,NSTOP,NES,TEST,NP
      DIMENSION RHOA(MW),RHOP(MW),RHOPP(MW),CHDV(MW),FINT(MW),WC(MW),
     1          RHOD(2),RHODD(2)
      DATA CHDV/MW*0.D0/
  101 FORMAT(/,' TOTALE:   (Rydberg)',
     2      //,11X,'EONE =',F14.6,' VINT =',F14.6,' EKIN =',F14.6,
     1       /,11X,'ENUC =',F14.6,' ECOR =',F14.6,
     2       /,11X,'EVAL =',F14.6,' EXCT =',F14.6,' EXCC =',F14.6,
     3       /,11X,'ETOT =',F14.6)
  102 FORMAT(/,11X,'Valence energies in the frozen core approximation',
     1       ' (Rydbergs)')
C
C     Calculate the core one-electron energies for soft-core
C
      eonec(ita,it)=0.d0
      DO 19 J=1,NORB
      IF(NCORB(J).EQ.0) THEN
         eonec(ita,it)=eonec(ita,it)+NEL(J)*DEN(J)
      ENDIF
   19 CONTINUE
C
C     Transform from Hartree to Rydberg
C
      eonec(ita,it)=2.d0*eonec(ita,it)
C
C     Sum of one-electron energies
C
      EONE=0.D0
      NELV=0
      DO 22 J=1,NORB
      IF(NCORB(J).EQ.1) THEN
         NELV=NELV+NEL(J)
         EONE=EONE+NEL(J)*DEN(J)
      ENDIF
   22 CONTINUE
      EONE=2.D0*EONE
C
      DO 24 IR=1,NP
   24 CHDV(IR)=CHDT(IR)-CHDC(IR)
C
C     Integral of total potential
C
      DO 25 IR=1,NP
   25 FINT(IR)=DV(IR)*CHDV(IR)*DR(IR)
      CALL SIMPN(FINT,DX,NP,VINT)
      VINT=2.D0*VINT
C
C     Integral of nuclear potential
C
      CALL SIMPN(CHDV,DX,NP,FI)
      ENUC=-2.D0*Z*FI
C
C     Integral of electrostatic core potential
C
      POTB=2.D0*(Z-NELV)/DR(NP)
      CALL POISON(CHDC,DR,WC,POTB,NP,DX)
      DO 33 IR=1,NP
   33 FINT(IR)=WC(IR)*CHDV(IR)*DR(IR)
      CALL SIMPN(FINT,DX,NP,ECOR)
C
C     Integral of electrostatic valence potential
C
      POTB=2.D0*NELV/DR(NP)
      CALL POISON(CHDV,DR,WC,POTB,NP,DX)
      DO 35 IR=1,NP
   35 FINT(IR)=WC(IR)*CHDT(IR)*DR(IR)
      CALL SIMPN(FINT,DX,NP,EVAL)
C
C     Total exchange-correlation potential
C
      IXC=IEX+1
      DO 36 IR=1,NP
      FPRR=FOURPI*DR(IR)**2
      CHRT=CHDT(IR)
      IF(CHRT.LT.1.E-8) CHRT=1.D-8
   36 RHOA(IR)=CHRT/FPRR
      IF(IXC.EQ.5.OR.IXC.GE.8) CALL DIFFN(RHOA,RHOP,RHOPP,NP,DX)
      DO 37 IR=1,NP
      R=DR(IR)
      RCE=R*R
      RHO=RHOA(IR)
      RHO1=0.5D0*RHO
      RHO2=RHO1
      IF(IXC.EQ.5.OR.IXC.GE.8) THEN
         RHOD(1)=RHOP(IR)/R
         RHODD(1)=(RHOPP(IR)-RHOP(IR))/RCE
         RHOD(1)=RHOD(1)/2.D0
         RHOD(2)=RHOD(1)
         RHODD(1)=RHODD(1)/2.D0
         RHODD(2)=RHODD(1)
      ENDIF
      CALL XCPOT(IXC,RHO1,RHO2,RHO,RHOD,RHODD,R,VXC1,VXC2,EXC)
      FINT(IR)=EXC*CHDT(IR)*DR(IR)
   37 CONTINUE
      CALL SIMPN(FINT,DX,NP,EXCT)
C
C     Core exchange-correlation potential
C
      DO 38 IR=1,NP
      FPRR=FOURPI*DR(IR)**2
      CORT=CHDC(IR)
      IF(CORT.LT.1.E-8) CORT=1.D-8
   38 RHOA(IR)=CORT/FPRR
      IF(IXC.EQ.5.OR.IXC.GE.8) CALL DIFFN(RHOA,RHOP,RHOPP,NP,DX)
      DO 39 IR=1,NP
      R=DR(IR)
      RCE=R*R
      RHO=RHOA(IR)
      RHO1=0.5D0*RHO
      RHO2=RHO1
      IF(IXC.EQ.5.OR.IXC.GE.8) THEN
         RHOD(1)=RHOP(IR)/R
         RHODD(1)=(RHOPP(IR)-RHOP(IR))/RCE
         RHOD(1)=RHOD(1)/2.D0
         RHOD(2)=RHOD(1)
         RHODD(1)=RHODD(1)/2.D0
         RHODD(2)=RHODD(1)
      ENDIF
      CALL XCPOT(IXC,RHO1,RHO2,RHO,RHOD,RHODD,R,VXC1,VXC2,EXC)
      FINT(IR)=EXC*CHDC(IR)*DR(IR)
   39 CONTINUE
      CALL SIMPN(FINT,DX,NP,EXCC)
   30 CONTINUE
      EKIN=EONE-VINT
      ETOT=EKIN+ENUC+0.5D0*ECOR+0.5D0*EVAL+EXCT-EXCC
      etotv=etot
      WRITE(M6,101) EONE,VINT,EKIN,ENUC,ECOR,EVAL,EXCT,EXCC,ETOT
      WRITE(M6,102)
      WRITE(M6,'(/,11X,A,F14.6)') 'Kinetic energy   ',EKIN
      WRITE(M6,'(11X,A,F14.6)') 'El-ion energy    ',ENUC
      WRITE(M6,'(11X,A,F14.6)') 'El-el energy     ',0.5*(ECOR+EVAL)
      WRITE(M6,'(11X,A,F14.6)') 'Exc energy       ',EXCT-EXCC
      WRITE(M6,'(/,11X,A,F14.6)') 'Total energy     ',ETOT
      RETURN
      END
      SUBROUTINE ATOMt(it,ita,etota)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the total atomic energy.                          *
C   *                                                                *
C   ******************************************************************
      USE softcore
      USE message
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MW=999,MO=30,ML=4)
      COMMON/ATMM/DEN(MO),DQ1(MO),DFL(MO),DGC(MW,MO),DPC(MW,MO),RMAX,
     2            DR1,FOURPI,RWAT,CHDT(MW),CHDC(MW),DEXV,DEXE,VMIX,
     3            TESTE,TESTY,TESTV,TETS,NODEA(ML),NCORB(MO),NQN(MO),
     4            NQL(MO),NK(MO),NMAX(MO),NEL(MO),NORB,IWAT,IEX,JWSA,
     5            JWOA,JRIA,NITER,IZ,ION,ICEL,NPRNA
      COMMON/DIRA/DV(MW),DR(MW),DP(MW),DQ(MW),DX,Z,NSTOP,NES,TEST,NP
      DIMENSION RHOA(MW),RHOP(MW),RHOPP(MW),CHDV(MW),FINT(MW),WC(MW),
     1          RHOD(2),RHODD(2)
      DATA CHDV/MW*0.D0/
  101 FORMAT(/,' TOTALT:   (Rydberg)',
     2      //,11X,'EONE =',F14.6,' VINT =',F14.6,' EKIN =',F14.6,
     1       /,11X,'ENUC =',F14.6,' ECOR =',F14.6,
     2       /,11X,'EVAL =',F14.6,' EXCT =',F14.6,' EXCC =',F14.6,
     3       /,11X,'ETOT =',F14.6)
  102 FORMAT(/,11X,'Valence energies in the frozen core approximation',
     1       ' (Rydbergs)')
C
C     Calculate the core one-electron energies for soft-core
C
      eonec(ita,it)=0.d0
      DO 19 J=1,NORB
      eonec(ita,it)=eonec(ita,it)+NEL(J)*DEN(J)
   19 CONTINUE
C
C     Transform from Hartree to Rydberg
C
      eonec(ita,it)=2.d0*eonec(ita,it)
C
C     Sum of one-electron energies
C
      EONE=0.D0
      NELV=0
      DO 22 J=1,NORB
      NELV=NELV+NEL(J)
      EONE=EONE+NEL(J)*DEN(J)
   22 CONTINUE
      EONE=2.D0*EONE
C
      DO 24 IR=1,NP
   24 CHDV(IR)=CHDT(IR)-CHDC(IR)
C
C     Integral of total potential
C
      DO 25 IR=1,NP
   25 FINT(IR)=DV(IR)*CHDT(IR)*DR(IR)
      CALL SIMPN(FINT,DX,NP,VINT)
      VINT=2.D0*VINT
C
C     Integral of nuclear potential
C
      CALL SIMPN(CHDT,DX,NP,FI)
      ENUC=-2.D0*Z*FI
C
C     Integral of electrostatic core potential
C
      POTB=2.D0*(Z-NELV)/DR(NP)
      CALL POISON(CHDC,DR,WC,POTB,NP,DX)
      DO 33 IR=1,NP
   33 FINT(IR)=WC(IR)*CHDT(IR)*DR(IR)
      CALL SIMPN(FINT,DX,NP,ECOR)
C
C     Integral of electrostatic valence potential
C
      POTB=2.D0*NELV/DR(NP)
      CALL POISON(CHDV,DR,WC,POTB,NP,DX)
      DO 35 IR=1,NP
   35 FINT(IR)=WC(IR)*CHDT(IR)*DR(IR)
      CALL SIMPN(FINT,DX,NP,EVAL)
C
C     Total exchange-correlation potential
C
      IXC=IEX+1
      DO 36 IR=1,NP
      FPRR=FOURPI*DR(IR)**2
      CHRT=CHDT(IR)
      IF(CHRT.LT.1.E-8) CHRT=1.D-8
   36 RHOA(IR)=CHRT/FPRR
      IF(IXC.EQ.5.OR.IXC.GE.8) CALL DIFFN(RHOA,RHOP,RHOPP,NP,DX)
      DO 37 IR=1,NP
      R=DR(IR)
      RCE=R*R
      RHO=RHOA(IR)
      RHO1=0.5D0*RHO
      RHO2=RHO1
      IF(IXC.EQ.5.OR.IXC.GE.8) THEN
         RHOD(1)=RHOP(IR)/R
         RHODD(1)=(RHOPP(IR)-RHOP(IR))/RCE
         RHOD(1)=RHOD(1)/2.D0
         RHOD(2)=RHOD(1)
         RHODD(1)=RHODD(1)/2.D0
         RHODD(2)=RHODD(1)
      ENDIF
      CALL XCPOT(IXC,RHO1,RHO2,RHO,RHOD,RHODD,R,VXC1,VXC2,EXC)
      FINT(IR)=EXC*CHDT(IR)*DR(IR)
   37 CONTINUE
      CALL SIMPN(FINT,DX,NP,EXCT)
      EXCC=0.d0
C
      EKIN=EONE-VINT
      ETOT=EKIN+ENUC+0.5D0*ECOR+0.5D0*EVAL+EXCT-EXCC
      etota=etot
      WRITE(M6,101) EONE,VINT,EKIN,ENUC,ECOR,EVAL,EXCT,EXCC,ETOT
      WRITE(M6,102)
      WRITE(M6,'(/,11X,A,F14.6)') 'Kinetic energy   ',EKIN
      WRITE(M6,'(11X,A,F14.6)') 'El-ion energy    ',ENUC
      WRITE(M6,'(11X,A,F14.6)') 'El-el energy     ',0.5*(ECOR+EVAL)
      WRITE(M6,'(11X,A,F14.6)') 'Exc energy       ',EXCT-EXCC
      WRITE(M6,'(/,11X,A,F14.6)') 'Total energy     ',ETOT
      RETURN
      END
