      SUBROUTINE ovrlps(npr)
C   ******************************************************************
C   *                                                                *
C   *    Set up the distances to the overlapping potential spheres   *
C   *    and the corresponding weights and types.                    *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE csts    ; USE control_data ; USE fandgeq
      USE lattice ; USE message      ; USE radialmesh
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tov0
      INTEGER,      DIMENSION(:), ALLOCATABLE :: type0
      REAL(KIND=8), PARAMETER :: tol = 1.d-6
      REAL(KIND=8) :: s, sjqq, vx, vy, vz, dr, omega
      INTEGER :: npr, iq, it, ita, iv, jqq, jt, iqq, iov, jov
C
C     First establishe the maximum size of the arrays
C
      ALLOCATE(type0(nv),tov0(nv),nov(mnta,nt),ovt(mnta,nt,nt))
C
      nov=0
      DO 20 iq=1,nq
      it=itq(iq)
      DO 21 ita=1,nta(it)
      IF(nov(ita,it).GT.0) GO TO 21
C
      s=hsr(ita,it)
      iov=0
      type0=0
      tov0=0.d0
      DO 22 iv=1,nv
      jqq=jqbas(iv)
      iqq=iqbas(iv)
      IF(iq.EQ.iqq) THEN
         jt=itq(jqq)
         sjqq=SUM(conc(1:nta(jt),jt)*hsr(1:nta(jt),jt))
         vx=tx(iv)+qx(jqq)-qx(iq)
         vy=ty(iv)+qy(jqq)-qy(iq)
         vz=tz(iv)+qz(jqq)-qz(iq)
         dr=SQRT(vx*vx+vy*vy+vz*vz)*alat
         IF(dr.LT.tol) GO TO 22
         IF((s+sjqq-dr).GT.tol) THEN
            jov=0
   23       jov=jov+1
            IF(jov.GT.iov) GO TO 24
            IF(type0(jov).EQ.jt.AND.ABS(tov0(jov)-dr).LT.tol) THEN
               GO TO 22
            ENDIF
            GO TO 23
   24       iov=iov+1
            type0(iov)=jt
            tov0(iov)=dr
         ENDIF
      ENDIF
   22 CONTINUE
      nov(ita,it)=iov
   21 CONTINUE
   20 CONTINUE
C
      dimov=MAXVAL(nov)
      IF(dimov.EQ.0) dimov=1
      DEALLOCATE(type0,tov0)
      ALLOCATE(type(mnta,nt,dimov),wov(mnta,nt,dimov))
      ALLOCATE(tov(mnta,nt,dimov),jrov(mnta,nt,dimov))
      ALLOCATE(okac(mnta,nt,dimov))
C
C     Loop for the atoms from the unit cell
C
      nov=0
      vols=0.d0
      DO 30 iq=1,nq
      WRITE(m6,100) iq
      it=itq(iq)
      DO 31 ita=1,nta(it)
      IF(nov(ita,it).GT.0) GO TO 31
      s=hsr(ita,it)
      iov=0
C
C     Loop for the sites from the cluster
C
      DO 32 iv=1,nv
      jqq=jqbas(iv)
      iqq=iqbas(iv)
      IF(iq.EQ.iqq) THEN
         jt=itq(jqq)
         sjqq=SUM(conc(1:nta(jt),jt)*hsr(1:nta(jt),jt))
C
         vx=tx(iv)+qx(jqq)-qx(iq)
         vy=ty(iv)+qy(jqq)-qy(iq)
         vz=tz(iv)+qz(jqq)-qz(iq)
         dr=SQRT(vx*vx+vy*vy+vz*vz)*alat
         IF(dr.LT.tol) GO TO 32
C
C        Check if the spheres overlap
C
         IF((s+sjqq-dr).GT.tol) THEN
C
C           Check if there is another overlapping sphere of 
C           the same kind and situated at the same distance
C
            jov=0
   33       jov=jov+1
            IF(jov.GT.iov) GO TO 34
            IF(type(ita,it,jov).EQ.jt.AND.
     .         ABS(tov(ita,it,jov)-dr).LT.tol) THEN
               wov(ita,it,jov)=wov(ita,it,jov)+1
               GO TO 32
            ENDIF
            GO TO 33
C
   34       iov=iov+1
            type(ita,it,iov)=jt
            wov(ita,it,iov)=1
            tov(ita,it,iov)=dr
            omega=(s+sjqq)/dr-1.d0
            okac(ita,it,iov)=(dr**5.d0)*(omega**4.d0)
         ENDIF
      ENDIF
   32 CONTINUE
      nov(ita,it)=iov
C
C     Set up the last mesh point before tov-s
C
      IF(nov(ita,it).GT.0) THEN
         DO 35 iov=1,nov(ita,it)
         jt=type(ita,it,iov)
         sjqq=SUM(conc(1:nta(jt),jt)*hsr(1:nta(jt),jt))
         jrov(ita,it,iov)=
     .   DLOG((tov(ita,it,iov)-sjqq+tol)/r1(ita,it))/dx+2
         IF(ri(jrov(ita,it,iov)-1,ita,it).GT.(tov(ita,it,iov)-sjqq))THEN
            WRITE(msgio,*) ri(jrov(ita,it,iov)-1,ita,it),
     .                     tov(ita,it,iov)-sjqq
            STOP
         ENDIF
         IF(ri(jrov(ita,it,iov),ita,it).LT.(tov(ita,it,iov)-sjqq))THEN
            WRITE(msgio,*) ri(jrov(ita,it,iov)-1,ita,it),
     .                     tov(ita,it,iov)-sjqq
            STOP
         ENDIF
   35    CONTINUE
      ELSE
         jrov(ita,it,1)=jsrs(ita,it)
      ENDIF
C
C     Set up the correspondence between neighbour and type
C
      DO 36 iov=1,nov(ita,it)
      jt=type(ita,it,iov)
   36 ovt(ita,it,jt)=iov
C
C     Set up the volume of the potential spheres
C
      vols=vols+mmt(it)*conc(ita,it)*fourpi*s*s*s/3.d0
C
C     Scale the overlap correction to the one electron energies
C
      iov=nov(ita,it)
      IF(iov.GE.1) THEN
         okac(ita,it,1:iov)=pi*okac(ita,it,1:iov)/24.d0
      ENDIF
C
      IF(npr.EQ.1) THEN
         WRITE(m6,110) nov(ita,it)
         DO 37 iov=1,nov(ita,it)
         jt=type(ita,it,iov)
         sjqq=SUM(conc(1:nta(jt),jt)*hsr(1:nta(jt),jt))
         omega=(s+sjqq)/tov(ita,it,iov)-1.d0
   37    WRITE(m6,111) type(ita,it,iov),wov(ita,it,iov),
     .     ovt(ita,it,jt),jrov(ita,it,iov),tov(ita,it,iov),
     .     tov(ita,it,iov)/alat,okac(ita,it,iov),omega*100.d0
      ENDIF
   31 CONTINUE
   30 CONTINUE
C
C     Set up the difference between the total volume and the 
C     volume of the potential spheres
C
      voli=vol-vols
      WRITE(m6,120) vol,vols,voli
C
  100 FORMAT(/,' OVRLPS:   Set up the overlapping spheres for IQ =',i3)
  110 FORMAT(/,11x,'Number of different overlapping spheres NOV =',i3,
     . //,11x,' jt  wov   ovt jrov   tov(Bohr)   tov(a)      ',
     .        'okac     omega',/)
  111 FORMAT(9x,4(2x,i3),1x,f10.6,2x,f10.6,2x,f10.6,2x,f5.2,' %')
  120 FORMAT(/,' OVRLPS:   VOL =',f12.6,' VOLS =',f12.6,' VOLI =',f12.6)
      RETURN
      END
