      SUBROUTINE mbares(iq,ir1,ir2,nrlm)
C   ******************************************************************
C   *                                                                *
C   *  Evaluate the bare slope matrix S0 and energy derivatives .    *
C   *  The output is :  B = t3*t4 + t3*S0*t3 where t3 and t4 are     *
C   *  the screening parameters calculated in "screen".              *
C   *                                                                *
C   *  The convention is that the first index of the matrix is       *
C   *  primed and the second unprimed: i.e. we evaluate              *
C   *                                                                *
C   *                  S0(R',L';R,L)                                 *
C   *                                                                *
C   *  for the vectors R'-R = Q'-Q-T shorter than dmax. Additionally *
C   *  a Watson sphere with radius rwats = dmax + ws*sigma + dwats   *
C   *  is introduced. The new S0 has the structure:                  *
C   *                                                                *
C   *        S0(R',R)   S0(R',W)                                     *
C   *                                                                *
C   *        S0(W ,R)   S0(W ,W)                                     *
C   *                                                                *
C   *  where : N_R = - J_R' * S0(R',R)                               *
C   *          J_W = - J_R' * S0(R',W)                               *
C   *          N_R = - N_W  * S0(W ,R)                               *
C   *          S0(W W) = zero.                                       *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE factorial
      USE control_data
      USE control_text
      USE lattice
      USE lmtolm
      USE screening
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: rxp, ryp, rzp, dx, dy, dz, dt, dow, cl, cm,
     .                cn, t3p, t3, sum
      INTEGER :: ir1, ir2, irl, irp, jqp, irlm, lmp, lp,
     .           jrl, ir, jq, jrlm, lm, l, jd, id, nrlm, iq
C
      irl=-nlm
      DO 20 irp=ir1,ir2
      jqp=jqbas(irp)
      irl=irl+nlm
      rxp=rx(irp)
      ryp=ry(irp)
      rzp=rz(irp)
C
C     On-site term  (t3*t4)
C
      DO 21 lm=1,nlm
      l=llx(lm)
      irlm=irl+lm
      dett(irlm)=dfac(l,jqp)
      sb(irlm,irlm)=cd(1,l,jqp,0)
   21 sbd(irlm,irlm,1:nder)=cd(1,l,jqp,1:nder)
C
C     Off-site terms (t3*S0*t3)
C
      DO 22 ir=irp+1,ir2
      jq=jqbas(ir)
      jrl=(ir-ir1)*nlm
      dx=rx(ir)-rxp
      dy=ry(ir)-ryp
      dz=rz(ir)-rzp
      dt=dsqrt(dx*dx+dy*dy+dz*dz)
      dow=dt/ws
      cl=dx/dt
      cm=dy/dt
      cn=dz/dt
C
      CALL gtneum(lmax+lmax,dow,cl,cm,cn)
      CALL s0lplk(nlm)
C
      DO 23 lmp=1,nlm
      lp=llx(lmp)
      t3p=tmat(3,lp,jqp,0)
      irlm=irl+lmp
      DO 23 lm=1,nlm
      l=llx(lm)
      t3=tmat(3,l,jq,0)
      jrlm=jrl+lm
      sb(irlm,jrlm)=t3p*sbare(lmp,lm,0)*t3
      sc(0:nder)=0.d0
      DO 24 jd=0,nder
      DO 24 id=0,jd
   24 sc(jd)=sc(jd)+ifib(jd,id)*sbare(lmp,lm,id)*
     1                tmat(3,l,jq,jd-id)
      DO 25 jd=1,nder
      sum=0.d0
      DO 26 id=0,jd
   26 sum=sum+ifib(jd,id)*tmat(3,lp,jqp,id)*sc(jd-id)
   25 sbd(irlm,jrlm,jd)=sum
   23 CONTINUE
   22 CONTINUE
      IF(wats.NE.'Y') GO TO 20
C
C     Off-site terms for the Watson sphere (r(W) = r(ir1))
C
      dx=rx(ir1)-rxp
      dy=ry(ir1)-ryp
      dz=rz(ir1)-rzp
      dt=dsqrt(dx*dx+dy*dy+dz*dz)
      dow=dt/ws
      IF(DABS(dt).GT.err) THEN
         cl=dx/dt
         cm=dy/dt
         cn=dz/dt
      ELSE
         cl=0.d0
         cm=0.d0
         cn=0.d0
      ENDIF
C
      CALL gtbess(lmax+lmaxw,dow,cl,cm,cn)
      CALL s0lplj(dow)
C
      jrl=nrlm
      DO 27 lmp=1,nlm
      lp=llx(lmp)
      t3p=tmat(3,lp,jqp,0)
      irlm=irl+lmp
      DO 27 lm=1,nlmw
      l=llx(lm)
      t3=-alwats(l,iq,0)
      jrlm=jrl+lm
      sb(irlm,jrlm)=t3p*sbare(lmp,lm,0)*t3
      sc(0:nder)=0.d0
      DO 28 jd=0,nder
      DO 28 id=0,jd
   28 sc(jd)=sc(jd)-ifib(jd,id)*sbare(lmp,lm,id)*
     .                alwats(l,iq,jd-id)
      DO 29 jd=1,nder
      sum=0.d0
      DO 30 id=0,jd
   30 sum=sum+ifib(jd,id)*tmat(3,lp,jqp,id)*sc(jd-id)
   29 sbd(irlm,jrlm,jd)=sum
   27 CONTINUE
   20 CONTINUE
C
C     On-site term for the Watson sphere
C
      irl=nrlm
      DO 31 lm=1,nlmw
      l=llx(lm)
      irlm=irl+lm
      dett(irlm)=1.d0
      sb(irlm,irlm)=-alwats(l,iq,0)
      sbd(irlm,irlm,1:nder)=-alwats(l,iq,1:nder)
   31 CONTINUE
C
      RETURN
      END
