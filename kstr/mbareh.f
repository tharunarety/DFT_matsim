      SUBROUTINE mbareh(iq,ir1,ir2,nrlmw)
C   ******************************************************************
C   *                                                                *
C   *  Evaluate the bare slope matrix S0 and energy derivatives .    *
C   *  The convention is that from "mbares", but here we evaluate :  *
C   *                                                                *
C   *    S0(R'H',RL) and  S0(R'H',WL'')                              *
C   *                                                                *
C   *  where : H' = from (0,0) to (lmaxh,lmaxh)                      *
C   *          L  = from (0,0) to (lmax ,lmax )                      *
C   *          L''= from (0,0) to (lmaxw,lmaxw)                      *
C   *                                                                *
C   *  Note :  Only R' = Q' are calculated.                          *
C   ******************************************************************
      USE basis
      USE control_data
      USE control_text
      USE lattice
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: rxp, ryp, rzp, dx, dy, dz, dt, dow, cl, cm,
     .                cn
      INTEGER :: iq, ir1, ir2, nrlmw, irp, jrl, ir, lmp, nlmwh
C
      IF(wats.EQ.'Y') THEN
        ALLOCATE(sbare(nlmh,nlmw,0:nder))
      ELSE
        ALLOCATE(sbare(nlmh,nlm,0:nder))
      ENDIF
      sbd(1:nlmh,1:nrlmw,0:nder)=0.d0
C
      irp=ir1
      rxp=rx(irp)
      ryp=ry(irp)
      rzp=rz(irp)
C
C     Off-site terms
C
      DO 20 ir=irp+1,ir2
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
      CALL gtneum(lmaxh+lmax,dow,cl,cm,cn)
      CALL s0lplk(nlmh)
C
      sbd(1:nlmh,jrl+1:jrl+nlm,0:nder)=sbare(1:nlmh,1:nlm,0:nder)
   20 CONTINUE
      IF(wats.NE.'Y') GO TO 30
C
C     Off-site terms for the Watson sphere (r(W) = r(ir1))
C
      jrl=(ir2-ir1)*nlm+nlm
      nlmwh=MIN0(nlmh,nlmw)
      DO 21 lmp=1,nlmwh
   21 sbd(lmp,jrl+lmp,0)=-1.d0
C
   30 DEALLOCATE(sbare)
      RETURN
      END
