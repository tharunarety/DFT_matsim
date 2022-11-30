      SUBROUTINE slfcls(iq,nvsl,zsl,rsl)
C   ******************************************************************
C   *                                                                *
C   *    Set up the latice vectors arround the site IQ.              *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE control_data
      USE lattice
      USE message
      USE radialmesh
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: tol = 1.d-6
      REAL(KIND=8), DIMENSION(nv,4) :: rsl
      REAL(KIND=8), DIMENSION(nv)   :: zsl
      REAL(KIND=8) :: vx, vy, vz, dr
      INTEGER :: iq, jv, iv, jqq, iqq, itp, nvsl
C
      jv=0
      DO 21 iv=1,nv
        jqq=jqbas(iv)
        iqq=iqbas(iv)
        IF(iq.EQ.iqq) THEN
          vx=tx(iv)+qx(jqq)-qx(iq)
          vy=ty(iv)+qy(jqq)-qy(iq)
          vz=tz(iv)+qz(jqq)-qz(iq)
          dr=DSQRT(vx*vx+vy*vy+vz*vz)
          IF(dr.GT.tol) THEN
            jv=jv+1
            rsl(jv,1)=vx/dr
            rsl(jv,2)=vy/dr
            rsl(jv,3)=vz/dr
            rsl(jv,4)=dr*alat
            itp=itq(jqq)
            zsl(jv)=SUM(conc(1:nta(itp),itp)*nz(1:nta(itp),itp))
          ENDIF
        ENDIF
   21 CONTINUE
      nvsl=jv
C
      RETURN
      END
