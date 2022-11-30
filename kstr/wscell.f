      SUBROUTINE wscell(si,sc,nv)
C   ******************************************************************
C   *                                                                *
C   *    Set up the Wigner-Seitz cell.                               *
C   *    Calculate the inscribed and circumscribed spheres radii.    *
C   *                                                                *
C   *    The surrounding sites are given by:                         *
C   *                                                                *
C   *     (xn,yn,zn) perpendiculars to the bisector planes;          *
C   *                (halves of the bonds)                           *
C   *     dn         distances to the bisector planes;               *
C   *                                                                *
C   ******************************************************************
      USE message
      USE voronoi
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: tol = 1.d-6, erro = 1.d-10
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: dscp
      INTEGER,      DIMENSION(:), ALLOCATABLE :: iisc
      REAL(KIND=8) :: si, sc
      REAL(KIND=8) :: ri2, rj2, rk2, delta, deltx, delty, deltz
      REAL(KIND=8) :: xsc1, xsc2, xsc3, dr, rxr, div, xxt
      REAL(KIND=8) :: nx, ny, nz, mx, my, mz, mpx, mpy, mpz
      REAL(KIND=8) :: mppx, mppy, mppz, vx, vy, vz, vp
      INTEGER :: nv, iv, jv, kv, ivn, i
      INTEGER :: ntr, itr, isc, jtr, jsc, ktr, ksc
C
      ALLOCATE(dscp(nv),iisc(nv))
      dscp=0.d0
C
C     Calculate the inscribed sphere radius Si as the distance to
C     the closest lying bisector plane 
C
      si=dn(1)
C
C     Find the intersections of three bisector planes
C
      isc=0
      DO 20 iv=1,nvn
      ri2=dn(iv)*dn(iv)
      DO 20 jv=1,nvn
      IF(jv.NE.iv) THEN
         rj2=dn(jv)*dn(jv)
         DO 21 kv=1,nvn
         IF(kv.NE.iv.AND.kv.NE.jv) THEN
            rk2=dn(kv)*dn(kv)
            delta=xn(iv)*yn(jv)*zn(kv)+yn(iv)*zn(jv)*xn(kv)+
     .            zn(iv)*xn(jv)*yn(kv)-zn(iv)*yn(jv)*xn(kv)-
     .            yn(iv)*xn(jv)*zn(kv)-xn(iv)*zn(jv)*yn(kv)
            IF(DABS(delta).GT.tol) THEN
               deltx=ri2   *yn(jv)*zn(kv)+yn(iv)*zn(jv)*rk2+
     .               zn(iv)*rj2   *yn(kv)-zn(iv)*yn(jv)*rk2-
     .               yn(iv)*rj2   *zn(kv)-ri2   *zn(jv)*yn(kv)
               xsc1=deltx/delta
               delty=xn(iv)*rj2   *zn(kv)+ri2   *zn(jv)*xn(kv)+
     .               zn(iv)*xn(jv)*rk2   -zn(iv)*rj2   *xn(kv)-
     .               ri2   *xn(jv)*zn(kv)-xn(iv)*zn(jv)*rk2
               xsc2=delty/delta
               deltz=xn(iv)*yn(jv)*rk2   +yn(iv)*rj2   *xn(kv)+
     .               ri2   *xn(jv)*yn(kv)-ri2   *yn(jv)*xn(kv)-
     .               yn(iv)*xn(jv)*rk2   -xn(iv)*rj2   *yn(kv)
               xsc3=deltz/delta
C
C     Drop those intersections which are farther than any of the planes
C
               i=0
   22          i=i+1
               dr=dn(i)*dn(i)
               rxr=xsc1*xn(i)+xsc2*yn(i)+xsc3*zn(i)
               IF((rxr-dr).GT.erro) GO TO 21
               IF(i.LT.nvn) GO TO 22
C
C     Intersection (xsc1,xsc2,xsc3) is a vertex, save it
C
               IF(isc.EQ.0) THEN
                  isc=isc+1
                  xsc(1,isc)=xsc1
                  xsc(2,isc)=xsc2
                  xsc(3,isc)=xsc3
                  dsc(isc)=DSQRT(xsc1*xsc1+xsc2*xsc2+xsc3*xsc3)
                  dscp(isc)=dsc(isc)
               ELSE
                  jsc=0
   23             jsc=jsc+1
                  div=DABS(xsc1-xsc(1,jsc))+
     .                DABS(xsc2-xsc(2,jsc))+
     .                DABS(xsc3-xsc(3,jsc))
                  IF(div.LT.tol) GO TO 21
                  IF(jsc.LT.isc) GO TO 23
                  isc=isc+1
                  IF(isc.GT.nv) THEN
                     WRITE(m6,101) isc,nv
                     STOP
                  ENDIF
                  xsc(1,isc)=xsc1
                  xsc(2,isc)=xsc2
                  xsc(3,isc)=xsc3
                  dsc(isc)=DSQRT(xsc1*xsc1+xsc2*xsc2+xsc3*xsc3)
                  dscp(isc)=dsc(isc)
               ENDIF
            ENDIF
         ENDIF
   21    CONTINUE
      ENDIF
   20 CONTINUE
C
C     NSC total vertex was found
C
      nsc=isc
C
C     Select those bisector planes which generate these vertexes
C
      ivn=0
      DO 30 iv=1,nvn
      div=dn(iv)*dn(iv)
      ntr=0
      DO 31 isc=1,nsc
      xxt=xsc(1,isc)*xn(iv)+xsc(2,isc)*yn(iv)+xsc(3,isc)*zn(iv)
      IF(DABS(xxt-div).LE.tol) THEN
         ntr=ntr+1
         IF(ntr.GT.nv) THEN
            WRITE(m6,102) ntr,nv
            STOP
         ENDIF
         iisc(ntr)=isc
      ENDIF
   31 CONTINUE
      IF(ntr.GE.3) THEN
         ivn=ivn+1
         nn(ivn)=ntr
         inn(ivn,1:ntr)=iisc(1:ntr)
         xn(ivn)=xn(iv)
         yn(ivn)=yn(iv)
         zn(ivn)=zn(iv)
         dn(ivn)=dn(iv)
         iqn(ivn)=iqn(iv)
      ENDIF
   30 CONTINUE
C
C     NVN bisector planes generate the NSC vertexes
C
      nvn=ivn
C
C     Order (trigonometrical) the vertexes on each facet of the polyhedron
C
      DO 40 iv=1,nvn
      ntr=nn(iv)
      iisc(1:ntr)=inn(iv,1:ntr)
C
C     Normal to the facet
C
      nx=xn(iv)/dn(iv)
      ny=yn(iv)/dn(iv)
      nz=zn(iv)/dn(iv)
C
C     First point M: ITR = 1, ISC = IISC(1)
C
      jsc=iisc(1)
      DO 41 itr=1,ntr-1
      isc=jsc
      mx=xsc(1,isc)
      my=xsc(2,isc)
      mz=xsc(3,isc)
C
C     Take any other vertex but M as the second point M'
C
      DO 42 jtr=1,ntr
      jsc=iisc(jtr)
      IF(jsc.NE.isc) THEN
         mpx=xsc(1,jsc)-mx
         mpy=xsc(2,jsc)-my
         mpz=xsc(3,jsc)-mz
C
C     Test if MM' is an edge and if it is in the trigonometrical direction or not
C
         DO 43 ktr=1,ntr
         ksc=iisc(ktr)
         IF(ksc.NE.isc.AND.ksc.NE.jsc) THEN
            mppx=xsc(1,ksc)-mx
            mppy=xsc(2,ksc)-my
            mppz=xsc(3,ksc)-mz
C
            CALL cross(mpx,mpy,mpz,mppx,mppy,mppz,vx,vy,vz)
            vp=vx*nx+vy*ny+vz*nz
C
C      M' is not the right point, continue searching
C
            IF(vp.LT.erro) GO TO 42
         ENDIF
   43    CONTINUE
C
C      M' is the right point, save and continue with M := M'
C
         GO TO 44
      ENDIF
   42 CONTINUE
C
   44 inn(iv,itr+1)=jsc
   41 CONTINUE
   40 CONTINUE
C
      sc=MAXVAL(dscp)
C
      DEALLOCATE(dscp,iisc)
C
  101 FORMAT(/,' WSCELL:** ISC =',i5,' greater than NV  =',i5)
  102 FORMAT(/,' WSCELL:** NTR =',i5,' greater than NV  =',i5)
  103 FORMAT(/,' WSCELL:** IVP =',i5,' greater than NVN =',i5)
      RETURN
      END
