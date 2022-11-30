      SUBROUTINE blatts
C   ******************************************************************
C   *                                                                *
C   *     Set up the Voronoi polyhedron for each lattice site and    *
C   *     determine the radius of the inscribed and circumscribed    *
C   *     spheres.                                                   *
C   *                                                                *
C   *     Sc :  circumscribed sphere radius,                         *
C   *     Si :  inscribed sphere radius,                             *
C   *     WSa:  average Wigner-Seitz radius in units a,              *
C   *     WST:  Voronoi polyhedra radius in units a.                 *
C   *                                                                *
C   ******************************************************************
      USE basis ; USE csts ; USE lattice ; USE message ; USE voronoi
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: csx, csy, csz, dc
      REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: rrp
      REAL(KIND=8), PARAMETER :: tol = 1.d-8, omlim = 0.15d0
      REAL(KIND=8) :: volt, d0, d1, twod1, div, si, sc
      REAL(KIND=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
      REAL(KIND=8) :: sumv, vitr, omeg, xx, yy, zz
      INTEGER,      DIMENSION(:),   ALLOCATABLE :: index, jqn, nnb
      INTEGER,      DIMENSION(:,:), ALLOCATABLE :: iqjq
      INTEGER :: iqq, ivv, ir, iq, jq, nvv, i, indx, ivn, ntr
      INTEGER :: itr, itrb, itrc, iscb, iscc, iv, isc, iqp, frst
C
      ALLOCATE(csx(nr),csy(nr),csz(nr),dc(nr))
      ALLOCATE(jqn(nr),index(nr))
      ALLOCATE(xn(nr),yn(nr),zn(nr),dn(nr))
      ALLOCATE(xsc(3,nr),dsc(nr))
      ALLOCATE(nn(nr),inn(nr,nr),iqn(nr))
      ALLOCATE(wi(nq),wc(nq),wst(nq),rrp(nq,nr))
      ALLOCATE(iqjq(nq,nr),nnb(nq))
C
      frst=1
      volt=0.d0
      DO 20 iqq=1,nq
      ivv=0
      DO 21 ir=1,nr
      jq=jqbas(ir)
      iq=iqbas(ir)
      IF(iq.EQ.iqq) THEN
         ivv=ivv+1
         csx(ivv)=tx(ir)+qx(jq)-qx(iqq)
         csy(ivv)=ty(ir)+qy(jq)-qy(iqq)
         csz(ivv)=tz(ir)+qz(jq)-qz(iqq)
         dc(ivv)=
     .   DSQRT(csx(ivv)*csx(ivv)+csy(ivv)*csy(ivv)+csz(ivv)*csz(ivv))
         jqn(ivv)=jq
      ENDIF
  21  CONTINUE
      nvv=ivv
C
      CALL qsort(dc,index,nvv)
C
      d0=dc(1)
      IF(d0.GT.tol) THEN
         WRITE(m6,100) d0
         STOP
      ENDIF
      d1=dc(2)
      twod1=10.d0*d1
C
      nvn=0
      DO 22 i=2,nvv
      indx=index(i)
      div=dc(i)
      IF(div.LE.twod1) THEN
         nvn=nvn+1
         dn(nvn)=div/2.d0
         xn(nvn)=csx(indx)/2.d0
         yn(nvn)=csy(indx)/2.d0
         zn(nvn)=csz(indx)/2.d0
         iqn(nvn)=jqn(indx)
      ENDIF
   22 CONTINUE
C
C     Set up the Voronoi polyhedra
C
      CALL wscell(si,sc,nr)
C
C     Test the volume of the tetrahedra and the Wigner-Seitz cell
C
      x1=0.d0
      y1=0.d0
      z1=0.d0
      sumv=0.d0
      nnb(iqq)=nvn
      DO 23 ivn=1,nvn
      x2=xn(ivn)
      y2=yn(ivn)
      z2=zn(ivn)
      rrp(iqq,ivn)=2.d0*dn(ivn)
      iqjq(iqq,ivn)=iqn(ivn)
      ntr=nn(ivn)
      DO 23 itr=1,ntr
      itrb=itr-1
      IF(itr.EQ.1) itrb=ntr
      iscb=inn(ivn,itrb)
      itrc=itr
      iscc=inn(ivn,itrc)
      x3=xsc(1,iscb)
      y3=xsc(2,iscb)
      z3=xsc(3,iscb)
      x4=xsc(1,iscc)
      y4=xsc(2,iscc)
      z4=xsc(3,iscc)
      CALL volumt(xx,yy,zz,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,vitr)
      IF(vitr.LE.-tol) WRITE(m6,101) ivn,itr,vitr
   23 sumv=sumv+vitr
      volt=volt+sumv
      IF(sumv.LE.0.d0) THEN
         WRITE(m6,102) sumv
         STOP
      ENDIF
      wst(iqq)=(3.d0*sumv/4.d0/pi)**(1.d0/3.d0)
      wst(iqq)=wst(iqq)/ws
      WRITE(m6,103) iqq,nsc,nvn,sumv,vol/nq,si,sc,
     .              wst(iqq)*ws,wst(iqq)
      IF(nprn.EQ.2) THEN
         DO 24 iv=1,nvn
         WRITE(m6,104) iv,2.*xn(iv),2.*yn(iv),2.*zn(iv),2.*dn(iv),nn(iv)
         DO 25 itr=1,nn(iv)
         isc=inn(iv,itr)
   25    WRITE(m6,105) xsc(1:3,isc),dsc(isc)
   24    CONTINUE
      ENDIF
      wi(iqq)=si
      wc(iqq)=sc
   20 CONTINUE
C
      IF(DABS(volt-vol).GT.tol) THEN
         WRITE(m6,110) volt,vol
         IF(msgl.NE.0) WRITE(msgio,110) volt,vol
      ENDIF
C
      DO 30 iq=1,nq
      nvn=nnb(iq)
      DO 31 ivn=1,nvn
      iqp=iqjq(iq,ivn)
      omeg=ws*(wst(iq)+wst(iqp))/rrp(iq,ivn)-1.d0
      IF(omeg.GT.omlim) THEN
         IF(frst.EQ.1) THEN
         WRITE(m6,120)
            frst=2
         ENDIF
         WRITE(m6,121) iq,ivn,iqp,omeg*100.d0
      ENDIF
   31 CONTINUE
   30 CONTINUE
C
      DEALLOCATE(xn,yn,zn,dn,xsc,dsc,nn,inn,iqn)
      DEALLOCATE(csx,csy,csz,dc,index,rrp,jqn,iqjq,nnb)
C
  100 FORMAT(/,' BLATTS:** D0 =',e15.6,' greater than 0')
  101 FORMAT(/,' Overhanging tetrahedra IVN =',i3,' ITR =',i3,
     .         ' Vol =',f10.6)
  102 FORMAT(/,' BLATTS:** The volume of a tetrahedra is',
     .         ' negative SUMV =',f10.6,/)
  103 FORMAT(/,' BLATTS:   IQ  =',i3,' NSC =',i3,' NVN =',i3,
     .       ' V(t)=',f9.6,' V(a)=',f9.6,//,11x,'Si =',f9.6,
     .       ' Sc =',f9.6,' WST=',f9.6,' WST/WS=',f9.6)
  104 FORMAT(/,5x,i4,4f12.8,i4,/)
  105 FORMAT(9x,4f12.8)
  110 FORMAT(/,' BLATTS:** The volume of the tetrahedra =',f10.6,
     .       /,'           does not add up to the WS volume =',f10.6)
  120 FORMAT(/,' WARNING:  ASA overlap may be too large')
  121 FORMAT(/,' WARNING:  IQ = ',i3,' IVN = ',i3,' IQP = ',i3,
     .         ' omega =',f6.2,'%')
      RETURN
      END
      SUBROUTINE volumt(x,y,z,x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,v)
C   ******************************************************************
C   *                                                                *
C   *    Calculate the volume and centre of gravity of a tetrahedron *
C   *    spanned by vectors:                                         *
C   *                                                                *
C   *        (X1,Y1,Z1);(X2,Y2,Z2);(X3,Y3,Z3);(X4,Y4,Z4).            *
C   *                                                                *
C   *    The volume is negative when (X1,Y1,Z1) falls out of the     *
C   *    triangle defined by (X2,Y2,Z2);(X3,Y3,Z3);(X4,Y4,Z4).       *
C   *                                                                *
C   ******************************************************************
      IMPLICIT NONE
      REAL(KIND=8) :: x, y, z, v, p1, p2, p3
      REAL(KIND=8) :: x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4
      REAL(KIND=8) :: a2, ab, b2, m1, m2, m3, aa, bb
C
      x=0.25d0*(x1+x2+x3+x4)
      y=0.25d0*(y1+y2+y3+y4)
      z=0.25d0*(z1+z2+z3+z4)
C
      a2=(x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2) 
      b2=(x4-x2)*(x4-x2)+(y4-y2)*(y4-y2)+(z4-z2)*(z4-z2) 
      ab=(x3-x2)*(x4-x2)+(y3-y2)*(y4-y2)+(z3-z2)*(z4-z2)
      aa=(b2+ab)/(a2+b2)
      bb=(a2-ab)/(a2+b2)
      m1=aa*(x3-x2)+bb*(x4-x2)
      m2=aa*(y3-y2)+bb*(y4-y2)
      m3=aa*(z3-z2)+bb*(z4-z2)
C
      CALL cross(m1,m2,m3,x4-x3,y4-y3,z4-z3,p1,p2,p3)
      v=(x2-x1)*p1+(y2-y1)*p2+(z2-z1)*p3
      v=v/6.d0
C
      RETURN
      END
