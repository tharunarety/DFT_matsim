      SUBROUTINE screen(ta,ts,zm,nzm)
C   ******************************************************************
C   *                                                                *
C   * Set up the screening matrices for the energies zm.             *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE logderivative
      USE potential
      USE slope
      USE radialmesh
      IMPLICIT NONE
      INTEGER :: nzm, lz, is, it, ita, l
      COMPLEX(KIND=8), DIMENSION(nzm,4,0:lmax,mnta,nt,ns,0:1) :: ts
      COMPLEX(KIND=8), DIMENSION(nzm,4,0:lmax,nt,ns)     :: ta
      COMPLEX(KIND=8), DIMENSION(4,0:lmax,nt,0:1)      :: ka
      COMPLEX(KIND=8), DIMENSION(4,0:lmax,mnta,nt,0:1) :: ks
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8) :: x1, x2, x3, x4, x1d, x2d, x3d, x4d, z
      REAL(KIND=8)    :: w2
C
C     Loop for complex energy
C
      w2=sws*sws
      DO 20 is=1,ns
      DO 20 lz=1,nzm
      z=(zm(lz)-vmtz(is))*w2
      DO 21 it=1,nt
   21 sgm(0:lmax,1,it)=sigmt(0:lmax,it)/sws
      CALL trmtrz(z,sgm,ka,itrans,1,0)
      DO 22 it=1,nt
      DO 22 ita=1,nta(it)
      sgm(0:lmax,ita,it)=hsr(ita,it)/sws
      IF(localmt(ita,it).EQ.1) THEN
         sgm(0:lmax,ita,it)=wsm(ita,it)/sws
      ENDIF
   22 CONTINUE
      CALL trmtrz(z,sgm,ks,itrans,mnta,0)
      ka(1:4,0:lmax,1:nt,1)=w2*ka(1:4,0:lmax,1:nt,1)
      DO 23 it=1,nt
      DO 23 ita=1,nta(it)
   23 ks(1:4,0:lmax,ita,it,1)=w2*ks(1:4,0:lmax,ita,it,1)
C
      DO 20 it=1,nt
      ta(lz,1:4,0:lmax,it,is)=ka(1:4,0:lmax,it,0)
      DO 20 ita=1,nta(it)
      DO 20 l=0,lmax
C
      x1= ka(1,l,it,0)*ks(4,l,ita,it,0)-ka(2,l,it,0)*ks(3,l,ita,it,0)
      x2=-ka(1,l,it,0)*ks(2,l,ita,it,0)+ka(2,l,it,0)*ks(1,l,ita,it,0)
      x3= ka(3,l,it,0)*ks(4,l,ita,it,0)-ka(4,l,it,0)*ks(3,l,ita,it,0)
      x4=-ka(3,l,it,0)*ks(2,l,ita,it,0)+ka(4,l,it,0)*ks(1,l,ita,it,0)
      ts(lz,1,l,ita,it,is,0)=x1
      ts(lz,2,l,ita,it,is,0)=x2
      ts(lz,3,l,ita,it,is,0)=x3
      ts(lz,4,l,ita,it,is,0)=x4
C
      x1d= ka(1,l,it,0)*ks(4,l,ita,it,1)-ka(2,l,it,0)*ks(3,l,ita,it,1)
     .    +ka(1,l,it,1)*ks(4,l,ita,it,0)-ka(2,l,it,1)*ks(3,l,ita,it,0)
      x2d=-ka(1,l,it,0)*ks(2,l,ita,it,1)+ka(2,l,it,0)*ks(1,l,ita,it,1)
     .    -ka(1,l,it,1)*ks(2,l,ita,it,0)+ka(2,l,it,1)*ks(1,l,ita,it,0)
      x3d= ka(3,l,it,0)*ks(4,l,ita,it,1)-ka(4,l,it,0)*ks(3,l,ita,it,1)
     .    +ka(3,l,it,1)*ks(4,l,ita,it,0)-ka(4,l,it,1)*ks(3,l,ita,it,0)
      x4d=-ka(3,l,it,0)*ks(2,l,ita,it,1)+ka(4,l,it,0)*ks(1,l,ita,it,1)
     .    -ka(3,l,it,1)*ks(2,l,ita,it,0)+ka(4,l,it,1)*ks(1,l,ita,it,0)
      ts(lz,1,l,ita,it,is,1)=x1d
      ts(lz,2,l,ita,it,is,1)=x2d
      ts(lz,3,l,ita,it,is,1)=x3d
      ts(lz,4,l,ita,it,is,1)=x4d
C
   20 CONTINUE
C
      RETURN
      END
