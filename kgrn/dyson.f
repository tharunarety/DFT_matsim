      SUBROUTINE dyson(is,g0,g,dtil,nzm,nzlin,dy)
C   ******************************************************************
C   *                                                                *
C   *   Solve Dyson equation for each sort:                          *
C   *                                                                *
C   *      g(z) = g0(z) + g0(z)*(D - D0)*g(z).                       *
C   *                                                                *
C   *   Note: D is a matrix D(l'm';lm) and D0 is a vector D0(lm).    *
C   *                                                                *
C   *   For ordered sublattice D = D0 and, therefore, g = g0.        *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE kinkmatrix 
      USE logderivative
      IMPLICIT NONE
      INTEGER :: is, lz, nzm, nzlin, dy, iq, it, ntait, ita, l, m, lm
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: g0
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: dtil
      COMPLEX(KIND=8), DIMENSION(nzlin,mnta,nq,nlm,nlm,ns) :: g
      COMPLEX(KIND=8), DIMENSION(nlm,nlm) :: ma, mb, mc
      COMPLEX(KIND=8) :: dfd
C
C     Loop for the sites and sorts
C
      DO 20 lz=1,nzm
      DO 20 iq=1,nq
      it=itq(iq)
      ntait=nta(it)
C
      IF(ntait.EQ.1.AND.dy.EQ.0) THEN
C
C        Ordered sublattice
C
         g(lz,1,iq,1:nlm,1:nlm,is)=g0(lz,iq,1:nlm,1:nlm,is)
      ELSE
C
C        Disordered sublattice
C
         mc(1:nlm,1:nlm)=g0(lz,iq,1:nlm,1:nlm,is)
C
         DO 21 ita=1,ntait
C
         mb(1:nlm,1:nlm)=dtil(lz,iq,1:nlm,1:nlm,is)
C
         DO 22 l=0,lmax
         dfd=dfi(lz,ita,it,l,is,1)
         DO 22 m=-l,l
         lm=l*l+l+m+1
   22    mb(lm,lm)=mb(lm,lm)-dfd
C
         ma=MATMUL(mc,mb)
C
         DO 23 lm=1,nlm
   23    ma(lm,lm)=ma(lm,lm)+zone
C
         CALL mtrxinv(nlm,ma,mc,mb,info)
         IF(info.NE.0) STOP 'DYSON: INFO <> 0'
C
C        New path operator
C
         g(lz,ita,iq,1:nlm,1:nlm,is)=mb(1:nlm,1:nlm)
C
   21    CONTINUE
      ENDIF
   20 CONTINUE
C
      RETURN
      END
