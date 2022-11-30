      SUBROUTINE pathop0(is,g0,gi,bg0,bgdos,dtil,nzm,nzlin)
C   ******************************************************************
C   *                                                                *
C   * 1. Calculate the Green's functions for the alloy components:   *
C   *                                                                *
C   *     gi(z) = g0(z) + g0(z) * [Di(z) - D0(z)] * gi(z).           *
C   *                                                                *
C   * 2. Find the new coherent potential function D(z):              *
C   *                                                                *
C   *  D(z) = D0(z) + [sum(sort) c(sort) gi(z)(sort)]^-1 - g0(z)^-1. *
C   *                                                                *
C   * 3. Solve the approximate CPA equation and update coherent      *
C   *    Green's functions:                                          *
C   *                                                                *
C   *      g(z) = g0(z) + g0(z) * [D(z) - D0(z)] * g(z).             *
C   *                                                                *
C   *      or:                                                       *
C   *                                                                *
C   *      g(z) = sum(sort) c(sort) gi(z)(sort).                     *
C   *                                                                *
C   *   Note: D and D0 are matrices of (l'm';lm).                    *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE control_data
      USE control_text
      USE kinkmatrix 
      USE logderivative
      USE partialwaves
      IMPLICIT NONE
      INTEGER :: is, nzm, nzlin, lz, iq, it, ita, ntait, l, m, lm
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: g0
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: dtil
      COMPLEX(KIND=8), DIMENSION(nzlin,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,ns) :: bg0
      COMPLEX(KIND=8), DIMENSION(nzlin,nq,nlm) :: bgdos
      COMPLEX(KIND=8), DIMENSION(nlm,nlm) :: ma, mb, mc, md
      COMPLEX(KIND=8) :: dfd
C
C     Eq. 1 Green's functions for the alloy components
C     a) For ordered alloy gi is calculated from g0 by
C        solving the single-site Dyson equation.
C     b) For disordered alloy gi is the alloy component Green's
C        function and g0 is the coherent Green's function
C
      CALL dyson(is,g0,gi,dtil,nzm,nzlin,1)
C
C     Find the new bg0 = g * S^dot
C     Update coherent potential and Green's functions
C
      DO 21 lz=1,nzm
      DO 21 iq=1,nq
      it=itq(iq)
      ntait=nta(it)
C
      IF(ntait.GT.1) THEN
C
C        Eq. 2 New coherent potential function for sublattice iq
C
         ma=zero
         DO 22 ita=1,ntait
C
         ma(1:nlm,1:nlm)=ma(1:nlm,1:nlm)+
     .                   conc(ita,it)*gi(lz,ita,iq,1:nlm,1:nlm,is)
   22    CONTINUE
C
C        [sum(sort) c(sort) gi(z)(sort)]^-1
C
         CALL mtrxinv(nlm,ma,unil,mb,info)
         IF(info.NE.0) STOP 'PATHOP0-1: INFO <> 0'
C
C        g0(z)^-1
C
         mc(1:nlm,1:nlm)=g0(lz,iq,1:nlm,1:nlm,is)
         CALL mtrxinv(nlm,mc,unil,md,info)
         IF(info.NE.0) STOP 'PATHOP0-2: INFO <> 0'
C
         mc(1:nlm,1:nlm)=mb(1:nlm,1:nlm)-md(1:nlm,1:nlm)
         dtil(lz,iq,1:nlm,1:nlm,is)=dtil(lz,iq,1:nlm,1:nlm,is)-
     .                              mc(1:nlm,1:nlm)
C
C        Update the coherent Green's functions g0
C
         mb(1:nlm,1:nlm)=g0(lz,iq,1:nlm,1:nlm,is)
         md=MATMUL(mb,mc)
         md=md+unil
         CALL mtrxinv(nlm,md,mb,mc,info)
         IF(info.NE.0) STOP 'PATHOP0-3: INFO <> 0'
C
         g0(lz,iq,1:nlm,1:nlm,is)=mc(1:nlm,1:nlm)
C
C        Update bg0
C
         DO 23 lm=1,nlm
         dfd=mc(lm,lm)/mb(lm,lm)
         bg0(lz,iq,lm,is)=bg0(lz,iq,lm,is)*dfd
   23    bgdos(lz,iq,lm)=bg0(lz,iq,lm,is)
C
      ELSE
C
C        For ordered case bg0 and g0 are not updated
C
         DO 24 lm=1,nlm
         dfd=gi(lz,1,iq,lm,lm,is)/g0(lz,iq,lm,lm,is)
   24    bgdos(lz,iq,lm)=bg0(lz,iq,lm,is)*dfd
      ENDIF
C
   21 CONTINUE
C
      RETURN
      END
