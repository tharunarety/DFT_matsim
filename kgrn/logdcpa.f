      SUBROUTINE logdcpa(dtil,nzm)
C   ******************************************************************
C   *                                                                *
C   *    Set up the coherent potential function:                     *
C   *                                                                *
C   *        D^-1 = sum(sort) c(sort) * D^-1(sort).                  *
C   *                                                                *
C   *   For ordered sublattice this function is the logarithmic      *
C   *   derivative (D^a) calculated at the hard sphere (see logder0) *
C   *   and it is (lm,l'm') diagonal.                                *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE control_data
      USE logderivative
      IMPLICIT NONE
      INTEGER :: nzm, is, iq, it, ita, l, m, lm, lz
      COMPLEX(KIND=8), DIMENSION(nzm,nq,nlm,nlm,ns) :: dtil
C
      dtil=zero
      DO 20 is=1,ns
      DO 20 iq=1,nq
      it=itq(iq)
C
      DO 20 l=0,lmax
      DO 20 m=-l,l
      lm=l*l+l+m+1
C
      DO 21 lz=1,nzm
      DO 22 ita=1,nta(it)
      dtil(lz,iq,lm,lm,is)=dtil(lz,iq,lm,lm,is)+
     .                     conc(ita,it)/dfi(lz,ita,it,l,is,1)
   22 CONTINUE
      dtil(lz,iq,lm,lm,is)=zone/dtil(lz,iq,lm,lm,is)
   21 CONTINUE
C
   20 CONTINUE
C
      RETURN
      END
