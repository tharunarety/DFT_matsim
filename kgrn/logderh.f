      SUBROUTINE logderh(zm,nzm,lmaxact)
C   ******************************************************************
C   *                                                                *
C   * Calculate the free electron solutions for complex energies ZM  *
C   * from the higher l from a_H = 0 to Sc.                          *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE message
      USE partialwaves
      USE potential
      USE radialmesh
      IMPLICIT NONE
      INTEGER         :: lmaxact, nzm, lz, it, ita, l, jrn, is, ir
      COMPLEX(KIND=8), DIMENSION(0:lmaxact) :: jl, hl
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8) :: kap2
      REAL(KIND=8)    :: r, rw, rw0
C
C     Loop for complex energy
C
      DO 20 is=1,ns
      DO 20 lz=1,nzm
      kap2=zm(lz)-vmtz(is)
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jrn=jrsm(ita,it)
      IF(fcd.EQ.'Y') jrn=MAX0(jrn,jwsc(ita,it))
C
C     Free electron solutions
C
      DO 21 ir=1,jrn
      r=ri(ir,ita,it)
C
      CALL besslz(kap2*r*r,0,lmaxact,jl(0),hl(0))
      rw0=r/sws
      rw=rw0**lmax
      DO 21 l=lmax+1,lmaxact
      rw=rw*rw0
   21 cf(ir,l,ita,it,is,lz)=-jl(l)*rw*r
C
   20 CONTINUE
C
      RETURN
      END
