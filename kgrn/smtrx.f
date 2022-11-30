      SUBROUTINE smtrx(nlmact,pkx,pky,pkz)
C   ******************************************************************
C   *                                                                *
C   *    Evaluate slope matrices for a crystal by means              *
C   *    of the Bloch sum.                                           *
C   *                                                                *
C   *   *On entry "PKX,PKY,PKZ" contain the rectangular coordinates  *
C   *    of the Bloch vector in units of 1/A and "SALPL" contain     *
C   *    the transfer matrices to be summed over the translations    *
C   *    "TX,TY,TZ".                                                 *
C   *                                                                *
C   *   *On exit "SAA" contain the complex structure matrix :        *
C   *                                                                *
C   *    S(K;QP,LP;Q,L) = sum(T) EXP(I*K*T)*S(QP,LP;Q+T,L)           *
C   *                                                                *
C   *    where "Q" is the position of the atom in the primitive      *
C   *    cell and "T" is a translation vector.                       *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE control_text ; USE lattice
      USE message      ; USE radialmesh   ; USE slope
      USE taylor
      IMPLICIT NONE
      COMPLEX(KIND=8) :: phs
      REAL(KIND=8)    :: pkx, pky, pkz, dot
      INTEGER         :: jl, jq, il, iq, iv1, iv2, iv
      INTEGER         :: nlmact, lm, lmp, jd
C
      sa=zero
      IF(expan.EQ.'D') sap=zero
C
C     Bloch sum
C
      jl=-nlm
      DO 20 jq=1,nq
      jl=jl+nlm
      DO 20 iq=1,nq
      il=(iq-1)*nlmact
      iv1=nviq(iq-1)+1
      iv2=nviq(iq)
      DO 20 iv=iv1,iv2
      IF(jqbas(iv).NE.jq) GO TO 20
C
      dot=pkx*tx(iv)+pky*ty(iv)+pkz*tz(iv)
      phs=CMPLX(COS(dot),SIN(dot),8)
      DO 21 jd=0,nder
      DO 21 lm=1,nlm
      DO 21 lmp=1,nlmact
      sa(il+lmp,jl+lm,jd)=sa(il+lmp,jl+lm,jd)+
     .                    phs*salpl(lmp,lm,iv,jd)
   21 CONTINUE
      IF(expan.EQ.'D') THEN
         DO 22 jd=0,nderp
         DO 22 lm=1,nlm
         DO 22 lmp=1,nlmact
         sap(il+lmp,jl+lm,jd)=sap(il+lmp,jl+lm,jd)+
     .                        phs*salplp(lmp,lm,iv,jd)
   22    CONTINUE
      ENDIF
   20 CONTINUE
C
      RETURN
      END
