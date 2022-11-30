      SUBROUTINE logder(zm,nzm)
C   ******************************************************************
C   *                                                                *
C   * Calculate the partial waves for complex energies ZM.           *
C   * Logarithmic derivative of the partial waves:                   *
C   *                                                                *
C   *  [r*fi(z,r)]'' = [l*(l+1)/r^2 + v(r) - z]*r*fi(z,r)            *
C   *                                                                *
C   * defined by:                                                    *
C   *                         fi'(s)                                 *
C   *           D{fi(s)} = s* ------                                 *
C   *                         fi(s)                                  *
C   *                                                                *
C   * Note: the partial waves are multiplied by r.                   *
C   *                                                                *
C   ******************************************************************
      USE diracparam   ; USE control_data ; USE logderivative
      USE partialwaves ; USE potential    ; USE radialmesh
      USE temporary
      IMPLICIT NONE
      INTEGER :: nzm, is, it, ita, il, l, k, lz, ir, jsr
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8) :: z, pu, qu, cfint, smm, smmq
      REAL(KIND=8)    :: s, s2, vu, smm1, smm2, beta
C
      DO 20 is=1,ns
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      s=hsr(ita,it)
      s2=s*s
      jsr=jsrs(ita,it)
C
      vu=v(jsr,ita,it,is)/s2
      DO 20 l=0,lmax
      il=l+1
      k=-l-1
C
C     Solve the Dirac equation for complex energies
C
      DO 21 lz=1,nzm
      z=zm(lz)
C
      CALL cdirac(z,k,it,ita,il,is,0,beta)
C
C     Logarithmic derivative at potential sphere for z
C
      pu=p(jsr)
      qu=q(jsr)
      dfi(lz,ita,it,l,is,1)=s*(1.d0-(vu-z)/csq)*qu/pu+l
C
C     Energy derivative of the logarithmic derivative for z
C
      DO 22 ir=1,jsr
      cfint=p(ir)*p(ir)*ri(ir,ita,it)
      fi(ir)=REAL(cfint,8)
   22 fip(ir)=AIMAG(cfint)
      CALL simpn(fi,dx,jsr,smm1)
      CALL simpn(fip,dx,jsr,smm2)
      smm=CMPLX(smm1,smm2,8)+p(1)*p(1)*ri(1,ita,it)/(1.d0+2.d0*beta)
C
C     Include the scalar relativistic correction to D_dot
C
      IF(srde.EQ.1) THEN
         DO 23 ir=1,jsr
         cfint=q(ir)*q(ir)*ri(ir,ita,it)
         fi(ir)=REAL(cfint,8)
   23    fip(ir)=AIMAG(cfint)
         CALL simpn(fi,dx,jsr,smm1)
         CALL simpn(fip,dx,jsr,smm2)
         smmq=CMPLX(smm1,smm2,8)+q(1)*q(1)*ri(1,ita,it)/(1.d0+2.d0*beta)
         dfi(lz,ita,it,l,is,2)=-s*(smm+smmq/csq)/p(jsr)/p(jsr)
      ELSE
         dfi(lz,ita,it,l,is,2)=-s*smm/p(jsr)/p(jsr)
      ENDIF
C
C     Save radial solutions (need not be normalized !!)
C
      cf(1:jsr,l,ita,it,is,lz)=p(1:jsr)
C
   21 CONTINUE
   20 CONTINUE
C
      RETURN
      END
