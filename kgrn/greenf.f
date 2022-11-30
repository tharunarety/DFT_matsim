      SUBROUTINE greenf(lin,nzm1,nzm2,gi,hgh,zm,wgm,nzm,nzlin)
C   ******************************************************************
C   *                                                                *
C   * Calculate the imaginar part of the r-diagonal Green's          *
C   * function :                                                     *
C   *                                                                *
C   * Imag [ phi^a_Rl'(r,z) * g^a_RL'RL(z) * phi^a_Rl(r,z) ]/pi      *
C   *                                                                *
C   * where g^a_RL'RL(z) is the k-integrated path operator,          *
C   *                                                                *
C   * and phi^a_Rl(r,z) = phi_Rl(r,z)/varphi_Rl(a,z)                 *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data  ; USE control_text
      USE density ; USE logderivative ; USE partialwaves
      USE radialmesh ; USE slope
      IMPLICIT NONE
      INTEGER         :: lin, nzm, nzm1, nzm2, nzlin
      INTEGER         :: is, iq, it, ita, ir, twol
      INTEGER         :: lp, mp, lmp, l, m, lm, lpp, jrn, lz, i
      COMPLEX(KIND=8), DIMENSION(nzlin,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(nzm,nq,ns,lmax+1:lmaxt) :: hgh
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm, wgm
      COMPLEX(KIND=8) :: path, path0, grr, rad
C
      DO 20 is=1,ns
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      jrn=jrsm(ita,it)
C
      DO 21 l=0,lmax
C
      DO 22 m=-l,l
      lm=l*l+l+m+1
      DO 22 lz=nzm1,nzm2
C
C     g+1/D
C
      path=gi(lz,ita,iq,lm,lm,is)+zone/dfi(lz,ita,it,l,is,1)
      path=wgm(lz)*path
C
      DO 24 ir=1,jrn
      rad=cf(ir,l,ita,it,is,lz)
      grr=rad*rad*path
C
C     sum_i -1/(z-e_i)
C
      DO 23 i=1,necr(l,ita,it,is)
      path0=nocr(l,ita,it,is,i)/(zm(lz)-ecr(l,ita,it,is,i))
   23 grr=grr-wgm(lz)*path0*cfr(ir,l,ita,it,is,i)
C
   24 chdn(ita,it,is,ir)=chdn(ita,it,is,ir)+AIMAG(grr)
   22 CONTINUE
C
   21 CONTINUE
C
C     Higher tails
C
      IF(lin.EQ.0) THEN
         DO 26 l=lmax+1,lmaxt
         DO 26 lz=nzm1,nzm2
         path=wgm(lz)*hgh(lz,iq,is,l)
         DO 26 ir=1,jrn
         rad=cf(ir,l,ita,it,is,lz)
         grr=rad*rad*path
         chdh(ita,it,is,ir)=chdh(ita,it,is,ir)+AIMAG(grr)
   26    CONTINUE
      ENDIF
C
C     Calculate the r-integrated off-diagonal Green's function 
C
      IF(func.NE.'ASA'.AND.lin.LE.1) THEN
         DO 30 lp=0,lmax
         DO 30 mp=-lp,lp
         lmp=lp*lp+lp+mp+1
         DO 30 l=0,lmax
         DO 30 m=-l,l
         lm=l*l+l+m+1
         DO 30 lz=nzm1,nzm2
         path=gi(lz,ita,iq,lmp,lm,is)
         IF(lmp.EQ.lm) THEN
C
C           g+1/D
C
            path=path+zone/dfi(lz,ita,it,l,is,1)
         ENDIF
         path=wgm(lz)*path
C
         DO 30 lpp=0,lmax2
         grr=cfm(lp,l,lpp,ita,it,is,lz)*path
         IF(lmp.EQ.lm) THEN
C
C           sum_i -1/(z-e_i)
C
            DO 32 i=1,necr(l,ita,it,is)
            path0=nocr(l,ita,it,is,i)/(zm(lz)-ecr(l,ita,it,is,i))
            grr=grr-wgm(lz)*path0*cfrm(lp,lpp,ita,it,is,i)
   32       CONTINUE
         ENDIF
         grnfm(iq,is,lmp,lm,lpp)=grnfm(iq,is,lmp,lm,lpp)+
     .                           conc(ita,it)*AIMAG(grr)
   30    CONTINUE

C
      ENDIF
C
   20 CONTINUE
C
      RETURN
      END
