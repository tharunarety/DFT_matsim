      SUBROUTINE grnfcd(efg,nzm1,nzm2,gah,gi,zm,wgm,nzm)
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
      USE atomicdens ; USE control_data ; USE control_text
      USE density ; USE logderivative ; USE message
      USE partialwaves ; USE realgaunt ; USE radialmesh
      USE slope
      IMPLICIT NONE
      INTEGER         :: nzm, nzm1, nzm2
      INTEGER         :: is, iq, it, ita, ir, ign, j, lz
      INTEGER         :: lp, lmp, mp, l, lm, m, lmpp, jrc, i
      COMPLEX(KIND=8), DIMENSION(nzm,nq,ns,nlmf,nlmf) :: gah
      COMPLEX(KIND=8), DIMENSION(nzm,mnta,nq,nlm,nlm,ns) :: gi
      COMPLEX(KIND=8), DIMENSION(nzm)       :: zm, wgm
      COMPLEX(KIND=8), DIMENSION(nzm1:nzm2) :: path
      COMPLEX(KIND=8) :: path0, grr
      REAL(KIND=8)    :: efg, gnti
C
      WRITE(m6,100)
      IF(msgl.NE.0) WRITE(msgio,100)
C
      DO 20 is=1,ns
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      jrc=jwsc(ita,it)
C
      DO 21 lp=0,lmaxf
      DO 21 mp=-lp,lp
      lmp=lp*lp+lp+mp+1
      DO 21 l=0,lmaxf
      DO 21 m=-l,l
      lm=l*l+l+m+1
C
      IF(lm.LE.nlm.AND.lmp.LE.nlm) THEN
         path(nzm1:nzm2)=gi(nzm1:nzm2,ita,iq,lmp,lm,is)
      ELSE
         path(nzm1:nzm2)=gah(nzm1:nzm2,iq,is,lmp,lm)
      ENDIF
      IF(lm.EQ.lmp.AND.l.LE.lmax) THEN
         DO 22 lz=nzm1,nzm2
C
C        g+1/D
C
         path(lz)=path(lz)+zone/dfi(lz,ita,it,l,is,1)
   22    CONTINUE
      ENDIF
      path(nzm1:nzm2)=wgm(nzm1:nzm2)*path(nzm1:nzm2)
C
      DO 24 lz=nzm1,nzm2
      DO 24 ir=1,jrc
      grr=cf(ir,lp,ita,it,is,lz)*cf(ir,l,ita,it,is,lz)*path(lz)
C
C     sum_i -1/(z-e_i)
C
      IF(lm.EQ.lmp.AND.l.LE.lmax) THEN
         DO 23 i=1,necr(l,ita,it,is)
         path0=nocr(l,ita,it,is,i)/(zm(lz)-ecr(l,ita,it,is,i))
   23    grr=grr-wgm(lz)*path0*cfr(ir,l,ita,it,is,i)
      ENDIF
C
      DO 25 j=1,ngnt(lmp,lm)
      ign=ignt(lmp,lm,j)
      gnti=gnt(ign)
      lmpp=lmppg(ign)
      chdl(ita,iq,is,lmpp,ir)=chdl(ita,iq,is,lmpp,ir)+gnti*AIMAG(grr)
   25 CONTINUE
   24 CONTINUE
C
   21 CONTINUE
   20 CONTINUE
C
  100 FORMAT(/,' GRNFCD: Construct the full charge density')
      RETURN
      END
