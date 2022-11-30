      SUBROUTINE kkrfcd(efg,tfm)
C   ******************************************************************
C   *                                                                *
C   * Solve the kink equation and construct the total charge density.*
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE control_data ; USE control_text ; USE csts
      USE density      ; USE dosmom       ; USE energymesh
      USE force        ; USE greenfunc    ; USE logderivative
      USE message      ; USE realgaunt    ; USE moments
      USE radialmesh   ; USE potparam     ; USE slope
      IMPLICIT NONE
      INTEGER, PARAMETER :: prnt = 0
      REAL(KIND=8) :: efg, dn
      INTEGER      :: tfm, is, it, ita, l, ir
      INTEGER      :: lz, iq, lm, lmp, lmaxf2
C
      IF(tfm.EQ.1) THEN
C
C        Store the full charge density chdl in chd0
C
         ALLOCATE(chd0(mnta,nq,ns,nlmh,dimr))
         chd0=chdl
C
C        Do not compute SCA density if lmaxh=lmax
C
         IF(lmaxh.EQ.lmax) RETURN
         lmaxf=lmax
         nlmf=(lmaxf+1)*(lmaxf+1)
         nlmqf=nlmf*nq
         DEALLOCATE(chdl,qlmn)
         IF(func.NE.'ASA') THEN
            DEALLOCATE(gnt,lmpg,lmg,lmppg,lppg)
            DEALLOCATE(ngnt,ignt)
         ENDIF
      ELSE
         lmaxf=lmaxh
         nlmf=(lmaxf+1)*(lmaxf+1)
         nlmqf=nlmf*nq
      ENDIF
C
      ALLOCATE(rp(dimr),rq(dimr))
      ALLOCATE(gah(nzm,nq,ns,nlmf,nlmf))
      ALLOCATE(bg0(nzm,nq,nlm,ns),dtilz(nzm,nq,nlm,nlm,ns))
      ALLOCATE(slop(nlmqf,nlmqf),tmp(nlmqf,nlmq))
      CALL kkrarr(nzm,lmaxf,1)
      gah=zero
C
C     1. Solve kink equation for the complex contour ZM
C
      CALL setexp(zm,nzm,nz2)
      CALL expans(zm,nzm)
C
C     Solve the Dirac equation and make the logarithmic derivatives
C
      CALL logder(zm,nzm)
C
C     Set up the logarithmic derivative at s^m if requisted
C
      IF(lclmff.NE.0) CALL lclmtz(tmz,zm,nzm)
C
C     Integrate the single-site Green's function on the real axis
C
      CALL singlg(efg,prnt)
C
C     Make the logarithmic derivative of the backward extrapolated
C     free electron solution : fi^a = f^a + g^a * D{fi0^a(a)}/(-a*d)
C
      CALL screen(taz,tsz,zm,nzm)
      CALL logder0(taz,tsz,tmz,zm,nzm)
C
C     Initial guess for the coherent potential function
C
      CALL logdcpa(dtilz,nzm)
C
C     Make the higher basis functions
C
      CALL logderh(zm,nzm,lmaxf)
C
C     Find the k-integrated Green's function 1/(D - S) on the contour,
C     and construct the g*S, S*g*S matrices
C
      dn=0.d0
      DO 20 is=1,ns
      CALL fcdpth(is,gah,bg0,gi,dtilz,nzm)
C
      CALL doscpa(is,bg0,gi,dos,zm,nzm)
C
C     Find the number of electrons on the contour
C
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
C
      DO 22 l=0,lmax
   22 dos(1:nzm,ita,it,l)=wgm(1:nzm)*dos(1:nzm,ita,it,l)
      DO 23 l=0,lmax
      IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
         tnos0(ita,it,l,is)=SUM(AIMAG(dos(1:nz2,ita,it,l)))
      ENDIF
   23 tnos(ita,it,l,is)=SUM(AIMAG(dos(1:nzm,ita,it,l)))
C
      DO 24 l=0,lmax
   24 dos(1:nzm,ita,it,l)=zm(1:nzm)*dos(1:nzm,ita,it,l)
      DO 25 l=0,lmax
      IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
         emom0(ita,it,l,is)=SUM(AIMAG(dos(1:nz2,ita,it,l)))
      ENDIF
   25 emom(ita,it,l,is)=SUM(AIMAG(dos(1:nzm,ita,it,l)))
      dn=dn+conc(ita,it)*mmt(it)*SUM(tnos(ita,it,0:lmax,is))
   21 CONTINUE
   20 CONTINUE
      tnos=spinfc*tnos/pi
      emom=spinfc*emom/pi
      dn=spinfc*dn/pi
      WRITE(m6,100) dn,elt
      IF(msgl.NE.0) WRITE(msgio,100) dn,elt
      IF(zmsh.EQ.'M'.OR.zmsh.EQ.'m'.OR.zmsh.EQ.'f') THEN
         tnos0=spinfc*tnos0/pi
         emom0=spinfc*emom0/pi
      ENDIF
      IF(prnt.EQ.1) CALL prinfo('KKRFCD')
C
C     Calculate the imaginar part of the r-diagonal Green's function 
C
      DEALLOCATE(slop,tmp,bg0,dtilz)
      IF(tfm.EQ.1) THEN
         DEALLOCATE(salpl)
         IF(expan.EQ.'D') DEALLOCATE(salplp)
      ENDIF
C
      IF(tfm.EQ.1) THEN
         lmaxf2=2*lmaxf
         nlmf=(lmaxf2+1)*(lmaxf2+1)
         CALL mgaunt(lmaxf,lmaxf,lmaxf2)
      ELSE
         CALL mgaunt(lmaxf,lmaxf,lmaxf)
      ENDIF
      CALL ordgnt(nlmf)
C
      ALLOCATE(chdl(mnta,nq,ns,nlmf,dimr))
C
      chdl=0.d0
      CALL grnfcd(efg,1,nzm,gah,gi,zm,wgm,nzm)
C
C     Set up the kinetic force components
C
      IF(tfm.EQ.0) THEN
         ALLOCATE(fsy(nq),fsz(nq),fsx(nq))
         fsy=0.d0
         fsz=0.d0
         fsx=0.d0
         IF(frc.EQ.'Y') THEN
            ALLOCATE(fny(nq),fnz(nq),fnx(nq))
            fny=0.d0
            fnz=0.d0
            fnx=0.d0
            CALL kinstr(1,nzm,gah,zm,wgm,nzm)
         ENDIF
      ENDIF
C
C     Add the higher tails contribution to the density
C
      IF(tfm.EQ.1.AND.lmaxt.GT.lmaxf) THEN
         DO 30 is=1,ns
         DO 31 iq=1,nq
         it=itq(iq)
         DO 32 ita=1,nta(it)
         DO 33 ir=1,dimr
         chdl(ita,iq,is,1,ir)=chdl(ita,iq,is,1,ir)+
     .                        chdh(ita,it,is,ir)/sqfpi/mmt(it)
   33    CONTINUE
   32    CONTINUE
   31    CONTINUE
   30    CONTINUE
      ENDIF
C
      DEALLOCATE(rp,rq,gah)
      CALL kkrarr(nzm,lmaxf,0)
C
      IF(tfm.EQ.1) lmaxf=lmaxf2
C
  100 FORMAT(/,' KKRFCD: NOS(Ef) =',f10.6,' ELT =',f10.6)
      RETURN
      END
