      SUBROUTINE logdrr0(e,il,it,ita,is,dr,drd,rsl,zdzd,zdz)
C   ******************************************************************
C   *                                                                *
C   * Calculate the logarithmic derivative of the backward           *
C   * extrapolated free electron solution for real energy e:         *
C   *                                                                *
C   *   fi^a = f^a - g^a * D{fi0^a(a)}/(a*d),                        *
C   *                                                                *
C   * where :                                                        *
C   *                                                                *
C   *   f^a = t^a(1)*n + t^a(2)*j and g^a = - t^a(3)*n - t^a(4)*j    *
C   *                                                                *
C   * t^a is the screening matrix and                                *
C   *                                                                *
C   *                      D{fi(s)} - D{f^a(s)}    f^a(s)            *
C   * D{fi0^a(a)} = a*d * ---------------------- * ------            *
C   *                      D{fi(s)} - D{g^a(s)}    g^a(s)            *
C   *                                                                *
C   * Renormalize the partial wave : phi^a(r) = phi(r)/fi(a) .       *
C   *                                                                *
C   * Calculate the free electron solutions for real energy e        *
C   * from the potential sphere S to the circumscribed sphere Sc.    *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE diracparam
      USE lattice
      USE potential
      USE potparam
      USE radialmesh
      USE slope
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(0:lmax) :: jl, hl
      REAL(KIND=8), DIMENSION(4,0:1)  :: ta, ts, tm
      REAL(KIND=8) :: z, drs, dr, sgma, kap2, chi
      REAL(KIND=8) :: x1, x2, x3, x4, xn, xs, rw, bs, nm, ga, fa
      REAL(KIND=8) :: fm, gm, dm
      REAL(KIND=8) :: drsd, drd, x1d, x2d, x3d, x4d, zdzd, zdz
      REAL(KIND=8) :: e, s, w2, beta, vu, r, fi0a, fint, fintq
      INTEGER      :: rsl, jsr, jrn, jsi, ir, it, ita,il, l, k, is
C
      w2=sws*sws
      s=hsr(ita,it)
      jsr=jsrs(ita,it)
      jrn=jrsm(ita,it)
      IF(fcd.EQ.'Y') jrn=MAX0(jrn,jwsc(ita,it))
      vu=v(jsr,ita,it,is)/s/s
C
      l=il-1
      k=-l-1
      CALL dirac(e,k,it,ita,il,is,jrn,0,beta,0)
C
C     Logarithmic derivative and its energy derivative for e
C
      drs=s*(1.d0-(vu-e)/csq)*rq(jsr)/rp(jsr)+l
C
      fi(1:jsr)=rp(1:jsr)*rp(1:jsr)*ri(1:jsr,ita,it)
      CALL simpn(fi,dx,jsr,fint)
      fint=fint+rp(1)*rp(1)*ri(1,ita,it)/(1.d0+2.d0*beta)
C
C     Include the scalar relativistic correction to D_dot
C
      IF(srde.EQ.1) THEN
         fi(1:jsr)=rq(1:jsr)*rq(1:jsr)*ri(1:jsr,ita,it)
         CALL simpn(fi,dx,jsr,fintq)
         fintq=fintq+rq(1)*rq(1)*ri(1,ita,it)/(1.d0+2.d0*beta)
         drsd=-s*(fint+fintq/csq)/rp(jsr)/rp(jsr)
      ELSE
         drsd=-s*fint/rp(jsr)/rp(jsr)
      ENDIF
C
C     Set up the logarithmic derivative at s^m if requisted
C
      chi=rp(jsr)
      IF(localmt(ita,it).EQ.1) CALL 
     .                         lclmtx(drs,drsd,chi,e,l,ita,it,is,tm)
      dm=drs
C
C     Makes the screening matrices t^a and t^s and energy derivatives
C
      z=(e-vmtz(is))*w2
      sgma=sigmt(l,it)/sws
      CALL trmtrx(z,sgma,ta,l,itrans,0)
      sgma=s/sws
      IF(localmt(ita,it).EQ.1) sgma=wsm(ita,it)/sws
      CALL trmtrx(z,sgma,ts,l,itrans,0)
      ta(1:4,1)=w2*ta(1:4,1)
      ts(1:4,1)=w2*ts(1:4,1)
C
C     Makes logarithmic derivatives at the hard sphere
C
      x1= ta(1,0)*ts(4,0)-ta(2,0)*ts(3,0)
      x2=-ta(1,0)*ts(2,0)+ta(2,0)*ts(1,0)
      x3= ta(3,0)*ts(4,0)-ta(4,0)*ts(3,0)
      x4=-ta(3,0)*ts(2,0)+ta(4,0)*ts(1,0)
C
      xs=drs*x1+x2
      xn=drs*x3+x4
      fi0a=-0.5*sws/(sigmt(l,it)*s)*chi*xn
      dr=xs/xn
C
      IF(rsl.EQ.0) RETURN
C
      x1d= ta(1,0)*ts(4,1)-ta(2,0)*ts(3,1)
     .    +ta(1,1)*ts(4,0)-ta(2,1)*ts(3,0)
      x2d=-ta(1,0)*ts(2,1)+ta(2,0)*ts(1,1)
     .    -ta(1,1)*ts(2,0)+ta(2,1)*ts(1,0)
      x3d= ta(3,0)*ts(4,1)-ta(4,0)*ts(3,1)
     .    +ta(3,1)*ts(4,0)-ta(4,1)*ts(3,0)
      x4d=-ta(3,0)*ts(2,1)+ta(4,0)*ts(1,1)
     .    -ta(3,1)*ts(2,0)+ta(4,1)*ts(1,0)
      drd=(drsd*x1+drs*x1d+x2d)/xn-
     . dr*(drsd*x3+drs*x3d+x4d)/xn
      drd=sigmt(l,it)*drd
C
C     Free electron solutions for r > S,...,RI(JRN)
C
      rp(1:jsr)=rp(1:jsr)/fi0a
      DO 20 ir=jsr+1,jrn
      r=ri(ir,ita,it)
      rw=r/sws
C
      IF(localmt(ita,it).EQ.1.AND.r.LE.wsm(ita,it)) THEN
         kap2=e-vmtzr(it,is)
         CALL bessl(kap2*r*r,0,l,jl(0),hl(0))
         bs=jl(l)*rw**l
         nm=hl(l)*rw**(-l-1)
         fm= tm(1,0)*nm+tm(2,0)*bs
         gm=-tm(3,0)*nm-tm(4,0)*bs
         rp(ir)=(fm+gm*dm)*r*chi/fi0a/s
      ELSE
         kap2=e-vmtz(is)
         CALL bessl(kap2*r*r,0,l,jl(0),hl(0))
         bs=jl(l)*rw**l
         nm=hl(l)*rw**(-l-1)
         fa= ta(1,0)*nm+ta(2,0)*bs
         ga=-ta(3,0)*nm-ta(4,0)*bs
         rp(ir)=(fa+ga*dr)*r
      ENDIF
   20 CONTINUE
      IF(fcd.EQ.'Y') THEN
         jsi=jwsi(ita,it)-shf
         fi(1:jrn)=rp(1:jrn)/ri(1:jrn,ita,it)
         CALL diffn(fi,rhop,rhopp,jrn,dx)
         zdzd=-rhop(jsi)*rhop(jsi)/drd
         zdz =-rhop(jsi)*fi(jsi)/drd
      ENDIF
C
      RETURN
      END
