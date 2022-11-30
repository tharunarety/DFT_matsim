      SUBROUTINE logder0(ta,ts,tm,zm,nzm)
C   ******************************************************************
C   *                                                                *
C   * Calculate the logarithmic derivative of the backward           *
C   * extrapolated free electron solution for complex energies ZM:   *
C   *                                                                *
C   *   fi^a = f^a + g^a * D{fi0^a(a)},                              *
C   *                                                                *
C   * where :                                                        *
C   *                                                                *
C   *   f^a = t^a(1)*n + t^a(2)*j and g^a = - t^a(3)*n - t^a(4)*j    *
C   *                                                                *
C   * t^a is the screening matrix and                                *
C   *                                                                *
C   *                      D{fi(s)} - D{f^a(s)}    f^a(s)            *
C   *   D{fi0^a(a)} =  -  ---------------------- * ------            *
C   *                      D{fi(s)} - D{g^a(s)}    g^a(s)            *
C   *                                                                *
C   * Renormalize the partial wave : phi^a(r) = phi(r)/fi(a) .       *
C   *                                                                *
C   * Calculate the free electron solutions for complex energies ZM  *
C   * from the potential sphere S to the circumscribed sphere Sc.    *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE logderivative
      USE message
      USE partialwaves
      USE potential
      USE radialmesh
      USE temporary
      IMPLICIT NONE
      INTEGER :: nzm, lz, iq, it, ita, l, lp, lpp, is, i, m, ir
      INTEGER :: jrn, jsr, jws
      COMPLEX(KIND=8), DIMENSION(nzm,4,0:lmax,mnta,nt,ns,0:1) :: ts
      COMPLEX(KIND=8), DIMENSION(nzm,4,0:lmax,nt,ns) :: ta
      COMPLEX(KIND=8), DIMENSION(nzm,4,0:lmax,mnta,nt,ns) :: tm
      COMPLEX(KIND=8), DIMENSION(0:lmax) :: jl, hl
      COMPLEX(KIND=8), DIMENSION(nzm) :: zm
      COMPLEX(KIND=8) :: kap2, fi0a, d, dd, bs, nm, ga, fa, fm, gm, dm
      COMPLEX(KIND=8) :: x1, x2, x3, x4, x1d, x2d, x3d, x4d, xn, xs, chi
      REAL(KIND=8)    :: r, rw, rephi, imphi
C
C     Loop for complex energy
C
      DO 20 is=1,ns
      DO 20 lz=1,nzm
C
C     Makes logarithmic derivatives at the hard sphere
C
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
C
      jsr=jsrs(ita,it)
      DO 20 l=0,lmax
      x1=ts(lz,1,l,ita,it,is,0)
      x2=ts(lz,2,l,ita,it,is,0)
      x3=ts(lz,3,l,ita,it,is,0)
      x4=ts(lz,4,l,ita,it,is,0)
      x1d=ts(lz,1,l,ita,it,is,1)
      x2d=ts(lz,2,l,ita,it,is,1)
      x3d=ts(lz,3,l,ita,it,is,1)
      x4d=ts(lz,4,l,ita,it,is,1)
C
      dm=dfi(lz,ita,it,l,is,1)
      xs=dfi(lz,ita,it,l,is,1)*x1+x2
      xn=dfi(lz,ita,it,l,is,1)*x3+x4
      d=xs/xn
      dd=(dfi(lz,ita,it,l,is,2)*x1+dfi(lz,ita,it,l,is,1)*x1d+x2d)/xn-
     . d*(dfi(lz,ita,it,l,is,2)*x3+dfi(lz,ita,it,l,is,1)*x3d+x4d)/xn
C
      dfi(lz,ita,it,l,is,1)=d*sigmt(l,it)
      dfi(lz,ita,it,l,is,2)=dd*sigmt(l,it)
      gsng(lz,ita,it,l,is)=dd/d
C
C     Get the norm function 1/N
C
      IF(localmt(ita,it).NE.0) THEN
         chi=fi0m(lz,ita,it,l,is)
      ELSE
         chi=cf(jsr,l,ita,it,is,lz)
      ENDIF
      fi0a=-0.5*sws/(sigmt(l,it)*hsr(ita,it))*chi*xn
      cf(1:jsr,l,ita,it,is,lz)=cf(1:jsr,l,ita,it,is,lz)/fi0a
C
C     Free electron solutions for r > S,...,RI(JRN)
C
      jrn=jrsm(ita,it)
      IF(fcd.EQ.'Y') jrn=MAX0(jrn,jwsc(ita,it))
      DO 24 ir=jsr+1,jrn
      r=ri(ir,ita,it)
      rw=r/sws
C
      IF(localmt(ita,it).EQ.1.AND.r.LE.wsm(ita,it)) THEN
         kap2=zm(lz)-vmtzr(it,is)
         CALL besslz(kap2*r*r,0,l,jl(0),hl(0))
         bs=jl(l)*rw**l
         nm=hl(l)*rw**(-l-1)
         fm= tm(lz,1,l,ita,it,is)*nm+tm(lz,2,l,ita,it,is)*bs
         gm=-tm(lz,3,l,ita,it,is)*nm-tm(lz,4,l,ita,it,is)*bs
         cf(ir,l,ita,it,is,lz)=(fm+gm*dm)*r*chi/fi0a/hsr(ita,it)
      ELSE
         kap2=zm(lz)-vmtz(is)
         CALL besslz(kap2*r*r,0,l,jl(0),hl(0))
         bs=jl(l)*rw**l
         nm=hl(l)*rw**(-l-1)
         fa= ta(lz,1,l,it,is)*nm+ta(lz,2,l,it,is)*bs
         ga=-ta(lz,3,l,it,is)*nm-ta(lz,4,l,it,is)*bs
         cf(ir,l,ita,it,is,lz)=(fa+ga*d)*r
      ENDIF
   24 CONTINUE
C
   20 CONTINUE
C
      IF(func.EQ.'ASA') RETURN
C
C     Integrate the partial waves for the multipole moments
C
      DO 40 is=1,ns
      DO 40 it=1,nt
      DO 40 ita=1,nta(it)
      jws=jwss(ita,it)
      DO 40 lp=0,lmax
      DO 40 l=lp,lmax
      DO 41 lz=1,nzm
      DO 41 lpp=0,lmax2
      DO 42 ir=1,jws
      r=ri(ir,ita,it)
      rw=(r/sws)**lpp
      xs=cf(ir,lp,ita,it,is,lz)*cf(ir,l,ita,it,is,lz) 
      fi(ir) =REAL(xs,8)*rw*r
   42 fip(ir)=AIMAG(xs)*rw*r
      CALL simpn(fi,dx,jws,rephi)
      CALL simpn(fip,dx,jws,imphi)
   41 cfm(lp,l,lpp,ita,it,is,lz)=CMPLX(rephi,imphi,8)
      IF(l.GT.lp) THEN
         cfm(l,lp,0:lmax2,ita,it,is,1:nzm)=
     .   cfm(lp,l,0:lmax2,ita,it,is,1:nzm)
      ENDIF
   40 CONTINUE
C
      RETURN
      END
