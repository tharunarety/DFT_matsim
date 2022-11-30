      SUBROUTINE slfcfg(fixg,prnt,p,f,d,fbar,g,lg,vmtz0)
C   ******************************************************************
C   *                                                                *
C   *  Set up the optimized overlapping muffin-tin potential by      *
C   *  solving the "f" and "g" equations self-consistently.          *
C   *  Input:                                                        *
C   *    p   : the spherical part of the full-potential;             *
C   *    fbar: the average of the full-potential over the unit cell; *
C   *    d   : valence density*r^2.                                  *
C   *  Output:                                                       *
C   *    f   : spherical potential [ v(r) = f(r) + g ];              *
C   *    g   : the muffin-tin zero.                                  *
C   *                                                                *
C   *  Four cases:                                                   *
C   *    1) fixg = 0 : f+g equations are solved self-consistently.   *
C   *    2) fixg = 1 : f equation is solved self-consistently;       *
C   *                  g is fixed at fbari + vmtz0.                  *
C   *    3) fixg = 2 : f is fixed at the spherical part of p;        *
C   *                  g is fixed at fbari.                          *
C   *    3) fixg = 3 : f is fixed at the spherical part of p;        *
C   *                  g is fixed at fbari + vmtz0.                  *
C   *  Note:                                                         *
C   *  The subroutine is set up for logarithmic radial mesh!!        *
C   *                                                                *
C   ******************************************************************
      USE atomicdens ; USE csts ; USE control_data ; USE fandgeq
      USE message ; USE radialmesh ; USE temporary
      IMPLICIT NONE
      INTEGER :: fixg, prnt
      REAL(KIND=8), DIMENSION(mnta,nt,dimr) :: p, f, fn, d
      REAL(KIND=8), DIMENSION(nt) :: lg, dlg
      REAL(KIND=8) :: fbar, g, gn, fbari, vmtz0
C
      REAL(KIND=8), PARAMETER ::  mixfg = 0.9d0
      REAL(KIND=8), PARAMETER ::  tol = 1.d-6, toll = 1.d-6
      REAL(KIND=8) :: smm, smmp, fint, fit, r, rp, x, t, s, sj, err
      INTEGER, PARAMETER :: nloop = 200
      INTEGER :: it, ita, jsr, jli, jt, jta, ir, jr, ntait
      INTEGER :: kr, n, iov, jov, loop, iso, jso, novit
C
C     Initial guess for g, the non-overlapping case, Eq. (9)
C
      lg=0.d0
      dlg=0.d0
      IF(ABS(voli).GT.1.d-8) THEN
         smm=0.d0
         DO 20 it=1,nt
         ntait=nta(it)
         DO 21 ita=1,ntait
         jsr=jsrs(ita,it)
         fi(1:jsr)=d(ita,it,1:jsr)*p(ita,it,1:jsr)*
     .             ri2(1:jsr,ita,it)/ri(1:jsr,ita,it)
         CALL simpn(fi,dx,jsr,fint)
         dlg(it)=dlg(it)+conc(ita,it)*d(ita,it,jsr)
         lg(it)=lg(it)+conc(ita,it)*d(ita,it,jsr)*p(ita,it,jsr)
         smm=smm+conc(ita,it)*mmt(it)*fint
   21    CONTINUE
   20    CONTINUE
         fbari=vol*fbar/voli-fourpi*smm/voli
      ELSE
         smm=0.d0
         smmp=0.d0
         DO 22 it=1,nt
         ntait=nta(it)
         DO 23 ita=1,ntait
         jsr=jsrs(ita,it)
         dlg(it)=dlg(it)+conc(ita,it)*d(ita,it,jsr)
         lg(it)=lg(it)+conc(ita,it)*d(ita,it,jsr)*p(ita,it,jsr)
         smm=smm+conc(ita,it)*mmt(it)*d(ita,it,jsr)*p(ita,it,jsr)
         smmp=smmp+conc(ita,it)*mmt(it)*d(ita,it,jsr)
   23    CONTINUE
   22    CONTINUE
         fbari=smm/smmp
      ENDIF
C
      DO 24 it=1,nt
      ntait=nta(it)
      IF(dlg(it).LT.1.d-8) THEN
         lg(it)=0.d0
         DO 25 ita=1,ntait
         jsr=jsrs(ita,it)
         lg(it)=lg(it)+conc(ita,it)*p(ita,it,jsr)
   25    CONTINUE
      ELSE
         lg(it)=lg(it)/dlg(it)
      ENDIF
   24 CONTINUE
C
      IF(fixg.EQ.1.OR.fixg.EQ.3) THEN
         g=fbari+vmtz0
      ELSE
         g=fbari
      ENDIF
C
C     Initial guess for f the non-overlapping case, Eq. (6)
C
      f=0.d0
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      jsr=jsrs(ita,it)
      DO 30 ir=1,jsr
   30 f(ita,it,ir)=ri(ir,ita,it)*(p(ita,it,ir)-g)
C
      IF(fixg.GE.2) GO TO 60
C
C     Self-consistent solution of f and g, Eqs. (3) and (5)
C
      loop=1
   40 CONTINUE
C
C     New f(r), Eq. (3)
C
      smm=0.d0
      DO 50 it=1,nt
      DO 50 ita=1,nta(it)
      jsr=jsrs(ita,it)
C
C     Minimal overlap radius for IT
C
      novit=nov(ita,it)
      iso=MINVAL(jrov(ita,it,1:novit))
      DO 51 ir=1,jsr
      r=ri(ir,ita,it)
C
      fn(ita,it,ir)=r*(p(ita,it,ir)-g)
C
      IF(ir.GE.iso) THEN
C
         DO 52 iov=1,novit
C
         IF(ir.LT.jrov(ita,it,iov)) GO TO 52
C
         t=tov(ita,it,iov)/2.d0
         n=wov(ita,it,iov)
         jt=type(ita,it,iov)
C
         fint=0.d0
         DO 53 jta=1,nta(jt)
         sj=hsr(jta,jt)
         jov=ovt(jta,jt,it)
         jso=jrov(jta,jt,jov)-1
C
C        Overlapping spheres
C
         kr=0
         IF((r-2.d0*t+sj).GT.tol) THEN
            DO 54 jr=jso,jsrs(jta,jt)
            rp=ri(jr,jta,jt)
            IF((rp-2.d0*t+r).GE.-tol) THEN
               kr=kr+1
               fi(kr)=f(jta,jt,jr)*rp
            ELSE
               jli=jr
            ENDIF
   54       CONTINUE
C
            CALL gensim(fi,dx,kr,fint)
C
C           Integral from x = 2t-s to r(jli+1), the first mesh point
C           inside the overlap region, using a linear interpolation
C           between r(jli) (< x) and r(jli+1)
C
            x=2.d0*t-r
            fit=0.5d0*(ri(jli+1,jta,jt)-x)*(f(jta,jt,jli+1)+
     .          f(jta,jt,jli)+
     .          (f(jta,jt,jli+1)-f(jta,jt,jli))*(x-ri(jli,jta,jt))/
     .          (ri(jli+1,jta,jt)-ri(jli,jta,jt)))
            fint=fint+conc(jta,jt)*fit
C
         ENDIF
   53    CONTINUE
         fn(ita,it,ir)=fn(ita,it,ir)-n*fint/(4.d0*t)
   52    CONTINUE
      ENDIF
C
   51 CONTINUE
      fn(ita,it,jsr+1:dimr)=0.d0
C
C     New g, Eq. (5)
C
      IF(fixg.EQ.0) THEN
         fi(1:jsr)=ri2(1:jsr,ita,it)*fn(ita,it,1:jsr)
         CALL simpn(fi,dx,jsr,fint)
         smm=smm+conc(ita,it)*mmt(it)*fint
      ENDIF
   50 CONTINUE
C
      IF(fixg.EQ.0) THEN
         gn=fbar-4.d0*pi*smm/vol
      ELSE
         gn=g
      ENDIF
C
      err=ABS(SUM(f(1:mnta,1:nt,1:dimr)-fn(1:mnta,1:nt,1:dimr)))/dimr+
     .    ABS(g-gn)
      IF(err.GT.toll) THEN
         DO 41 it=1,nt
         DO 41 ita=1,nta(it)
         jsr=jsrs(ita,it)
   41    f(ita,it,1:jsr)=(1.d0-mixfg)*f(ita,it,1:jsr)+
     .                   mixfg*fn(ita,it,1:jsr)
         g=(1.d0-mixfg)*g+mixfg*gn
         loop=loop+1
         IF(loop.LT.nloop) GO TO 40
         IF(prnt.EQ.1) THEN
            WRITE(m6,100) nloop
            IF(msgl.NE.0) WRITE(msgio,100) nloop
         ENDIF
C
C        Use the non-overlapping guess if it is not converged, Eqs. (6) and (9)
C
         f=0.d0
         IF(fixg.EQ.1.OR.fixg.EQ.3) THEN
            g=fbari+vmtz0
         ELSE
            g=fbari
         ENDIF
         DO 42 it=1,nt
         DO 42 ita=1,nta(it)
         jsr=jsrs(ita,it)
         DO 42 ir=1,jsr
   42    f(ita,it,ir)=ri(ir,ita,it)*(p(ita,it,ir)-g)
      ELSE
         IF(prnt.EQ.1) THEN
            WRITE(m6,110) loop,g,fbari
            IF(msgl.NE.0) WRITE(msgio,110) loop,g,fbari
         ENDIF
      ENDIF
C
   60 CONTINUE
C
C     Prepare output: v(r) = f(r) + g and g
C
      DO 61 it=1,nt
      DO 61 ita=1,nta(it)
      jsr=jsrs(ita,it)
   61 f(ita,it,1:jsr)=f(ita,it,1:jsr)/ri(1:jsr,ita,it)+g
C
  100 FORMAT(/,' SLFCFG:  Not converged NLOOP =',i4)
  110 FORMAT(/,' SLFCFG:  Converged in ',i3,' iterations g =',f10.6,
     .         ' Fbari =',f10.6)
      RETURN
      END
