      SUBROUTINE rootf(e0,ebefx,mbrack,froot,il,it,ita,is,prnt)
C   ******************************************************************
C   *                                                                *
C   *    Find the root of the D{fi0^a}                               *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE message
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: tol = 1.d-12, stol = 1.d-12
      REAL(KIND=8)  :: tol1
      REAL(KIND=8)  :: ebefx, e0, fa, fb, fc, step, drd, zdzd, zdz
      REAL(KIND=8)  :: a, b, c, d, e, s, p, q, r, xm
      INTEGER, PARAMETER :: mroot = 50
      INTEGER :: mbrack, froot, il, it, ita, is, prnt
      INTEGER :: ipass, ibrack, iroot
C
C     Test for the first point
C
      ipass=1
      CALL logdrr0(e0,il,it,ita,is,fa,drd,0,zdzd,zdz)
      IF(ABS(fa).LT.tol) THEN
         froot=1
         IF(prnt.EQ.1) THEN
            IF(msgl.EQ.1) WRITE(msgio,100) 
     .                    ipass,il-1,it,ita,is,e0,fa
            WRITE(m6,100) ipass,il-1,it,ita,is,e0,fa
         ENDIF
         RETURN
      ENDIF
C
C     Bracket root
C
      a=e0
      b=a
      step=ebefx/mbrack
      fc=fa
      DO 10 ibrack=1,mbrack
      ipass=ipass+1
      b=b+step
      CALL logdrr0(b,il,it,ita,is,fb,drd,0,zdzd,zdz)
      IF(ABS(fb).LT.tol) THEN
         e0=b
         froot=1
         IF(prnt.EQ.1) THEN
            IF(msgl.EQ.1) WRITE(msgio,100) 
     .                    ipass,il-1,it,ita,is,e0,fb
            WRITE(m6,100) ipass,il-1,it,ita,is,e0,fb
         ENDIF
         RETURN
      ENDIF
      IF(fa*fb.LT.0.d0) THEN
         IF((fa-fc)*(fb-fa).GE.stol) GO TO 20
      ENDIF
      a=b
      fc=fa
      fa=fb
   10 CONTINUE
C
C     No root found
C
      froot=0
      RETURN
C
C     Root has been bracket
C
   20 c = a
      fc = fa
      d = b - a
      e = d
   30 IF(ABS(fc).GE.ABS(fb)) GO TO 40
      a = b
      b = c
      c = a
      fa = fb
      fb = fc
      fc = fa
C
C     Convergence test
C
   40 tol1 = 2.0d0*meps*ABS(b) + 0.5d0*tol
      xm = .5*(c - b)
      IF(ABS(xm).LE.tol1) GO TO 90
      IF(ABS(fb).LE.tol) GO TO 90
C
C     is bisection necessary
C
      IF(ABS(e).LT.tol1) GO TO 70
      IF(ABS(fa).LE.ABS(fb)) GO TO 70
C
C     is quadratic interpolation possible
C
      IF(ABS(a-c).GT.tol) GO TO 50
C
C     linear interpolation
C
      s = fb/fa
      p = 2.0d0*xm*s
      q = 1.0d0 - s
      GO TO 60
C
C     inverse quadratic interpolation
C
   50 q = fa/fc
      r = fb/fc
      s = fb/fa
      p = s*(2.0d0*xm*q*(q - r) - (b - a)*(r - 1.0d0))
      q = (q - 1.0d0)*(r - 1.0d0)*(s - 1.0d0)
C
C     adjust signs
C
   60 IF(p.GT.0.0d0) q = -q
      p = ABS(p)
C
C     is interpolation acceptable
C
      IF((2.0d0*p).GE.(3.0d0*xm*q - ABS(tol1*q))) GO TO 70
      IF(p.GE.ABS(0.5d0*e*q)) GO TO 70
      e = d
      d = p/q
      GO TO 80
C
C     bisection
C
   70 d = xm
      e = d
C
C     complete step
C
   80 a = b
      fa = fb
      IF(ABS(d).GT.tol1) THEN
         b = b + d
      ELSE
         b = b + DSIGN(tol1, xm)
      ENDIF
      ipass=ipass+1
      CALL logdrr0(b,il,it,ita,is,fb,drd,0,zdzd,zdz)
      IF((fb*(fc/ABS(fc))).GT.0.0d0) GO TO 20
      GO TO 30
C
   90 e0=b
      froot=1
      IF(prnt.EQ.1) THEN
         IF(msgl.EQ.1) WRITE(msgio,100) 
     .                 ipass,il-1,it,ita,is,e0,fb
         WRITE(m6,100) ipass,il-1,it,ita,is,e0,fb
      ENDIF
C
  100 FORMAT(/,' ROOTF: IPASS =',i3,' L,IT,ITA,IS =',4i3,
     .         ' E0= ',f10.6,' D(e0)= ',1pe10.2)
      RETURN
      END
