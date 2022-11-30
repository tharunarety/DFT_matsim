      SUBROUTINE singlg(ef,prnt)
C   ******************************************************************
C   *                                                                *
C   * Find the poles of the regular part of the single-site Green's  *
C   * function on the real axis from EB to EF+ZX(NX):                *
C   *                                                                *
C   *     G(r,r,e) = fi(e)*D_dot(e)/D(e)*fi(e)                       *
C   *                                                                *
C   * Only the roots of the D(e)=D{fi0^a(a)} are considered !!!      *
C   *                                                                *
C   ******************************************************************
      USE energymesh ; USE control_data ; USE control_text
      USE message    ; USE radialmesh   ; USE partialwaves
      USE potparam   ; USE temporary; USE csts
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: tol = 3.d-8, step = 0.01d0
      REAL(KIND=8), PARAMETER :: stp = 0.0002d0
      REAL(KIND=8), DIMENSION(5) :: dr5, dr5d, dr5dd
      REAL(KIND=8) :: ef, efx, ebefx, e, edr, dr, drd, drdecr
      REAL(KIND=8) :: zdzd, zdz, e1, e2, zm1, fint, bandw, r, rw
      INTEGER      :: prnt, mbrack, ip, froot, i, j
      INTEGER      :: it, ita, il, l, lpp, is, jws, jrn, ir
C
C     Top of the complex contour zm(nzm)+(nx-1)/2*hx = zx(nx)
C
      zm1=ABS(10.d0*AIMAG(zx(nx)))
      efx=REAL(zx(nx),8)+MAX(zm1,0.5d0)
C
C     Botom of the contour: zm(1) = eb
C
      necr=0
      DO 20 is=1,ns
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jrn=jrsm(ita,it)
      IF(fcd.EQ.'Y') jrn=MAX0(jrn,jwsc(ita,it))
      jws=jwss(ita,it)
      DO 20 l=0,lmax
      i=0
      il=l+1
      froot=0
      e=eb
   21 ebefx=efx-e
C
      mbrack=MAX(INT(ebefx/step),3)
      IF(iter.GT.1) THEN
         IF(pan.EQ.1) THEN
            ip=pan
         ELSE
            e1=(efx+e)/2.d0-eny(l,ita,it,is,1)
            e2=(efx+e)/2.d0-eny(l,ita,it,is,2)
            IF(ABS(e1).LT.ABS(e2)) THEN
               ip=1
            ELSE
               ip=2
            ENDIF
         ENDIF
         bandw=ABS(top(l,ita,it,is,ip)-bot(l,ita,it,is,ip))/5.d0
         bandw=MAX(bandw,0.01d0*step)
         IF(step.GE.bandw) mbrack=INT(ebefx/bandw)
      ENDIF
      mbrack=MAX(mbrack,1)
C
      CALL rootf(e,ebefx,mbrack,froot,il,it,ita,is,prnt)
C
      IF(froot.EQ.1) THEN
         CALL logdrr0(e,il,it,ita,is,dr,drd,1,zdzd,zdz)
         IF(ABS(dr).GT.1.d-6) THEN
            WRITE(m6,100) e,dr
            STOP
         ENDIF
         i=i+1
         IF(i.GT.dimecr) THEN
            WRITE(m6,101) i,dimecr
            STOP
         ENDIF
         necr(l,ita,it,is)=i
         ecr(l,ita,it,is,i)=e
         cfrr(l,ita,it,is,i)=zdzd
         cfrt(l,ita,it,is,i)=zdz
         cfr(1:jrn,l,ita,it,is,i)=rp(1:jrn)*rp(1:jrn)/drd
         DO 41 lpp=0,lmax2
         DO 42 ir=1,jws
         r=ri(ir,ita,it)
         rw=(r/sws)**lpp
   42    fi(ir)=rp(ir)*rp(ir)*rw*r
         CALL simpn(fi,dx,jws,fint)
   41    cfrm(l,lpp,ita,it,is,i)=fint/drd
C
         drdecr=drd
         DO 22 j=1,5
         edr=e+stp*(j-3)
         CALL logdrr0(edr,il,it,ita,is,dr,drd,0,zdzd,zdz)
   22    dr5(j)=sigmt(l,it)*dr
         CALL diffn(dr5,dr5d,dr5dd,5,stp)
         nocr(l,ita,it,is,i)=drdecr/dr5d(3)
         IF(prnt.NE.0) THEN
            WRITE(m6,102) nocr(l,ita,it,is,i)
            IF(msgl.NE.0) WRITE(msgio,102) nocr(l,ita,it,is,i)
         ENDIF
C
         e=e+tol
         GO TO 21
      ENDIF
   20 CONTINUE
C
  100 FORMAT(/,' SINGLG:   error in root finding, E =',f10.6,
     .         ' D(e) =',f10.6)
  101 FORMAT(/,' SINGLG:   too many roots i =',i3,' dimecr =',i3)
  102 FORMAT(/,' SINGLG:   NOCR =',f10.6)
      RETURN
      END
