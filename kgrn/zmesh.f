      SUBROUTINE zmesh(ef,npr)
C   ******************************************************************
C   *                                                                *
C   *    Establish the complex energy mesh on which the Green's      *
C   *    function is to be calculated.                               *
C   *                                                                *
C   *    NPR   = 1 : Print the energy contour.                       *
C   *          = 0 : No print.                                       *
C   *                                                                *
C   *    ZMSH  ='C': Generate the semi-circular contour along which  *
C   *                the Green's function is integrated by means of  *
C   *                an n-point Gauss expression.                    *
C   *          ='M': Generate two semi-circular contours along       *
C   *                which the Green's function is integrated.       *
C   *          ='E': Generate the semi-elliptic contour along which  *
C   *                the Green's function is integrated.             *
C   *                The small axis (b) along the imaginar axis is   *
C   *                smaller or equal with 0.5 Ry.                   *
C   *          ='c': Equivalent with C but with dense mesh in both   *
C   *                ends of the contour.                            *
C   *          ='m': Equivalent with M but with dense mesh in both   *
C   *                ends of the upper contour.                      *
C   *          ='e': Equivalent with E but with dense mesh in both   *
C   *                ends of the contour.                            *
C   *          ='D': Generate 'horizontal' mesh to examine the       *
C   *                state density.                                  *
C   *          ='F': Generate contour and residues used to include   *
C   *                the Fermi function.                             *
C   *          ='f': Generate contour and residues used to include   *
C   *                the Fermi function plus a lower contour.        *
C   *          ='S': Fermi surface calculation.                      *
C   *          ='T': Small eliptic contour for electronic entropy.   *
C   *                                                                *
C   ******************************************************************
      USE botomtop     ; USE control_data ; USE control_text
      USE csts         ; USE energymesh   ; USE message
      USE text
      IMPLICIT NONE
      COMPLEX(KIND=8) :: zx0
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: ags, wgs
      REAL(KIND=8) :: ef, de, epsx, yf, x0, dyf, emlz, fib, fia
      REAL(KIND=8) :: zm1, zmi, zmr, y0, fif, z0, zr, length, hd
      REAL(KIND=8) :: hxmax, deptha, cags, sags, drda, a, b, ba
      REAL(KIND=8) :: boltzk = 6.336d-6
      INTEGER      :: lz, iz, ires, npr, lx, nzmh, nzl
C
      IF(zmsh.EQ.'F') THEN
         nzm=nz1+nz2+nz3+nres
         ALLOCATE(ags(MAX0(nz1,nz2)),wgs(MAX0(nz1,nz2)))
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
C        Contour with Fermi function of temperature TFERMI (Ry)
C
         zm1=ef-depth
         IF(zm1.GT.ebt) THEN
            WRITE(m6,120) zm1,ebt,depth
         ENDIF

C        Gausian mesh on lower semi-circle
C
         eb=zm1
         y0=-2.d0*delta
         fib=-pi
         fia=ATAN(y0/depth)
         fif=2.d0*fia
         zr=SQRT(y0*y0+depth*depth)-delta
         zr=0.5d0*zr/COS(fia)
         z0=eb+zr
         CALL wagaus(fib,fif,ags,wgs,nz1)
         DO 20 lz=1,nz1
         zmr=z0+zr*COS(ags(lz))
         zmi=zr*SIN(ags(lz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=CMPLX(0.d0,wgs(lz),8)*(zm(lz)-z0)
     .            /(1.d0+EXP((zm(lz)-ef)/tfermi))
   20    CONTINUE
C
C        Gausian mesh on intermediate semi-circle
C
         zr=delta
         x0=ef
         fib=pi+fia
         fif=0.5d0*pi
         x0=ef
         CALL wagaus(fib,fif,ags,wgs,nz2)
         lz=nz1
         DO 21 iz=1,nz2
         lz=lz+1
         zmr=x0+zr*COS(ags(iz))
         zmi=y0+zr*SIN(ags(iz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=CMPLX(0.d0,wgs(iz),8)*(zm(lz)-CMPLX(x0,y0,8))
     .             /(1.d0+EXP((zm(lz)-ef)/tfermi))
   21    CONTINUE
C
C        Gaussian mesh on semi-circle from (ef,-delta) to (ef+delta,0)
C
         zr=delta
         fib=-pi/2.d0
         fif=0.d0
         epsx=1.01*fib
         CALL lgaus(fib,fif,nz3,epsx,ags,wgs)
         DO 22 iz=1,nz3
         lz=lz+1
         zmr=ef+zr*COS(ags(iz))
         zmi=zr*SIN(ags(iz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=CMPLX(0.d0,wgs(iz),8)*(zm(lz)-ef)
     .            /(1.d0+EXP((zm(lz)-ef)/tfermi))
   22    CONTINUE
C
C        Get the poles of the Fermi function
C
         lz=nzm+1
         yf=pi*tfermi
         dyf=2.d0*pi*tfermi
         DO 23 ires=1,nres
         lz=lz-1
         wgm(lz)=CMPLX(0.d0,dyf,8)
         yf=yf-dyf
   23    zm(lz)=CMPLX(ef,yf,8)
C
C        The small linear contour
C
         nx0=(nx+1)/2
         zx0=zm(nzm)
         hxmax=REAL(2.d0*(zx0-zm(+1))/(nx-1),8)
         IF(hx.GT.hxmax) THEN
            WRITE(m6,121) hx,hxmax
            hx=hxmax
         ENDIF
         DO 24 lx=1,nx
   24    zx(lx)=zx0+hx*(lx-nx0)
         IF(npr.EQ.1) THEN
            WRITE(m6,100) 
     .      'Contour with Fermi function, T =',tfermi/boltzk
            WRITE(m6,110) elt,eb,ef,depth,delta,epsx
            DO 25 lz=1,nzm
            WRITE(m6,113) lz,zm(lz),wgm(lz)
            IF(lz.EQ.nz1) WRITE(m6,'()')
            IF(lz.EQ.nz1+nz2) WRITE(m6,'()')
            IF(lz.EQ.nz1+nz2+nz3) WRITE(m6,'()')
   25       CONTINUE
            WRITE(m6,114)
            DO 26 lx=1,nx
   26       WRITE(m6,113) lx,zx(lx)
         ENDIF
         DEALLOCATE(ags,wgs)
      ELSEIF(zmsh.EQ.'T') THEN
C
C        Elliptic contour for electronic entropy
C
         nzm=MIN0(2*nz3,nz1)
         ALLOCATE(ags(nzm),wgs(nzm))
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
         zm1=ef-MIN(delta,depth)
         eb=zm1
         z0=ef
         a=MIN(delta,depth)
         b=-pi*tfermi/2.d0
         ba=b/a
         fib=-pi/2.d0
         fif=0.d0
         nzmh=nzm/2
C
C        Generate logarithmic Gaussian mesh on the semi-elliptic contour
C
         CALL lgaus(fib,fif,nzmh,eps,ags,wgs)
C
         DO 90 lz=1,nzmh
         cags=DCOS(ags(lz))
         sags=DSIN(ags(lz))
         zr=a/DSQRT(cags*cags+sags*sags/ba/ba)
         zmr=z0+zr*DCOS(ags(lz))
         zmi=zr*DSIN(ags(lz))
C
C        dl = i * (z-z0) * d alpha + (z-z0) * d r/d alpha 1/r d alpha
C
         drda=-(zmr-z0)*zmi*(1.d0/b/b-1.d0/a/a)
         zm(lz+nzmh)=CMPLX(zmr,zmi,8)
         wgm(lz+nzmh)=wgs(lz)*((0.d0,1.d0)+drda)*(zm(lz+nzmh)-z0)
         iz=nzmh-lz+1
         zmr=z0-zr*COS(ags(lz))
         zm(iz)=CMPLX(zmr,zmi,8)
         drda=-(zmr-z0)*zmi*(1.d0/b/b-1.d0/a/a)
         wgm(iz)=wgs(lz)*((0.d0,1.d0)+drda)*(zm(iz)-z0)
   90    CONTINUE
         zx(1:nx)=zm(nzm)
C
         IF(npr.EQ.1) THEN
            WRITE(m6,101) 'Semi-elliptic contour'
            WRITE(m6,115) elt,eb,ef,MIN(delta,depth),ba,eps
            DO 91 lz=1,nzm
   91       WRITE(m6,113) lz,zm(lz),wgm(lz)
         ENDIF
         DEALLOCATE(ags,wgs)
      ELSEIF(zmsh.EQ.'C'.OR.zmsh.EQ.'c') THEN
C
C        Contour with step function at the Fermi level
C
         nzm=nz1
         ALLOCATE(ags(nzm),wgs(nzm))
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
         zm1=ef-depth
         IF(zm1.GT.ebt) THEN
            WRITE(m6,120) zm1,ebt,depth
         ENDIF
         eb=zm1
         z0=0.5d0*(ef+eb)
         zr=0.5d0*(ef-eb)
C
         nzmh=nzm
         fib=-pi
         fif=0.d0
C
         IF(zmsh.EQ.'c') THEN
            nzmh=nzm/2
            fib=-pi/2.d0
         ENDIF
C
C        Generate logarithmic Gaussian mesh on the semi-circle
C
         CALL lgaus(fib,fif,nzmh,eps,ags,wgs)
         DO 30 lz=1,nzmh
         zmr=z0+zr*COS(ags(lz))
         zmi=zr*SIN(ags(lz))
         IF(zmsh.EQ.'C') THEN
            zm(lz)=CMPLX(zmr,zmi,8)
            wgm(lz)=wgs(lz)*(0.d0,1.d0)*(zm(lz)-z0)
         ELSE
            zm(lz+nzmh)=CMPLX(zmr,zmi,8)
            wgm(lz+nzmh)=wgs(lz)*(0.d0,1.d0)*(zm(lz+nzmh)-z0)
            iz=nzmh-lz+1
            zmr=z0-zr*COS(ags(lz))
            zm(iz)=CMPLX(zmr,zmi,8)
            wgm(iz)=wgs(lz)*(0.d0,1.d0)*(zm(iz)-z0)
         ENDIF
   30    CONTINUE
         zm1=ABS(100.d0*AIMAG(zm(1)))
         hd=depth/2.d0
         eb=eb-MIN(zm1,hd)
C
C        The small linear contour
C
         nx0=(nx+1)/2
         zx0=zm(nzm)
         hxmax=REAL(2.d0*(zx0-zm(+1))/(nx-1),8)
         IF(hx.GT.hxmax) THEN
            WRITE(m6,121) hx,hxmax
            hx=hxmax
         ENDIF
         DO 31 lx=1,nx
   31    zx(lx)=zx0+hx*(lx-nx0)
C
         IF(npr.EQ.1) THEN
            WRITE(m6,101) 'Semi-circular contour'
            WRITE(m6,110) elt,eb,ef,depth,delta,eps
            DO 32 lz=1,nzm
   32       WRITE(m6,113) lz,zm(lz),wgm(lz)
            WRITE(m6,114)
            DO 33 lx=1,nx
   33       WRITE(m6,113) lx,zx(lx)
         ENDIF
         DEALLOCATE(ags,wgs)
      ELSEIF(zmsh.EQ.'f') THEN
         nzl=nz2
         nzm=nz1+nz2+nz3+nres+nzl
         ALLOCATE(ags(MAX0(nz1,nz2)),wgs(MAX0(nz1,nz2)))
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
C        Lower panel starting from Ef-ELIM down
C
         zm1=ef-elim-depth
         IF(zm1.GT.ebt) THEN
            WRITE(m6,120) zm1,ebt,depth
         ENDIF
         eb=zm1
         z0=0.5d0*(ef-elim+eb)
         zr=0.5d0*(ef-elim-eb)
         fib=-pi
         fif=0.d0
C
         CALL wagaus(fib,fif,ags,wgs,nzl)
         DO 10 lz=1,nzl
         zmr=z0+zr*COS(ags(lz))
         zmi=zr*SIN(ags(lz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=wgs(lz)*(0.d0,1.d0)*(zm(lz)-z0)
   10    CONTINUE
C
C        Contour with Fermi function of temperature TFERMI (Ry)
C
         zm1=ef-depth
         IF(zm1.GT.ebt) THEN
            WRITE(m6,120) zm1,ebt,depth
         ENDIF

C        Gausian mesh on lower semi-circle
C
         eb=zm1
         y0=-2.d0*delta
         fib=-pi
         fia=ATAN(y0/depth)
         fif=2.d0*fia
         zr=SQRT(y0*y0+depth*depth)-delta
         zr=0.5d0*zr/COS(fia)
         z0=eb+zr
         CALL wagaus(fib,fif,ags,wgs,nz1)
         lz=nzl
         DO 11 iz=1,nz1
         lz=lz+1
         zmr=z0+zr*COS(ags(iz))
         zmi=zr*SIN(ags(iz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=CMPLX(0.d0,wgs(iz),8)*(zm(lz)-z0)
     .            /(1.d0+EXP((zm(lz)-ef)/tfermi))
   11    CONTINUE
C
C        Gausian mesh on intermediate semi-circle
C
         zr=delta
         x0=ef
         fib=pi+fia
         fif=0.5d0*pi
         x0=ef
         CALL wagaus(fib,fif,ags,wgs,nz2)
         lz=nz1+nzl
         DO 12 iz=1,nz2
         lz=lz+1
         zmr=x0+zr*COS(ags(iz))
         zmi=y0+zr*SIN(ags(iz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=CMPLX(0.d0,wgs(iz),8)*(zm(lz)-CMPLX(x0,y0,8))
     .             /(1.d0+EXP((zm(lz)-ef)/tfermi))
   12    CONTINUE
C
C        Gaussian mesh on semi-circle from (ef,-delta) to (ef+delta,0)
C
         zr=delta
         fib=-pi/2.d0
         fif=0.d0
         epsx=1.01*fib
         CALL lgaus(fib,fif,nz3,epsx,ags,wgs)
         DO 13 iz=1,nz3
         lz=lz+1
         zmr=ef+zr*COS(ags(iz))
         zmi=zr*SIN(ags(iz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=CMPLX(0.d0,wgs(iz),8)*(zm(lz)-ef)
     .            /(1.d0+EXP((zm(lz)-ef)/tfermi))
   13    CONTINUE
C
C        Get the poles of the Fermi function
C
         lz=nzm+1
         yf=pi*tfermi
         dyf=2.d0*pi*tfermi
         DO 14 ires=1,nres
         lz=lz-1
         wgm(lz)=CMPLX(0.d0,dyf,8)
         yf=yf-dyf
         zm(lz)=CMPLX(ef,yf,8)
   14    CONTINUE
C
C        The small linear contour
C
         nx0=(nx+1)/2
         zx0=zm(nzm)
         hxmax=REAL(2.d0*(zx0-zm(+1))/(nx-1),8)
         IF(hx.GT.hxmax) THEN
            WRITE(m6,121) hx,hxmax
            hx=hxmax
         ENDIF
         DO 15 lx=1,nx
   15    zx(lx)=zx0+hx*(lx-nx0)
C
         eb=elim-depth
         IF(npr.EQ.1) THEN
            WRITE(m6,100) 
     .      'Contour with Fermi function, T =',tfermi/boltzk
            WRITE(m6,110) elt,eb,ef,depth,delta,epsx
            DO 16 lz=1,nzm
            WRITE(m6,113) lz,zm(lz),wgm(lz)
            IF(lz.EQ.nzl) WRITE(m6,'()')
            IF(lz.EQ.nzl+nz1) WRITE(m6,'()')
            IF(lz.EQ.nzl+nz1+nz2) WRITE(m6,'()')
            IF(lz.EQ.nzl+nz1+nz2+nz3) WRITE(m6,'()')
   16       CONTINUE
            WRITE(m6,114)
            DO 17 lx=1,nx
   17       WRITE(m6,113) lx,zx(lx)
         ENDIF
         DEALLOCATE(ags,wgs)
      ELSEIF(zmsh.EQ.'M'.OR.zmsh.EQ.'m') THEN
C
C        Double contour with step function at the Fermi level
C
         nzm=nz1+nz2
         ALLOCATE(ags(nzm),wgs(nzm))
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
C        Lower panel starting from Ef-ELIM down
C
         zm1=ef-elim-depth
         IF(zm1.GT.ebt) THEN
            WRITE(m6,120) zm1,ebt,depth
         ENDIF
         eb=zm1
         z0=0.5d0*(ef-elim+eb)
         zr=0.5d0*(ef-elim-eb)
         fib=-pi
         fif=0.d0
C
         CALL wagaus(fib,fif,ags,wgs,nz2)
         DO 40 lz=1,nz2
         zmr=z0+zr*COS(ags(lz))
         zmi=zr*SIN(ags(lz))
         zm(lz)=CMPLX(zmr,zmi,8)
         wgm(lz)=wgs(lz)*(0.d0,1.d0)*(zm(lz)-z0)
   40    CONTINUE
C
C        Upper panel starting from Ef down
C        The upper contour should not overlap with the lower one !!
C
         deptha=MIN(elim,depth)
C
         zm1=ef-deptha
         z0=0.5d0*(ef+zm1)
         zr=0.5d0*(ef-zm1)
C
         nzmh=nz1
         fib=-pi
         fif=0.d0
C
         IF(zmsh.EQ.'m') THEN
            nzmh=nz1/2
            fib=-pi/2.d0
         ENDIF
C
         CALL lgaus(fib,fif,nzmh,eps,ags,wgs)
         DO 41 lz=1,nzmh
         zmr=z0+zr*COS(ags(lz))
         zmi=zr*SIN(ags(lz))
         IF(zmsh.EQ.'M') THEN
            zm(lz+nz2)=CMPLX(zmr,zmi,8)
            wgm(lz+nz2)=wgs(lz)*(0.d0,1.d0)*(zm(lz+nz2)-z0)
         ELSE
            zm(lz+nzmh+nz2)=CMPLX(zmr,zmi,8)
            wgm(lz+nzmh+nz2)=wgs(lz)*(0.d0,1.d0)*(zm(lz+nzmh+nz2)-z0)
            iz=nzmh-lz+1
            zmr=z0-zr*COS(ags(lz))
            zm(iz+nz2)=CMPLX(zmr,zmi,8)
            wgm(iz+nz2)=wgs(lz)*(0.d0,1.d0)*(zm(iz+nz2)-z0)
         ENDIF
   41    CONTINUE
         zm1=ABS(10.d0*AIMAG(zm(1)))
         hd=depth/2.d0
         eb=eb-MIN(zm1,hd)
C
C        The small linear contour
C
         nx0=(nx+1)/2
         zx0=zm(nzm)
         hxmax=REAL(2.d0*(zx0-zm(nz2+1))/(nx-1),8)
         IF(hx.GT.hxmax) THEN
            WRITE(m6,121) hx,hxmax
            hx=hxmax
         ENDIF
         DO 42 lx=1,nx
   42    zx(lx)=zx0+hx*(lx-nx0)
C
         IF(npr.EQ.1) THEN
            WRITE(m6,101) 'Two semi-circular contours'
            WRITE(m6,111) elt,eb,ef,deptha,elim,eps
            DO 43 lz=1,nzm
   43       WRITE(m6,113) lz,zm(lz),wgm(lz)
            WRITE(m6,114)
            DO 44 lx=1,nx
   44       WRITE(m6,113) lx,zx(lx)
         ENDIF
         DEALLOCATE(ags,wgs)
      ELSEIF(zmsh.EQ.'S') THEN
         nzm=1
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
         eb=ef-depth
         wgm(1)=-CMPLX(imagz/ABS(imagz),0.d0)
         zmi=imagz
         zm(1)=CMPLX(ef,zmi,8)
         IF(npr.EQ.1) THEN
            WRITE(m6,113) 1,zm(1),wgm(1)
            WRITE(m6,101) 'Fermi surface: z= (ef,imagz)'
         ENDIF
         zx=zm(1)
      ELSEIF(zmsh.EQ.'D') THEN
         WRITE(m6,101) 'Horizontal mesh for the DOS'
         nzm=nzd
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
C        Linear contour below the real axis
C
         wgm(1:nzm)=0.d0
         IF(nzm.GT.1) THEN
            length=1.2d0*depth
            zm1=ef-length
            de=length/(nzm-50)
            DO 51 lz=1,nzm
            emlz=zm1+(lz-1)*de
   51       zm(lz)=CMPLX(emlz,imagz,8)
            eb=zm1
            zm1=ABS(100.d0*AIMAG(zm(1)))
            hd=depth/2.d0
            eb=eb-MIN(zm1,hd)
         ELSE
            zm(1)=CMPLX(ebt,imagz,8)
         ENDIF
         zx=zm(nzm)
      ELSEIF(zmsh.EQ.'E'.OR.zmsh.EQ.'e') THEN
C
C        Elliptic contour with step function at the Fermi level
C
         nzm=nz1
         ALLOCATE(ags(nzm),wgs(nzm))
         ALLOCATE(zm(nzm),wgm(nzm),zx(nx))
C
         zm1=ef-depth
         IF(zm1.GT.ebt) THEN
            WRITE(m6,120) zm1,ebt,depth
         ENDIF
         eb=zm1
         z0=0.5d0*(ef+eb)
         a=0.5d0*(ef-eb)
         ba=1.d0
         b=ba*a
         IF(b.GT.0.5d0) THEN
            b=0.5d0
            ba=b/a
         ENDIF
         fib=-pi
         fif=0.d0
         nzmh=nzm
C
         IF(zmsh.EQ.'e') THEN
            fib=-pi/2.d0
            nzmh=nzm/2
         ENDIF
C
C        Generate logarithmic Gaussian mesh on the semi-elliptic contour
C
         CALL lgaus(fib,fif,nzmh,eps,ags,wgs)
C
         DO 60 lz=1,nzmh
         cags=DCOS(ags(lz))
         sags=DSIN(ags(lz))
         zr=a/DSQRT(cags*cags+sags*sags/ba/ba)
         zmr=z0+zr*DCOS(ags(lz))
         zmi=zr*DSIN(ags(lz))
C
C        dl = i * (z-z0) * d alpha + (z-z0) * d r/d alpha 1/r d alpha
C
         drda=-(zmr-z0)*zmi*(1.d0/b/b-1.d0/a/a)
         IF(zmsh.EQ.'E') THEN
            zm(lz)=CMPLX(zmr,zmi,8)
            wgm(lz)=wgs(lz)*((0.d0,1.d0)+drda)*(zm(lz)-z0)
         ELSE
            zm(lz+nzmh)=CMPLX(zmr,zmi,8)
            wgm(lz+nzmh)=wgs(lz)*((0.d0,1.d0)+drda)*(zm(lz+nzmh)-z0)
            iz=nzmh-lz+1
            zmr=z0-zr*COS(ags(lz))
            zm(iz)=CMPLX(zmr,zmi,8)
            drda=-(zmr-z0)*zmi*(1.d0/b/b-1.d0/a/a)
            wgm(iz)=wgs(lz)*((0.d0,1.d0)+drda)*(zm(iz)-z0)
         ENDIF
   60    CONTINUE
         zm1=ABS(10.d0*AIMAG(zm(1)))
         hd=depth/2.d0
         eb=eb-MIN(zm1,hd)
C
C        The small linear contour
C
         nx0=(nx+1)/2
         zx0=zm(nzm)
         hxmax=REAL(2.d0*(zx0-zm(+1))/(nx-1),8)
         IF(hx.GT.hxmax) THEN
            WRITE(m6,121) hx,hxmax
            hx=hxmax
         ENDIF
         DO 61 lx=1,nx
   61    zx(lx)=zx0+hx*(lx-nx0)
C
         IF(npr.EQ.1) THEN
            WRITE(m6,101) 'Semi-elliptic contour'
            WRITE(m6,112) elt,eb,ef,depth,ba,eps
            DO 62 lz=1,nzm
   62       WRITE(m6,113) lz,zm(lz),wgm(lz)
            WRITE(m6,114)
            DO 63 lx=1,nx
   63       WRITE(m6,113) lx,zx(lx)
         ENDIF
         DEALLOCATE(ags,wgs)
      ENDIF
C
  100 FORMAT(/,' ZMESH:    ',a,f10.2,' K')
  101 FORMAT(/,' ZMESH:    ',a,/)
  110 FORMAT(/,11x,'ELT  =',f10.6,
     .       ' EB   =',f10.6,' EFX=',f10.6,/,11x,'DEPTH=',f10.6,
     .       ' DELTA=',f10.6,' EPS=',f10.6,//,11x,'LZ',5x,'Re(z)',
     .       5x,'Im(z)',4x,'Weight',/)
  111 FORMAT(/,11x,'ELT  =',f10.6,
     .       ' EB   =',f10.6,' EFX=',f10.6,/,11x,'DEPTHA=',f10.6,
     .       ' ELIM =',f10.6,' EPS=',f10.6,//,11x,'LZ',5x,'Re(z)',
     .       5x,'Im(z)',4x,'Weight',/)
  112 FORMAT(/,' ZMESH:    ELT  =',f10.6,
     .       ' EB   =',f10.6,' EFX=',f10.6,/,11x,'DEPTH=',f10.6,
     .       ' b/a  =',f10.6,' EPS=',f10.6,//,11x,'LZ',5x,'Re(z)',
     .       5X,'Im(z)',4X,'Weight',/)
  113 FORMAT(10x,i3,1x,4f10.6)
  114 FORMAT(/,11x,'LX',5x,'Re(z)',5x,'Im(z)',/)
  115 FORMAT(' ZMESH:    ELT  =',f10.6,
     .       ' EB   =',f10.6,' EFX=',f10.6,/,11x,'DELTA=',f10.6,
     .       ' b/a  =',f10.6,' EPS=',f10.6,//,11x,'LZ',5x,'Re(z)',
     .       5x,'Im(z)',4x,'Weight',/)
  120 FORMAT(/,' ZMESH:**  Start of contour, ZM1 =',f10.6,', above ',
     .       'bottom, EB =',f10.6,/,11X,'Increase DEPTH =',f10.6)
  121 FORMAT(/,' ZMESH: hx=',f10.6,' is changed to hxmax=',f10.6)
      RETURN
      END
