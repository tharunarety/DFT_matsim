      SUBROUTINE wavefc(numb,ef)
C   ******************************************************************
C   *                                                                *
C   *    Find potential parameters.                                  *
C   *                                                                *
C   ******************************************************************
      USE control_data ; USE diracparam ; USE dosmom
      USE message      ; USE potential  ; USE potparam
      USE radialmesh   ; USE temporary  ; USE totalenergy
      USE text
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(4) :: pws, d
      REAL(KIND=8) :: ef, s, s2, vw, ded, e, pu, qu, fint
      REAL(KIND=8) :: m0, a0, b0, om, om0, oml, omx, bm, beta
      REAL(KIND=8) :: sfisql, tnosa, emoma
      INTEGER      :: jws, it, ita, is, ip, il, l, lp1, k
      INTEGER      :: numb, ie, nob = 0
C
      IF(numb.EQ.0) THEN
        ALLOCATE(omm(0:lmax,mnta,nt,ns,pan),sfim(0:lmax,mnta,nt,ns,pan))
        ALLOCATE(fmofp(0:lmax,mnta,nt,ns,pan))
        ALLOCATE(signfi(0:lmax,mnta,nt,ns,pan))
        ALLOCATE(amy(0:lmax,mnta,nt,ns,pan),dny(0:lmax,mnta,nt,ns,pan))
        ALLOCATE(tl(0:lmax,mnta,nt,ns,pan),d3(0:lmax,mnta,nt,ns,pan))
        ALLOCATE(wndw(0:lmax,mnta,nt,ns,pan),vl(0:lmax,mnta,nt,ns,pan))
        ALLOCATE(d1(0:lmax,mnta,nt,ns,pan),d2(0:lmax,mnta,nt,ns,pan))
        ALLOCATE(phil(0:lmax,mnta,nt,ns),logdl(0:lmax,mnta,nt,ns))
        ALLOCATE(epmatrix(0:lmax,mnta,nt,ns))
      ENDIF
C
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      s=ws(ita,it)
      s2=s*s
      jws=jwss(ita,it)
C
      DO 20 is=1,ns
      vw=v(jws,ita,it,is)/s2
      DO 20 ip=1,pan
      DO 20 l=0,lmax
      il=l+1
      IF(numb.EQ.0) THEN
         IF(tpot.NE.'Y') THEN
            IF(pan.EQ.1) THEN
               eny(l,ita,it,is,ip)=ef-depth/2.d0
            ELSE
               IF(ip.EQ.pan) THEN
                  eny(l,ita,it,is,ip)=ef-depth/2.d0
               ELSE
                  eny(l,ita,it,is,ip)=elim-depth/2.d0
               ENDIF
            ENDIF
         ENDIF
      ELSE
         IF(pan.EQ.1) THEN
            tnosa=tnos(ita,it,l,is)
            emoma=emom(ita,it,l,is)
         ELSEIF(pan.EQ.2) THEN
            IF(ip.EQ.pan) THEN
               tnosa=tnos(ita,it,l,is)-tnos0(ita,it,l,is)
               emoma=emom(ita,it,l,is)-emom0(ita,it,l,is)
            ELSE
               tnosa=tnos0(ita,it,l,is)
               emoma=emom0(ita,it,l,is)
            ENDIF
         ENDIF
         IF(tnosa.GT.1.d-6) eny(l,ita,it,is,ip)=emoma/tnosa
      ENDIF
      lp1=l+1
      k=-l-1
      IF(l.le.1) THEN
         ded=0.02d0
      ELSEIF(l.GT.1.AND.l.LE.3) THEN
         ded=0.005d0
      ELSEIF(l.GT.3) THEN
         ded=0.001d0
      ENDIF
C
C     Solve the Dirac equation for ENY-DED, ENY, and ENY+DED
C
      DO 21 ie=1,4
      IF(ie.LE.3) THEN
         e=eny(l,ita,it,is,ip)+(ie-2)*ded
      ELSE
C
C        Get radial solutions at Ef for Hopfield parameter
C
         e=ef
      ENDIF
C
      CALL dirac(e,k,it,ita,il,is,jws,0,beta,0)
C
      pu=rp(jws)
      qu=rq(jws)
      d(ie)=s*(1.d0-(vw-e)/csq)*qu/pu+l
C
      fi(1:jws)=rp(1:jws)*rp(1:jws)*ri(1:jws,ita,it)
      CALL simpn(fi,dx,jws,fint)
      fint=fint+rp(1)**2*ri(1,ita,it)/(1.d0+2.d0*beta)
      pws(ie)=pu/DSQRT(fint)
   21 CONTINUE
C
C     Laurent parameters
C
      m0=1.d0/pws(2)**2/s         !! m = (s*phi^2)^(-1)/s^2
      a0=-pws(2)/2.d0/ded*(pws(3)-pws(1))/s
      b0=-pws(2)*(pws(3)-2.d0*pws(2)+pws(1))/s/3.d0/ded/ded/s2
      IF(nob.EQ.0) THEN
         b0=0.d0
      ELSE
         IF(b0.LT.1.D-06) THEN
            WRITE(m6,101) it,il,is,b0,a0,m0,ded
            b0=0.d0
         ENDIF
      ENDIF
      dny(l,ita,it,is,ip)=d(2)             !! d(eny)  logarithmic derivative
C
C     Second order energies
C
      om=(d(2)+lp1)/(1.d0+(d(2)+lp1)*a0)/m0
      om0=d(2)/(1.d0+d(2)*a0)/m0
      oml=(d(2)-l)/(1.d0+(d(2)-l)*a0)/m0
      omx=1.d0/m0/a0
      omm(l,ita,it,is,ip)=om/s2
      bm=b0*m0
C
C     Standard parameters
C
C     SFIM   : WS*Phi^2(-l-1,WS)
C     FMOFP  : Phi(-)/Phi(+)
C     WNDW   : SQRT(< Phi-dot^2 >) = SQRT(small p)
C
      fmofp(l,ita,it,is,ip)=1.d0-(2.d0*lp1-1.d0)/(lp1+d(2)+1.d0/a0)
      sfim(l,ita,it,is,ip)=1.d0/m0/s2/(1.d0+a0*(d(2)+lp1))**2
      sfisql=1.d0/m0/s2/(1.d0+a0*(d(2)-l))**2
      wndw(l,ita,it,is,ip)=s2*SQRT(bm)
      signfi(l,ita,it,is,ip)=SIGN(1.d0,pws(2))
C
C     Third order estimates
C
      cc(l,ita,it,is,ip)=eny(l,ita,it,is,ip)+om/s2/(1.d0+bm*om**2)
      bot(l,ita,it,is,ip)=eny(l,ita,it,is,ip)+om0/s2/(1.d0+bm*om0**2)
      top(l,ita,it,is,ip)=eny(l,ita,it,is,ip)+omx/s2/(1.d0+bm*omx**2)
      vl(l,ita,it,is,ip)=eny(l,ita,it,is,ip)+oml/s2/(1.d0+bm*oml**2)
      amy(l,ita,it,is,ip)=2.d0/s2/sfim(l,ita,it,is,ip)*
     .             (1.d0+bm*om**2)**2/(1.d0-bm*om**2)
      tl(l,ita,it,is,ip)=(2*l+3)/s2/sfisql*
     .            (1.d0+bm*oml**2)**2/(1.d0-bm*oml**2)
      d1(l,ita,it,is,ip)=dny(l,ita,it,is,ip)+1.d0/(a0+DSQRT(b0/m0))
      d2(l,ita,it,is,ip)=dny(l,ita,it,is,ip)+1.d0/a0
      d3(l,ita,it,is,ip)=dny(l,ita,it,is,ip)+1.d0/(a0-DSQRT(b0/m0))
C
C     Save parameters for the Hopfield parameter
C
      phil(l,ita,it,is)=pws(4)/s
      logdl(l,ita,it,is)=d(4)
   20 CONTINUE
C
C     Compute the electron-phonon matrix elements
C
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      DO 30 is=1,ns
      epmatrix(0:lmax,ita,it,is)=0.d0
      DO 30 l=0,lmax-1
      lp1=l+1
C
      epmatrix(l,ita,it,is)=-phil(l,ita,it,is)*phil(lp1,ita,it,is)*
     .     ((logdl(l,ita,it,is)-l)*(logdl(lp1,ita,it,is)+l+2.d0)+
     .     (ef-vw)*s2)
   30 CONTINUE
C
  101 FORMAT(' WAVEFC:** B0 small or negative',/,' IT,IL,IS=',3i2,
     .       ' B=',f10.6,' A=',f10.6,' M=',f10.6,' DED=',f10.6)
      RETURN
      END
