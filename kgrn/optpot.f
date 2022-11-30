      SUBROUTINE optpot(fixg,prnt,ef,vmix)
C   ******************************************************************
C   *                                                                *
C   *  Set up the optimized overlapping muffin-tin potential by      *
C   *  solving the "f" and "g" equations self-consistently.          *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE csts ; USE control_data ; USE control_text ; USE density
      USE fandgeq ; USE message ; USE pota ; USE potential
      USE radialmesh ; USE temporary
      IMPLICIT NONE
      INTEGER :: fixg, prnt
      REAL(KIND=8), DIMENSION(ns) :: vmtzn
      REAL(KIND=8), DIMENSION(nt,ns) :: vmtzrn
      REAL(KIND=8), DIMENSION(nt) :: lg
      REAL(KIND=8), DIMENSION(mnta,nt,dimr) :: fp, f, dens
      REAL(KIND=8) :: ef, fspin, fint, fbar, g
      REAL(KIND=8) :: vlms, vlmw, volit, vmix
      INTEGER :: is, it, ita, jsr, jws, jrn
C
      IF(ns.EQ.1) THEN
         fspin=1.d0
      ELSE
         fspin=.5d0
      ENDIF
C
      DO 20 is=1,ns
      fbar=0.d0
      volit=0.d0
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      jsr=jsrs(ita,it)
      jws=jwss(ita,it)
      jrn=jrsm(ita,it)
C
      IF(fixg.GE.2) THEN
C
C        The "f" and "g" equations are weighted by the density
C
         dens(ita,it,1:jrn)=chdo(ita,it,is,1:jrn)-
     .                      fspin*cor(ita,it,1:jrn)
C
C        Read the spherical part of the full-potential and
C        weight by the density
C
         fp(ita,it,1:jsr)=fullp(1:jsr,ita,it,is)/ri2(1:jsr,ita,it)
C
C        Find the density weighted VolI
C
         IF(ABS(vol-vols).GE.1.d-8) THEN
C
C           Non ASA case (jsr<>jws)
C
            fi(1:jrn)=dens(ita,it,1:jrn)*ri(1:jrn,ita,it)
            CALL simpn(fi,dx,jsr,vlms)
            CALL simpn(fi,dx,jws,vlmw)
            volit=volit+conc(ita,it)*fourpi*(vlmw-vlms)*mmt(it)
C
C           Get the average full-potential within the SCA
C
            fi(1:jws)=dens(ita,it,1:jws)*
     .                fullp(1:jws,ita,it,is)/ri(1:jws,ita,it)
            CALL simpn(fi,dx,jws,fint)
            fbar=fbar+conc(ita,it)*fourpi*fint*mmt(it)
         ELSE
C
C           ASA case (jws=jsr) Fbar is not needed !!
C
         ENDIF
      ELSE
C
C        The "f" and "g" equations are not weighted by the density
C
         dens(ita,it,1:jsr)=ri2(1:jsr,ita,it)
C
C        Read the spherical part of the full-potential
C
         fp(ita,it,1:jsr)=fullp(1:jsr,ita,it,is)/ri2(1:jsr,ita,it) 
C
C        Get the average full-potential within the SCA
C
         fi(1:jws)=fullp(1:jws,ita,it,is)*ri(1:jws,ita,it)
         CALL simpn(fi,dx,jws,fint)
         fbar=fbar+conc(ita,it)*fourpi*fint*mmt(it)
      ENDIF
   21 CONTINUE
      fbar=fbar/vol
      IF(fixg.EQ.2) voli=volit
C
C     Solve the f+g equations
C
      CALL slfcfg(fixg,prnt,fp,f,dens,fbar,g,lg,vmtz0)
C
      IF(prnt.EQ.1) THEN
         WRITE(m6,100) fbar,g,voli
         IF(msgl.NE.0) WRITE(msgio,100) fbar,g,voli
      ENDIF
C
C     Save the spherical potential for the Schrodinger equation
C
      DO 22 it=1,nt
      DO 23 ita=1,nta(it)
      jsr=jsrs(ita,it)
C
C     Local muffin-tin zero (does not depend on ita)
C
      IF(localmt(ita,it).EQ.0) THEN
         vmtzrn(it,is)=g
      ELSE
         vmtzrn(it,is)=lg(it)
      ENDIF
C
      v(1:jsr,ita,it,is)=f(ita,it,1:jsr)*ri2(1:jsr,ita,it) 
      CALL diff(dx,jsr,it,ita,is)
   23 CONTINUE
   22 CONTINUE
C
C     Global muffin-tin zero
C
      vmtzn(is)=g
C
   20 CONTINUE
C
C     Mix the muffin-tin zero
C
      vmtz=vmix*vmtzn+(1.d0-vmix)*vmtz
      vmtzr=vmix*vmtzrn+(1.d0-vmix)*vmtzr
      IF(prnt.EQ.1) THEN
         WRITE(m6,110) vmtz
         IF(msgl.EQ.1) WRITE(msgio,110) vmtz
         IF(lclmff.NE.0) THEN
            DO 24 it=1,nt
            WRITE(m6,120) it,vmtzr(it,1:ns)
            IF(msgl.EQ.1) WRITE(msgio,120) it,vmtzr(it,1:ns)
   24       CONTINUE
         ENDIF
      ENDIF
C
  100 FORMAT(/,' OPTPOT:   Fbar =',f10.6,' g =',f10.6,' VolI =',f10.6)
  110 FORMAT(/,11x,'VMTZ(Up,Dwn) = ',2f14.6)
  120 FORMAT(/,11x,'Local muffin-tin zero for IT = ',i3,' is ',2f14.6)
      RETURN
      END
