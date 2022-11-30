      SUBROUTINE setrms
C   ******************************************************************
C   *                                                                *
C   *  Set up the radial logarithmic mesh as:                        *
C   *                                                                *
C   *      r(i) = r1 * exp[(i-1) * dx]   i = 1, MAX(jris,jwsc),      *
C   *                                                                *
C   *      r(i) = ws   for  i = jwss      ASA (=WS) sphere,          *
C   *                                                                *
C   *      r(i) >~ wsc  for  i = jwsc     circumscribed sphere,      *
C   *                                                                *
C   *      r(i) ~< wsi  for  i = jwsi     inscribed sphere,          *
C   *                                                                *
C   *      r(i) = s    for  i = jsrs      potential sphere,          *
C   *                                                                *
C   *      jris = jwss + 2.                                          *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE message
      USE radialmesh
      IMPLICIT NONE
      REAL(KIND=8) :: sc, si, sw, s
      INTEGER :: iq, it, ita, js, ji, dimrs, dimrc
C
      ALLOCATE(jsrs(mnta,nt),jwsi(mnta,nt),jwsc(mnta,nt),jrsm(mnta,nt))
      jsrs=0
      jwsi=0
      jwsc=0
      jrsm=0
C
C     Set the maximum size af the radial arrays
C
      jwss=jris-2
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      sc=wsc(iq)
      si=wsi(iq)
      sw=ws(ita,it)
      s=hsr(ita,it)
C
C     The inscribed and circumscribed spheres do not depend on ita
C     but the radial meshes do !!
C
      jwsc(ita,it)=LOG(sc/sw)/dx+1+jwss(ita,it)
      ji=LOG(si/sw)/dx+jwss(ita,it)
      IF(2*(ji/2).EQ.ji) ji=ji-1
      jwsi(ita,it)=ji
C
      js=LOG(s/sw)/dx+jwss(ita,it)
      IF(2*(js/2).EQ.js) js=js+1
      jsrs(ita,it)=js
   20 CONTINUE
      DO 21 it=1,nt
      DO 21 ita=1,nta(it)
      jrsm(ita,it)=MAX0(jsrs(ita,it),jwss(ita,it))+2
   21 hsr(ita,it)=ws(ita,it)*dexp(dx*(jsrs(ita,it)-jwss(ita,it)))
C
      dimrs=MAXVAL(jrsm)+2
      dimrc=MAXVAL(jwsc)
C
      dimr=MAX0(dimrc,dimrs)
C
      WRITE(m6,100) dimrs,dimrc
      WRITE(m6,110)
      DO 22 iq=1,nq
      it=itq(iq)
      DO 22 ita=1,nta(it)
      WRITE(m6,120) it,ita,jwsi(ita,it),jwss(ita,it),jwsc(ita,it),
     .              jsrs(ita,it)
   22 WRITE(m6,130) wsi(iq),ws(ita,it),wsc(iq),hsr(ita,it)
C
  100 FORMAT(/,' SETRMS:   dimension for potential spheres     : ',i3,
     .       //,11x,'dimension for circumscribed spheres : ',i3,/)
  110 FORMAT(11x,'IT ITA    Si      w       Sc      S',/)
  120 FORMAT(11x,i2,i3,5i8,/)
  130 FORMAT(17x,6f8.4)
      RETURN
      END
