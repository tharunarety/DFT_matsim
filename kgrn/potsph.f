      SUBROUTINE potsph
C   ******************************************************************
C   *                                                                *
C   *    Establish the optimal potential sphere radii according to:  *
C   *                                                                *
C   *      1. v_R(S_R) ~ v_const for each site R;                    *
C   *                                                                *
C   *      2. v_const has the maximal possible value;                *
C   *                                                                *
C   *      3. S_R < Si_R * maximal linear overlap (~50%).            *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE message
      USE potential
      USE radialmesh
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) , DIMENSION(mnta,nt) :: pmax
      REAL(KIND=8)   :: potir, maxpot, vconst, potdif, potdif0
      INTEGER        :: it, ita, ir, jsr, jsrit
C
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jsr=jsrs(ita,it)
      maxpot=-200.d0
      DO 21 ir=1,jsr
      potir=SUM(v(ir,ita,it,1:ns))/ns/ri2(ir,ita,it)
      IF(potir.GT.maxpot) maxpot=potir
   21 CONTINUE
   20 pmax(ita,it)=maxpot
C
      vconst=200.d0
      DO 30 it=1,nt
      DO 30 ita=1,nta(it)
      IF(pmax(ita,it).LT.vconst) vconst=pmax(ita,it)
   30 CONTINUE
C
      DO 40 it=1,nt
      DO 40 ita=1,nta(it)
      potdif0=100.d0
      jsr=jsrs(ita,it)
      DO 41 ir=1,jsr
      potir=SUM(v(ir,ita,it,1:ns))/ns/ri2(ir,ita,it)
      potdif=ABS(potir-vconst)
      IF(potdif.LT.potdif0) THEN
         potdif0=potdif
         jsrit=ir
      ENDIF
   41 CONTINUE
      IF(2*(jsrit/2).EQ.jsrit) THEN
         potir=SUM(v(jsrit-1,ita,it,1:ns))/ns/ri2(jsrit-1,ita,it)
         potdif=ABS(potir-vconst)
         potir=SUM(v(jsrit+1,ita,it,1:ns))/ns/ri2(jsrit+1,ita,it)
         potdif0=ABS(potir-vconst)
         IF(potdif.GT.potdif0) THEN
            jsrit=jsrit+1
         ELSE
            jsrit=jsrit-1
         ENDIF
      ENDIF
      IF(jsrit.NE.jsrs(ita,it)) THEN
         potir=SUM(v(jsr,ita,it,1:ns))/ns/ri2(jsr,ita,it)
         WRITE(m6,100) it,ita,hsr(ita,it),jsrs(ita,it),ri(jsrit,ita,it),
     .                 jsrit,potir,vconst
         IF(msgl.NE.0) WRITE(msgio,110) 
     .                 it,ita,hsr(ita,it),jsrs(ita,it),ri(jsrit,ita,it),
     .                 jsrit,potir,vconst
      ELSE
         WRITE(m6,120) it,ita
         IF(msgl.NE.0) WRITE(msgio,130) it,ita
      ENDIF
      hsr(ita,it)=ri(jsrit,ita,it)
      jsrs(ita,it)=jsrit
      jrsm(ita,it)=MAX0(jsrit,jwss(ita,it))+2
   40 CONTINUE
C
  100 FORMAT(/,' POTSPH:   potential sphere radius changed for',
     .         ' IT=',i3,' ITA=',i3,//,11x,'S_old = ',f10.6,' (',i3,
     .         ') --> S_new = ',f10.6,' (',i3,')',//,11x,
     .         'V(S_old) = ',f10.6,' --> V(S_new) = ',f10.6)
  110 FORMAT(/,' PTSPH: Potential sphere radius changed for',
     .         ' IT=',i3,' ITA=',i3,//,8x,'S_old = ',f10.6,' (',i3,
     .         ') --> S_new = ',f10.6,' (',i3,')',//,8x,
     .         'V(S_old) = ',f10.6,' --> V(S_new) = ',f10.6)
  120 FORMAT(/,' POTSPH:   potential sphere radius unchanged for',
     .         ' IT=',i3,' ITA=',i3)
  130 FORMAT(/,' PTSPH: Potential sphere radius unchanged for',
     .         ' IT=',i3,' ITA=',i3)
      RETURN
      END
