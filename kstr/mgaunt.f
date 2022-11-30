      SUBROUTINE mgaunt
C   ******************************************************************
C   *                                                                *
C   *    Set up the non-zero real Gaunt numbers I(l'm';lm;l''m'').   *
C   *    The convention is:   l'  = 0,lmaxh                          *
C   *                         l   = 0,lmax (< or =lmaxh)             *
C   *                         l'' = |l'-l|,l'+l .                    *
C   *    For l',l < lmax only the upper triangle is calculated.      *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE lmtolm
      USE realgaunt
      USE message
      IMPLICIT NONE
      REAL(KIND=8) :: rgnt
      INTEGER :: lmp, lp, myp, lm, l, my, lpp, mypp, ign,
     .           maxlpp, minlpp, lm0, nlm0
C
C
C     Initialize Clebsch-Gordan parameters
C
      CALL setclb
C
      IF(msgl.EQ.1) WRITE(msgio,'(a,i5)')
     .          ' Generate real Gaunt numbers'
      ign=0
      DO 10 lmp=1,nlmh
      lp=llx(lmp)
      myp=mmx(lmp)
      lm0=lmp
      IF(lmp.GT.nlm) lm0=1
      DO 10 lm=lm0,nlm
      l=llx(lm)
      my=mmx(lm)
      minlpp=IABS(lp-l)
      maxlpp=lp+l
      DO 10 lpp=minlpp,maxlpp
      DO 10 mypp=-lpp,lpp
      CALL gaunty(l,my,lp,myp,lpp,mypp,rgnt)
      IF(DABS(rgnt).GT.1.d-8) ign=ign+1
   10 CONTINUE
      ngaunt=ign
C
      ALLOCATE(gnt(ngaunt))
      ALLOCATE(lmg(ngaunt),lmpg(ngaunt),lmppg(ngaunt))
      ALLOCATE(kpow(ngaunt))
C
      IGN=0
      DO 11 lmp=1,nlmh
      lp=llx(lmp)
      myp=mmx(lmp)
      lm0=lmp
      IF(lmp.GT.nlm) lm0=1
      DO 11 lm=lm0,nlm
      l=llx(lm)
      my=mmx(lm)
      minlpp=IABS(lp-l)
      maxlpp=lp+l
      DO 11 lpp=minlpp,maxlpp
      DO 11 mypp=-lpp,lpp
      CALL gaunty(l,my,lp,myp,lpp,mypp,rgnt)
      IF(DABS(rgnt).GT.1.d-8) THEN
        ign=ign+1
        gnt (ign)=rgnt
        lmpg(ign)=lmp
        lmg (ign)=lm
        lmppg(ign)=lpp*lpp+lpp+mypp+1
        kpow(ign)=(lp+l-lpp)/2
      ENDIF
   11 CONTINUE
C
      nlm0=MAX0(nlm,nlmw)
      ign=0
      DO 20 lmp=1,nlm
      lp=llx(lmp)
      myp=mmx(lmp)
      DO 20 lm=1,nlm0
      l=llx(lm)
      my=mmx(lm)
      minlpp=IABS(lp-l)
      maxlpp=lp+l
      DO 20 lpp=minlpp,maxlpp
      DO 20 mypp=-lpp,lpp
      CALL gaunty(l,my,lp,myp,lpp,mypp,rgnt)
      IF(DABS(rgnt).GT.1.d-8) ign=ign+1
   20 CONTINUE
      ngauntw=ign
C
      ALLOCATE(gntw(ngauntw))
      ALLOCATE(lmgw(ngauntw),lmpgw(ngauntw),lmppgw(ngauntw))
      ALLOCATE(jpow(ngauntw))
C
      IGN=0
      DO 21 lmp=1,nlm
      lp=llx(lmp)
      myp=mmx(lmp)
      DO 21 lm=1,nlm0
      l=llx(lm)
      my=mmx(lm)
      minlpp=IABS(lp-l)
      maxlpp=lp+l
      DO 21 lpp=minlpp,maxlpp
      DO 21 mypp=-lpp,lpp
      CALL gaunty(l,my,lp,myp,lpp,mypp,rgnt)
      IF(DABS(rgnt).GT.1.d-8) THEN
        ign=ign+1
        gntw(ign)=rgnt
        lmpgw(ign)=lmp
        lmgw(ign)=lm
        lmppgw(ign)=lpp*lpp+lpp+mypp+1
        jpow(ign)=(lp-l+lpp)/2
      ENDIF
   21 CONTINUE
C
      WRITE(m6,100) ngaunt,ngauntw
      IF(nprn.EQ.2) THEN
         WRITE(m6,110)
         DO 30 ign=1,ngaunt
         WRITE(m6,130) ign,lmpg(ign),lmg(ign),lmppg(ign),gnt(ign)
   30    CONTINUE
      ENDIF
C
  100 FORMAT(/,' MGAUNT:   NGAUNT =',i6,' NGAUNTW =',i6)
  110 FORMAT(/,11X,'Gaunt coefficients in the real harmonic',
     1       ' representation'//,11x,'IG  LMP LM  LMPP',
     2       6X,'GAUNT',/)
  130 FORMAT(10x,4i4,3x,2f10.6)
      RETURN
      END
