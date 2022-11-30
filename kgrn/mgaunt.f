      SUBROUTINE mgaunt(lmaxp,lmax,lmaxpp)
C   ******************************************************************
C   *                                                                *
C   *   Set up the non-zero real Gaunt numbers I(l'm';lm;l''m'').    *
C   *                                                                *
C   ******************************************************************
      USE realgaunt
      USE message
      IMPLICIT NONE
      REAL(KIND=8) :: rgnt
      INTEGER :: lmaxp, lmax, lmaxpp
      INTEGER :: lp, mp, l, m, lpp, mpp, ign, maxlpp, minlpp
C
C     Initialize Clebsch-Gordan parameters
C
      IF(ngaunt.EQ.0) CALL setclb
C
      WRITE(m6,'(/,a,i5)') ' MGAUN:    Generate real Gaunt numbers'
      ign=0
      DO 20 lp   =0,lmaxp
      DO 20 mp   =-lp,lp
      DO 20 l    =0,lmax
      DO 20 m    =-l,l
      minlpp     =IABS(lp-l)
      maxlpp     =MIN0(lp+l,lmaxpp)
      DO 20 lpp  =minlpp,maxlpp
      DO 20 mpp  =-lpp,lpp
      CALL gaunty(l,m,lp,mp,lpp,mpp,rgnt)
      IF(ABS(rgnt).GT.1.d-8) ign=ign+1
   20 CONTINUE
      ngaunt=ign
C
      ALLOCATE(gnt(ngaunt))
      ALLOCATE(lmpg(ngaunt),lmg(ngaunt),lmppg(ngaunt),lppg(ngaunt))
C
      IGN=0
      DO 21 lp   =0,lmaxp
      DO 21 mp   =-lp,lp
      DO 21 l    =0,lmax
      DO 21 m    =-l,l
      minlpp     =IABS(lp-l)
      maxlpp     =MIN0(lp+l,lmaxpp)
      DO 21 lpp  =minlpp,maxlpp
      DO 21 mpp  =-lpp,lpp
      CALL gaunty(l,m,lp,mp,lpp,mpp,rgnt)
      IF(ABS(rgnt).GT.1.d-8) THEN
        ign       =ign+1
        gnt  (ign)=rgnt
        lmpg (ign)=lp*lp+lp+mp+1
        lmg  (ign)=l*l+l+m+1
        lppg (ign)=lpp
        lmppg(ign)=lpp*lpp+lpp+mpp+1
      ENDIF
   21 CONTINUE
C
      WRITE(m6,100) ngaunt
      IF(nprn.EQ.2) THEN
         WRITE(m6,110)
         DO 30 ign=1,ngaunt
         WRITE(m6,120) ign,lmpg(ign),lmg(ign),lmppg(ign),gnt(ign)
   30    CONTINUE
      ENDIF
C
  100 FORMAT(/,' MGAUNT:   NGAUNT =',i6)
  110 FORMAT(/,11X,'Gaunt coefficients in the real harmonic',
     1       ' representation'//,11x,' IG LMP  LM  LMPP',
     2       6X,'GAUNT',/)
  120 FORMAT(10x,4i4,3x,f10.6)
C
      RETURN
      END
