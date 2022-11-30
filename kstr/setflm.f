      SUBROUTINE setflm
C   ******************************************************************
C   *                                                                *
C   *    Set factorials and l,m sequences and initialize harmonics.  *
C   *                                                                *
C   ******************************************************************
      USE factorial
      USE control_data
      USE csts
      USE lmtolm
      USE screening
      IMPLICIT NONE
      INTEGER :: lmaxhw, nlmhw, nfac2
      INTEGER :: i, j, lp, l, m, lm
C
      nfac2=nlh+MAX0(nlw,nl)-2
      lmaxhw=MAX0(lmaxh,MAX0(lmaxw,lmax))
      nlmhw=MAX0(nlmh,MAX0(nlmw,nlm))
      ALLOCATE(fac2(0:nfac2),efac(0:nfac2),sig(0:nfac2))
      ALLOCATE(edot(0:nfac2,0:nder))
      ALLOCATE(fack(0:lmaxh,0:lmax),facj(0:lmax,0:lmaxw))
      ALLOCATE(llx(nlmhw),mmx(nlmhw))
C
C     Factorial (2l-1)!!
C
      fac2(0)=1.d0
      DO 20 i=1,nfac2
   20 fac2(i)=fac2(i-1)*(2*i-1)
C
C     sig(l)=(-1)**l
C     efac(l)=(-kap2)**l
C
      efac(0)=1.d0
      edot(0,0)=1.d0
      DO 21 j=1,nder
   21 edot(0,j)=0.d0
      sig(0)=1.d0
      DO 22 l=1,nfac2
      efac(l)=-kap2*efac(l-1)
      edot(l,0)=efac(l)
      sig(l)=-sig(l-1)
      DO 22 j=1,nder
   22 edot(l,j)=-l*edot(l-1,j-1)
C
C     fack(lp,l)=-8*pi/(2l-1)!!/(2lp-1)!!*sig(l)
C     facj(lp,l)=-8*pi*(2l-1)!!/(2lp-1)!!
C
      DO 23 lp=0,lmaxh
      DO 23 l=0,lmax
   23 fack(lp,l)=-8.d0*pi*sig(l)/fac2(lp)/fac2(l)
      DO 24 lp=0,lmax
      DO 24 l=0,lmaxw
   24 facj(lp,l)=-8.d0*pi*fac2(l)/fac2(lp)
C
C     Real harmonic sequence for l=s,p,d, ...
C
      DO 30 l=0,lmaxhw
      DO 30 m=-l,l
      lm=l*l+l+m+1
      llx(lm)=l
   30 mmx(lm)=m
C
C     Initialize real harmonics
C
      CALL realhinit(nfac2)
C
C     Initialize Gaunt numbers
C
      CALL mgaunt
C
      RETURN
      END
