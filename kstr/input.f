      SUBROUTINE input(minl,maxl)
C   ******************************************************************
C   *                                                                *
C   *   Read control data etc.                                       *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE message
      IMPLICIT NONE
      INTEGER :: minl, maxl
C
      READ(5,'(6(8x,i2))') nl,nlh,nlw,nder,itrans,nprn
      READ(5,'(3(10x,f10.6))') kap2,dmax,dwats
      IF(nl.LT.1) nl=1
      nlm=nl*nl
      lmax=nl-1
      IF(nlh.LT.nl.OR.high.NE.'Y') THEN
         nlh=nl
      ENDIF
      nlmh=nlh*nlh
      lmaxh=nlh-1
      wats='Y'
      IF(nlw.LT.1) wats='N'
      nlmw=nlw*nlw
      lmaxw=nlw-1
      minl=-MAX0(nder,1)
      maxl=lmaxh+MAX0(lmax,lmaxw)+MAX0(nder,1)
C
      IF(msgl.EQ.1) WRITE(msgio,'(a,i5)')
     1           ' Message level (MSGL) =',msgl
      WRITE(m6,101) mode,store,nprn,nl,nlh,nlw,nder,dmax,dwats
      WRITE(m6,102) kap2,itrans
C
  101 FORMAT(/,' INPUT:',4x,'MODE  =',a3,' STORE =',a3,' NPRN  =',i3,
     1     //,11x,'NL   =',i3,' NLH  =',i3,' NLW  =',i3,' NDER =',i3,
     2      //,11x,'DMAX    =',f10.6,' DWATS  =',f10.6)
  102 FORMAT(/,11x,'(K*W)^2 =',f10.6,' ITRANS =',i2)
      RETURN
      END
