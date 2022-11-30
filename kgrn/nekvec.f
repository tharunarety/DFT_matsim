      SUBROUTINE nekvec(nk)
C   ******************************************************************
C   *                                                                *
C   *  Select the non-equivalent k-vecrots.                          *
C   *                                                                *
C   ******************************************************************
      USE bzmesh
      USE message
      IMPLICIT NONE
      REAL(KIND=8), PARAMETER :: err = 1.d-6
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xf, yf, zf, wf
      REAL(KIND=8)  :: xr, yr, zr, wr, diffk, sumw
      INTEGER       :: nk, nf, lk, i, itr
C
      ALLOCATE(xf(nk),yf(nk),zf(nk),wf(nk))
C
C     Save the first vector
C
      xf(1)=fkx(1)
      yf(1)=fky(1)
      zf(1)=fkz(1)
      wf(1)=fkw(1)
C
      nf=1
C
C     Loop for the k-vectors
C
      DO 20 lk=2,nk
C
C     Choose a vectors (2nd,3rd,4th, etc.)
C
      xr=fkx(lk)
      yr=fky(lk)
      zr=fkz(lk)
      wr=fkw(lk)
      itr=0
C
C     Check if this vector was already saved
C
      DO 21 i=1,nf
      diffk=ABS(xr-xf(i))+ABS(yr-yf(i))+ABS(zr-zf(i))
      IF(diffk.LT.err) THEN
C
C        Old vector f(i) is the same as the new fk(lk)
C        Neglect the new and increase the weight for the old
C
         wf(i)=wf(i)+wr
      ELSE
         itr=itr+1
      ENDIF
   21 CONTINUE
C
C     The fk(lk) was not saved yet, save it
C
      IF(itr.EQ.nf) THEN
         nf=nf+1
         xf(nf)=xr
         yf(nf)=yr
         zf(nf)=zr
         wf(nf)=wr
      ENDIF
   20 CONTINUE
C
C     Save the non-equivalent k-vectors in fk array
C
      IF(nk.NE.nf) THEN
         nk=nf
         DEALLOCATE(fkx,fky,fkz,fkw)
         ALLOCATE(fkx(nk),fky(nk),fkz(nk),fkw(nk))
C
         fkx(1:nk)=xf(1:nk)
         fky(1:nk)=yf(1:nk)
         fkz(1:nk)=zf(1:nk)
         fkw(1:nk)=wf(1:nk)
         sumw=SUM(fkw(1:nk))
         IF(ABS(sumw-1.d0).GT.1.d-8) THEN
            WRITE(m6,100) sumw
            STOP
         ENDIF
      ENDIF
C
      DEALLOCATE(xf,yf,zf,wf)
C
  100 FORMAT(/,' NEKVEC:**  SUMW = ',f12.8,' it should be 1.')
      RETURN
      END
