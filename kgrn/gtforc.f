      SUBROUTINE gtforce
C   ******************************************************************
C   *                                                                *
C   *  Force calculation.                                            *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE force
      IMPLICIT NONE
C
      ALLOCATE(fxy(nq),fxz(nq),fxx(nq))
      fxy=0.d0
      fxz=0.d0
      fxx=0.d0
      CALL xchstr
      ALLOCATE(fey(nq),fez(nq),fex(nq))
      fey=0.d0
      fez=0.d0
      fex=0.d0
      CALL elestr
      ALLOCATE(fiy(nq),fiz(nq),fix(nq))
      fiy=0.d0
      fiz=0.d0
      fix=0.d0
      CALL slfstr
      ALLOCATE(fcy(nq),fcz(nq),fcx(nq))
      fcy=0.d0
      fcz=0.d0
      fcx=0.d0
      ALLOCATE(fqy(nq),fqz(nq),fqx(nq))
      fqy=0.d0
      fqz=0.d0
      fqx=0.d0
      CALL nspstr
      CALL forcep
C
      RETURN
      END
