C     Last change:  HLS  25 Feb 1999    1:06 pm
      SUBROUTINE SETCLB
C   ******************************************************************
C   *                                                                *
C   *    Initialize factorials used in CLBSHG                        *
C   *                                                                *
C   ******************************************************************
      USE clebsch_gordon
      IMPLICIT NONE
      REAL(KIND=8) :: X
      INTEGER :: I
      ALLOCATE(H(MCLB),J(MCLB))
      H(1)=1.0D0
      J(1)=0
      X=0.D0
      DO 20 I=2,MCLB
      X=X+1.0D0
      H(I)=H(I-1)*X
      J(I)=J(I-1)
   21 IF(H(I).LT.10.0D0) GO TO 20
      H(I)=0.01D0*H(I)
      J(I)=J(I)+2
      GO TO 21
   20 CONTINUE
      RETURN
      END
