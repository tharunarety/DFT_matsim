      SUBROUTINE JOBNAM
C   ******************************************************************
C   *                                                                *
C   *    Generate standard file names based on the JOB.              *
C   *                                                                *
C   ******************************************************************
      USE control_text
      IMPLICIT NONE
      CHARACTER(LEN=11) :: JOBT
      INTEGER :: IJ, I1, I, I6
      JOBT(1:11)=JOB(1:10)//' '
      IJ=INDEX(JOBT,' ')-1
      I1=INDEX(FOR001,' ')-1
      I=I1+IJ+4
      FOR001(1:I)=FOR001(1:I1)//JOBT(1:IJ)//'.tfm'
      FOR001(I+1:60)='                                        '
      FOR002(1:I)=FOR001(1:I1)//JOBT(1:IJ)//'.tfh'
      FOR002(I+1:60)='                                        '
      I6=INDEX(FOR006,' ')-1
      I=I6+IJ+4
      FOR006(1:I)=FOR006(1:I6)//JOBT(1:IJ)//'.prn'
      FOR006(I+1:60)='                                        '
      RETURN
      END
