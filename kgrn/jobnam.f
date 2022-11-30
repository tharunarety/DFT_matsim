      SUBROUTINE jobnam
C   ******************************************************************
C   *                                                                *
C   *    Generate standard file names based on the JOB.              *
C   *                                                                *
C   ******************************************************************
      USE control_text
      USE message
      IMPLICIT NONE
      INTEGER :: ij, i, i2, i3, i6, i9, i10, i11, ispace
C
      IJ=INDEX(JOB,' ')-1
      I2=INDEX(FOR002,' ')-1
      I=I2+IJ+4
      FOR002(1:I)=FOR002(1:I2)//JOB(1:IJ)//'.pot'
      DO ispace=I+1,80
         FOR002(ispace:ispace)=' '
      ENDDO
C
      I3=INDEX(FOR003,' ')-1
      I=I3+IJ+4
      FOR003(1:I)=FOR003(1:I3)//JOB(1:IJ)//'.atm'
      DO ispace=I+1,80
         FOR003(ispace:ispace)=' '
      ENDDO
C
      I6=INDEX(FOR006,' ')-1
      I=I6+IJ+4
      FOR006(1:I)=FOR006(1:I6)//JOB(1:IJ)//'.prn'
      DO ispace=I+1,80
         FOR006(ispace:ispace)=' '
      ENDDO
C
      I9=INDEX(FOR009,' ')-1
      I=I9+IJ+4
      FOR009(1:I)=FOR009(1:I9)//JOB(1:IJ)//'.zms'
      DO ispace=I+1,80
         FOR009(ispace:ispace)=' '
      ENDDO
C
      I10=INDEX(FOR010,' ')-1
      I=I10+IJ+4
      FOR010(1:I)=FOR010(1:I10)//JOB(1:IJ)//'.chd'
      DO ispace=I+1,80
         FOR010(ispace:ispace)=' '
      ENDDO
C
      I11=INDEX(FOR011,' ')-1
      I=I11+IJ+4
      FOR011(1:I)=FOR011(1:I11)//JOB(1:IJ)//'.tmp'
      DO ispace=I+1,80
         FOR011(ispace:ispace)=' '
      ENDDO
C
      RETURN
      END
