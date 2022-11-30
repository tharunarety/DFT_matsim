      SUBROUTINE FHNDLR(DATO,CLOCK)
C   ******************************************************************
C   *                                                                *
C   *    Reads file names, opens for I/O, and prints program header  *
C   *    and file names.                                             *
C   *                                                                *
C   *   *At run time HP must be set according to:                    *
C   *                                                                *
C   *    HP  System                                                  *
C   *                                                                *
C   *    Y   HP                                                      *
C   *    N   All others                                              *
C   *                                                                *
C   ******************************************************************
C   *                                                                *
C   *    Files and their content:                                    *
C   *                                                                *
C   *    FOR001 (Out)   : Madelung potential matrix.                 *
C   *    FOR005 (In)    : Control data.                              *
C   *    FOR006 (Out)   : Printer.                                   *
C   *                                                                *
C   ******************************************************************
      USE control_data
      USE control_text
      USE message
      IMPLICIT NONE
      CHARACTER(LEN=10) :: IDDATA
      CHARACTER(LEN=9)  :: DATO
      CHARACTER(LEN=5)  :: CLOCK
      CHARACTER(LEN=1)  :: HP
  100 FORMAT(7X,A)
  101 FORMAT(' ',10X,2A)
  102 FORMAT(/,' FHNDLR: ** Wrong input file, header = ',A)
C
C     Read system ID, file names and control data from standard input
C
      READ(5,'(A,9X,A)') IDDATA,HP
      IF(IDDATA.NE.'BMDL') THEN
         WRITE(6,102) IDDATA
         STOP
      ENDIF
      MSGIO=6
      M6=7
C
C     Open standard output for HP730
C
      IF(HP.EQ.'Y') THEN
         M6=6
         MSGIO=7
         OPEN(MSGIO)
      ENDIF
      READ(5,'(10X,A,2(7X,I3))') JOB,MSGL,NPRN
      READ(5,100) FOR001
      READ(5,100) FOR006
C
      CALL JOBNAM
C
C     Open for printer output, and print header and file names
C
      OPEN(M6,FILE=FOR006,STATUS='UNKNOWN',FORM='FORMATTED')
      CALL HTIMER(DATO,CLOCK)
      WRITE(M6,'(2A,16X,2A)') ' FHNDLR:   Program: BMDL, Version: ',
     1      VERSION,'   TM: ',CLOCK
      READ(5,'(A)') TXT
      TXT(59:67)=DATO
      WRITE(M6,'(11X,A,/)') TXT
      WRITE(M6,101) 'JOB   : ',JOB
      WRITE(M6,101) 'FOR001: ',FOR001
      WRITE(M6,101) 'FOR006: ',FOR006
C
C     Open structure-constant output file
C
      OPEN(1,FILE=FOR001,STATUS='UNKNOWN',FORM='UNFORMATTED')
      RETURN
      END
