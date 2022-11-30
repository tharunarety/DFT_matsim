      SUBROUTINE fhndlr
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
C   *    Files and their content                                     *
C   *                                                                *
C   *    FOR001 (Out)   : Most localized slope matrices for the      *
C   *                     lattice generated in LATT3D.               *
C   *    FOR005 (In)    : Control data.                              *
C   *    FOR006 (Out)   : Printer.                                   *
C   *                                                                *
C   ******************************************************************
      USE control_text
      USE message
      IMPLICIT NONE
      CHARACTER(LEN=10) :: iddata
      CHARACTER(LEN=9)  :: dato
      CHARACTER(LEN=5)  :: clock
      CHARACTER(LEN=1)  :: hp
C
C     Read file names and control data from standard input
C
      READ(5,'(a,9x,a)') iddata,hp
      IF(iddata.NE.'KSTR') THEN
         WRITE(6,102) iddata
         STOP
      ENDIF
C
C     Open standard output for HP730
C
      msgio=6
      m6=7
      IF(hp.EQ.'HP') THEN
         msgio=7
         m6=6
         OPEN(msgio)
      ENDIF
      READ(5,'(10x,a,7x,i3,4(9x,a))') job,msgl,mode,store,high
      READ(5,100) for001
      READ(5,100) for006
C
      CALL jobnam
C
C     Open for printer output, and print header and file names
C
      OPEN(m6,FILE=for006,STATUS='UNKNOWN',FORM='FORMATTED')
      READ(5,'(a)') txt
C
      CALL htimer(dato,clock)
C
      WRITE(m6,'(2a,16x,2a)') ' FHNDLR:   Program: KSTR, Version: ',
     .           version,'   TM: ',clock
      txt(59:67)=dato
      WRITE(m6,'(11x,a,/)') txt
      WRITE(m6,101) 'JOB   : ',job
      WRITE(m6,101) 'FOR001: ',for001
      IF(high.EQ.'Y') WRITE(m6,101) 'For002: ',FOR002
      WRITE(M6,101) 'FOR006: ',for006
C
C     Open output files if requested
C
      IF(store.EQ.'Y') THEN
         OPEN(1,FILE=for001,STATUS='UNKNOWN',FORM='UNFORMATTED')
         IF(high.EQ.'Y') THEN
            OPEN(2,FILE=for002,STATUS='UNKNOWN',FORM='UNFORMATTED')
         ENDIF
      ENDIF
C
  100 FORMAT(7X,a)
  101 FORMAT(11x,10a,2a)
  102 FORMAT(/,' FHNDLR: **Unidentified input, IDDATA =',a)
      RETURN
      END
