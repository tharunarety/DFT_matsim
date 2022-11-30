      SUBROUTINE fhndlr
C   ******************************************************************
C   *                                                                *
C   *    Reads file names, opens for I/O, and prints program header  *
C   *    and file names.                                             *
C   *                                                                *
C   ******************************************************************
C   *                                                                *
C   *    Fortran units, files and their content:                     *
C   *                                                                *
C   *    FOR001 (IN)   : Slope matrices generated by KSTR,           *
C   *                    one or two different kappa values.          *
C   *    FOR002 (I/O)  : Control data, charge densities, potentials, *
C   *                    etc. generated in an earlier run of KGRN    *
C   *                    with STRT = 'B' or 'C'. This is the file    *
C   *                    used as input to PGRN.                      *
C   *    FOR003 (I/O)  : Control data, charge densities, potentials, *
C   *                    etc. generated in an earlier run of KGRN    *
C   *                    with STRT = 'A'.                            *
C   *    FOR004 (I/O)  : Madelung matrix.                            *
C   *    FOR005 (IN)   : Specifies K-mesh, etc.                      *
C   *    FOR006 (OUT)  : Printer.                                    *
C   *    FOR008 (I/O)  : Temporary storage for ATOMA.                *
C   *    FOR009 (OUT)  : File for z-mesh data.                       *
C   *    FOR010 (OUT)  : Full-charge density for FCD calculation.    *
C   *                                                                *
C   ******************************************************************
      USE control_text ; USE message ; USE text
      IMPLICIT NONE
      CHARACTER(LEN=4) :: iddata
      INTEGER :: ifor
C
C     Read file names and control data from standard input
C
      READ(5,'(a)') iddata
      IF(iddata.NE.'KGRN') THEN
         WRITE(6,100) iddata
         STOP
      ENDIF
      msgio=6
      m6=7
C
      READ(5,'(7x,a)') job
      READ(5,'(9x,a,7x,i3,3(9x,a))') strt,msgl,expan,fcd,func
      READ(5,101) for001
      READ(5,101) for001p
      READ(5,101) for002
      READ(5,101) for003
      READ(5,101) for004
      READ(5,101) for006
      READ(5,101) for009
      READ(5,101) for010
      READ(5,101) for011
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
      txt(51:55)=clock
      txt(56:58)=' / '
      txt(59:67)=dato
      WRITE(m6,'(2a,/)') ' TASK:     ',txt
      WRITE(m6,'(2a,/)') ' FHNDLR:   Program: KGRN, Version: ',
     .           version
      WRITE(m6,102) 'JOB:    ',job
      ifor=INDEX(for001,' ')-1
      WRITE(m6,102) 'FOR001: ',for001(1:ifor)
      ifor=INDEX(for001p,' ')-1
      IF(expan.EQ.'D') WRITE(m6,102) 'FOR001: ',for001p(1:ifor)
      ifor=INDEX(for002,' ')-1
      WRITE(m6,102) 'FOR002: ',for002(1:ifor)
      ifor=INDEX(for003,' ')-1
      WRITE(m6,102) 'FOR003: ',for003(1:ifor)
      ifor=INDEX(for004,' ')-1
      WRITE(m6,102) 'FOR004: ',for004(1:ifor)
      ifor=INDEX(for006,' ')-1
      WRITE(m6,102) 'FOR006: ',for006(1:ifor)
      ifor=INDEX(for009,' ')-1
      WRITE(m6,102) 'FOR009: ',for009(1:ifor)
      ifor=INDEX(for010,' ')-1
      WRITE(m6,102) 'FOR010: ',for010(1:ifor)
      ifor=INDEX(for011,' ')-1
      WRITE(m6,102) 'FOR011: ',for011(1:ifor)
C
C     Open input files
C
      OPEN(1,FILE=for001,STATUS='OLD',FORM='UNFORMATTED')
      OPEN(4,FILE=for004,STATUS='OLD',FORM='UNFORMATTED')
C
C     Open output file
C
      OPEN(2,FILE=for002,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(3,FILE=for003,STATUS='UNKNOWN',FORM='UNFORMATTED')
      OPEN(9,FILE=for009,STATUS='UNKNOWN',FORM='UNFORMATTED')
      IF(fcd.EQ.'Y') THEN
      OPEN(10,FILE=for010,STATUS='UNKNOWN',FORM='UNFORMATTED')
      ENDIF
C
      WRITE(m6,103) msgl,strt,func,expan,fcd
      IF(strt.NE.'A'.AND.strt.NE.'B'.AND.strt.NE.'N') THEN
         WRITE(m6,104)
         STOP
      ENDIF
C
  100 FORMAT(/,' FHNDLR: **Unidentified input, IDDATA =',a)
  101 FORMAT(7X,a)
  102 FORMAT(' ',10x,2a)
  103 FORMAT(/,11x,'MSGL = ',i1,' STRT =',a,' FUNC = ',a,
     .             ' EXPAN = ',a,' FCD = ',a)
  104 FORMAT(/,' FHNDLR:  **STRT must be one of the following:',/,
     .       /,11x,'A     for the first start',
     .       /,11x,'B     for restart an old calculation with',
     .             ' old k, z mesh parameters',
     .       /,11x,'N     for restart an old calculation with',
     .             ' new k, z mesh parameters')
      RETURN
      END