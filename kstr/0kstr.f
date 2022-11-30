C   ******************************************************************
C   *                                                                *
C   *                       *** KSTR ***                             *
C   *                                                                *
C   *  Real space calculation of the energy dependent slope matrix   *
C   *                                                                *
C   *                          L. Vitos                              *
C   *                 Condensed Matter Theory Group                  *
C   *            Physics Department, Uppsala University              *
C   *                   S-75121 Uppsala, Sweden                      *
C   *                                                                *
C   *                        H. L. Skriver                           *
C   *                     Department of Physics                      *
C   *                Technical University of Denmark                 *
C   *                        DK-2800 Lyngby                          *
C   *                                                                *
C   ******************************************************************
      PROGRAM kstr
      USE control_data
      USE control_text
      USE message
      IMPLICIT NONE
      CHARACTER(LEN=9) :: dato
      CHARACTER(LEN=5) :: clock
      INTEGER :: npa, maxl, minl
C
      version='LV 23 Jan 2008'
C
      CALL fhndlr
      CALL setcst
      CALL input(minl,maxl)
      CALL primv
      CALL primkr
      CALL set3d
      CALL latt3d
      CALL blatts
C
C     Initialize Bessel functions
C
      CALL bessinit(minl,maxl)
      CALL screen
      CALL setflm
      CALL madl3
      IF(mode.EQ.'S') THEN
         CALL layer(npa)
         IF(store.EQ.'Y') CALL storel(npa)
         CALL salpl
      ELSEIF(mode.EQ.'B') THEN
         IF(store.EQ.'Y') CALL storel(npa)
         CALL salpl
      ENDIF
      IF(high.EQ.'Y') CLOSE(1,STATUS='DELETE')
      CALL htimer(dato,clock)
C
      WRITE(m6,'(/,2a)') ' KSTR:     Finished at : ',clock
      END
