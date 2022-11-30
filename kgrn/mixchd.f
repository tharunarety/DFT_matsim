      SUBROUTINE mixchd(mix,mixkey,lin)
C   ******************************************************************
C   *                                                                *
C   *    Mix charge densities for use in next iteration.             *
C   *                                                                *
C   *   *On entry:                                                   *
C   *                                                                *
C   *    CHDN   : Output electron density from present iteration.    *
C   *    CHDO   : Input electron density to present iteration.       *
C   *                                                                *
C   *    MIXKEY: = 0: No mixing                                      *
C   *              1: Simple mix of charge density                   *
C   *              2: Mix charge- and spin densities separately      *
C   *                                                                *
C   *   *On exit:                                                    *
C   *                                                                *
C   *    CHDO   : Input electron density for next iteration.         *
C   *                                                                *
C   ******************************************************************
      USE density
      USE control_data
      USE control_text
      USE message
      USE moments
      USE radialmesh
      USE temporary
      IMPLICIT NONE
      REAL(KIND=8) :: mix, f1, f2, sf1, sf2, smix
      INTEGER      :: mixkey, it, ita, jrn, lin
C
      f1=1.d0-mix
      f2=mix
      IF(mixkey.EQ.0) THEN
         chdo=chdn
      ELSEIF(mixkey.EQ.1) THEN
         IF(lin.EQ.1) qlmo=f1*qlmo+f2*qlmn
         DO 20 it=1,nt
         DO 20 ita=1,nta(it)
         jrn=jrsm(ita,it)
   20    chdo(ita,it,1,1:jrn)=f1*chdo(ita,it,1,1:jrn)+
     .                    f2*chdn(ita,it,1,1:jrn)
      ELSEIF(mixkey.EQ.2) THEN
         IF(lin.EQ.1) qlmo=f1*qlmo+f2*qlmn
         smix=0.1d0
         sf1=1.d0-smix
         sf2=smix
         DO 30 it=1,nt
         DO 30 ita=1,nta(it)
         jrn=jrsm(ita,it)
C
         fi(1:jrn)=f1*(chdo(ita,it,1,1:jrn)+chdo(ita,it,2,1:jrn))+
     .             f2*(chdn(ita,it,1,1:jrn)+chdn(ita,it,2,1:jrn))
         fip(1:jrn)=sf1*(chdo(ita,it,1,1:jrn)-chdo(ita,it,2,1:jrn))+
     .              sf2*(chdn(ita,it,1,1:jrn)-chdn(ita,it,2,1:jrn))
         chdo(ita,it,1,1:jrn)=(fi(1:jrn)+fip(1:jrn))/2.d0
         chdo(ita,it,2,1:jrn)=(fi(1:jrn)-fip(1:jrn))/2.d0
   30    CONTINUE
C
      ELSE
         WRITE(m6,'(/,a,i3,a)') ' MIXCHD:** MIXKEY =',mixkey,
     1           ' not implemented'
         STOP
      ENDIF
C
      RETURN
      END
