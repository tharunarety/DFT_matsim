      SUBROUTINE renorm
C   ******************************************************************
C   *                                                                *
C   *    Renormalize core and atomic electron densities.             *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE control_data
      USE density
      USE message
      USE moments
      USE pota
      USE potential
      USE radialmesh
      USE temporary
      USE text
      IMPLICIT NONE
      INTEGER :: it, ita, jrn, jws, jsr
      REAL(KIND=8) :: zws, cws, dc, dz, facz, facc
C
      ALLOCATE(cor(mnta,nt,dimr))
C
      WRITE(m6,101)
      DO 20 it=1,nt
      DO 20 ita=1,nta(it)
      jws=jwss(ita,it)
      jrn=jrsm(ita,it)
      jsr=jsrs(ita,it)
C
C     Integrate charge densities up to WS
C
      rhoa(1:jws)=atmc(ita,it,1:jws)*ri(1:jws,ita,it)
      CALL simpn(rhoa,dx,jws,zws)
      rhoa(1:jws)=core(ita,it,1:jws)*ri(1:jws,ita,it)
      CALL simpn(rhoa,dx,jws,cws)
C
      WRITE(m6,'(11x,2i3,2x,2e15.8)') it,ita,cws,zws
C
      dc=nz(ita,it)-eln(ita,it)-cws
      dz=nz(ita,it)-zws+qtr(ita,it)
      IF(ABS(dc).GT.0.5) THEN
         WRITE(m6,103)
         STOP
      ENDIF
C
C     Renormalize charge densities
C
      facc=dc*3.d0/ws(ita,it)**3
      cor(ita,it,1:jws)=core(ita,it,1:jws)+facc*ri2(1:jws,ita,it)
      cor(ita,it,jws+1:jrn)=0.d0
C
      IF(tpot.NE.'N') GO TO 20
C
      facz=dz*3.d0/ws(ita,it)**3
      IF(ns.EQ.1) THEN
         chdo(ita,it,1,1:jrn)=chdo(ita,it,1,1:jrn)+
     .                        facz*ri2(1:jrn,ita,it)
      ELSE
         facz=0.5d0*facz
         chdo(ita,it,1,1:jrn)=chdo(ita,it,1,1:jrn)+
     .                        facz*ri2(1:jrn,ita,it)
         chdo(ita,it,2,1:jrn)=chdo(ita,it,2,1:jrn)+
     .                        facz*ri2(1:jrn,ita,it)
      ENDIF
C
      IF(ns.EQ.1) THEN
         rhoa(1:jrn)=chdo(ita,it,1,1:jrn)*ri(1:jrn,ita,it)
      ELSE
         rhoa(1:jrn)=(chdo(ita,it,1,1:jrn)+chdo(ita,it,2,1:jrn))*
     .                ri(1:jrn,ita,it)
      ENDIF
      CALL simpn(rhoa,dx,jws,zws)
      qs(ita,it)=zws-nz(ita,it)
   20 CONTINUE
C
  101 FORMAT(/,' RENORM:   IT',9x,'CWS',13x,'ZWS')
  102 FORMAT(11x,i2,2x,2e15.8)
  103 FORMAT(/,' RENORM:** Core contains wrong number of electrons')
      RETURN
      END
