      SUBROUTINE scaws
C   ******************************************************************
C   *                                                                *
C   *    Scale atomic radii to the actual volume given by the        *
C   *    average Wigner-Seitz raius.                                 *
C   *                                                                *
C   ******************************************************************
      USE atomicdens
      USE control_data
      USE control_text
      USE lattice
      USE message
      USE radialmesh
      USE slope
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(mnta,nt) :: ws0
      REAL(KIND=8) :: sumws, sumwst
      INTEGER :: iq, it, ita, i, l, m
C
      ALLOCATE(sigmt(0:lmax,nt))
      ws0=ws                  !! Wigner-Seitz radii in units of wst
C
      sumws=0.d0
      sumwst=0.d0
      alat=sws/wsa            !! lattice constant
      wst=wst*alat            !! Wigner-Seitz radius
      wsi=wsi*alat            !! inscribed sphere radius
      wsc=wsc*alat            !! circumscribed sphere radius
      DO iq=1,nq
         it=itq(iq)
         DO ita=1,nta(it)
         ws(ita,it)=ws0(ita,it)*wst(iq)
         sumws =sumws+conc(ita,it)*ws(ita,it)**3
         ENDDO
         sumwst=sumwst+wst(iq)**3
         sigmt(0:lmax,it)=sigma(0:lmax,iq)*alat !! hard sphere radius
      ENDDO
C
C     Renormalize ws(it)
C
      sumws=(sumws/sumwst)**(1.d0/3.d0)
      ws=ws/sumws
C
C     Potential sphere
C
      hsr=ws*hsr           !! potential sphere radius
C
C     Local muffin-tin zero radius
C
      wsm=wsm*hsr
C
      i=0
      DO iq=1,nq
      DO l=0,lmax
      DO m=-l,l
         i=i+1
         asc(i)=sigma(l,iq)*alat !! hard sphere radius
      ENDDO
      ENDDO
      ENDDO
C
      WRITE(m6,100)
      DO iq=1,nq
         it=itq(iq)
         DO ita=1,nta(it)
         IF(lclmff.EQ.1.AND.wsm(ita,it).GT.hsr(ita,it)) 
     .      localmt(ita,it)=1
         WRITE(m6,101) iq,it,ita,ws(ita,it),hsr(ita,it),wsm(ita,it),
     .                 wsi(iq),wsc(iq),sigmt(0:lmax,it)
         ENDDO
      ENDDO
C
  100 FORMAT(/,' SCAWS:    IQ IT ITA  WS',5x,'S',6x,'Sm',5x,'Si',5x,
     .         'Sc',5x,'hard spheres (s,p,...)',/)
  101 FORMAT(10x,3i3,11f7.4)
      RETURN
      END
