      SUBROUTINE fullbz(x,y,z,h,xs,ys,zs,ws,nf)
C   ******************************************************************
C   *                                                                *
C   * 1.Increase the irreducibile Brillouin zone if it is needed.    *
C   *                                                                *
C   * 2.Set up the star of the k=(x,y,z) vector.                     *
C   *                                                                *
C   *   This subroutine uses the symmetry operators generated        *
C   *   in ROTM3D.                                                   *
C   *                                                                *
C   ******************************************************************
      USE bzmesh
      USE control_text
      USE message
      USE symmetry
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(48) :: xs, ys, zs, ws
      REAL(KIND=8), PARAMETER :: err = 1.d-6
      REAL(KIND=8)  :: x, y, z, h
      INTEGER       :: nf, i, irot, itr
C
      IF(fllbz.NE.'Y') THEN
C
C        Irreducible Brillouin zone
C
         i=1
         ys(i)=y
         zs(i)=z
         xs(i)=x
C
         DO 20 irot=2,nrot
         IF(imsymm(irot).NE.0) THEN
            i=i+1
            ys(i)=y*ugam(irot,1,-1,-1)+z*ugam(irot,1,-1,0)+
     .            x*ugam(irot,1,-1, 1)
            zs(i)=y*ugam(irot,1, 0,-1)+z*ugam(irot,1, 0,0)+
     .            x*ugam(irot,1, 0, 1)
            xs(i)=y*ugam(irot,1, 1,-1)+z*ugam(irot,1, 1,0)+
     .            x*ugam(irot,1, 1, 1)
         ENDIF
   20    CONTINUE
         nf=i
         ws(1:nf)=h/nf
      ELSEIF(fllbz.EQ.'Y') THEN
C
C        Full Brillouin zone
C
         i=0
         DO 21 irot=1,nrot
         IF(ibzrot(irot).NE.0) THEN
            i=i+1
            ys(i)=y*ugam(irot,1,-1,-1)+z*ugam(irot,1,-1,0)+
     .            x*ugam(irot,1,-1, 1)
            zs(i)=y*ugam(irot,1, 0,-1)+z*ugam(irot,1, 0,0)+
     .            x*ugam(irot,1, 0, 1)
            xs(i)=y*ugam(irot,1, 1,-1)+z*ugam(irot,1, 1,0)+
     .            x*ugam(irot,1, 1, 1)
         ENDIF
   21    CONTINUE
         nf=i
         ws(1:nf)=h/nf
      ENDIF
C
      RETURN
      END
