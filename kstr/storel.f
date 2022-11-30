      SUBROUTINE storel(npa)
C   ******************************************************************
C   *                                                                *
C   *   Store lattice informations.                                  *
C   *                                                                *
C   ******************************************************************
      USE basis
      USE control_data
      USE control_text
      USE lattice
      USE madelung_matrix
      USE screening
      USE message
      IMPLICIT NONE
      INTEGER :: npa
C
      wst=wst*ws
      WRITE(1) txt,mode,'N'
      WRITE(1) nq,nl,nlm,nlmq,nr
      WRITE(1) nriq(0:nq),iqbas(1:nr),jqbas(1:nr)
      WRITE(1) npa,npw,incl(1:npw)
      WRITE(1) ws,wst(1:nq),wi(1:nq),wc(1:nq)
      WRITE(1) bsx,bsy,bsz,qx(1:nq),qy(1:nq),qz(1:nq)
      WRITE(1) bkx,bky,bkz,boa,coa,alf,bet,gam
      WRITE(1) tx(1:nr),ty(1:nr),tz(1:nr)
      WRITE(1) kap2,nder,itrans
      WRITE(1) sigma(0:lmax,1:nq)
      WRITE(1) tmat(1:4,0:lmax,1:nq,0:nder)
      WRITE(1) cmdl,vmadl(1:nq,1:nq)
      IF(HIGH.EQ.'Y') THEN
         WRITE(2) txt,mode,'Y'
         WRITE(2) nq,nl,nlm,nlmq,nr,nlh,nlmh,lat
         WRITE(2) nriq(0:nq),iqbas(1:nr),jqbas(1:nr)
         WRITE(2) npa,npw,incl(1:npw)
         WRITE(2) ws,wst(1:nq),wi(1:nq),wc(1:nq)
         WRITE(2) bsx,bsy,bsz,qx(1:nq),qy(1:nq),qz(1:nq)
         WRITE(2) bkx,bky,bkz,boa,coa,alf,bet,gam
         WRITE(2) tx(1:nr),ty(1:nr),tz(1:nr)
         WRITE(2) kap2,nder,itrans
         WRITE(2) sigma(0:lmax,1:nq)
         WRITE(2) tmat(1:4,0:lmax,1:nq,0:nder)
         WRITE(2) cmdl,vmadl(1:nq,1:nq)
      ENDIF
C
      DEALLOCATE(tx,ty,tz)
      DEALLOCATE(vmadl)
C
      RETURN
      END
