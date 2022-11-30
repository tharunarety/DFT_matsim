      SUBROUTINE rotchd
C   ********************************************************************
C   *                                                                  *
C   *  Calculate the total charge density from the n^IBZ(r) as         *
C   *     n_L(r) = U(T)_LL'*n^IBZ_L'(r) , where  U(T)_LL' is the       *
C   *  rotation unitary matrix correspondent to symmerty operator T.   *
C   *                                                                  *
C   ********************************************************************
      USE atomicdens
      USE control_data
      USE density
      USE message
      USE radialmesh
      USE symmetry
      IMPLICIT NONE
      REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: chdl0
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: smm
      REAL(KIND=8) :: wq
      INTEGER :: iq, is, it, ita, jq, jt
      INTEGER :: l, m, mp, lm, lmp, irot, iqrot, jrc
C
      ALLOCATE(chdl0(mnta,nq,ns,nlmf,dimr),smm(dimr))
C
      WRITE(m6,100)
      IF(msgl.NE.0) WRITE(msgio,100)
C
C     Initialize a temporary block
C
      chdl0=0.d0
C
C     Loop for spin and site
C
      DO 20 is=1,ns
      DO 20 iq=1,nq
      it=itq(iq)
      DO 20 ita=1,nta(it)
      wq=wqst(iq)
      jrc=jwsc(ita,it)
C
C     Loop for radial mesh and lm
C
      DO 21 l=0,lmaxf
      DO 21 m=-l,l
      lm=l*l+l+m+1
C
C     Sum for symmetry operators
C
      DO 22 irot=1,nrot
      IF(ibzr(irot,iq).NE.0) THEN
         smm(1:jrc)=0.d0
C
C        Sum for l'm'
C
         DO 23 mp=-l,l
         lmp=l*l+l+mp+1
   23    smm(1:jrc)=smm(1:jrc)+
     .              ugam(irot,l,m,mp)*chdl(ita,iq,is,lmp,1:jrc)
C
C        If symmetry element IROT transform IQ to IQ' rearrange
C        correspondent partial densities
C
         iqrot=iprmt(irot,iq)
         chdl0(ita,iqrot,is,lm,1:jrc)=chdl0(ita,iqrot,is,lm,1:jrc)+
     .                                smm(1:jrc)*wq
      ENDIF
   22 CONTINUE
   21 CONTINUE
C
   20 CONTINUE
C
      chdl=chdl0
      DEALLOCATE(chdl0,smm)
C
  100 FORMAT(/,' ROTCHD:   Performe the 3D rotations')
      RETURN
      END
