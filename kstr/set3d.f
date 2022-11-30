      SUBROUTINE SET3D
C   ******************************************************************
C   *                                                                *
C   *    Initialize the basis vectors (QX,QY,QZ) and the number NQ   *
C   *    of atoms in the primitive cell.                             *
C   *                                                                *
C   *   *On entry (QX3,QY3,QZ3) contain the NQ3 basis vectors read   *
C   *    from input.                                                 *
C   *                                                                *
C   *   *On exit (QX,QY,QZ) contain the NQ=NQ3 basis vectors.        *
C   *                                                                *
C   ******************************************************************
      use basis
      use control_data
      use message
      IMPLICIT NONE
      INTEGER :: iq
c
      ALLOCATE(qx(nq3),qy(nq3),qz(nq3))
c
      DO 20 IQ=1,NQ3
      QX(IQ)=QX3(IQ)
      QY(IQ)=QY3(IQ)
      QZ(IQ)=QZ3(IQ)
   20 CONTINUE
      NQ=NQ3
      NLMQ=NLM*NQ
      WRITE(M6,101)
c
      DEALLOCATE(qx3,qy3,qz3)
c
  101 FORMAT(/,' SET3D:',4X,'Basis vectors in 3D cell are set.')
      RETURN
      END
