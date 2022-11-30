C   ******************************************************************
C   *                                                                *
C   *    Modules for BMDL.                                           *
C   *                                                                *
C   ******************************************************************
      module basis
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: qx3, qy3, qz3
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: qx, qy, qz
        REAL(KIND=8), DIMENSION(3) :: bsx, bsy, bsz
        REAL(KIND=8), DIMENSION(3) :: bkx, bky, bkz
        REAL(KIND=8) :: sws, vol, boa, coa, alf, bet, gam
        INTEGER :: nq, nq3
      end module
      module control_data
        INTEGER :: nl2, nl22, nlm2, nl, nlm, nlmq
      end module
      module control_text
        CHARACTER(LEN=69) :: txt
        CHARACTER(LEN=60) :: for001, for006
        CHARACTER(LEN=15) :: version
        CHARACTER(LEN=10) :: job
        CHARACTER(LEN=2) :: mode
      end module
      module csts
        INTEGER, PARAMETER :: float = SELECTED_REAL_KIND(15,307)
        REAL(float) :: one = 1.0_float
        REAL(KIND=8) :: pi, twopi, fourpi, sqrt2, sqrtpi, sqfpi, r2h
      end module
      module factorial
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fac, dac
        INTEGER :: nfctrl
      end module
      module gaunt_coeff
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: c
        INTEGER, DIMENSION(:), ALLOCATABLE :: ll, mm
      end module
      module madelung_lattice
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: asx, asy, asz, dr
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: akx, aky, akz, dg
        REAL(KIND=8) :: alamda, amax, bmax, rmax, gmax
        INTEGER :: numvr, nr0, numvg
      end module
      module madelung_matrix
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: vmd
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vmad
      end module
      module message
        INTEGER :: msgl, msgio, m6, nprn
      end module
      module s_matrix
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: srlmq, silmq
      end module
