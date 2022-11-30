C   ******************************************************************
C   *                                                                *
C   *    Modules for LIBRARY.                                        *
C   *                                                                *
C   ******************************************************************
      MODULE bessl_fact
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: facbl
        REAL(KIND=8), PARAMETER :: tol = 1.d-15
        INTEGER :: lmaxi,lmini
      END MODULE
      MODULE clebsch_gordon
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: h
        INTEGER,      PARAMETER ::  mclb=399
        INTEGER,      DIMENSION(:), ALLOCATABLE :: j
      END MODULE
      MODULE realhr_norm
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ylnorm
        REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: al, alinv, tlp1
        REAL(KIND=8), PARAMETER :: small=1.d-12
        REAL(KIND=8), PARAMETER :: sqr2 = 1.41421356237309515d0
        INTEGER :: lmaxi
      END MODULE
