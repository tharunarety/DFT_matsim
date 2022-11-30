C   ******************************************************************
C   *                                                                *
C   *    Modules for KSTR.                                           *
C   *                                                                *
C   ******************************************************************
      MODULE basis
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: qx3, qy3, qz3
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: qx, qy, qz
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: wi, wc, wst
        REAL(KIND=8), DIMENSION(3)              :: bsx, bsy, bsz
        REAL(KIND=8), DIMENSION(3)              :: bkx, bky, bkz
        REAL(KIND=8) :: ws, vol, boa, coa, alf, bet, gam
        INTEGER     , DIMENSION(:), ALLOCATABLE :: incl
        INTEGER      :: nq, nq3, npw, lat
      END MODULE
      MODULE bessel
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  :: hank, bess
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: fi, gi, gdot
      END MODULE
      MODULE control_data
        REAL(KIND=8)  :: dmax, kap2, dwats
        REAL(KIND=8)  :: err=1.d-10, errp=1.d-6
        INTEGER  :: nl, nlh, nlw, nder, itrans, nghbp
        INTEGER  :: lmax, nlm, lmaxh, nlmh, lmaxw, nlmw, nlmq
      END MODULE
      MODULE control_text
        CHARACTER*1  :: mode,store,high,wats
        CHARACTER*10 :: job
        CHARACTER*15 :: version
        CHARACTER*69 :: txt
        CHARACTER*60 :: for001, for002, for006
      END MODULE
      MODULE csts
        INTEGER, PARAMETER :: float = SELECTED_REAL_KIND(15,307)
        REAL(float) :: one = 1.0_float
        REAL(KIND=8) :: pi, twopi, fourpi, sqrtpi, sqrt2, sqfpi, r2h
      END MODULE
      MODULE factorial
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: fac2, efac, sig
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE  :: fack, facj, edot
        INTEGER, DIMENSION(:,:), ALLOCATABLE  :: ifib
      END MODULE
      MODULE lattice
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: rx, ry, rz
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: tx, ty, tz, rwats
        INTEGER,      DIMENSION(:), ALLOCATABLE  :: nriq, iqbas, jqbas
        INTEGER,      DIMENSION(:), ALLOCATABLE  :: nsh
        INTEGER                                  :: nr
      END MODULE
      MODULE lmtolm
        INTEGER,      DIMENSION(:), ALLOCATABLE  :: llx,mmx
      END MODULE
      MODULE madelung_matrix
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: vmadl
        REAL(KIND=8)                              :: cmdl
      END MODULE
      MODULE madelung_lattice
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: asx, asy, asz, dr
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: akx, aky, akz, dg
        REAL(KIND=8)  :: alamda, amax, bmax, rmax, gmax
        INTEGER       :: numvr, nr0, numvg
      END MODULE
      MODULE message
        INTEGER :: msgl, msgio, m6, nprn
      END MODULE
      MODULE realgaunt
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: gnt, gntw, gnth
        INTEGER     , DIMENSION(:), ALLOCATABLE :: lmg, lmpg, lmppg
        INTEGER     , DIMENSION(:), ALLOCATABLE :: lmgw, lmpgw, lmppgw
        INTEGER     , DIMENSION(:), ALLOCATABLE :: lmgh, lmpgh, lmppgh
        INTEGER     , DIMENSION(:), ALLOCATABLE :: kpow, jpow
        INTEGER     :: ngaunt, ngauntw, ngaunth
      END MODULE
      MODULE realharmonic
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE  :: ylm
      END MODULE
      MODULE screening
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: alwats
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: tmat, bigd, cd
        REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: dfac, sigma
      END MODULE
      MODULE temporary
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: sa, siq, sah
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: sbare, sbd
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: sb
        REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: work, sc, dett
        INTEGER, DIMENSION(:),   ALLOCATABLE :: ipvt
      END MODULE
      MODULE voronoi
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: xn,yn,zn,dn,dsc
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: xsc
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: inn
        INTEGER, DIMENSION(:), ALLOCATABLE :: iqn,nn
        INTEGER   :: nsc, nvn
      END MODULE
