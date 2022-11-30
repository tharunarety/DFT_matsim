C   ******************************************************************
C   *                                                                *
C   *    Modules for KGRN.                                           *
C   *                                                                *
C   ******************************************************************
      MODULE atomb
        REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: hwst, hwsi, hwsc
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: hws, hhsr, hqtr
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: hdexch
        REAL(KIND=8)  :: hsws, hefgs
        INTEGER,      DIMENSION(:),     ALLOCATABLE :: kitq, knta
        INTEGER       :: knt, knq, knl, kns
      END MODULE atomb
      MODULE atomdata
        REAL(KIND=8)  :: vmix, rwat, rmax, dx, dr1, test
        REAL(KIND=8)  :: teste, testy, testv
        INTEGER       :: iex, np, nes, niter, iwat, nprna
      END MODULE atomdata
      MODULE atomicdens
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: atmc, core
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: eln, conc
        INTEGER,      DIMENSION(:,:), ALLOCATABLE :: nodes
        INTEGER,      DIMENSION(:,:), ALLOCATABLE :: nz, ion
      END MODULE atomicdens
      MODULE botomtop
        REAL(KIND=8)  :: ebt, etp
      END MODULE botomtop
      MODULE bzmesh
        REAL(KIND=8), DIMENSION(3)              :: tkx, tky, tkz
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fkx, fky, fkz
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fsx, fsy, fsz
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fkw, ww, fsw
        REAL(KIND=8)  :: dkx, dky, dkz, dhx, boa, coa, alf, bet, gam
        INTEGER,      DIMENSION(:), ALLOCATABLE :: kx, ky, kz
        INTEGER       :: ibz, nkx, nky, nkz, ibz2, nkx2, nky2, nkz2
        INTEGER       :: nfs, fsts, fs
      END MODULE bzmesh
      MODULE bzmesh2
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: wk, akx, aky, akz
        REAL(KIND=8)  :: dkz2, dkp 
        INTEGER       :: npar, nper
      END MODULE bzmesh2
      MODULE density
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: grnfm
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: chdl, chd0
        REAL(KIND=8), DIMENSION(:,:,:,:),   ALLOCATABLE :: chdo, chdn
        REAL(KIND=8), DIMENSION(:,:,:,:),   ALLOCATABLE :: chdh, chde
        REAL(KIND=8), DIMENSION(:,:,:),     ALLOCATABLE :: chdr
        REAL(KIND=8), DIMENSION(:,:),       ALLOCATABLE :: qti
        REAL(KIND=8), DIMENSION(:,:),       ALLOCATABLE :: qtr, qtro
        REAL(KIND=8), DIMENSION(:,:),       ALLOCATABLE :: zn0
        REAL(KIND=8)  :: znt
      END MODULE density
      MODULE control_data
        COMPLEX(KIND=8) :: zero = (0.d0,0.d0), zone = (1.d0,0.d0)
        REAL(KIND=8)    :: delta, depth, imagz, elim, tole, tolef
        REAL(KIND=8)    :: eb, elt, exchf, amix, efmix, spinfc, meps
        REAL(KIND=8)    :: alphmd, tolcpa
        INTEGER,        DIMENSION(:), ALLOCATABLE :: itq, ityp, mmt, nta
        INTEGER,        DIMENSION(:,:), ALLOCATABLE :: vac
        INTEGER,        DIMENSION(:,:), ALLOCATABLE :: localmt
        INTEGER         :: iter, niter, nkvec, ixc, idebug
        INTEGER         :: ns, nt, mnta, nq, pan, fixvmtz
        INTEGER         :: nl, lmax, nlm, nl2, lmax2, nlm2, nlmq, lat
        INTEGER         :: lmaxh, nlmh, nlmqh, lmaxt, nlmt, nlmqt
        INTEGER         :: lmaxf, nlmf, nlmqf
        INTEGER         :: diml, dimecr = 4
        INTEGER         :: ncpa, icpa, lclmff
      END MODULE control_data
      MODULE control_text
        CHARACTER(LEN=1)   :: strt, afm, zmsh, pole, ops, dosc, fllbz
        CHARACTER(LEN=1)   :: fcd, frc, crt, softc, expan, stmp, empt
        CHARACTER(LEN=1)   :: cpa, oexpan
        CHARACTER(LEN=2)   :: mode
        CHARACTER(LEN=3)   :: func
        CHARACTER(LEN=40)  :: job
        CHARACTER(LEN=15)  :: version
        CHARACTER(LEN=80)  :: for001, for001p, for002, for003, for004
        CHARACTER(LEN=80)  :: for006, for009, for010, for011
      END MODULE control_text
      MODULE csts
        INTEGER, PARAMETER :: float = SELECTED_REAL_KIND(15,307)
        REAL(float) :: one = 1.0_float
        REAL(KIND=8) :: pi, twopi, fourpi, sqrtpi, sqrt2, sqfpi, r2h
      END MODULE csts
      MODULE diracparam
        INTEGER :: srde
        REAL(KIND=8)  :: test, clight, csq, expdxh, dxd8
        REAL(KIND=8)  :: a1, a2, a3, a4, a5, a6, fac3, fac4
      END MODULE diracparam
      MODULE dosmom
        REAL(KIND=8), DIMENSION(:,:,:,:),  ALLOCATABLE :: tnos, emom
        REAL(KIND=8), DIMENSION(:,:,:,:),  ALLOCATABLE :: tdos, entr
        REAL(KIND=8), DIMENSION(:,:,:,:),  ALLOCATABLE :: tnos0, emom0
        REAL(KIND=8), DIMENSION(:,:),    ALLOCATABLE :: amag
        REAL(KIND=8)  :: tmag, tmago, tmagoo
      END MODULE dosmom
      MODULE energymesh
        COMPLEX(KIND=8), DIMENSION(:),      ALLOCATABLE :: zm, wgm, zx
        COMPLEX(KIND=8), DIMENSION(:),      ALLOCATABLE :: zfcd, wfcd
        COMPLEX(KIND=8), DIMENSION(:),      ALLOCATABLE :: zm0, wgm0
        REAL(KIND=8)  :: eps, tfermi, hx
        INTEGER       :: nzm, nz1, nz2, nz3, nres, nx, nx0, nz0, nzd
        INTEGER       :: nfcd
      END MODULE energymesh
      MODULE fandgeq
        REAL(KIND=8),    DIMENSION(:,:,:), ALLOCATABLE :: tov, okac
        INTEGER,         DIMENSION(:,:),   ALLOCATABLE :: nov
        INTEGER,         DIMENSION(:,:,:), ALLOCATABLE :: type, wov
        INTEGER,         DIMENSION(:,:,:), ALLOCATABLE :: jrov, ovt
        REAL(KIND=8)   :: vols, voli
        INTEGER        :: dimov
      END MODULE fandgeq
      MODULE force
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: ylm, gt, gp
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: fsy, fsz, fsx
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: fny, fnz, fnx
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: fey, fez, fex
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: fiy, fiz, fix
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: fxy, fxz, fxx
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: fcy, fcz, fcx
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: fqy, fqz, fqx
      END MODULE force
      MODULE gaussi
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: thvec, fivec
        REAL(KIND=8),    DIMENSION(:),      ALLOCATABLE :: wth, wfi
        INTEGER       :: nth, nfi
      END MODULE gaussi
      MODULE greenfunc
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: gah, gd
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: ga, gax
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: tnosx
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: gi
        COMPLEX(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: hgh, hghx
        COMPLEX(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: bg0, bg0x
        COMPLEX(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: dos
        COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: bgdos
      END MODULE greenfunc
      MODULE kinkmatrix
        COMPLEX(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: kinkm, kinkmd
        COMPLEX(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: gak, unit
        COMPLEX(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: unil
        REAL(KIND=8),    DIMENSION(:),     ALLOCATABLE :: work, worl
        INTEGER :: info = 0
      END MODULE kinkmatrix
      MODULE lattice
        REAL(KIND=8), DIMENSION(3) :: bsx, bsy, bsz
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: qx, qy, qz
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: tx, ty, tz
        INTEGER,      DIMENSION(:), ALLOCATABLE :: nviq, iqbas, jqbas
        INTEGER       :: nv
      END MODULE lattice
      MODULE logderivative
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: dfi
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: dtilz
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: dtilx
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: gsng
        REAL(KIND=8),    DIMENSION(:,:,:),     ALLOCATABLE :: sgm
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: fi0m
      END MODULE logderivative
      MODULE message
        INTEGER       :: idsyst, iasc, msgl, msgio, m6, nprn
      END MODULE message
      MODULE moments
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: vmad
        REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: qlmn, qlmo
        REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: qlm0, qs, potmc
        REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: qcpa
        REAL(KIND=8), DIMENSION(:),       ALLOCATABLE :: qsca
        REAL(KIND=8), DIMENSION(:,:),     ALLOCATABLE :: madc
        INTEGER :: nlmmad
      END MODULE moments
      MODULE partialwaves
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: cfm
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: cf
        COMPLEX(KIND=8), DIMENSION(:),         ALLOCATABLE :: p, q
        REAL(KIND=8),  DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: cfr, cfrm
        REAL(KIND=8),    DIMENSION(:,:,:,:,:), ALLOCATABLE :: cfrr, cfrt
        REAL(KIND=8),    DIMENSION(:,:,:,:,:), ALLOCATABLE :: ecr, nocr
        INTEGER,         DIMENSION(:,:,:,:),   ALLOCATABLE :: necr
      END MODULE partialwaves
      MODULE poissonparam
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: e, f
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: a, c2, f1
        REAL(KIND=8)  :: c, edl
      END MODULE poissonparam
      MODULE pota
        CHARACTER(LEN=1), DIMENSION(:,:), ALLOCATABLE :: fixst
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: cor
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: split
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: splito, splitoo
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: fixsjta
      END MODULE pota
      MODULE potential
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: v, vp
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: fullp
        REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: vmdl, potm
        REAL(KIND=8), DIMENSION(:,:,:),   ALLOCATABLE :: pots, potw
        REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: vmtz
        REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: vmtzr
        REAL(KIND=8), DIMENSION(:,:),   ALLOCATABLE :: dexch, dexcho
        REAL(KIND=8)  :: mmom, cmdl, vmtz0
      END MODULE potential
      MODULE potparam
        REAL(KIND=8), DIMENSION(:),     ALLOCATABLE :: rp, rq
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: eny, omm
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: wndw, cc
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: vl, tl, bot
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: sfim, signfi
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: top, dny
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: d1, d2, d3
        REAL(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: fmofp, amy
        REAL(KIND=8), DIMENSION(:,:,:,:),   ALLOCATABLE :: phil, logdl
        REAL(KIND=8), DIMENSION(:,:,:,:),   ALLOCATABLE :: epmatrix
      END MODULE potparam
      MODULE radialmesh
        REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: wst, wsi, wsc
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: wsm
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: hsr, ws
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: r1, sigmt
        REAL(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: ri, ri2
        REAL(KIND=8)  :: sws, wsa, alat, dx, vol
        INTEGER,      DIMENSION(:,:), ALLOCATABLE :: jwss, jris
        INTEGER,      DIMENSION(:,:), ALLOCATABLE :: jrsm, jsrs
        INTEGER,      DIMENSION(:,:), ALLOCATABLE :: jwsi, jwsc
        INTEGER       :: dimr, shf
      END MODULE radialmesh
      MODULE realgaunt
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: gnt
        INTEGER,  DIMENSION(:,:,:), ALLOCATABLE :: ignt
        INTEGER,  DIMENSION(:,:),   ALLOCATABLE :: ngnt
        INTEGER     , DIMENSION(:), ALLOCATABLE :: lmpg, lmg, lmppg
        INTEGER     , DIMENSION(:), ALLOCATABLE :: lppg
        INTEGER       :: ngaunt
      END MODULE realgaunt
      MODULE slope
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: tsz
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:,:), ALLOCATABLE :: tsx
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: tmz
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: tmx
        COMPLEX(KIND=8), DIMENSION(:,:,:,:,:), ALLOCATABLE :: taz, tax
        COMPLEX(KIND=8),  DIMENSION(:,:,:),   ALLOCATABLE :: sa, sal
        COMPLEX(KIND=8),  DIMENSION(:,:,:),   ALLOCATABLE :: sap, sapl
        COMPLEX(KIND=8),  DIMENSION(:,:),     ALLOCATABLE :: slop, tmp
        COMPLEX(KIND=8),  DIMENSION(:),       ALLOCATABLE :: sa0, sap0
        REAL(KIND=8),     DIMENSION(:,:,:,:), ALLOCATABLE :: salpl, tmat
        REAL(KIND=8),     DIMENSION(:,:,:,:), ALLOCATABLE :: salplp
        REAL(KIND=8),     DIMENSION(:,:),     ALLOCATABLE :: sigma
        REAL(KIND=8),     DIMENSION(:),       ALLOCATABLE :: asc
        INTEGER           :: itrans
      END MODULE slope
      MODULE softcore
        CHARACTER(LEN=4),  DIMENSION(:,:), ALLOCATABLE :: symbols
        CHARACTER(LEN=24), DIMENSION(:,:), ALLOCATABLE :: configs
        REAL(KIND=8),      DIMENSION(:,:,:,:), ALLOCATABLE :: dens, dq1s
        REAL(KIND=8),      DIMENSION(:,:),   ALLOCATABLE :: eonec
        INTEGER,      DIMENSION(:,:,:), ALLOCATABLE :: nqns, nks, nels
        INTEGER,      DIMENSION(:,:,:), ALLOCATABLE :: ncorbs
        INTEGER,      DIMENSION(:,:), ALLOCATABLE :: izs, norbs, ions
        INTEGER :: softz
      END MODULE softcore
      MODULE symmetry
        REAL(KIND=8), DIMENSION(:,:,:,:), ALLOCATABLE :: ugam
        REAL(KIND=8), DIMENSION(:),   ALLOCATABLE :: wqst
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: ibzr, iprmt
        INTEGER, DIMENSION(:),   ALLOCATABLE :: iwrot, ibzrot, imsymm
        INTEGER       :: nrot, nprprt, incribz
      END MODULE symmetry
      MODULE text
        CHARACTER(LEN=63), DIMENSION(:), ALLOCATABLE :: tatmc
        CHARACTER(LEN=69), DIMENSION(:,:), ALLOCATABLE :: txtp
        CHARACTER(LEN=4),  DIMENSION(:,:), ALLOCATABLE :: ttxt
        CHARACTER*1,  DIMENSION(6) :: txtl = (/'s','p','d','f','g','h'/)
        CHARACTER*6,  DIMENSION(0:5,2) :: thead
        CHARACTER(LEN=1)   :: tpot, conv
        CHARACTER(LEN=3)   :: txch
        CHARACTER(LEN=5)   :: clock
        CHARACTER(LEN=9)   :: dato
        CHARACTER(LEN=69)  :: txts, txt
      END MODULE text
      MODULE taylor
        COMPLEX(KIND=8), DIMENSION(:,:,:), ALLOCATABLE :: tayl, tayld
        REAL(KIND=8)  :: kw20, kw20p
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: facd
        INTEGER,      DIMENSION(:), ALLOCATABLE :: clkw0
        INTEGER       :: nder, nderp, mder
      END MODULE taylor
      MODULE temporary
        REAL(KIND=8), DIMENSION(2) :: rhod, rhodd
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhosp1, rhosp2
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhxcd1, rhxcdd1
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhxcd2, rhxcdd2
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: fi, fip, wc, rhoa
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: rhov, rhop, rhopp
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: chdt
      END MODULE temporary
      MODULE totalenergy
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: eone, emadl, vint
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: enuc, ecor, eval
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: ents
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: exct, excc, ekin
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: etot, okae, vintc
        REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: etotcore
      END MODULE totalenergy
      MODULE volume
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: swsl, efl, etol
        INTEGER       :: nsws
      END MODULE volume
