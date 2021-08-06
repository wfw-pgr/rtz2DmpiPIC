module constants
  implicit none
  ! ---   Job    --- !
  character(20)               :: job           = 'sFRC02_'
  character(2)                :: coordinate    = 'RZ'
  character(3)                :: bctype        = 'pcw'
  character(3)                :: solver        = 'BCG'
  character(3)                :: equtype       = 'FRC'   ! ---  Sph or STk or FRC --- !
  integer                     :: psisw         = 3
  character(5)                :: CaseIO        = 'CaseO' ! CaseI / CaseO
  character(3)                :: CHSM          = 'ctr'   !  ctr  / co-
  character(3)                :: normType      = 'PIC'   ! 'MHD' / 'PIC' !
  character(3)                :: normInGS      = 'PIC'   ! 'MHD' / 'PIC' !
  character(100)              :: jobdir
  logical         , parameter :: normSW        = .false. ! Normalization !
  logical         , parameter :: RgridMode     = .false.
  logical         , parameter :: BgridMode     = .true.
  
  ! --- Parallel --- !
  integer         , parameter :: OMPNumThreads = 4
  integer         , parameter :: PEtot         = 16
  integer         , parameter :: PEpic         = 16
  integer                     :: myRank
  
  ! ---   Grid   --- !
  integer         , parameter :: Nr            = 512
  integer         , parameter :: Nz            = 1024     ! - Nr, Nz must be PEtot*Segment - !
  integer         , parameter :: Nmid          = 50
  double precision            :: rMin          = + 0.d0
  double precision            :: rMax          = + 2.d0
  double precision            :: zMin          = - 2.d0
  double precision            :: zMax          = + 2.d0
  
  ! ---  PIC Parameters          --- !
  double precision, parameter :: vthcv         = 0.08d0
  double precision, parameter :: wpewce        = 5.0d0
  double precision, parameter :: TiTe          = 1.0d0
  double precision, parameter :: mr            = 5.0d1
  double precision, parameter :: dr_Debye      = 1.2d0
  double precision, parameter :: dz_Debye      = 1.2d0
  double precision, parameter :: valfe         = 1.0d0 / wpewce

  ! ---  Iteration Settings      --- !
  integer         , parameter :: itermax1      = 3
  integer         , parameter :: itermax2      = 3
  double precision, parameter :: convergence1  = 1.d-4
  double precision, parameter :: convergence2  = 1.d-4
  double precision, parameter :: picard_alphaB = 0.5d0
  
  ! ---  Configuration Settings  --- !
  double precision, parameter :: pTop          = 1.0d-2
  double precision, parameter :: Betap         = 1.00d0
  double precision, parameter :: Betat         = 0.00d0
  double precision, parameter :: AspectR       = 0.0d0
  double precision, parameter :: Bmax          = 1.0d0
  double precision, parameter :: rhofloor      = 0.10d0
  
  ! ---  FluxFunc. g(psi),p(psi) --- !
  double precision, parameter :: epsilon       = 2.0d0
  double precision, parameter :: alpha         = 2.0d0
  double precision, parameter :: eta           = 0.15d0
  double precision, parameter :: c1            = 0.8d0 ! c1 + c2 = 1.d0 ( for p=p0 )
  double precision, parameter :: c2            = 0.2d0
  double precision, parameter :: beta_sep      = 0.2d0
  double precision, parameter :: Djh           = 0.20d0
  double precision, parameter :: lambda_       = 1.20d0
  
  ! ---  Limiter Settings        --- !
  !   :: relative value ( e.g. limposr * rleng -> limposr )   --- !
  integer         , parameter :: Nlim          = 2
  double precision            :: limposr(Nlim)
  double precision            :: limposz(Nlim)
  integer                     :: limidx        = 1
  data limposr/ 0.800d0, 0.15d0 /
  data limposz/ 0.500d0, 0.50d0 /
  ! data limposr/ 0.70d0, 0.50d0, 0.50d0 /
  ! data limposz/ 0.50d0, 0.85d0, 0.15d0 /
  
  ! ---  EF Feedback Settings    --- !
  double precision            :: Bv0
  double precision, parameter :: alpha_star    = 0.40d0
  double precision            :: Bsensor_r1    = 0.35d0
  double precision            :: Bsensor_r2    = 0.65d0

end module constants


module variables
  
  use constants, only : Nr, Nz
  implicit none
  ! ---    Grid   --- !
  double precision          :: r(Nr), z(Nz), dr, dz, drinv, dzinv, rinv(Nr)
  double precision          :: raxis, zaxis, aMinor
  
  ! ---    Equilibrium Quantities  --- !
  double precision          :: psi(Nz,Nr)
  double precision          :: psipast(Nz,Nr,2)
  double precision          :: rhs(Nz,Nr), Jphi(Nz,Nr)
  double precision          :: prs(Nz,Nr), gfunc(Nz,Nr), rho(Nz,Nr)
  double precision          :: psilim, psimin
  double precision          :: Itotal, g0, alpha_g
  
  ! ---    Boundary Conditions     --- !
  double precision          :: BdrcLR(Nr,2), BdrcTB(Nz,2)
 
end module variables
