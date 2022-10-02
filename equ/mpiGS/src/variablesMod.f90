module variablesMod
  implicit none

  ! ====================================================== !
  ! ===  constants                                     === !
  ! ====================================================== !
  
  ! ------------------------------------------------------ !
  ! --- [1] system parameters                          --- !
  ! ------------------------------------------------------ !
  integer         , parameter   :: cLen          = 300
  integer         , parameter   :: lun           = 50
  integer         , parameter   :: PEtot         = 2
  character(cLen) , parameter   :: parameterFile = "dat/mpiGS.lst"
  double precision, parameter   :: pi = 4.d0 * atan( 1.d0 )

  ! ------------------------------------------------------ !
  ! --- [2] abbreviations                              --- !
  ! ------------------------------------------------------ !
  integer         , parameter   :: z_  = 1
  integer         , parameter   :: r_  = 2
  integer         , parameter   :: lft_= 1, rgt_= 2
  integer         , parameter   :: top_= 1, bot_= 2


  ! ------------------------------------------------------ !
  ! --- [2] simulation settings                        --- !
  ! ------------------------------------------------------ !
  character(cLen)               :: job
  character(cLen)               :: solverType, equiliType
  character(cLen)               :: parityType, caseIOType
  character(cLen)               :: normalizationType
  character(cLen)               :: normalizationInGS
  character(cLen)               :: coordinateType, boundaryType, gridType
  logical                       :: normalizeField

  ! ------------------------------------------------------ !
  ! --- [3] coordinate settings                        --- !
  ! ------------------------------------------------------ !
  integer                       :: LI, LJ, Lmid
  double precision              :: zMin_mhd, zMax_mhd, rMin_mhd, rMax_mhd
  
  ! ------------------------------------------------------ !
  ! --- [4] PIC simulation parameters                  --- !
  ! ------------------------------------------------------ !
  double precision              :: vthcv, wpewce, TiTe
  double precision              :: mr, dr_Debye, dz_Debye
  double precision              :: valfe

  ! ------------------------------------------------------ !
  ! --- [5] GS-solver control parameters               --- !
  ! ------------------------------------------------------ !
  integer                       :: iterMax1    , iterMax2
  double precision              :: convergence1, convergence2
  double precision              :: Picard_alphaB

  ! ------------------------------------------------------ !
  ! --- [6] Equilibrium parameters                     --- !
  ! ------------------------------------------------------ !
  integer                       :: psisw
  double precision              :: pTop, beta_pol, beta_tor, aspectRatio, Bmax, rhoFloor
  double precision              :: epsilon, alpha, eta, c1, c2, beta_sep, Djh, lambda
  double precision              :: alpha_star, Bsensor_r1, Bsensor_r2

  integer         , parameter   :: nLimiter_max  = 5
  double precision              :: limiter_zPos(nLimiter_max), limiter_rPos(nLimiter_max)

  ! ====================================================== !
  ! ===  variables                                     === !
  ! ====================================================== !

  ! ------------------------------------------------------ !
  ! --- [1] system variables                           --- !
  ! ------------------------------------------------------ !
  character(cLen)               :: jobdir
  integer                       :: OMPNumThreads, myRank, PEpic

  ! ------------------------------------------------------ !
  ! --- [2] Grid variables                             --- !
  ! ------------------------------------------------------ !
  double precision, allocatable :: zAxis(:), rAxis(:), rInv(:)
  double precision              :: dz, dr, dzInv, drInv
  double precision              :: rLeng, zLeng, zMin, zMax, rMin, rMax
  double precision              :: aMinor, magAxis(2)

  ! ------------------------------------------------------ !
  ! --- [3] Equilibrium variables                      --- !
  ! ------------------------------------------------------ !
  double precision, allocatable :: psi(:,:), psi_past(:,:,:), psis(:,:)
  double precision, allocatable :: rhs(:,:), current (:,:), pressure(:,:)
  double precision, allocatable :: rho(:,:), gfunc   (:,:)
  double precision, allocatable :: bdrcLR(:,:), bdrcTB(:,:)
  double precision              :: psiLim, psiMin, Itot, g0, alpha_g

  double precision              :: Bv0
  integer                       :: nLimiter, limiter_index


  
end module variablesMod
