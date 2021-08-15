module constants
  implicit none
  ! ------------------------------------- !
  ! --- [1] Job / Directory Settings  --- !
  ! ------------------------------------- !
  character(20)                 :: job             = 'test01_01'
  character(20)                 :: savedjob        = 'test01_01'
  character(20)                 :: equJob          = 'test01_'
  character(20)                 :: equType         = 'stanierTube' !  "stanierTube" / "mpiGS"
  character(20)                 :: equGrid         = 'BgD'
  character(25)                 :: jobDir
  character(30)                 :: binDir, datDir, savDir
  character(50)                 :: equDir
  ! ------------------------------------- !
  ! --- [2] Parallelization Settings  --- !
  ! ------------------------------------- !
  integer, parameter            :: OMPNumThreads   = 2
  integer, parameter            :: PEtot           = 4
  integer                       :: myRank
  character(6)                  ::  cRank
  ! ------------------------------------- !
  ! --- [3] Sim.Size & #.Particles    --- !
  ! ------------------------------------- !
  !  -- [3-1] Sim.Size                --  !
  integer                       :: LIs, LJs
  integer, parameter            :: LI              = 256 + 1
  integer, parameter            :: LJ              = 128 + 1
  !  -- [3-2] #.of particles          --  !
  integer, parameter            :: ns              = 2
  ! integer, parameter            :: ppcMax          = 256
  integer, parameter            :: ppcMax          = 16
  ! integer, parameter            :: ppcMin          = 16
  integer, parameter            :: ppcMin          = 4
  integer, parameter            :: ppcMargin       = 1
  ! integer, parameter            :: ppcInit         = 128
  integer, parameter            :: ppcInit         = 8
  integer, parameter            :: nptTot          = (LI-2)*(LJ-2)*ppcMax
  integer                       :: nptMax          = nptTot / PEtot
  integer                       :: npt             = (LI-1)*(LJ-1)*ppcInit*ns/PEtot
  integer, parameter            :: APRperCell      = 1
  integer, parameter            :: APRvdfResl      = 1
  real   , parameter            :: ppcMinAlpha     = 1.0
  integer                       :: ppcMinThr(LJ)
  ! ------------------------------------- !
  ! --- [4] Physical  Parameters      --- !
  ! ------------------------------------- !
  double precision, parameter   :: mr              = 25.d0         !    mi / me
  double precision, parameter   :: isoThr          = 1.41421356d0  ! Tperp / Tpara  
  double precision, parameter   :: nuwce           = 1.d-6         ! Collision Frequency
  double precision, parameter   :: prtrb_coef      = 0.0d0
  double precision, parameter   :: InitJiJe        = 0.0d0
  ! ------------------------------------- !
  ! --- [5] Numerical Parameters      --- !
  ! ------------------------------------- !
  double precision, parameter   :: SmCoef          = 1.d-5
  double precision, parameter   :: ptExpandRatio   = 1.1d0
  double precision, parameter   :: alpha_wce       = 0.025d0
  double precision, parameter   :: alpha_wpe       = 0.025d0
  double precision, parameter   :: alpha_CFL       = 0.9d0
  ! ------------------------------------- !
  ! --- [6] kstep Intervals           --- !
  ! ------------------------------------- !
  integer, parameter            :: k_Enrgy         = 50
  integer, parameter            :: k_FlxIp         = 50
  integer, parameter            :: k_Field         = 1000
  integer, parameter            :: k_Momnt         = 1000
  integer, parameter            :: k_tShow         = 500
  integer, parameter            :: k_tChck         = 100
  integer, parameter            :: k_Ecrct         = 500
  integer, parameter            :: k_NumSP         = 10000
  integer, parameter            :: k_prtcl         = 10000000
  integer, parameter            :: k_SaveP         = 50000
  integer, parameter            :: k_ppcCt         = 10000000
  integer, parameter            :: k_vdfsm         = 20
  integer, parameter            :: k_pSort         = 200
  integer, parameter            :: k_chDDD         = 500
  ! ------------------------------------- !
  ! --- [7] Start & End Settings      --- !
  ! ------------------------------------- !
  integer                       :: kstart          =    0
  integer, parameter            :: itermax         = 1000
  logical, parameter            :: Flag__LoadSave  = .false.
  character(9)                  :: wallLimit       = "200:00:00"
  logical                       :: Flag__wallElapsed
  logical, parameter            :: Flag__InitialdivEfix = .false.
  ! ------------------------------------- !
  ! --- [8] DDD parameters            --- !
  ! ------------------------------------- !
  integer         , parameter   :: DDDnLine_Max    = 5
  double precision, parameter   :: DDDcritStep     = 0.05d0
  double precision, parameter   :: DDDcolThreshold = 0.2d0
  ! ------------------------------------- !
  ! --- [9] Boundary Conditions       --- !
  ! ------------------------------------- !
  character(10), parameter      :: Boundary1__EB   = 'cWall'    ! | 'cWall' | 'periodic' |
  character(10), parameter      :: Boundary2__EB   = 'cWall'
  character(10), parameter      :: Boundary1__Jc   = 'cWall'
  character(10), parameter      :: Boundary2__Jc   = 'cWall'
  character(10), parameter      :: Boundary1__pt   = 'reflect'  ! | 'reflect' | 'periodic' |
  character(10), parameter      :: Boundary2__pt   = 'reflect'
  ! ------------------------------------- !
  ! --- [10] Normalization            --- !
  ! ------------------------------------- !
  double precision, parameter   :: wce             = 1.d0
  double precision, parameter   :: cv              = 1.d0
  double precision, parameter   :: qeme            = 1.d0
  ! ------------------------------------- !
  ! --- [11] Delta                    --- !
  ! ------------------------------------- !
  double precision              :: dt, dr, dz, drinv, dzinv
  double precision              :: dtdr, dtdz, drdt, dzdt
  double precision              :: rMin, rMax, zMin, zMax
  ! ------------------------------------- !
  ! --- [12] Particle's species       --- !
  ! ------------------------------------- !
  integer                       :: np(ns), npc(ns)
  double precision              ::  q(ns),  qm(ns)
  double precision              :: wp(ns),  rm(ns)
  ! ------------------------------------- !
  ! --- [13] Load Variables           --- !
  ! ------------------------------------- !
  double precision              :: wpewce, vthcv, TiTe, dr_Debye, dz_Debye
  double precision              :: lDebye, vthi,  vthe, vAlfe,    rhofloor
  ! ------------------------------------- !
  ! --- [14] other constants          --- !
  ! ------------------------------------- !
  double precision, parameter   :: pi              =  4.d0 * atan( 1.d0 )
  ! ------------------------------------- !
  ! --- [15] Additional Mode          --- !
  ! ------------------------------------- !
  integer                       :: SmCnt           =  0
  logical, parameter            :: Filter__E       = .true.
  logical, parameter            :: Filter__B       = .true.
  logical, parameter            :: Filter__J       = .false.
  logical, parameter            :: divE__Fix       = .false.
  logical, parameter            :: Pcollider       = .false.
  logical, parameter            :: TrackMode       = .false.
  ! ------------------------------------- !
  ! --- [16] variable's indicator     --- !
  ! ------------------------------------- !
  integer, parameter            :: rp_ = 1, zp_ = 2
  integer, parameter            :: vr_ = 3, vt_ = 4, vz_ = 5
  integer, parameter            :: ro_ = 6, zo_ = 7, wp_ = 8
  integer, parameter            :: el_ = 1, io_ = 2
  integer, parameter            :: er_ = 1, et_ = 2, ez_ = 3
  integer, parameter            :: br_ = 4, bt_ = 5, bz_ = 6
  integer, parameter            :: jr_ = 1, jt_ = 2, jz_ = 3
  integer, parameter            :: rnk_= 1, frm_= 2, to_ = 3
  integer, parameter            :: LIs_= 4, LJs_= 5
  integer, parameter            :: iFr_= 6, iTo_= 7, inn_= 8
  integer, parameter            :: isl_= 9, iel_=10

end module constants


module variables
  use constants, only : LI, LJ, ns, PEtot
  implicit none
  ! ------------------------------------- !
  ! --- [1] General variables         --- !
  ! ------------------------------------- !
  integer                       :: kstep
  ! ------------------------------------- !
  ! --- [2] Geometry variables        --- !
  ! ------------------------------------- !
  double precision              :: x1Leng   , x1Lengloc, x2Leng, x2Lengloc, x1g(0:PEtot-1,3)
  double precision              :: zf(LI)   , zh(LI)
  double precision              :: rf(LJ+1) , rh(LJ+1) , rhinv(LJ+1), rfinv(LJ+1)
  double precision              :: gVf(LJ)  , rgVf(LJ) , rgVfInv(0:LJ+2)
  double precision              :: gVh(LJ)  , rgVh(LJ) , rgVhInv(0:LJ+2)
  double precision              :: wVf(LJ)  , rwVf(LJ) , rwVfInv(0:LJ+2)
  double precision              :: wVh(LJ)  , rwVh(LJ) , rwVhInv(0:LJ+2)
  double precision, allocatable :: x1s(:)
  ! ------------------------------------- !
  ! --- [3] System variables          --- !
  ! ------------------------------------- !
  integer                       :: pCtrlDij(2,2,PEtot)
  integer                       :: ijDomain(0:PEtot-1,10), npgPE(0:PEtot-1,5)
  double precision              :: neMax, neMin
  double precision              :: RhoPerNsp, NspRepres, volRepres
  double precision, allocatable :: sumupW(:,:,:,:)  , sumupE(:,:,:)
  double precision, allocatable :: sumupP(:,:,:,:,:)
  ! ------------------------------------- !
  ! --- [4] Physical variables        --- !
  ! ------------------------------------- !  
  ! ---  [4-1] Particle  --- !
  integer         , allocatable :: goban(:)      , pCellIdx(:,:)
  integer         , allocatable :: cellTop(:,:,:), perCell(:,:,:), perCellW(:,:,:,:)
  double precision, allocatable :: pxv(:,:,:)
  ! ---  [4-2] Field     --- !
  double precision, allocatable :: w_p(:,:,:)
  double precision, allocatable :: EBf(:,:,:), EBr(:,:,:), EBo(:,:,:)
  double precision, allocatable :: Jcr(:,:,:,:), JcrW(:,:,:,:,:)
  
end module variables
