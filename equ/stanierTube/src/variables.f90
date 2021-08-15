! =================================================================== !
! ===  constants Module  :: Constants & universal Settings        === !
! =================================================================== !
module constants
  implicit none
  ! ------------------------------------- !
  ! --- [1] Configuration             --- !
  ! ------------------------------------- !
  character(20)               :: job            =   'test01_'
  character(2)                :: coordinate     =   'RZ'
  character(4)                :: MergingType    =   'CoHx'  ! -- 'CoHx'/'CtrI'/'CtrO' -- !
  character(3)                :: normType       =   'PIC'
  character(3)                :: sftype         =   'CIC'
  logical                     :: normSW         =   .false.
  logical         , parameter :: RgridMode      =   .false.
  logical         , parameter :: BgridMode      =   .true.
  logical         , parameter :: Flag__ReverseBt=   .false.
  integer                     :: myRank
  integer                     :: PEtot
  integer         , parameter :: PEpic          =    4

  ! ------------------------------------- !
  ! --- [2] PIC Parameters            --- !
  ! ------------------------------------- !
  double precision, parameter :: vthcv          =   0.12d0
  double precision, parameter :: wpewce         =   5.0d0
  double precision, parameter :: TiTe           =   1.0d0
  double precision, parameter :: mr             =   25.0d0
  double precision, parameter :: dx1_Debye      =   1.3d0
  double precision, parameter :: dx2_Debye      =   1.3d0
  double precision, parameter :: valfe          =   1.0d0 / wpewce

  ! ------------------------------------- !
  ! --- [3] Equilibrium Settings      --- !
  ! ------------------------------------- !
  !  -- [3-1] Simulation Box size     --  !
  integer         , parameter :: N1             =    256
  integer         , parameter :: N2             =    128
  ! integer         , parameter :: N1             =    512
  ! integer         , parameter :: N2             =    256
  double precision            :: x1Min          =   -20.0d0
  double precision            :: x1Max          =   +20.0d0
  double precision            :: x2Min          =   + 0.0d0
  double precision            :: x2Max          =   +20.0d0
  !  -- [3-2] B-Field Strength        --  !
  double precision, parameter :: Bmax           =   0.9d0
  double precision, parameter :: Bt0            =   0.0d0
  double precision, parameter :: Bv0            = - 1.22d0
  !  -- (ref.val.) [LI,LJ]=[256,128] & jCenter=0.18 => Bv0=-0.16169
  
  ! ------------------------------------- !
  ! --- [4] Current ( FluxTube )      --- !
  ! ------------------------------------- !
  !  -- [4-1] Current Position rc     --  !
  logical         , parameter :: fixedRange     =   .false.
  double precision, parameter :: x1cnt          =   0.50d0
  double precision, parameter :: x2cnt          =   0.45d0
  double precision, parameter :: r_omega        =   0.30d0
  !  -- [4-2] Current Strength |j|    --  !
  double precision, parameter :: jCenter1       =   0.18d0
  double precision, parameter :: jCenter2       =   0.18d0
  double precision, parameter :: symmBreak      =   +0.0d0
  
  ! ------------------------------------- !
  ! --- [5] Pressure & rho            --- !
  ! ------------------------------------- !
  character(9)                :: BetaMode       =  "EnergyOfB"
  double precision, parameter :: desiredBetapol =   5.0d-2
  double precision, parameter :: desiredBeta    =   1.5d-1
  double precision, parameter :: rho0           =   1.0d0
  double precision, parameter :: rhofloor       =   0.2d0
  double precision, parameter :: pfloor         =   0.2d0

end module constants


! =================================================================== !
! ===  variables Module  :: Common variables                      === !
! =================================================================== !
module variables
  use constants, only :  N1, N2
  implicit none
  double precision   :: dx1, dx2
  double precision   :: x1(N1), x2(N2)
  double precision   :: Bfd(3,0:N1,0:N2), Jcr(3,0:N1,0:N2), uvc(3,0:N1,0:N2,2)
  double precision   :: rhs(N1,N2), Avp(N1,N2), pcf(N1,N2)
  double precision   :: prs(0:N1,0:N2), rho(0:N1,0:N2)
  double precision   :: psi(0:N1,0:N2), frp(0:N1,0:N2,3)
end module variables
