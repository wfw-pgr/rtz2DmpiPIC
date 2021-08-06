program myGSchsmMPI
  ! ------------------------------------- !
  ! --- [1] Variables Loading         --- !
  ! ------------------------------------- !
  use constants,  only : itermax  , myRank   , kstart   , Flag__wallElapsed
  use constants,  only : k_prtcl  , k_Field  , k_SaveP  , k_Ecrct  , k_Momnt, k_chDDD
  use constants,  only : k_FlxIp  , k_tShow  , k_Enrgy  , k_ppcCt  , k_vdfsm, k_tChck, k_pSort
  use constants,  only : Filter__E, Filter__B, Filter__J, divE__Fix, Pcollider, Flag__InitialdivEfix
  use constants,  only : Boundary1__pt
  use variables,  only : kstep
  ! ------------------------------------- !
  ! --- [2] Routines  Loading         --- !
  ! ------------------------------------- !
  use ptAccelMod, only : ptPushMove
  use ptcBdrcMod, only : particleBoundary
  use ptOrderMod, only : ptReOrdering
  use collideMod, only : collision
  use ppcCtrlMod, only : ppcCtrl
  use currentMod, only : DensityDecomposition
  use maxwellMod, only : Bfield          , Efield
  use qchargeMod, only : charge          , Ecorrect
  use smoothfMod, only : smooth_E_Linear , smooth_B_Linear, smooth_J_Linear
  use writingMod, only : writeConst      , writePrtcl     , writeField
  use utilityMod, only : InitiateMPI     , MakeJobFile    , TimeMeasure
  use displayMod, only : OpeningRemarks  , EndingRemarks  , MainStartRemarks
  use displayMod, only : DisplayParameter, PE_Attendance
  use sEnergyMod, only : CalcEnergy      , Flow_Thermal_Energy
  use sFluxIpMod, only : FluxCurrentAnalysis
  use diffDDDMod, only : DiffuseDecomposition
  use makeVDFMod, only : VDFSample
  use momentsMod, only : MomentAnalysis
  use InitialMod, only : Initiator
  use saveRunMod, only : MakeSavePoint
  implicit none
  include 'mpif.h'
  integer :: ierr

  ! ------------------------------------- !
  ! --- [3]       PIC START           --- !
  ! ------------------------------------- !
  call InitiateMPI
  call OpeningRemarks
  call PE_Attendance
  call TimeMeasure( 0 )

  ! ------------------------------------- !
  ! --- [4]  Initialize Simulation    --- !
  ! ------------------------------------- !
  call MakeJobFile
  kstep = kstart
  call Initiator
  call writeConst
  call TimeMeasure( 1 )
  if ( Flag__InitialdivEfix ) then
     call charge
     if ( divE__Fix ) call Ecorrect
  endif
  call TimeMeasure( 2 )

  call DisplayParameter
  if ( kstep .ne.0 ) call DensityDecomposition
  call writeField
  ! call ppcCtrl( 'COM' )
  if ( myRank.eq.0 ) write(6,'(2x,a)') '[Pre-Loop Diagnosis]'
  call FluxCurrentAnalysis
  call MomentAnalysis
  call CalcEnergy
  call Flow_Thermal_Energy
  call TimeMeasure( 4 )

  ! ------------------------------------- !
  ! --- [5]      Main Loop Start      --- !
  ! ------------------------------------- !
  if ( myRank.eq.0 ) call MainStartRemarks
  kstart= kstep + 1
  do kstep=kstart, itermax
     
     ! --- [5-1]  B-Field  --- !
     call Bfield

     ! --- [5-2] Filtering --- !
     if ( Filter__E ) call smooth_E_Linear
     if ( Filter__B ) call smooth_B_Linear
     call TimeMeasure( 6 )

     ! --- [5-3] particle  --- !
     call ptPushMove
     call TimeMeasure( 5 )
     if ( Pcollider ) call collision
     call TimeMeasure( 7 )

     ! --- [5-4]  current  --- !
     call DensityDecomposition
     if ( Filter__J ) call smooth_J_Linear
     call TimeMeasure( 3 )

     ! --- [5-5] ptcl. B.C.--- !
     call particleBoundary
     call TimeMeasure( 5 )
     
     ! --- [5-6]  E-Field  --- !
     call Efield
     call TimeMeasure( 6 )

     ! --- [5-7]  Charge   --- !
     if ( mod( kstep,k_Ecrct ).eq.0 ) then
        call charge
        if ( divE__Fix ) call Ecorrect
     endif
     call TimeMeasure( 2 )

     ! --- [5-8] ppc Ctrl. --- !
     if ( mod( kstep,k_ppcCt ).eq.0 ) call ppcCtrl( 'COM' )
     call TimeMeasure( 8 )
     
     ! --- [5-9] Diagnostics --- !
     if ( mod( kstep,k_prtcl ).eq.0 ) call writePrtcl
     if ( mod( kstep,k_Field ).eq.0 ) call writeField
     if ( mod( kstep,k_vdfsm ).eq.0 ) call VDFsample
     if ( mod( kstep,k_Momnt ).eq.0 ) call MomentAnalysis
     if ( mod( kstep,k_Momnt ).eq.0 ) call Flow_Thermal_Energy
     if ( mod( kstep,k_FlxIp ).eq.0 ) call FluxCurrentAnalysis
     if ( mod( kstep,k_Enrgy ).eq.0 ) call CalcEnergy
     if ( mod( kstep,k_pSort ).eq.0 ) call ptReOrdering( 'jCell' )
     if ( mod( kstep,k_chDDD ).eq.0 ) call DiffuseDecomposition
     call TimeMeasure( 4 )
     ! --- [5-10] Time Out / Save --- !
     if ( mod( kstep,k_tShow ).eq.0 ) call TimeMeasure( -1 )
     if ( mod( kstep,k_tChck ).eq.0 ) call TimeMeasure( -2 )
     if ( mod( kstep,k_SaveP ).eq.0 ) call MakeSavePoint( 'auto' )
     if ( Flag__wallElapsed ) exit
     call TimeMeasure( 1 )
     
  enddo

  ! ------------------------------------- !
  ! --- [6]       Post Process        --- !
  ! ------------------------------------- !
  ! ---   Time Write Out    --- !
  call TimeMeasure( -1 )
  ! ---  Save and Greetings --- !
  if ( .not.( Flag__wallElapsed ) ) kstep = kstep - 1
  call MakeSavePoint( 'finn' )
  call EndingRemarks
  call MPI_FINALIZE ( ierr )

end program myGSchsmMPI
