program EquMain

  ! ------------------------------------------------------ !
  ! --- [1] constants & variables / Routines           --- !
  ! ------------------------------------------------------ !
  use constants , only : myRank, PEtot, coordinate, RgridMode, BgridMode
  use myutil    , only : MakeJobFile, TimeMeasure
  use initialMod, only : InitCond
  use gs_solver , only : SolveGS
  use Field     , only : calc_Bfield_Egrid, calc_Bfield_Bgrid
  use RgridMod  , only : RgridField
  use BgridMod  , only : BgridField
  use diagnosis , only : EquDiagnosis
  use equset    , only : Set2GS_Egrid, Set2GS_Bgrid
  use writingMod, only : WriteField, WriteParam, WriteF_MPI
  use psifunc   , only : pgWrite
  implicit none
  include 'mpif.h'
  integer             :: ierr, PEtoth

  ! ------------------------------------------------------ !
  ! --- [2] Initial Preparation                        --- !
  ! ------------------------------------------------------ !
  !  -- [2-1] MPI settings                             --  !
  call MPI_Init( ierr )
  call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
  call MPI_Comm_Size( MPI_COMM_WORLD, PEtoth, ierr )
  if ( PEtot.ne.PEtoth ) then
     write(6,*) '[ERROR] PEtot is incompatible ', PEtot, &
          &     ' : ', PEtoth, '[ERROR]'
     stop
  endif
  
  !  -- [2-2] Opening Remarks                          --  !
  if ( myRank.eq.0 ) then
     write(6,*)
     write(6,*) '  -------------------------------------------  '
     write(6,*) '  --------- Begining of GS-Solver -----------  '
     write(6,*) '  -------------------------------------------  '
     write(6,*)
     write(6,*)
  endif

  ! ------------------------------------------------------ !
  ! --- [3] Initialization                             --- !
  ! ------------------------------------------------------ !
  call TimeMeasure( 0 )
  call MakeJobFile
  call InitCond
  call TimeMeasure( 1 )
  
  ! ------------------------------------------------------ !
  ! --- [4] Grad-Shafranov Eq. Solving                 --- !
  ! ------------------------------------------------------ !
  call SolveGS
  call TimeMeasure( 3 )

  ! ------------------------------------------------------ !
  ! --- [5] Field Calculation                          --- !
  ! ------------------------------------------------------ !
  if ( myRank.eq.0 ) then
     
     ! -- (1) Common Part          -- !
     call calc_Bfield_Egrid
     call EquDiagnosis
     call pgWrite
     
     ! -- (2) R-grid Single torus  -- !
     if ( RgridMode ) then
        call RgridField( coordinate, 'R' )
        call WriteField( 'R' )
        ! call WriteF_MPI( 'R' )
        call WriteParam( 'R' )
     endif
     ! -- (3) B-grid Single torus  -- !
     if ( BgridMode ) then
        call BgridField( coordinate, 'B' )
        call WriteField( 'B' )
        call WriteF_MPI( 'B' )
        call WriteParam( 'B' )
     endif
     
  endif
  call TimeMeasure( 5 )

  ! ------------------------------------------------------ !
  ! --- [6] Closing Program                            --- !
  ! ------------------------------------------------------ !
  call TimeMeasure( -1 )
  call MPI_Finalize( ierr )

end program EquMain
