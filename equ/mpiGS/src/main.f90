! ====================================================== !
! === mpiGS.f90                                      === !
! ====================================================== !
program main

  ! ------------------------------------------------------ !
  ! --- [1] constants & variables / Routines           --- !
  ! ------------------------------------------------------ !
  use variablesMod
  use utilitiesMod, only : TimeMeasure_MPI, show__programLogo_MPI, show__endLogo_MPI
  use initiatorMod, only : initialize__mpiGS
  ! use initialMod  , only : InitCond
  ! use gs_solver   , only : SolveGS
  ! use Field       , only : calc_Bfield_Egrid, calc_Bfield_Bgrid
  ! use RgridMod    , only : RgridField
  ! use BgridMod    , only : BgridField
  ! use diagnosis   , only : EquDiagnosis
  ! use equset      , only : Set2GS_Egrid, Set2GS_Bgrid
  ! use writingMod  , only : WriteField, WriteParam, WriteF_MPI
  ! use psifunc     , only : pgWrite
  implicit none
  include 'mpif.h'
  integer             :: ierr, PEtoth

  ! ------------------------------------------------------ !
  ! --- [2] Initial Preparation                        --- !
  ! ------------------------------------------------------ !
  !  -- [2-1] MPI settings                             --  !
  call MPI_Init     ( ierr )
  call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
  ! call MPI_Comm_Size( MPI_COMM_WORLD, PEtot , ierr )
  call MPI_Comm_Size( MPI_COMM_WORLD, PEtoth, ierr )
  if ( PEtot.ne.PEtoth ) then
     write(6,*) '[ERROR] PEtot is incompatible ', PEtot, &
          &     ' : ', PEtoth, '[ERROR]'
     stop
  endif
  !  -- [2-2] Opening Remarks                          --  !
  call show__programLogo_MPI
  
  ! ------------------------------------------------------ !
  ! --- [3] Initialization                             --- !
  ! ------------------------------------------------------ !
  call TimeMeasure_MPI( 0 )
  call initialize__mpiGS
  call TimeMeasure_MPI( 1 )
  
  ! ------------------------------------------------------ !
  ! --- [4] Grad-Shafranov Eq. Solving                 --- !
  ! ------------------------------------------------------ !
  ! call SolveGS
  call TimeMeasure_MPI( 3 )

  ! ------------------------------------------------------ !
  ! --- [5] Field Calculation                          --- !
  ! ------------------------------------------------------ !
  ! if ( myRank.eq.0 ) then
     
  !    ! -- (1) Common Part          -- !
  !    call calc_Bfield_Egrid
  !    call EquDiagnosis
  !    call pgWrite
     
  !    ! -- (2) R-grid Single torus  -- !
  !    if ( RgridMode ) then
  !       call RgridField( coordinate, 'R' )
  !       call WriteField( 'R' )
  !       ! call WriteF_MPI( 'R' )
  !       call WriteParam( 'R' )
  !    endif
  !    ! -- (3) B-grid Single torus  -- !
  !    if ( BgridMode ) then
  !       call BgridField( coordinate, 'B' )
  !       call WriteField( 'B' )
  !       call WriteF_MPI( 'B' )
  !       call WriteParam( 'B' )
  !    endif
     
  ! endif
  call TimeMeasure_MPI( 5 )

  ! ------------------------------------------------------ !
  ! --- [6] Closing Program                            --- !
  ! ------------------------------------------------------ !
  call TimeMeasure_MPI( -1 )
  call show__endLogo_MPI
  call MPI_Finalize( ierr )

end program main
