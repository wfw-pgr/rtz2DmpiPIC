module setupWVMod
  implicit none
  character(3)     :: wavetype
  integer          :: wvMode
  double precision :: wvAmpl
  double precision :: wvAngl
contains

  subroutine LoadWaveConfig
    use constants , only : equDir, myRank
    implicit none
    include 'mpif.h'
    integer        :: ierr
    character(100) :: FileName, buff
    
    if ( myRank.eq.0 ) then
       !    -- [1-1] Read Parameters From 'wave.conf' --  !
       FileName = trim(equDir) // '/' // 'wave.conf'
       open (30, file=trim(FileName), status='old', form='formatted' )
       read (30,'(2(a12,2x),a15)'  ) buff, buff, wavetype
       read (30,'(2(a12,2x),i15)'  ) buff, buff, wvMode
       read (30,'(2(a12,2x),e15.8)') buff, buff, wvAmpl
       read (30,'(2(a12,2x),e15.8)') buff, buff, wvAngl
       close(30)
    endif
    !       -- [1-2] BroadCast parameters             --  !
    call MPI_Bcast( wavetype, 3, MPI_CHARACTER       , 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( wvMode  , 1, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( wvAmpl  , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'(2x,a)') '---------------         [    Wave Config    ]          ---------------'
       write(6,'(4x,a12,1x,a12  )') 'wavetype :: ', wavetype
       write(6,'(4x,a12,1x,i12  )') 'wvMode   :: ', wvMode
       write(6,'(4x,a12,1x,e12.5)') 'wvAmpl   :: ', wvAmpl
       write(6,'(4x,a12,1x,e12.5)') 'wvAngl   :: ', wvAngl
       write(6,*)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
    endif

    return
  end subroutine LoadWaveConfig


  subroutine setupWave
    use constants , only : myRank
    use waveEMOMod, only : setupEMO
    use waveEPWMod, only : setupEPW
    use displayMod, only : displaySectionTitle
    implicit none

    ! --- [1] Load Wave Config  --- !
    if ( myRank.eq.0 ) then
       call displaySectionTitle( 'setupWave'       , '-', 8, 0, 'section'    )
       call displaySectionTitle( 'LOAD WAVE-CONFIG', '-', 4, 4, 'subsection' )
    endif
    call LoadWaveConfig

    ! --- [2] setup Wave Field  --- !
    if ( wavetype.eq.'EMO' ) call setupEMO( wvMode, wvAmpl, wvAngl )
    if ( wavetype.eq.'EPW' ) call setupEPW( wvMode, wvAmpl, wvAngl )

    return
  end subroutine setupWave

  
end module setupWVMod
