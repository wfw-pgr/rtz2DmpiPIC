module utilityMod
  use constants, only : wallLimit
  implicit none
  integer            :: rtime
  double precision   :: ctime, accTime(2,9)
  double precision   :: wallLimit_sec
  character(100)     :: QuitOrderFile
contains

  ! =================================================================== !
  ! ===  InitiateMPI :: Check  ( PEtot & OMPNumThreads )            === !
  ! =================================================================== !
  subroutine InitiateMPI
    use constants, only : myRank, cRank
    use constants, only : PEtot , OMPNumThreads
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer            :: PE, NThread, ierr
    
    ! -------------------------------- !
    ! --- [1] Initialize MPI       --- !
    ! -------------------------------- !
    call MPI_INIT     ( ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PE    , ierr )
    write(cRank,'(i6.6)') myRank
    ! -------------------------------- !
    ! --- [2] Check PEtot          --- !
    ! -------------------------------- !
    if ( ( PE.ne.PEtot ).and.( myRank.eq.0 ) ) then
       write(6,*) ' [ERROR] Incompatible PEtot [ERROR] '
       write(6,*) '           PE == ', PE, ' <--->   PEtot == ', PEtot
       stop
    endif
    ! -------------------------------- !
    ! --- [3] Check OMPNumThreads  --- !
    ! -------------------------------- !
    !$omp parallel shared(NThread)
    !$omp single
    !$ NThread = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    if ( OMPNumThreads.ne.NThread ) then
       if ( myRank.eq.0 ) then
          write(6,*) ' [ERROR] Incompatible OMP_NUM_THREADS [ERROR] '
          write(6,'((a,i4))') '    OMPNumThreads       == ', OMPNumThreads
          write(6,'((a,i4))') '    omp_get_num_threads == ', NThread
          stop
       endif
    endif
    return
  end subroutine InitiateMPI


  ! =================================================================== !
  ! ===  MakeJobFile :: Prepare job Directories                     === !
  ! =================================================================== !
  subroutine MakeJobFile
    use constants, only : job, jobDir, equDir, Flag__LoadSave
    use constants, only : binDir, datDir, savDir
    use constants, only : savedjob, myRank
    use constants, only : equType, equJob, equGrid
    implicit none
    character(150)     :: filename
    character(100)     :: save_datDir

    ! ------------------------------- !
    ! --- [1] Job Main Directory  --- !
    ! ------------------------------- !
    jobDir        = 'job/' // trim(job)      // '/'
    save_datDir   = 'job/' // trim(savedjob) // '/'//'dat/'
    datDir        = trim(jobDir) // 'dat/'
    binDir        = trim(jobDir) // 'bin/'
    savDir        = trim(jobDir) // 'sav/'
    QuitOrderFile = trim(jobDir) // 'dat/'   // 'QuitOrder.dat'
    equDir        = 'equ/' // trim(equType) // '/job/' // trim(equJob) // '/' // trim(equGrid)

    if ( myRank.eq.0 ) then
       ! ---------------------------- !
       ! --- [2] Make Directory   --- !
       ! ---------------------------- !
       !  -- [2-1] Make subDir -- !
       call system( "mkdir -p " // trim(jobDir)                )
       call system( "mkdir -p " // trim(binDir)                )
       call system( "mkdir -p " // trim(datDir)                )
       call system( "mkdir -p " // trim(savDir)                )
       call system( "mkdir -p " // trim(jobDir) // 'pdf/'      )
       call system( "mkdir -p " // trim(jobDir) // 'png/'      )
       call system( "mkdir -p " // trim(jobDir) // 'vtk/'      )
       call system( "mkdir -p " // trim(datDir) // 'ijDomain/' )
       !  -- [2-2] Copy equDir -- !
       call system( 'if [ -d '  // trim(jobDir) // 'equ' // ' ]; then'   // &
            &       ' ( rm -r ' // trim(jobDir) // 'equ' // ' ); fi'              )
       call system( 'cp -r '    // trim(equDir) // ' '   // trim(jobDir) // 'equ' )
       !  -- [2-3] Copy src    -- !
       call system( 'cp -r '    // 'src '                // trim(jobDir) // 'src' )
       ! -------------------------------------- !
       ! --- [3] Copy From Saved Directory  --- !
       ! -------------------------------------- !
       if ( ( Flag__LoadSave ).and.( savedjob.ne.job ) ) then
          call system('cp '//trim(save_datDir)//'energy.dat '           //trim(datDir)//'energy.dat'           )
          call system('cp '//trim(save_datDir)//'energy_kinetic.dat '   //trim(datDir)//'energy_kinetic.dat'   )
          call system('cp '//trim(save_datDir)//'energy_field.dat '     //trim(datDir)//'energy_field.dat'     )
          call system('cp '//trim(save_datDir)//'charge_residual.dat '  //trim(datDir)//'charge_residual.dat'  )
          call system('cp '//trim(save_datDir)//'constants.dat '        //trim(datDir)//'constants.dat'        )
          call system('cp '//trim(save_datDir)//'Current.dat '          //trim(datDir)//'Current.dat'          )
          call system('cp '//trim(save_datDir)//'MagneticFlux.dat '     //trim(datDir)//'MagneticFlux.dat'     )
          call system('cp '//trim(save_datDir)//'XOpoint.dat '          //trim(datDir)//'XOpoint.dat'          )
          call system('cp '//trim(save_datDir)//'ConstantsOfMotion.dat '//trim(datDir)//'ConstantsOfMotion.dat')
       endif
       ! -------------------------------------- !
       ! --- [4] BackUp / Renew .dat Files  --- !
       ! -------------------------------------- !
       if ( Flag__LoadSave ) then
          ! -- BackUp -- !
          call system('cp '//trim(datDir)//'energy.dat '           //trim(datDir)//'energy.dat.org'            )
          call system('cp '//trim(datDir)//'energy_kinetic.dat '   //trim(datDir)//'energy_kinetic.dat.org'    )
          call system('cp '//trim(datDir)//'energy_field.dat '     //trim(datDir)//'energy_field.dat.org'      )
          call system('cp '//trim(datDir)//'charge_residual.dat '  //trim(datDir)//'charge_residual.dat.org'   )
          call system('cp '//trim(datDir)//'Current.dat '          //trim(datDir)//'Current.dat.org'           )
          call system('cp '//trim(datDir)//'XOpoint.dat '          //trim(datDir)//'XOpoint.dat.org'           )
          call system('cp '//trim(datDir)//'MagneticFlux.dat '     //trim(datDir)//'MagneticFlux.dat.org'      )
          call system('cp '//trim(datDir)//'ConstantsOfMotion.dat '//trim(datDir)//'ConstantsOfMotion.dat.org' )
       else
          ! -- Renew  -- !
          filename = trim(datDir) // 'energy.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'energy_kinetic.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'energy_field.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'charge_residual.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'Current.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'MagneticFlux.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'XOpoint.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'ConstantsOfMotion.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(datDir) // 'Nparticles.dat'
          open( 50, file=trim(filename), form='formatted', status='replace' )
          close(50)
          filename = trim(QuitOrderFile)
          open (50, file=trim(filename), form='formatted', status='replace' )
          write(50,*) 0
          close(50)
       end if
    endif
    return
  end subroutine MakeJobFile


  ! =================================================================== !
  ! ===  TimeMeasure ::  Measure Time ( cpu & real )                === !
  ! =================================================================== !
  subroutine TimeMeasure( subject )
    use constants, only  : wallLimit, Flag__wallElapsed
    implicit none
    include 'mpif.h'
    integer             :: i, k, rtime_, trate, t_max, myRank, ierr, ihms(3), quitorder
    double precision    :: ctime_, hms(3)
    character(41)       :: fmt
    integer, intent(in) :: subject
    
    ! ----------------------------------------- !
    ! --- [0] Initialization :: subject=  0 --- !
    ! ----------------------------------------- !
    if ( subject.eq.0 ) then ! initialize
       call system_clock( rtime )
       call cpu_time    ( ctime )
       accTime(:,:)      = 0.d0
       read(wallLimit(1:3),*) ihms(1)
       read(wallLimit(5:6),*) ihms(2)
       read(wallLimit(8:9),*) ihms(3)
       wallLimit_sec     = dble(ihms(1))*3600.d0 + dble(ihms(2))*60.d0 + dble(ihms(3))
       Flag__wallElapsed = .false.
    endif
    ! ----------------------------------------- !
    ! --- [-1] Display Time  :: subject= -1 --- !
    ! ----------------------------------------- !
    if ( subject.eq.-1 ) then
       call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
       if ( myRank.eq.0 ) then
          write(6,*)
          write(6,*) ' ---------------   Computation Time   --------------- '
          write(6,*)
          accTime(:,9) = 0.d0
          do k=1, 2
             do i=1, 8
                accTime(k,9) = accTime(k,9) + accTime(k,i)
             enddo
             
             if ( k.eq.1 ) write(6,*) ' ----       CPU    Time      ---- '
             if ( k.eq.2 ) write(6,*) ' ----       Real   Time      ---- '
             hms = sec2hms( accTime(k,1) )
             fmt = '(a,2x,e12.5,3x,2(i3,a),f7.3,a,3x,f7.3,a)'
             write(6,fmt) ' init+other   :: ', accTime(k,1), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,1)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,2) )
             write(6,fmt) ' charge       :: ', accTime(k,2), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,2)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,3) )
             write(6,fmt) ' current      :: ', accTime(k,3), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,3)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,4) )
             write(6,fmt) ' write&energy :: ', accTime(k,4), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,4)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,5) )
             write(6,fmt) ' particle     :: ', accTime(k,5), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,5)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,6) )
             write(6,fmt) ' field        :: ', accTime(k,6), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,6)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,7) )
             write(6,fmt) ' T/A coll     :: ', accTime(k,7), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,7)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,8) )
             write(6,fmt) ' ppc control  :: ', accTime(k,8), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,8)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,9) )
             write(6,fmt) ' ---- total ---- ', accTime(k,9), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,9)/accTime(k,9)*100.0, ' % '
             write(6,*) ' ----------------------------- '
             write(6,*)

          enddo
       endif
    endif
    ! ----------------------------------------- !
    ! --- [-2]  Check  Time  :: subject= -2 --- !
    ! ----------------------------------------- !
    if ( subject.eq.-2 ) then
       call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
       if ( myRank.eq.0 ) then
          ! -------------------------- !
          ! --  Time Elapsed Check  -- !
          ! -------------------------- !
          accTime(:,9) = 0.d0
          do i=1, 8
             accTime(2,9) = accTime(2,9) + accTime(2,i)
          enddo
          if ( accTime(2,9).gt.wallLimit_sec ) then
             write(6,*) accTime(2,9), wallLimit_sec
             write(6,*) " [END]  time elapsed wallLimit_sec...  [END] "
             Flag__wallElapsed = .true.
          endif
          ! -------------------------- !
          ! -- Quit Order From File -- !
          ! -------------------------- !
          open (50,file=trim(QuitOrderFile),form='formatted',status='old',err=100)
          read (50,*) quitorder
          close(50)
          if ( quitorder.eq.1 ) then
             write(6,*) " [END]  QuitOrder was detected...  Save and EXIT.  [END] "
             Flag__wallElapsed = .true.
          endif
100       continue
       endif
       call MPI_BCAST(Flag__wallElapsed,1,MPI_Logical,0,MPI_COMM_WORLD,ierr)
    endif
    ! ----------------------------------------- !
    ! --- [1~8]  Count Up   :: subject= 1~8 --- !
    ! ----------------------------------------- !
    if ( ( subject.gt.0 ).and.( subject.le.8 ) ) then
       ! ----------------- !
       ! ---  CPU time --- !
       ! ----------------- !
       ctime_                = ctime
       call cpu_time( ctime )
       accTime(1,subject)    = accTime(1,subject) + ( ctime - ctime_ )
       ! ----------------- !
       ! --- Real time --- !
       ! ----------------- !
       rtime_                = rtime
       call system_clock( rtime, trate, t_max )
       if ( rtime.lt.rtime_ ) then
          accTime(2,subject) = accTime(2,subject) + ( ( t_max - rtime_ ) + rtime + 1 ) / dble( trate )
       else
          accTime(2,subject) = accTime(2,subject) + ( rtime - rtime_ ) / dble( trate )
       endif
    endif
    
    return
  end subroutine TimeMeasure
  

  ! =================================================================== !
  ! ===  sec2hms :: 17 years old -> real age                        === !
  ! =================================================================== !
  function sec2hms( sec )
    implicit none
    double precision, intent(in) :: sec
    real                         :: sec2hms(3)
    sec2hms(1) = real( int( sec ) / 3600 )
    sec2hms(2) = real( int( sec - 3600.0*sec2hms(1) ) / 60 )
    sec2hms(3) = real( sec - 3600.0*sec2hms(1) - 60.0*sec2hms(2) )
    return
  end function sec2hms

  
end module utilityMod
