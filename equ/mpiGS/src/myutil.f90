module myutil
  implicit none
  integer          :: rtime
  double precision :: ctime, accTime(2,9)
contains

  subroutine MakeJobFile

    use constants, only : job, jobdir, myRank
    implicit none
    integer            :: ierr
    character(100)     :: command

    !  --- [1] Job Directory / MyRank Define ---  !
    jobdir  = 'job/' // trim(job) // '/'

    !  --- [2] Directory Making --- !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*) '  Job Name                            :: ', trim(job)
       write(6,*) '  Data will be saved in the Directory :: ', trim(jobdir)
       write(6,*)
       !     -- [2-0] make job dir           --  !
       command = 'if [ ! -d '//trim(jobdir)// ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// ' ); fi'
       call system( trim(command) )
       !     -- [2-1] Bgrid dual   dir       --  !
       command = 'if [ ! -d '//trim(jobdir)// 'BgD/' // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// 'BgD/' // ' ); fi'
       call system( trim(command) )
       !     -- [2-5] Rgrid single dir       --  !
       command = 'if [ ! -d '//trim(jobdir)// 'RgD/' // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// 'RgD/' // ' ); fi'
       call system( trim(command) )
       !     -- [2-5] PDF dir                --  !
       command = 'if [ ! -d '//trim(jobdir)// 'pdf/' // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// 'pdf/' // ' ); fi'
       call system( trim(command) )
       !     -- [2-6] src copy                --  !
       command = 'if [ ! -d '//trim(jobdir)// 'src/' // ' ]; then'
       command = trim(command) // ' ( cp -r src ' //trim(jobdir)// 'src' // ' ); fi'
       call system( trim(command) )

       !  --- [3] Renew Output Files  ---  !
       open(20,file=trim(jobdir)//'equilibrium.dat',form='formatted',status='replace')
       close(20)
    endif
    
    return
  end subroutine MakeJobFile
  
  
  subroutine TimeMeasure( subject )
    
    implicit none
    include 'mpif.h'
    integer             :: i, k, rtime_, trate, t_max, myRank, ierr
    double precision    :: ctime_, hms(3)
    character(41)       :: fmt
    integer, intent(in) :: subject
    
    if ( subject.eq.0 ) then ! initialize       
       call system_clock( rtime )
       call cpu_time( ctime )
       accTime(:,:) = 0.d0
       
    else if ( subject.eq.-1 ) then
       call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
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
             fmt = '(a,2x,e12.5,3x,2(i3,a),f7.3,a,3x,f7.3,a)'
             hms = sec2hms( accTime(k,1) )
             write(6,fmt) ' Init+Other   :: ', accTime(k,1) &
                  & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,1)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,2) )
             write(6,fmt) ' myBiCGSTAB   :: ', accTime(k,2) &
                  & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,2)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,3) )
             write(6,fmt) ' check/psiD   :: ', accTime(k,3) &
                  & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,3)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,4) )
             write(6,fmt) ' updateJphi   :: ', accTime(k,4) &
                  & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,4)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,5) )
             write(6,fmt) ' WriteField   :: ', accTime(k,5) &
                  & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,5)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,9) )  
             write(6,fmt) ' ---- total ---- ', accTime(k,9) &
                  & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,9)/accTime(k,9)*100.0, ' % '
             write(6,*) ' ----------------------------- '
             write(6,*)

          enddo
       endif

    else if ( ( subject.gt.0 ).and.( subject.le.8 ) ) then
       ! --  CPU time -- !
       ctime_                = ctime
       call cpu_time( ctime )
       accTime(1,subject)    = accTime(1,subject) + ( ctime - ctime_ )
       ! -- Real time -- !
       rtime_                = rtime
       call system_clock( rtime, trate, t_max )
       if ( rtime .lt. rtime_ ) then
          accTime(2,subject) = accTime(2,subject) + ( ( t_max - rtime_ ) + rtime + 1 ) / dble( trate )
       else
          accTime(2,subject) = accTime(2,subject) + ( rtime - rtime_ )                 / dble( trate )
       endif
       
    endif
       
    return
  end subroutine TimeMeasure
  

  function sec2hms( sec )
    implicit none
    double precision, intent(in) :: sec
    real                         :: sec2hms(3)
    sec2hms(1) = real( int( sec                     ) / 3600     )
    sec2hms(2) = real( int( sec - 3600.0*sec2hms(1) ) / 60       )
    sec2hms(3) = real( sec - 3600.0*sec2hms(1) - 60.0*sec2hms(2) )
    return
  end function sec2hms

  
end module myutil
