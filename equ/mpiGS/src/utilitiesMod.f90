module utilitiesMod
  implicit none
  integer            :: rtime
  double precision   :: ctime, accTime(2,9)
contains

  ! =================================================================== !
  ! ===  TimeMeasure_MPI ::  Measure Time ( cpu & real )            === !
  ! =================================================================== !
  subroutine TimeMeasure_MPI( subject )
    implicit none
    include 'mpif.h'
    integer             :: i, k, rtime_, trate, t_max, myRank, ierr
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
    endif
    ! ----------------------------------------- !
    ! --- [-1] Display Time  :: subject= -1 --- !
    ! ----------------------------------------- !
    if ( subject.eq.-1 ) then
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
             hms = sec2hms( accTime(k,1) )
             fmt = '(a,2x,e12.5,3x,2(i3,a),f7.3,a,3x,f7.3,a)'
             write(6,fmt) ' Initialize   :: ', accTime(k,1), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,1)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,2) )
             write(6,fmt) ' MainSolver   :: ', accTime(k,2), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,2)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,3) )
             write(6,fmt) ' update Vars  :: ', accTime(k,3), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,3)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,4) )
             write(6,fmt) ' write File   :: ', accTime(k,4), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,4)/accTime(k,9)*100.0, ' % '
             hms = sec2hms( accTime(k,5) )
             write(6,fmt) ' postProcess  :: ', accTime(k,5), &
                  &         int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', &
                  &         accTime(k,5)/accTime(k,9)*100.0, ' % '
             write(6,*) ' ----------------------------- '
             write(6,*)

          enddo
       endif
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
  end subroutine TimeMeasure_MPI


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


  ! ===================================================== !
  ! ===  DisplayClock                                 === !
  ! ===================================================== !
  subroutine DisplayClock
    implicit none
    integer           :: date_time(8)
    character(len=10) :: sys_time(3)
    character(len=34) :: fmt
    call date_and_time(sys_time(1), sys_time(2), sys_time(3), date_time)
    fmt = '(a,1x,i4.4,2(a,i2.2),3x,3(i2.2,a))'
    write(6,fmt) '[DisplayClocK] ', date_time(1), '-', date_time(2), '-', date_time(3), &
         &                          date_time(5), ':', date_time(6), ':', date_time(7), ' --'
    return
  end subroutine DisplayClock


  ! ====================================================== !
  ! === print section                                  === !
  ! ====================================================== !
  subroutine print__section( name, frame, nameLen, sideLen, frameLen )
    implicit none
    integer           , intent(in) :: nameLen, frameLen, sideLen
    character(nameLen), intent(in) :: name
    character(1)      , intent(in) :: frame
    character(2)                   :: num1, num2, num3, num4
    character(100)                 :: fmt1, fmt2
    integer                        :: spaceLen, remain

    ! ------------------------------------------------------ !
    ! --- [1] prepare                                    --- !
    ! ------------------------------------------------------ !
    spaceLen = ( frameLen - nameLen - 2*sideLen ) / 2
    remain   = mod( ( frameLen - nameLen - 2*sideLen ),  2 )
    write(num1,"(i2)")  nameLen + remain
    write(num2,"(i2)")  sideLen
    write(num3,"(i2)") frameLen
    write(num4,"(i2)") spaceLen
    num2     = adjustl(num2)
    fmt1     = "(a" // trim(num3) // ")"
    fmt2     = "(a" // trim(num2) // "," // trim(num4) // "x,a" // trim(num1) // ","
    fmt2     = trim(fmt2) // trim(num4) // "x,a" // trim(num2) // ")"

    ! ------------------------------------------------------ !
    ! --- [2] write section                              --- !
    ! ------------------------------------------------------ !
    write(6,trim(fmt1)) repeat( frame, frameLen )
    write(6,trim(fmt2)) repeat( frame, sideLen  ), trim(adjustl(name)), repeat( frame, sideLen  )
    write(6,trim(fmt1)) repeat( frame, frameLen )
    
    return
  end subroutine print__section


  ! ====================================================== !
  ! === print section (MPI ver.)                       === !
  ! ====================================================== !
  subroutine print__section_MPI( name, frame, nameLen, sideLen, frameLen )
    implicit none
    include 'mpif.h'
    integer           , intent(in) :: nameLen, frameLen, sideLen
    character(nameLen), intent(in) :: name
    character(1)      , intent(in) :: frame
    character(2)                   :: num1, num2, num3, num4
    character(100)                 :: fmt1, fmt2
    integer                        :: spaceLen, remain, myRank, ierr

    ! ------------------------------------------------------ !
    ! --- [1] prepare                                    --- !
    ! ------------------------------------------------------ !
    spaceLen = ( frameLen - nameLen - 2*sideLen ) / 2
    remain   = mod( ( frameLen - nameLen - 2*sideLen ),  2 )
    write(num1,"(i2)")  nameLen + remain
    write(num2,"(i2)")  sideLen
    write(num3,"(i2)") frameLen
    write(num4,"(i2)") spaceLen
    num2     = adjustl(num2)
    fmt1     = "(a" // trim(num3) // ")"
    fmt2     = "(a" // trim(num2) // "," // trim(num4) // "x,a" // trim(num1) // ","
    fmt2     = trim(fmt2) // trim(num4) // "x,a" // trim(num2) // ")"

    ! ------------------------------------------------------ !
    ! --- [2] write section                              --- !
    ! ------------------------------------------------------ !
    call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
    if ( myRank.eq.0 ) then
       write(6,trim(fmt1)) repeat( frame, frameLen )
       write(6,trim(fmt2)) repeat( frame, sideLen  ), &
            &              trim(adjustl(name)), repeat( frame, sideLen  )
       write(6,trim(fmt1)) repeat( frame, frameLen )
    endif
    return
  end subroutine print__section_MPI

  
  ! ====================================================== !
  ! === display program title logo                     === !
  ! ====================================================== !
  subroutine show__programLogo
    implicit none
    character(1)  :: frame
    character(4)  :: side
    character(56) :: line1, line2, line3, line4, line5, line6

    frame = "="
    side  = "===="
    line1 = " _                  _                _           _      "
    line2 = "| |___ _ __ ___    / \   _ __   __ _| |_   _ ___(_)___  "
    line3 = "| / __| '_ ` _ \  / _ \ | '_ \ / _` | | | | / __| / __| "
    line4 = "| \__ \ | | | | |/ ___ \| | | | (_| | | |_| \__ \ \__ \ "
    line5 = "|_|___/_| |_| |_/_/   \_\_| |_|\__,_|_|\__, |___/_|___/ "
    line6 = "                                       |___/            "

    write(6,*)
    write(6,*)
    write(6,'((a70))') repeat( frame, 70 )
    write(6,'(a4,3x,a56,3x,a4)')  side,  line1,  side
    write(6,'(a4,3x,a56,3x,a4)')  side,  line2,  side
    write(6,'(a4,3x,a56,3x,a4)')  side,  line3,  side
    write(6,'(a4,3x,a56,3x,a4)')  side,  line4,  side
    write(6,'(a4,3x,a56,3x,a4)')  side,  line5,  side
    write(6,'(a4,3x,a56,3x,a4)')  side,  line6,  side
    write(6,'((a70))') repeat( frame, 70 )
    write(6,*)
    write(6,*)

    return
  end subroutine show__programLogo


  ! ====================================================== !
  ! === display program end logo                       === !
  ! ====================================================== !
  subroutine show__endLogo
    implicit none
    character(1)  :: frame
    character(4)  :: side
    character(20) :: line1, line2, line3, line4, line5

    frame = "="
    side  = "===="
    line1 = " _____ _   _ ____   "
    line2 = "| ____| \ | |  _ \  "
    line3 = "|  _| |  \| | | | | "
    line4 = "| |___| |\  | |_| | "
    line5 = "|_____|_| \_|____/  "

    write(6,*)
    write(6,*)
    write(6,"((a70))") repeat( frame, 70 )
    write(6,"(a4,21x,a20,21x,a4)")  side,  line1,  side
    write(6,"(a4,21x,a20,21x,a4)")  side,  line2,  side
    write(6,"(a4,21x,a20,21x,a4)")  side,  line3,  side
    write(6,"(a4,21x,a20,21x,a4)")  side,  line4,  side
    write(6,"(a4,21x,a20,21x,a4)")  side,  line5,  side
    write(6,"((a70))") repeat(  frame, 70 )
    write(6,*)
    write(6,*)

    return
  end subroutine show__endLogo


  ! ====================================================== !
  ! === display program title logo                     === !
  ! ====================================================== !
  subroutine show__programLogo_MPI
    implicit none
    include "mpif.h"
    integer       :: ierr, myRank
    character(1)  :: frame
    character(4)  :: side
    character(56) :: line1, line2, line3, line4, line5, line6

    frame = "="
    side  = "===="

    line1 = "                  _  ____ ____  "
    line2 = "  _ __ ___  _ __ (_)/ ___/ ___| "
    line3 = " | '_ ` _ \| '_ \| | |  _\___ \ "
    line4 = " | | | | | | |_) | | |_| |___) |"
    line5 = " |_| |_| |_| .__/|_|\____|____/ "
    line6 = "           |_|                  "

    call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*)
       write(6,'((a70))') repeat( frame, 70 )
       write(6,'(a4,15x,a32,15x,a4)')  side,  line1,  side
       write(6,'(a4,15x,a32,15x,a4)')  side,  line2,  side
       write(6,'(a4,15x,a32,15x,a4)')  side,  line3,  side
       write(6,'(a4,15x,a32,15x,a4)')  side,  line4,  side
       write(6,'(a4,15x,a32,15x,a4)')  side,  line5,  side
       write(6,'(a4,15x,a32,15x,a4)')  side,  line6,  side
       write(6,'((a70))') repeat( frame, 70 )
       write(6,*)
       write(6,*)
    endif

    return
  end subroutine show__programLogo_MPI


  ! ====================================================== !
  ! === display program end logo                       === !
  ! ====================================================== !
  subroutine show__endLogo_MPI
    implicit none
    include "mpif.h"
    integer       :: ierr, myRank
    character(1)  :: frame
    character(4)  :: side
    character(20) :: line1, line2, line3, line4, line5

    frame = "="
    side  = "===="
    line1 = " _____ _   _ ____   "
    line2 = "| ____| \ | |  _ \  "
    line3 = "|  _| |  \| | | | | "
    line4 = "| |___| |\  | |_| | "
    line5 = "|_____|_| \_|____/  "

    call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*)
       write(6,"((a70))") repeat( frame, 70 )
       write(6,"(a4,21x,a20,21x,a4)")  side,  line1,  side
       write(6,"(a4,21x,a20,21x,a4)")  side,  line2,  side
       write(6,"(a4,21x,a20,21x,a4)")  side,  line3,  side
       write(6,"(a4,21x,a20,21x,a4)")  side,  line4,  side
       write(6,"(a4,21x,a20,21x,a4)")  side,  line5,  side
       write(6,"((a70))") repeat(  frame, 70 )
       write(6,*)
       write(6,*)
    endif

    return
  end subroutine show__endLogo_MPI

  
end module utilitiesMod
