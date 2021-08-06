module displayMod
contains
  
  subroutine OpeningRemarks
    use constants, only : myRank, job, Flag__LoadSave
    implicit none
    character(10)      :: msg_LoadGame, msg__NewGame
    ! ---  Massage  --- !
    if ( Flag__LoadSave ) then
       msg__NewGame = '  NEW GAME'
       msg_LoadGame = '* CONTINUE'
    else
       msg__NewGame = '* NEW GAME'
       msg_LoadGame = '  CONTINUE'
    endif
    ! --- [Comment] --- !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'( 2x,a)') '======================================================================'
       write(6,'( 2x,a)') '----------------------------------------------------------------------'
       write(6,'(20x,a)') '                 ____ ___   ____ '
       write(6,'(20x,a)') ' _ __ ___  _   _|  _ \_ _| / ___|'
       write(6,'(20x,a)') "| '_ ` _ \| | | | |_) | | | |    "
       write(6,'(20x,a)') '| | | | | | |_| |  __/| | | |___ '
       write(6,'(20x,a)') '|_| |_| |_|\__, |_|  |___| \____|'
       write(6,'(20x,a)') '           |___/                 '
       write(6,'( 2x,a)') '----------------------------------------------------------------------'
       write(6,'( 2x,a)') '======================================================================'
       write(6,*)
       write(6,'(32x,a)') '<  START  >'
       write(6,*)
       write(6,'(32x,a)') msg__NewGame
       write(6,'(32x,a)') msg_LoadGame
       write(6,*)
       write(6,'(32x,a)') '<   JOB   >'
       write(6,'(35x,a)') trim(job)
       write(6,*)
       write(6,*)
       write(6,*)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,*)
       write(6,*)
    endif
    return
  end subroutine OpeningRemarks

  subroutine MainStartRemarks
    use constants, only : myRank, Flag__LoadSave
    implicit none
    ! --- [Comment] --- !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*)
       write(6,*)
       write(6,'(2x,a)') '======================================================================'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(4x,a)') " ____ ___ ____   _                            _             _   "
       write(6,'(4x,a)') "|  _ \_ _/ ___| | |    ___   ___  _ __    ___| |_ __ _ _ __| |_ "
       write(6,'(4x,a)') "| |_) | | |     | |   / _ \ / _ \| '_ \  / __| __/ _` | '__| __|"
       write(6,'(4x,a)') "|  __/| | |___  | |__| (_) | (_) | |_) | \__ \ || (_| | |  | |_ "
       write(6,'(4x,a)') "|_|  |___\____| |_____\___/ \___/| .__/  |___/\__\__,_|_|   \__|"
       write(6,'(4x,a)') "                                 |_|                            "
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(2x,a)') '======================================================================'
       write(6,*)
       write(6,*)
       write(6,*)
    endif
    return
  end subroutine MainStartRemarks

  
  subroutine EndingRemarks
    use constants, only : myRank, job
    implicit none
    ! --- [Comment] --- !
    if ( myRank.eq.0 ) then
       write(6,*) ' ---- END OF MAIN PROGRAM ---- '
       write(6,*) '    --- JOB NAME ==>> ', trim(job), ' --- '
       write(6,*)
       write(6,*)
       write(6,'( 2x,a)') '======================================================================'
       write(6,'( 2x,a)') ' -------------------------------------------------------------------- '
       write(6,'(25x,a)') ' _____ _   _ ____  '
       write(6,'(25x,a)') '| ____| \ | |  _ \ '
       write(6,'(25x,a)') '|  _| |  \| | | | |'
       write(6,'(25x,a)') '| |___| |\  | |_| |'
       write(6,'(25x,a)') '|_____|_| \_|____/ '
       write(6,*)
       write(6,'( 2x,a)') ' -------------------------------------------------------------------- '
       write(6,'( 2x,a)') '======================================================================'
       write(6,*)
       write(6,*)
    endif
    return
  end subroutine EndingRemarks
  

  subroutine PE_Attendance
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer                               :: ompRank, NThread
    integer                               :: myRank , PEtot
    integer                               :: hostname_len, ierr
    character(len=MPI_MAX_PROCESSOR_NAME) :: hostname

    ! -------------------------------------------------- !
    ! --- [1] Take Attendance of PE & OpenMP threads --- !
    ! -------------------------------------------------- !
    call MPI_Barrier  ( MPI_COMM_WORLD,         ierr )
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )
    call mpi_get_processor_name(hostname,hostname_len,ierr)
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(2x,a)') ' [PE_Attendance] Display < MPI Node >  < HOST NAME >  < THREAD #. >'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
    endif
    call MPI_Barrier  ( MPI_COMM_WORLD, ierr )
    !$omp parallel private(ompRank,NThread)
    !$ NThread = omp_get_num_threads()
    !$ ompRank = omp_get_thread_num() + 1
    !$ write(6,'(3x,a,i4,a,i4,3a,i4,a,i4)') " I am Process-", myRank, "/", PEtot, &
    !$     &   " running on ",trim(hostname), " --- ", ompRank, "-th thread of ", NThread
    !$omp end parallel
    call MPI_Barrier  ( MPI_COMM_WORLD, ierr )
    if ( myRank.eq.0 ) write(6,*)
    return
  end subroutine PE_Attendance

  
  subroutine DisplayParameter
    use constants, only          : jobDir, myRank
    use constants, only          : q, mr, wpewce, vthcv, nuwce, TiTe, pi
    use variables, only          : x1Leng
    implicit none
    character(100)              :: FileName
    double precision, parameter :: B0MKS = 0.01d0   ! in Tesla   ( MKS unit )
    double precision, parameter :: cvMKS = 3.00d8   ! in m/s     ( MKS unit )
    double precision, parameter :: qeMKS = 1.60d-19 ! in coulomb ( MKS unit )
    double precision, parameter :: meMKS = 9.11d-31 ! in kg      ( MKS unit )
    double precision, parameter :: epsMKS= 8.85d-12
    double precision            :: wceMKS, wpeMKS, wciMKS, wpiMKS, tceMKS, tpeMKS, tciMKS, tpiMKS
    double precision            :: vteMKS, vtiMKS, rceMKS, rciMKS
    double precision            :: lDeMKS, lDiMKS, vAeMKS, U0MKS , P0MKS, UkMKS
    double precision            :: x0MKS , x0iMKS, LxMKS , n0MKS , diMKS , deMKS, VAMKS, tAMKS
    double precision            :: E0MKS , J0MKS , I0MKS , IbMKS
    double precision            :: nu_ei , tauei

    ! --------------------------------- !
    ! --- [1] Normalized Constants  --- !
    ! --------------------------------- !
    wceMKS = qeMKS  / meMKS  * B0MKS
    wpeMKS = wpewce * wceMKS
    wciMKS = abs( q(2) / q(1) ) / mr * wceMKS
    wpiMKS = wpeMKS / sqrt( mr )
    tceMKS = ( 2.d0 * pi ) / wceMKS
    tciMKS = ( 2.d0 * pi ) / wciMKS
    tpeMKS = ( 2.d0 * pi ) / wpeMKS
    tpiMKS = ( 2.d0 * pi ) / wpiMKS
    
    vteMKS = vthcv * cvMKS
    vtiMKS = sqrt( TiTe / mr  ) * vteMKS
    rceMKS = vteMKS / wceMKS
    rciMKS = vtiMKS / wciMKS
    lDeMKS = vteMKS / wpeMKS
    lDiMKS = vtiMKS / wpiMKS
    vAeMKS =  cvMKS / wpewce
    
    E0MKS  = cvMKS * B0MKS
    x0MKS  = cvMKS / wceMKS
    LxMKS  = x0MKS * x1Leng
    x0iMKS = cvMKS / wciMKS
    diMKS  = cvMKS / wpiMKS
    deMKS  = cvMKS / wpeMKS
    VAMKS  = cvMKS /(wpewce * sqrt( mr ))
    tAMKS  = LxMKS / VAMKS

    n0MKS  = meMKS * ( wpeMKS / qeMKS )**2 * epsMKS
    J0MKS  = qeMKS * n0MKS * cvMKS
    I0MKS  = J0MKS * x0MKS * x0MKS
    UkMKS  = 0.5d0 * meMKS * cvMKS**2
    U0MKS  = UkMKS * n0MKS * ( x0MKS )**3
    P0MKS  = U0MKS / tceMKS
    IbMKS  = 0.d0
    
    nu_ei  = nuwce * wceMKS
    tauei  = 1.d0  / nu_ei

    if ( myRank.eq.0 ) then
       ! --------------------------------- !
       ! --- [2]  Display  Parameters  --- !
       ! --------------------------------- !
       write(6,*)
       call displaySectionTitle( 'DisplayParameter', '=', 8, 0, 'section' )
       write(6,*)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(2x,a)') ' [DisplayParameter] Display  <<Real>>  Physical  Value for Ref.'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,*) 
       write(6,*) ' ---- [ Normalization quantities ( in MKS unit ) ] ---- '
       write(6,*) 
       write(6,'(2(4x,a,e12.5))') ' wce ( MKS  ) :: ', wceMKS       , ' wpe ( MKS  ) :: ', wpeMKS
       write(6,'(2(4x,a,e12.5))') ' wci ( MKS  ) :: ', wciMKS       , ' wpi ( MKS  ) :: ', wpiMKS
       write(6,'(2(4x,a,e12.5))') ' tce ( MKS  ) :: ', tceMKS       , ' tpe ( MKS  ) :: ', tpeMKS
       write(6,'(2(4x,a,e12.5))') ' tci ( MKS  ) :: ', tciMKS       , ' tpi ( MKS  ) :: ', tpiMKS
       write(6,'(2(4x,a,e12.5))') ' vte ( MKS  ) :: ', vteMKS       , ' vti ( MKS  ) :: ', vtiMKS
       write(6,'(2(4x,a,e12.5))') ' rce ( MKS  ) :: ', rceMKS       , ' rci ( MKS  ) :: ', rciMKS
       write(6,'(2(4x,a,e12.5))') 'c/wce( MKS  ) :: ', x0MKS        , 'c/wci( MKS  ) :: ', x0iMKS
       write(6,'(2(4x,a,e12.5))') ' lDe ( MKS  ) :: ', lDeMKS       , ' lDi ( MKS  ) :: ', lDiMKS
       write(6,'(2(4x,a,e12.5))') ' vAe ( MKS  ) :: ', vAeMKS       , ' vAi ( MKS  ) :: ',  VAMKS
       write(6,*) 
       write(6,'(2(4x,a,e12.5))') ' B0  ( MKS  ) :: ', B0MKS        , ' E0  ( MKS  ) :: ', E0MKS
       write(6,'(2(4x,a,e12.5))') ' x0  ( MKS  ) :: ', x0MKS        , ' Lx  ( MKS  ) :: ', LxMKS
       write(6,'(1(4x,a,e12.5))') ' n0  ( MKS  ) :: ', n0MKS
       write(6,'(2(4x,a,e12.5))') ' J0  ( MKS  ) :: ', J0MKS        , ' I0  ( MKS  ) :: ', I0MKS
       write(6,'(2(4x,a,e12.5))') ' Uk  ( MKS  ) :: ', UkMKS
       write(6,'(2(4x,a,e12.5))') ' U0  ( MKS  ) :: ', U0MKS        , ' P0  ( MKS  ) :: ', P0MKS
       write(6,'(1(4x,a,e12.5))') ' Ib  ( MKS  ) :: ', IbMKS
       write(6,'(2(4x,a,e12.5))') ' de  ( MKS  ) :: ', deMKS        , ' di  ( MKS  ) :: ', diMKS
       write(6,*) 
       write(6,'(2(4x,a,e12.5))') ' VA  ( MKS  ) :: ', VAMKS        , ' tA  ( MKS  ) :: ', tAMKS
       write(6,'(2(4x,a,e12.5))') ' nu_ei        :: ', nu_ei        , ' tauei        :: ', tauei
       write(6,*) 
       write(6,*) 
       write(6,*) ' ---- [ Normalization quantities (  with  unit ) ] ---- '
       write(6,*) 
       write(6,'(2(4x,a,f12.5))') ' wce ( GHz  ) :: ', wceMKS*1.d-9 , ' wpe ( GHz  ) :: ', wpeMKS*1.d-9
       write(6,'(2(4x,a,f12.5))') ' vte ( km/s ) :: ', vteMKS*1.d-3 , ' rce (  mm  ) :: ', rceMKS*1.d+3
       write(6,'(2(4x,a,f12.5))') ' lDe (  mm  ) :: ', lDeMKS*1.d+3 , ' vAe ( km/s ) :: ', vAeMKS*1.d-3
       write(6,'(2(4x,a,f12.5))') ' wci ( MHz  ) :: ', wciMKS*1.d-6 , ' wpi ( MHz  ) :: ', wpiMKS*1.d-6
       write(6,'(2(4x,a,f12.5))') ' vti ( km/s ) :: ', vtiMKS*1.d-3 , ' rci (  mm  ) :: ', rciMKS*1.d+3
       write(6,'(2(4x,a,f12.5))') ' lDi (  mm  ) :: ', lDiMKS*1.d+3 , ' vAi ( km/s ) :: ',  VAMKS*1.d-3
       write(6,*)
       write(6,'(2(4x,a,f12.5))') ' B0  (  mT  ) :: ',  B0MKS*1.d+3 , ' E0  ( MV/m ) :: ',  E0MKS*1.d-6
       write(6,'(4x,a,f12.5,a)' ) ' n0  ( m^-3 ) :: ',  n0MKS*1.d-18, ' x 10^18 '
       write(6,'(2(4x,a,f12.5))') ' J0  (MA/m^2) :: ',  J0MKS*1.d-6 , ' I0  (  kA  ) :: ',  I0MKS*1.d-3
       write(6,'(2(4x,a,e12.5))') ' Uk  (   J  ) :: ',  UkMKS
       write(6,'(2(4x,a,e12.5))') ' U0  (   J  ) :: ',  U0MKS       , ' P0  (  kW  ) :: ',  P0MKS*1.d-3
       write(6,'(2(4x,a,f12.5))') ' Ib  (   A  ) :: ',  IbMKS
       write(6,'(2(4x,a,f12.5))') ' de  (  mm  ) :: ',  deMKS*1.d+3 , ' di  (  mm  ) :: ',  diMKS*1.d+3
       write(6,'(1(4x,a,f12.5))') ' x0  (  mm  ) :: ',  x0MKS*1.d+3 , ' Lx  (  mm  ) :: ',  LxMKS*1.d+3
       write(6,*)
       write(6,'(2(4x,a,f12.5))') ' VA  ( km/s ) :: ',  VAMKS*1.d-3 , ' tA  (  us  ) :: ',  tAMKS*1.d+6
       write(6,'(2(4x,a,f12.5))') 'nu_ei( kHz  ) :: ',  nu_ei*1.d-3 , 'tauei(  us  ) :: ',  tauei*1.d+6
       write(6,*)
       write(6,'(2x,a)'        ) '----------------------------------------------------------------------'
       write(6,*)
       
       ! --------------------------------- !
       ! --- [3]  ConvertedParams.dat  --- !
       ! --------------------------------- !
       FileName   = trim(jobDir)//'dat/'//'ConvertedParams.dat'
       open( 50, file=trim(FileName), form='formatted' )
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'wce', 'GHz'   , wceMKS*1.d-9, 'wpe', 'GHz' , wpeMKS*1.d-9
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'vte', 'km/s'  , vteMKS*1.d-3, 'rce', 'mm'  , rceMKS*1.d+3
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'lDe', 'mm'    , lDeMKS*1.d+3, 'vAe', 'km/s', vAeMKS*1.d-3
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'wci', 'MHz'   , wciMKS*1.d-6, 'wpi', 'MHz' , wpiMKS*1.d-6
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'vti', 'km/s'  , vtiMKS*1.d-3, 'rci', 'mm'  , rciMKS*1.d+3
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'lDi', 'mm'    , lDiMKS*1.d+3, 'vAi', 'km/s',  VAMKS*1.d-3
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'B0' , 'mT'    ,  B0MKS*1.d+3, 'E0' , 'MV/m',  E0MKS*1.d-6
       write(50,'(1(a15,4x,a15,4x,e12.5))') 'n0' , '10^{18}m^{-3}'  ,  n0MKS*1.d-18
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'J0' , 'MA/m^2',  J0MKS*1.d-6, 'I0' , 'kA'  ,  I0MKS*1.d-3
       write(50,'(1(a15,4x,a15,4x,e12.5))') 'Uk' , 'J'     ,  UkMKS
       write(50,'(1(a15,4x,a15,4x,e12.5))') 'U0' , 'J'     ,  U0MKS,       'P0' , 'kW'  ,  P0MKS*1.d-3
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'Ib' , 'A'     ,  IbMKS
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'de' , 'mm'    ,  deMKS*1.d+3, 'di' , 'mm'  ,  diMKS*1.d+3
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'x0' , 'mm'    ,  x0MKS*1.d+3, 'Lx' , 'mm'  ,  LxMKS*1.d+3
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'VA' , 'km/s'  ,  VAMKS*1.d-3, 'tA' , 'us'  ,  tAMKS*1.d+6
       write(50,'(1(a15,4x,a15,4x,f12.5))') 'nu_ei', 'kHz' ,  nu_ei*1.d-3, 'tauei','us' ,  tauei*1.d+6
       close(50)
       write(6,'(2x,a10,a45,2x,a6)') "* SAVE :: ", trim(FileName), '[ OK ]'
       write(6,*)
       call displaySectionTitle( 'DisplayParameter -- END --', '=', 8, 0, 'section' )
       write(6,*)
       write(6,*)
    endif
    return
  end subroutine DisplayParameter


  subroutine displaySectionTitle( title, mark, nMark , nSpace, sectionType )
    implicit none
    character(*) , intent(in) :: title
    character    , intent(in) :: mark
    integer      , intent(in) :: nMark, nSpace
    character(*) , intent(in) :: sectionType
    integer      , parameter  :: windowLength = 70
    integer      , parameter  :: commonIndent =  2
    integer                   :: nchar, boxLength, nSpaceF, nSpaceB
    integer                   :: nFrontSpace
    character(2)              :: cFrontSpace
    character(7)              :: fmt
    character(100)            :: Frame, BoxedTitle
    
    nchar          = len( trim(title) )
    if ( trim(sectionType).eq.'subsection' ) then
       boxLength   = nMark*2 + nSpace*2 + nchar
       nFrontSpace = ( windowLength - boxLength ) / 2 + commonIndent
       nSpaceF     = nSpace
       nSpaceB     = nSpace
    endif
    if ( trim(sectionType).eq.'section'    ) then
       boxLength   = windowLength
       nFrontSpace = commonIndent
       nSpaceF     = ( windowLength - nchar - 2*nMark ) / 2
       nSpaceB     =   windowLength - nchar - 2*nMark - nSpaceF
    endif
    Frame          = repeat( mark,boxLength )
    BoxedTitle     = repeat( mark,nMark     ) // repeat( ' ' , nSpaceF ) // title // &
         &            repeat( ' ' ,nSpaceB   ) // repeat( mark, nMark   )
    write(cFrontSpace,'(i2)') nFrontSpace
    fmt            = '(' // cFrontSpace  // 'x,a)'
    write(6,*)
    write(6,fmt) trim(Frame)
    write(6,fmt) trim(BoxedTitle)
    write(6,fmt) trim(Frame)
    write(6,*)
    
    return
  end subroutine displaySectionTitle
  
  
end module displayMod
