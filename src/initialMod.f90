module InitialMod
  implicit none
  integer :: j0Repres = 0
contains

  ! =================================================================== !
  ! ===  Initiator  :: Initiate PIC simulation                      === !
  ! =================================================================== !
  subroutine Initiator
    use constants , only : LIs, LJs, Flag__LoadSave, myRank
    use constants , only : Boundary1__EB, Boundary1__jc
    use setupGSMod, only : setupGS
    use displayMod, only : displaySectionTitle
    use allocatMod, only : InitVariables
    use randGenMod, only : InitRandomSeed
    use saveRunMod, only : ReadSavePoint
    use fMPIComMod, only : allocExchangeTable_Field    , setupExchangeTable_Field
    use jMPIComMod, only : allocExchangeTable_Current  , setupExchangeTable_Current
    use rMPIComMod, only : allocExchangeTable_Relocated, setupExchangeTable_Relocated
    use pMPIComMod, only : setupExchangeTable_Particle
    use ptcBdrcMod, only : particleBoundary
    implicit none
    
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*)
       call displaySectionTitle( 'InitialMod', '=', 8, 0, 'section' )
    endif
    ! ------------------------------------------------------------ !
    ! --- [1]  Initialize ( Variables / Comm. / RandomSeed )   --- !
    ! ------------------------------------------------------------ !
    call LoadPEinfo
    if ( myRank.eq.0 ) call displaySectionTitle( 'INITIALIZATION', '-', 4, 4, 'subsection' )
    call InitVariables
    call InitRandomSeed

    call allocExchangeTable_Field    ( LJs )
    call allocExchangeTable_Current  ( LJs )
    call allocExchangeTable_Relocated( LJs )
    call setupExchangeTable_Field    ( LIs, LJs, Boundary1__EB )
    call setupExchangeTable_Current  ( LIs, LJs, Boundary1__jc )
    call setupExchangeTable_Relocated( LIs, LJs, Boundary1__EB )
    call setupExchangeTable_Particle
    
    ! ------------------------------------------------------------ !
    ! --- [2] ( Load SaveFile ) or ( LoadParameter & setupGS ) --- !
    ! ------------------------------------------------------------ !
    if ( Flag__LoadSave ) then
       ! -- [2-1] Continue From SaveFile -- !
       call ReadSavePoint
       call GridSetting
       call particleBoundary
    else
       ! -- [2-2] Load New Profile       -- !
       call LoadParameter
       call GridSetting
       call setupGS
    endif
    if ( myRank.eq.0 ) then
       call displaySectionTitle( 'InitialMod -- END --', '=', 8, 0, 'section' )
       write(6,*)
       write(6,*)
    endif
    return
  end subroutine Initiator


  ! =================================================================== !
  ! ===  LoadParamter  :: Load parameter.dat                        === !
  ! =================================================================== !
  subroutine LoadParameter
    use constants , only : equDir
    use constants , only : npt, ns, LI, LJ, myRank, PEtot
    use constants , only : qeme, mr, wce, cv, wpewce, vthcv, TiTe
    use constants , only : qm, q, wp, np, npc, rm
    use constants , only : dr, dz, dt, lDebye, dr_Debye, dz_Debye, drinv, dzinv
    use constants , only : rMin, rMax, zMin, zMax
    use constants , only : vthi, vthe, isoThr, rhofloor, vAlfe
    use constants , only : alpha_wce, alpha_wpe, alpha_CFL
    use constants , only : er_, et_, ez_, br_, bt_, bz_
    use displayMod, only : displaySectionTitle
    implicit none
    include 'mpif.h'
    integer             :: LIr, LJr, PEtoth, ierr
    double precision    :: dt_CFL, dt_wpe, dt_wce, ppcAve
    character(100)      :: FileName, buff

    ! ------------------------------ !
    ! --- [1] Load parameter.dat --- !
    ! ------------------------------ !
    if ( myRank.eq.0 ) then
       call displaySectionTitle( 'LOAD PARAMETER', '-', 4, 4, 'subsection' )
       !    -- [1-1] Read Parameters From 'parameter.dat' --  !
       FileName = trim(equDir) // '/parameter.dat'
       open (30, file=trim(FileName), status='old', form='formatted' )
       read (30,'(2(a12,2x),i15)'  ) buff, buff, LIr
       read (30,'(2(a12,2x),i15)'  ) buff, buff, LJr
       read (30,'(2(a12,2x),e15.8)') buff, buff, dz
       read (30,'(2(a12,2x),e15.8)') buff, buff, dr
       read (30,'(2(a12,2x),e15.8)') buff, buff, dz_Debye
       read (30,'(2(a12,2x),e15.8)') buff, buff, dr_Debye
       read (30,'(2(a12,2x),e15.8)') buff, buff, wpewce
       read (30,'(2(a12,2x),e15.8)') buff, buff, vthcv
       read (30,'(2(a12,2x),e15.8)') buff, buff, valfe
       read (30,'(2(a12,2x),e15.8)') buff, buff, TiTe
       read (30,'(2(a12,2x),e15.8)') buff, buff, rhofloor
       read (30,'(2(a12,2x),e15.8)') buff, buff, zMin
       read (30,'(2(a12,2x),e15.8)') buff, buff, zMax
       read (30,'(2(a12,2x),e15.8)') buff, buff, rMin
       read (30,'(2(a12,2x),e15.8)') buff, buff, rMax
       read (30,'(2(a12,2x),i15)'  ) buff, buff, PEtoth
       close(30)
       !    -- [1-2] Incompatible Grid Check    --  !
       if ( (LI.ne.LIr).or.(LJ.ne.LJr) ) then
          write(6,*) '  [ERROR] LI != LIr  or  LJ != LJr not [ERROR]'
          write(6,*) '         LI == ', LI, '<-> LIr == ', LIr
          write(6,*) '         LJ == ', LJ, '<-> LJr == ', LJr
          call MPI_Finalize( ierr )
          stop
       endif
       !    -- [1-3] Incompatible # of PE Check --  !
       if ( PEtot.ne.PEtoth ) then
          write(6,*) '  [ERROR] PEtot incompatible [ERROR]'
          write(6,*) '    PEtot<PIC> == ', PEtot, '<==>   PEtot<equ> == ', PEtoth
          call MPI_Finalize( ierr )
          stop
       endif
    endif
    !    -- [1-4] BroadCast parameters -- !
    call MPI_Bcast( dz      , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( dr      , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( dz_Debye, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( dr_Debye, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( wpewce  , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( vthcv   , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( valfe   , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( TiTe    , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( rhofloor, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( zMin    , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( zMax    , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( rMin    , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( rMax    , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

    ! -------------------------- !
    ! --- [2] Set Parameters --- !
    ! -------------------------- !
    q(1)    = -1.d0
    q(2)    = +1.d0
    qm(1)   = -qeme
    qm(2)   = +qeme / mr
    wp(1)   = wce   * wpewce
    wp(2)   = wp(1) / sqrt( mr )
    npc(1)  = npt   / 2
    npc(2)  = npt   / 2
    npt     = np(1) + np(2)
    np(1)   = npc(1)
    np(2)   = npc(2)
    rm(1)   = q(1)  / qm(1)
    rm(2)   = q(2)  / qm(2)
    vthe    = vthcv * cv
    vthi    = vthe  * sqrt( TiTe / mr )
    lDebye  = vthe  / wp(1)
    drinv   = 1.d0  / dr
    dzinv   = 1.d0  / dz
    
    ! ------------------------------------ !
    ! --- [3] Determin Time Step :: dt --- !
    ! ------------------------------------ !
    dt_wce  = alpha_wce / wce
    dt_wpe  = alpha_wpe / wp(1)
    dt_CFL  = alpha_CFL / ( cv * sqrt( drInv**2 + dzInv**2 ) )
    dt      = min( dt_CFL, dt_wpe, dt_wce )
    ppcAve  = dble(npt) / dble( LI*LJ )

    ! ------------------------------------ !
    ! --- [4] Display Parameters & dt  --- !
    ! ------------------------------------ !
    if ( myRank.eq.0 ) then
       ! -- [4-1] Show Parameters on Display -- !
       write(6,*)
       write(6,'(2x,a)') '---------------         [     Parameter     ]          ---------------'
       write(6,'(2x,a)') '-- constants -- '
       write(6,'(4x,2(a10  ,1x))'    ) 'LI', 'LJ'
       write(6,'(4x,2(i10  ,1x))'    )  LI ,  LJ
       write(6,'(4x,4(a10  ,1x))'    ) 'dr', 'dz', 'dr_Debye', 'dz_Debye'
       write(6,'(4x,4(f10.5,1x))'    )  dr ,  dz ,  dr_Debye ,  dz_Debye
       write(6,'(4x,3(a10  ,1x))'    ) '   zMin   ', '   ---   ', '   zMax   '
       write(6,'(4x,f10.5,10x,f10.5)')  zMin, zMax
       write(6,'(4x,7(a10  ,1x))'    ) '   rMin   ', '   ---   ', '   rMax   '
       write(6,'(4x,f10.5,10x,f10.5)')  rMin, rMax
       write(6,'(4x,5(a10  ,1x))'    ) 'wpewce', 'vthcv', 'vAlfe', 'TiTe', 'rhofloor'
       write(6,'(4x,5(f10.3,1x))'    )  wpewce,   vthcv,   vAlfe,   TiTe,   rhofloor
       write(6,*) 
       write(6,'(2x,a)') '-- electrons -- '
       write(6,'(4x,5(a10  ,1x))') 'qe', 'qeme', 'wpe', 'vthe', 'me'
       write(6,'(4x,5(f10.5,1x))')  q(1), qm(1), wp(1), vthe, rm(1)
       write(6,*) 
       write(6,'(2x,a)') '--   ions    -- '
       write(6,'(4x,5(a10  ,1x))') 'qi', 'qimi', 'wpi', 'vthi', 'mi'
       write(6,'(4x,5(f10.5,1x))')  q(2), qm(2), wp(2), vthi, rm(2)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,*)
       
       ! -- [4-2] Output Particle Parameters -- !
       write(6,*)
       write(6,'(2x,a)') '---------------         [     ppc &  DT     ]          ---------------'
       write(6,*)
       write(6,'(10x,a)'            ) '--- Particle / Cell  ---'
       write(6,'(10x,a,1(f12.6,1x))') ' p.p.c. ( Average )        :: ', ppcAve
       write(6,*)
       if ( ppcAve.lt.100.d0 ) then
          write(6,*)
          write(6,'(2x,a)') '!!! ========================================================== !!!'
          write(6,'(2x,a)') ' [ALLERT] Particle / cell ( averaged ) is less than 100 [ALLERT]'
          write(6,'(2x,a)') '!!! ========================================================== !!!'
       endif
       ! -- [4-3] DeterMin Time Step :: dt   -- !
       write(6,*)
       write(6,'(10x,a)'            ) '---    Time  Step    ---'
       write(6,'(10x,a,1(f12.6,1x))') '     dt (   CFL   )        :: ', dt_CFL
       write(6,'(10x,a,1(f12.6,1x))') '     dt (  1/wpe  )        :: ', dt_wpe
       write(6,'(10x,a,1(f12.6,1x))') '     dt (  1/wce  )        :: ', dt_wce
       write(6,'(10x,a)'            ) '---------------------------------------------'
       write(6,'(10x,a,1(f12.6,1x))') '       dt:PIC              :: ', dt
       write(6,*)
       write(6,'( 2x,a)') '----------------------------------------------------------------------'

    endif
    return
  end subroutine LoadParameter


  ! =================================================================== !
  ! ===  LoadPEinfo  ::  Load PEinfo.dat                            === !
  ! =================================================================== !
  subroutine LoadPEinfo
    use constants, only : PEtot , myRank, equDir, datDir, savedjob, Flag__LoadSave
    use constants, only : LI, LJ, LIs, LJs, ns, nptMax, ppcMax
    use constants, only : rnk_, frm_, to_, LIs_, LJs_, iFr_, iTo_, inn_, isl_, iel_
    use variables, only : ijDomain, kstep
    implicit none
    include 'mpif.h'
    integer            :: iPE, ierr
    character(  8)     :: cstep
    character(100)     :: FileName, buff
    
    ! ------------------------------------- !
    ! --- [1] Load PEinfo.dat           --- !
    ! ------------------------------------- !
    if ( myRank.eq.0 ) then
       if ( Flag__LoadSave ) then
          write(cstep,'(i8.8)') kstep
          FileName  = 'job/'//trim(savedjob)//'/'//'sav/const/ijDomain'//cstep//'.dat'
       else
          FileName  = trim(equDir)//'/'//'ijDomain.dat'
       endif
       open (50,file=trim(FileName),form='formatted')
       read (50,*) buff
       do iPE=0, PEtot-1
          read(50,*) ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
               &     ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
               &     ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
       enddo
       close(50)
       write(6,'(2x,a,60a)') ' [LoadPEinfo] Load File ::', FileName
       write(6,'(2x,10(a6,1x))') '# PE', 'From', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
    endif
    ! ------------------------------------- !
    ! --- [2] Set PE Info.              --- !
    ! ------------------------------------- !
    call MPI_Bcast( ijDomain, 10*PEtot, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    LIs    =  ijDomain(myRank,LIs_)
    LJs    =  ijDomain(myRank,LJs_)
    nptMax = (LIs-2) * (LJs-2) * ppcMax
    if ( .not.(Flag__LoadSave) ) then
       do iPE=0, PEtot-1
          ijDomain( iPE,iFr_) = ijDomain(    iPE,frm_) + 1
          ijDomain( iPE,iTo_) = ijDomain(    iPE,to_ ) - 1
          ijDomain( iPE,inn_) = ijDomain(    iPE,LIs_) - 2
          ijDomain( iPE,isl_) = 2
          ijDomain( iPE,iel_) = ijDomain(    iPE,LIs_) - 1
       enddo
       ijDomain(      0,iFr_) = ijDomain(      0,frm_)
       ijDomain(PEtot-1,iTo_) = ijDomain(PEtot-1,to_ )
       ijDomain(      0,inn_) = ijDomain(      0,LIs_) - 1
       ijDomain(PEtot-1,inn_) = ijDomain(PEtot-1,LIs_) - 1
       ijDomain(      0,isl_) = 1
       ijDomain(PEtot-1,iel_) = ijDomain(PEtot-1,LIs_)
    end if
    ! ------------------------------------- !
    ! --- [3] write PEinfo.dat          --- !
    ! ------------------------------------- !
    if ( myRank.eq.0 ) then
       !  -- [3-1] Display              --  !
       do iPE=0, PEtot-1
          write(6,'(2x,10(i6,1x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
               &                    ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
               &                    ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
       enddo
       !  -- [3-2] write File           --  !
       write(cstep,'(i8.8)') kstep
       FileName = trim(datDir) // 'ijDomain/' // 'ijDomain'//cstep//'.dat'
       open( 50,file=trim(FileName),form='formatted',status='replace' )
       write(50,'(10(a12,2x))') '# PE', 'From', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
       do iPE=0, PEtot-1
          write(50,'(10(i12,2x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
               &                   ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
               &                   ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
       enddo
       close(50)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(2x,a10,a45,2x,a6)') "* SAVE :: ", trim(FileName), '[ OK ]'
       write(6,*)
    endif
    
    return
  end subroutine LoadPEinfo


  ! =================================================================== !
  ! ===  GridSetting  ::  Grid Making                               === !
  ! =================================================================== !
  subroutine GridSetting
    use constants , only : LI, LJ, LIs, LJs, pi, myRank, PEtot, Flag__LoadSave
    use constants , only : dt, dr, dz, dtdr, drdt, dtdz, dzdt
    use constants , only : ppcMin, ppcMinAlpha, ppcMinThr, zMin, zMax, rMin, rMax
    use variables , only : x1Leng, x1Lengloc, x2Leng, x2Lengloc, volRepres, ijDomain, x1g, x1s
    use variables , only : rf,  rh,   zf, zh, rfInv, rhInv
    use variables , only : gVf, rgVf, rgVfInv, gVh, rgVh, rgVhInv
    use variables , only : wVf, rwVf, rwVfInv, wVh, rwVh, rwVhInv
    use displayMod, only : displaySectionTitle
    implicit none
    include 'mpif.h'
    integer            :: i, j, iPE, ierr
    double precision   :: ALJInv, coef1, coef2
    integer, parameter :: rnk_  =1, frm_  =2, to_  =3, LIs_=4, iFr_  =6
    integer, parameter :: x1Len_=1, x1Frm_=2, x1To_=3

    ! ----------------------- !
    ! --- [1] Define Grid --- !
    ! ----------------------- !
    !  -- [1-1] delta     --  !
    dtdr      =   dt / dr
    drdt      =   dr / dt
    dtdz      =   dt / dz
    dzdt      =   dz / dt
    x1Lengloc = dble(LIs-2) * dz
    x2Lengloc = dble(LJs-2) * dr
    x2Leng    = x2Lengloc
    !  -- [1-2] x1g / x1s --  !
    do iPE=0, PEtot-1
       x1g(iPE,x1Len_) = dble( ijDomain(iPE,LIs_) - 2 ) * dz
       x1g(iPE,x1Frm_) = dble( ijDomain(iPE,iFr_) - 2 ) * dz
       x1g(iPE,x1To_ ) = x1g(iPE,x1Frm_) + x1g(iPE,x1Len_)
    enddo
    x1g(0,x1Frm_) = 0.d0
    allocate( x1s(LIs) )
    do i=1, LIs
       x1s(i) = dz * dble(i-2) + x1g(myRank,x1Frm_) + 0.5d0*dz
    enddo
    !  -- [1-3] x1Leng    --  !
    call MPI_ALLREDUCE( x1Lengloc, x1Leng, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    if ( myRank.eq.0 ) then
       call displaySectionTitle( 'GRID  SETTINGS', '-', 4, 4, 'subsection' )
       write(6,'(2x,a)'        ) ' [GridSetting] '
       write(6,'(2x,3(a9,1x),3(a12,1x))') '< PE #.>','< iFrom >','< iTo >','< x1Leng >','< x1From >', '< x1To >'
       write(6,'(2x,a)'        ) '----------------------------------------------------------------------'
       do iPE=0, PEtot-1
          write(6,'(2x,3(i9,1x),3(e12.5,1x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
               &                                 x1g(iPE,x1Len_),  x1g(iPE,x1Frm_),  x1g(iPE,x1To_)
       enddo
       write(6,'(2x,a)'        ) '----------------------------------------------------------------------'
       write(6,*)
    endif

    !   -- [1-2] Set Grid --  !
    !    - (1) rf         -   !
    do j=1, LJs+1
       rf(j)       = dr * dble(j-2) + rMin
       rh(j)       = dr * dble(j-2) + rMin + 0.5d0*dr
    enddo
    !    - (2) rfInv      -   !
    rfInv(:) = 0.d0
    rhInv(:) = 0.d0
    do j=1, LJs+1
       if ( rf(j).ne.0.d0 ) rfInv(j) = 1.d0 / rf(j)
       if ( rh(j).ne.0.d0 ) rhInv(j) = 1.d0 / rh(j)
    enddo
    !    - (3) zf / zh    -   !
    do i=1, LI
       zf(i)       = dz * dble(i-2)
       zh(i)       = dz * dble(i-2) + 0.5d0*dz
    enddo
    
    ! ------------------------------------------ !
    ! --- [2] Calculate Volume of each cell  --- !
    ! ------------------------------------------ !
    !  -- [2-1] Volume for Full Grid : [gVf] --  !
    if ( rMin.eq.0.d0 ) then
       gVf(1)      = 0.d0
       gVf(2)      = pi * ( rh(  2)**2               ) * dz
    else
       gVf(1)      = 0.d0
       gVf(2)      = pi * ( rh(  2)**2 - rh(   2)**2 ) * dz
    endif
    do j=3, LJs-1
       gVf(j)      = pi * ( rh(  j)**2 - rh( j-1)**2 ) * dz
    enddo
    gVf(LJs)       = pi * ( rh( LJs)**2 - rh(LJs-1)**2 ) * dz
    ! gVf(LJs)       = pi * ( rf( LJs)**2 - rh(LJs-1)**2 ) * dz
    !  -- [2-2] Volume for half grid : [gVh] --  !
    if ( rMin.eq.0.d0 ) then
       gVh( 1)     = pi * ( rf(  1)**2 - rf(   2)**2 ) * dz
    else
       gVh( 1)     = pi * ( rf(  2)**2 - rf(   1)**2 ) * dz
    end if
    do j=2, LJs-1
       gVh(j)      = pi * ( rf(j+1)**2 - rf(   j)**2 ) * dz
    enddo
    gVh(LJs)       = pi * ( (rf(LJs)+dr)**2 - rf(LJs)**2 ) * dz
    !  -- [2-3] Weighted Volume at Half Grid For Cylindrical Geometry : [wVf,wVh] -- !
    !   -        (( See :: J.P.Verboncoeur (2001) J.Comput.Phys.   )) -  !
    if ( rMin.eq.0.d0 ) then
       coef1       = 4.d0 /  3.d0 * pi*dr*dz
       coef2       = 1.d0 / 12.d0 * pi*dr*dz
       do j=3, LJs-1
          wVf(j)   = coef1 * rf(j) + coef2 * ( 4.d0*rf(j-1) + dr ) &
               &                   + coef2 * ( 4.d0*rf(j+1) - dr )
          wVh(j)   = coef1 * rh(j) + coef2 * ( 4.d0*rh(j-1) + dr ) &
               &                   + coef2 * ( 4.d0*rh(j+1) - dr )
       enddo
       ! -- Boundary Region for Full Grid :: wVf -- !
       wVf(   2)   = + 2.d0*pi*dr*dz * ( 1.d0/3.0d0*rf(    2) +  5.d0/ 64.d0*dr ) &
            &        + coef2         * (      4.0d0*rf(    3) -              dr )
       wVf(   3)   =        pi*dr*dz * ( 7.d0/24.d0*rf(    2) + 17.d0/192.d0*dr ) &
            &        + coef1         * (            rf(    3)                   ) &
            &        + coef2         * (      4.0d0*rf(    4) -              dr )
       wVf(LJs-1)  = + coef2         * (      4.0d0*rf(LJs-2) +              dr ) &
            &        + coef1         * (            rf(LJs-1)                   ) &
            &        +      pi*dr*dz * ( 7.d0/24.d0*rf(LJs  ) - 17.d0/192.d0*dr ) 
       wVf(LJs  )  = + coef2         * (      4.0d0*rf(LJs-1) +              dr ) &
            &        + 2.d0*pi*dr*dz * ( 1.d0/3.0d0*rf(LJs  ) -  5.d0/ 64.d0*dr )
       wVf(   1)   = wVf(2)
       ! -- Boundary Region for Half Grid :: wVh -- !
       wVh(    2)  = coef1*rh(    2) + coef2*( 4.d0*rh(    3) - dr )
       wVh(LJs-1)  = coef1*rh(LJs-1) + coef2*( 4.d0*rh(LJs-2) + dr )
       wVh(    1)  = wVh(    2)
       wVh(LJs  )  = wVh(LJs-1)
    endif
    !   -- [2-4] relative Volume against [volRepres] ( == j0repres grid Volume ) -- !
    if ( .not.( Flag__LoadSave ) )  then
       if ( j0Repres.eq.0 ) j0Repres = LJs / 2
       volRepres   = gVh( j0Repres )
    endif
    do j=1, LJs
       rgVf(j)     = gVf(j) / volRepres
       rgVh(j)     = gVh(j) / volRepres
       rwVh(j)     = wVh(j) / volRepres
       rwVf(j)     = wVf(j) / volRepres
    enddo
    do j=2, LJs
       rgVfInv(j)  =   1.d0 / rgVf(j)
       rgVhInv(j)  =   1.d0 / rgVh(j)
       rwVfInv(j)  =   1.d0 / rwVf(j)
       rwVhInv(j)  =   1.d0 / rwVh(j)
    enddo
    if ( rMin.eq.0.d0 ) then
       rgVfInv(0)  = 0.d0
       rgVfInv(1)  = 0.d0        ! or relatively larger edge fluctuations !
       rgVhInv(0)  = 0.d0
       rgVhInv(1)  = 0.d0
       ! rgVfInv(0) = rgVfInv(4)
       ! rgVhInv(0) = rgVhInv(3)
       ! rgVfInv(1) = rgVfInv(3)        ! or relatively larger edge fluctuations !
    else
       rgVhInv(0)  = 0.d0
       rgVfInv(1)  = 1.d0 / rgVf(1)
       rgVfInv(0)  = 0.d0
       rgVfInv(1)  = 0.d0
    endif
    rgVfInv(LJs+1) = 0.d0
    rgVfInv(LJs+2) = 0.d0
    rgVhInv(LJs+1) = 0.d0
    rgVhInv(LJs+2) = 0.d0

    ! --------------------------------- !
    ! --- [3] ppcMinThr Definition  --- !
    ! --------------------------------- !
    if ( ppcMinAlpha.eq.1.d0 ) then
       do j=1, LJs
          ppcMinThr(j) = ppcMin
       enddo
    endif
    if ( ppcMinAlpha.ne.1.d0 ) then
       ALJInv = 1.d0 / ( 0.2d0 * dble( LJs ) )
       do j=1, LJs
          ppcMinThr(j) = ( ( ppcMinAlpha -1.d0 ) / ( dble(j)*ALJInv + 1.d0 )**2 + 1.d0 ) * dble( ppcMin )
       enddo
    endif

    ! --------------------------------- !
    ! --- [4] Display Grid Settings --- !
    ! --------------------------------- !
    if ( myRank.eq.0 ) then
       if ( ( ( rMax-rMin ) - x2Leng ).gt.1.d-8 ) then
          write(6,*) ' [CAUTION] rLeng != rMax - rMin [CAUTION] '
          write(6,'(10x,3(a,1x,e12.5))') 'rMax :: ', rMax, ' rMin :: ', rMin, ' rLeng :: ', x2Leng
          write(6,'(10x,1(a,1x,e12.5))') 'diff :: ', rMax-rMin-x2Leng
       endif
       write(6,'(10x,a)'            ) '---  Simulation Box  ---'
       write(6,'(10x,2(a,f10.6,1x))') 'x1Leng (z)  :: ', x1Leng, '--- x2Leng (r)  :: ', x2Leng
       write(6,'(10x,2(a,f10.6,1x))') 'x1Min(zMin) :: ', zMin  , '--- x1Max(zMax) :: ', zMax
       write(6,'(10x,2(a,f10.6,1x))') 'x2Min(rMin) :: ', rMin  , '--- x2Max(rMax) :: ', rMax
       write(6,*)
       write(6,'(10x,a)'            ) '--- rf & rh ( repr.) ---'
       write(6,'(10x,a,3(f10.6,1x))') 'rf(2), rf(LJs/2), rf(LJs)  :: ', rf(2), rf(LJs/2), rf(LJs)
       write(6,'(10x,a,3(f10.6,1x))') 'rh(2), rh(LJs/2), rh(LJs)  :: ', rh(2), rh(LJs/2), rh(LJs)
       write(6,*)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,*)
       write(6,*)
    endif
    
    return
  end subroutine GridSetting


end module InitialMod
