module saveRunMod
  implicit none
  logical       , parameter  :: Flag__SortParticle = .true.
contains
  
  subroutine MakeSavePoint( mode )
    use constants
    use variables
    use ptOrderMod, only      : ptReOrdering
    use displayMod, only      : displaySectionTitle
    implicit none
    include 'mpif.h'
    integer                  :: iPE, ierr
    character(8)             :: cstep
    character(100)           :: jobcgFile, paramFile, ptspcFile, fieldFile
    character(100)           :: prtclFile, oldEBFile, mpicgFile, PEinfFile
    character(37)            :: fmt
    character(4), intent(in) :: mode

    ! --------------------------- !
    ! --- [0]  Preparation    --- !
    ! --------------------------- !
    if ( mode.eq.'auto' ) cstep = '00000000'
    if ( mode.eq.'finn' ) write(cstep,'(i8.8)') kstep
    call system( 'mkdir -p '//trim(savDir)//'const/'  )
    call system( 'mkdir -p '//trim(savDir)//'field/'  )
    call system( 'mkdir -p '//trim(savDir)//'oldEB/'  )
    call system( 'mkdir -p '//trim(savDir)//'prtcl/'  )
    call system( 'mkdir -p '//trim(savDir)//'mpicg/'  )
    jobcgFile = trim(savDir)//'const/'//trim(job)//'_jobcg_'//cstep//'.dat'
    paramFile = trim(savDir)//'const/'//trim(job)//'_param_'//cstep//'.bin'
    ptspcFile = trim(savDir)//'const/'//trim(job)//'_ptspc_'//cstep//'.bin'
    fieldFile = trim(savDir)//'field/'//trim(job)//'_field_'//cstep//'_'//cRank//'.bin'
    oldEBFile = trim(savDir)//'oldEB/'//trim(job)//'_oldEB_'//cstep//'_'//cRank//'.bin'
    prtclFile = trim(savDir)//'prtcl/'//trim(job)//'_prtcl_'//cstep//'_'//cRank//'.bin'
    mpicgFile = trim(savDir)//'mpicg/'//trim(job)//'_mpicg_'//cstep//'_'//cRank//'.dat'
    PEinfFile = trim(savDir)//'const/'//          'ijDomain'//cstep//'.dat'
    fmt       = '(2x,"[MakeSavePoint] ",a60," [ OK ]")'
    if ( Flag__SortParticle ) call ptReOrdering( 'iCell' )
    
    ! -- CONFIG FILES -- !
    if ( myRank.eq.0 ) then
       write(6,*)
       if ( myRank.eq.0 ) call displaySectionTitle( 'Make SavePoint', '-', 4, 4, 'subsection' )
       write(6,'(2x,a)'      ) '--------------------------------------------------------------'
       write(6,'(2x,a,i10,a)') '--- MakeSavePoint routine is called at kstep == ', kstep, '---'
       write(6,'(2x,a)'      ) '--------------------------------------------------------------'
       ! --------------------------- !
       ! --- [1] Data Parameters --- !
       ! --------------------------- !
       open (70, file=jobcgFile, status='replace', form='formatted' )
       write(70,'(a)'  ) trim(job)
       write(70,'(a)'  ) trim(jobdir)
       write(70,'(i12)') kstep
       write(70,'(i12)') ns
       write(70,'(i12)') npt
       write(70,'(i12)') LI
       write(70,'(i12)') LJ
       write(70,'(i12)') ppcMax
       write(70,'(i12)') ppcMin
       write(70,'(i12)') ppcMargin
       close(70)
       write(6,fmt) trim(jobcgFile)
       
       ! ------------------------------- !
       ! --- [2] Physical Constants  --- !
       ! ------------------------------- !
       open (70, file=trim(paramFile), status ='replace',      &
            &    form='unformatted'  , convert='LITTLE_ENDIAN' )
       write(70) dt
       write(70) dr
       write(70) dz
       write(70) mr
       write(70) isoThr
       write(70) TiTe
       write(70) wpewce
       write(70) vthcv
       write(70) dr_Debye
       write(70) dz_Debye
       write(70) vthe
       write(70) vthi
       write(70) rMin
       write(70) rMax
       write(70) zMin
       write(70) zMax
       write(70) SmCoef
       write(70) NspRepres
       write(70) volRepres
       close(70)
       write(6,fmt) trim(paramFile)
       
       ! ------------------------------- !
       ! --- [3] Particle Parameters --- !
       ! ------------------------------- !
       open (70, file=trim(ptspcFile), status ='replace',      &
            &    form='unformatted'  , convert='LITTLE_ENDIAN' )
       write(70) dble(np)
       write(70) q
       write(70) qm
       write(70) wp
       close(70)
       write(6,fmt) trim(ptspcFile)
    endif
    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    
    ! ------------------------------- !
    ! --- [4]     Field Data      --- !
    ! ------------------------------- !
    !  -- [4-1]   EBf   --  !
    open (70, file=trim(fieldFile), status ='replace', &
         &    form='unformatted'  , convert='LITTLE_ENDIAN' )
    write(70) EBf
    close(70)
    write(6,fmt) trim(fieldFile)
    !  -- [4-2]   EBo   --  !
    open (70, file=trim(oldEBFile), status ='replace', &
         &    form='unformatted'  , convert='LITTLE_ENDIAN' )
    write(70) EBo
    close(70)
    write(6,fmt) trim(oldEBFile)
    
    ! ------------------------------- !
    ! --- [5]  Particles :: pxv   --- !
    ! ------------------------------- !
    open (70, file=trim(prtclFile), status ='replace' &
         &  , form='unformatted'  , convert='LITTLE_ENDIAN' )
    write(70) pxv
    close(70)
    write(6,fmt) trim(prtclFile)
    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    
    ! ------------------------------- !
    ! --- [6]  MPI Information    --- !
    ! ------------------------------- !
    open (70, file=trim(mpicgFile), status ='replace', form='formatted' )
    write(70,*) myRank
    write(70,*) PEtot
    write(70,*) np (1)
    write(70,*) np (2)
    write(70,*) npc(1)
    write(70,*) npc(2)
    write(70,*) npt
    write(70,*) nptMax
    write(70,*) ijDomain(myRank,frm_)
    write(70,*) ijDomain(myRank,to_ )
    write(70,*) ijDomain(myRank,LIs_)
    write(70,*) ijDomain(myRank,LJs_)
    close(70)
    write(6,fmt) trim(mpicgFile)
    call MPI_Barrier( MPI_COMM_WORLD, ierr )

    ! ------------------------------- !
    ! --- [7]   PE Information    --- !
    ! ------------------------------- !
    if ( myRank.eq.0 ) then
       open (70, file=trim(PEinfFile), status ='replace', form='formatted' )
       write(70,'(10(a12,2x))') '# PE', 'From', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
       do iPE=0, PEtot-1
          write(70,'(10(i12,2x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
               &                   ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
               &                   ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
       enddo
       close(70)
       write(6,fmt) trim(PEinfFile)
    end if
    call MPI_Barrier( MPI_COMM_WORLD, ierr )


    ! ------------------- !
    ! ---     END     --- !
    ! ------------------- !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*) ' ---------------------------      saved     --------------------------- '        
       write(6,*)
       write(6,*)
    endif

    return
  end subroutine MakeSavePoint

  
  subroutine ReadSavePoint
    use constants
    use variables
    implicit none
    include 'mpif.h'
    character(100)                :: jobcgFile, paramFile, ptspcFile, fieldFile
    character(100)                :: mpicgFile, prtclFile, oldEBFile
    character(100)                :: jobr, jobdirr, LoadDir
    character(37)                 :: fmt
    character(8)                  :: cstep
    integer                       :: k, kstepr, ns_r, npt_r, LI_r, LJ_r, LIs_r, LJs_r, iFrom, iTo
    integer                       :: ppcMaxr, ppcMinr, ppcMarginr, ierr, myRank_r, PEtot_r
    integer                       :: iPE, npInfo_s(5), npInfo(5,PEtot)
    double precision              :: mrr, isoThrr, SmCoefr
    double precision              :: npr(ns), qr(ns), qmr(ns), wpr(ns)

    ! -------------------------- !
    ! --- [0]  preparation   --- !
    ! -------------------------- !
    write(cstep,'(i8.8)') kstep
    LoadDir   = 'job/' // trim(savedjob) // '/' // 'sav/'
    jobcgFile = trim(LoadDir)//'const/'//trim(savedjob)//'_jobcg_'//cstep//'.dat'
    paramFile = trim(LoadDir)//'const/'//trim(savedjob)//'_param_'//cstep//'.bin'
    ptspcFile = trim(LoadDir)//'const/'//trim(savedjob)//'_ptspc_'//cstep//'.bin'
    fieldFile = trim(LoadDir)//'field/'//trim(savedjob)//'_field_'//cstep//'_'//cRank//'.bin'
    oldEBFile = trim(LoadDir)//'oldEB/'//trim(savedjob)//'_oldEB_'//cstep//'_'//cRank//'.bin'
    prtclFile = trim(LoadDir)//'prtcl/'//trim(savedjob)//'_prtcl_'//cstep//'_'//cRank//'.bin'
    mpicgFile = trim(LoadDir)//'mpicg/'//trim(savedjob)//'_mpicg_'//cstep//'_'//cRank//'.dat'
    fmt       = '(2x,"[ReadSavePoint] ",a60," [ OK ]")'
    if ( myRank.eq.0 ) write(6,'(2x,a,i10)') '[ReadSavePoint]   called at kstep == ', kstep
    
    ! --------------------------- !
    ! --- [1] Data Parameters --- !
    ! --------------------------- !
    if ( myRank.eq.0 ) then
       open (70, file=trim(jobcgFile), status='old', form='formatted', action='read' )
       read (70,'(a)'  ) jobr
       read (70,'(a)'  ) jobdirr
       read (70,'(i12)') kstepr
       read (70,'(i12)') ns_r
       read (70,'(i12)') npt_r
       read (70,'(i12)') LI_r
       read (70,'(i12)') LJ_r
       read (70,'(i12)') ppcMaxr
       read (70,'(i12)') ppcMinr
       read (70,'(i12)') ppcMarginr
       close(70)
       if (    kstep .ne. kstepr    ) write(6,*) '[CAUTION] Incompatible kstep     :: ', kstep    , kstepr
       if (      ns  .ne. ns_r      ) write(6,*) '[CAUTION] Incompatible ns        :: ', ns       , ns_r
       if (      LI  .ne. LI_r      ) write(6,*) '[CAUTION] Incompatible LI        :: ', LI       , LI_r
       if (      LJ  .ne. LJ_r      ) write(6,*) '[CAUTION] Incompatible LJ        :: ', LJ       , LJ_r
       if (   ppcMax .ne. ppcMaxr   ) write(6,*) '[CAUTION] Incompatible ppcMax    :: ', ppcMax   , ppcMaxr
       if (   ppcMin .ne. ppcMinr   ) write(6,*) '[CAUTION] Incompatible ppcMin    :: ', ppcMin   , ppcMinr
       if ( ppcMargin.ne.ppcMarginr ) write(6,*) '[CAUTION] Incompatible ppcMargin :: ', ppcMargin, ppcMarginr
       write(6,fmt) trim(jobcgFile)
    endif
    
    ! ---------------------------- !
    ! --- [2] Other Parameters --- !
    ! ---------------------------- !
    !  -- [2-1] Parameter Load --  !
    if ( myRank.eq.0 ) then
       open (70, file  =trim(paramFile), status ='old'          , form='unformatted', &
            &    action='read'         , convert='LITTLE_ENDIAN'                      )
       read (70) dt
       read (70) dr
       read (70) dz
       read (70) mrr
       read (70) isoThrr
       read (70) TiTe
       read (70) wpewce
       read (70) vthcv
       read (70) dr_Debye
       read (70) dz_Debye
       read (70) vthe
       read (70) vthi
       read (70) rMin
       read (70) rMax
       read (70) zMin
       read (70) zMax
       read (70) SmCoefr
       read (70) NspRepres
       read (70) volRepres
       close(70)
       if ( mrr    .ne. mr          ) write(6,*) '[CAUTION] Incompatible mr        :: ', mr    , mrr
       if ( isoThr .ne. isoThrr     ) write(6,*) '[CAUTION] Incompatible isoThr    :: ', isoThr, isoThrr
       if ( SmCoef .ne. SmCoefr     ) write(6,*) '[CAUTION] Incompatible SmCoef    :: ', SmCoef, SmCoefr
    endif
    call MPI_Bcast( dt       , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( dz       , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( dr       , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( TiTe     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( vthcv    , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( wpewce   , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( dr_Debye , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( dz_Debye , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( vthe     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( vthi     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( zMin     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( zMax     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( rMin     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( rMax     , 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( NspRepres, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( volRepres, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    !  -- [2-2] Calculation -- !
    drInv     = 1.d0 / dr
    drdt      =   dr / dt
    dtdr      =   dt / dr
    dzInv     = 1.d0 / dz
    dzdt      =   dz / dt
    dtdz      =   dt / dz
    valfe     = 1.d0 / wpewce
    RhoPerNsp = 1.d0 / NspRepres
    if ( myRank.eq.0 ) write(6,fmt) trim(paramFile)

    ! -------------------------------------------------- !
    ! --- [3] Particle Information for each Species  --- !
    ! -------------------------------------------------- !
    !  -- [3-1] Data Load -- !
    if ( myRank.eq.0 ) then
       open (70, file  =trim(ptspcFile), status='old', form   ='unformatted',  &
            &    action='read'         , recl  =8*ns , convert='LITTLE_ENDIAN' )
       read (70) npr
       read (70) qr
       read (70) qmr
       read (70) wpr
       close(70)
    endif
    call MPI_Bcast( qr , ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( qmr, ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( wpr, ns, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    !  -- [3-2] Data Substition -- !
    do k=1, ns
       q(k)  = qr(k)
       qm(k) = qmr(k)
       wp(k) = wpr(k)
       rm(k) = q(k) / qm(k)
    enddo
    lDebye   = vthe / wp(1)
    if ( myRank.eq.0 ) write(6,fmt) trim(ptspcFile)
    
    ! ----------------------- !
    ! --- [4]    Field    --- !
    ! ----------------------- !
    ! -- [4-1]  EBf  -- !
    open (70, file  =trim(fieldFile), status='old'      , form   ='unformatted' &
         &  , action='read'         , recl  =8*6*LIs*LJs, convert='LITTLE_ENDIAN' )
    read (70) EBf
    close(70)
    ! -- [4-2]  EBo  -- !
    open (70, file  =trim(oldEBFile), status='old'      , form   ='unformatted' &
         &  , action='read'         , recl  =8*6*LIs*LJs, convert='LITTLE_ENDIAN' )
    read (70) EBo
    close(70)
    if ( myRank.eq.0 ) then
       write(6,fmt) trim(oldEBFile)
       write(6,fmt) trim(fieldFile)
    endif
    
    ! ----------------------- !
    ! --- [5]  MPI Info.  --- !
    ! ----------------------- !
    open (70, file=trim(mpicgFile), status='old', form='formatted' )
    read (70,*) myRank_r
    read (70,*) PEtot_r
    read (70,*) np(1)
    read (70,*) np(2)
    read (70,*) npc(1)
    read (70,*) npc(2)
    read (70,*) npt
    read (70,*) nptMax
    read (70,*) iFrom
    read (70,*) iTo
    read (70,*) LIs_r
    read (70,*) LJs_r
    close(70)
    npt_r = npc(1) + npc(2)
    if ( myRank.ne.myRank_r ) write(6,*) '[CAUTION] Incompatible myRank    :: ', myRank, myRank_r
    if ( PEtot .ne.PEtot_r  ) write(6,*) '[CAUTION] Incompatible PEtot     :: ', PEtot , PEtot_r 
    if ( LIs   .ne.LIs_r    ) write(6,*) '[CAUTION] Incompatible LIs       :: ', LIs   , LIs_r   
    if ( LJs   .ne.LJs_r    ) write(6,*) '[CAUTION] Incompatible LJs       :: ', LJs   , LJs_r   
    if ( npt   .ne.npt_r    ) write(6,*) '[CAUTION] Incompatible npt       :: ', npt   , npt_r   
    if ( myRank.eq.0        ) write(6,fmt) trim(mpicgFile)
    ! -- Re:Allocate pxv -- !
    deallocate( pxv, goban, pCellIdx )
    allocate( pxv(8,nptMax,ns), goban(nptMax), pCellIdx(nptMax,ns) )

    ! ----------------------- !
    ! --- [6]  Particle   --- !
    ! ----------------------- !
    if ( myRank.eq.0 ) write(6,'(2x,a)') '[ReadSavePoint]   particle File  :: Now Loading....'
    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    open( 70, File  =trim(prtclFile), status='old'     , form   ='unformatted'  , &
         &    action='read'         , recl  =8*8*nptMax, convert='LITTLE_ENDIAN'  )
    read( 70) pxv
    close(70)
    if ( myRank.eq.0 ) write(6,fmt) trim(prtclFile)
    call MPI_Barrier( MPI_COMM_WORLD, ierr )

    ! ----------------------- !
    ! --- [7]   Message   --- !
    ! ----------------------- !
    npInfo_s(1) = np (1)
    npInfo_s(2) = np (2)
    npInfo_s(3) = npc(1)
    npInfo_s(4) = npc(2)
    npInfo_s(5) = nptMax
    call MPI_Gather( npInfo_s, 5             , MPI_INTEGER, &
         &           npInfo  , 5             , MPI_INTEGER, &
         &           0       , MPI_COMM_WORLD, ierr     )
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'           ) '[ReadSavePoint]'
       write(6,'(2x,a)'           ) '----------           Loaded particle Information.           ----------'
       write(6,'(10x,a)'          ) '-- electrons --'
       write(6,'(10x,5(a8   ,1x))') 'q', 'qm', 'wp', 'vthe', 'me'
       write(6,'(10x,5(f10.4,1x))')  q(1), qm(1), wp(1), vthe, rm(1)
       write(6,'(10x,a)'          ) '--   ions    --'
       write(6,'(10x,5(a8   ,1x))') 'q', 'qm', 'wp', 'vthi', 'mi'
       write(6,'(10x,5(f10.4,1x))')  q(2), qm(2), wp(2), vthi, rm(2)
       write(6,'(2x,a)'           ) '----------------------------------------------------------------------'
       write(6,*)
       write(6,'(2x,a)'           ) '[ReadSavePoint]'
       write(6,'(2x,a)'           ) '----------              particle #. on each PE              ----------'
       write(6,'(2x,6(a10,1x))'   ) 'myRank', 'np(1)', 'np(2)', 'npc(1)', 'npc(2)', 'nptMax'
       write(6,'(2x,a)'           ) '----------------------------------------------------------------------'
       do iPE=0, PEtot-1
          write   (6,'(2x,6(i10,1x))')  iPE, npInfo(1,iPE+1), npInfo(2,iPE+1), npInfo(3,iPE+1), &
               &                             npInfo(4,iPE+1), npInfo(5,iPE+1)
       enddo
    endif
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'           ) '----------------------------------------------------------------------'
       write(6,*)
    endif

    
    return
  end subroutine ReadSavePoint

end module saveRunMod
