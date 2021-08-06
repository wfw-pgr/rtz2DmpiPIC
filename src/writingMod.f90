module writingMod
  implicit none
  integer         , parameter   :: nWrite        = 20
  double precision, allocatable :: w_w(:,:,:)
contains
  
  subroutine writeField
    use constants , only : ns, LIs, LJs, q, binDir, k_NumSP, myRank, cRank, PEtot
    use constants , only : er_, et_, ez_, br_, bt_, bz_, jr_, jt_, jz_
    use constants, only  : rnk_, frm_, to_, LIs_, LJs_, iFr_, iTo_, inn_, isl_, iel_
    use variables , only : kstep, ijDomain
    use variables , only : EBf, EBo, Jcr, w_p, neMax, neMin
    use momentsMod, only : nCountUp_CC, nCountUp_SP
    implicit none
    include 'mpif.h'
    integer             :: i, j, iPE, NumSP(LIs,LJs,ns), ierr
    character(8)        :: cStep
    character(100)      :: FileName, binDirh
    double precision    :: neMins, neMaxs

    ! --------------------------------- !
    ! --- [0]    Preparation        --- !
    ! --------------------------------- !
    write(cStep,'(i8.8)') kstep
    binDirh = trim(binDir) // 'kstep' // cStep // '/'
    if ( myRank.eq.0 ) call system( 'mkdir -p ' // trim(binDirh) )
    call system( 'mkdir -p ' // trim(binDirh) )
    allocate( w_w(LIs,LJs,nWrite) )

    ! --------------------------------- !
    ! --- [1] Obtain Density / Flow --- !
    ! --------------------------------- !
    !  -- [1-1] Obtain Density      --  !
    call nCountUp_CC( w_p )
    !  -- [1-2] Density Max/Min     --  !
    neMaxs       = 0.d0
    neMins       = 1.d8
    !$omp parallel default(none) &
    !$omp shared(neMaxs,neMins,w_p,LIs,LJs) private(i,j)
    !$omp do reduction(max:neMaxs) reduction(min:neMins)
    do j=3, LJs-2
       do i=3, LIs-2
          neMaxs = max( neMaxs, w_p(i,j,1) )
          neMins = min( neMins, w_p(i,j,1) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    call MPI_AllReduce( neMins, neMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( neMaxs, neMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    
    ! --------------------------------- !
    ! --- [2]  Substituition        --- !
    ! --------------------------------- !
    !$omp parallel default(none) &
    !$omp shared(w_w,EBf,EBo,Jcr,w_p,q,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          w_w(i,j, 1) = EBf(br_,i,j)
          w_w(i,j, 2) = EBf(bt_,i,j)
          w_w(i,j, 3) = EBf(bz_,i,j)
          w_w(i,j, 4) = EBf(er_,i,j)
          w_w(i,j, 5) = EBf(et_,i,j)
          w_w(i,j, 6) = EBf(ez_,i,j)
          w_w(i,j, 7) = EBo(br_,i,j)
          w_w(i,j, 8) = EBo(bt_,i,j)
          w_w(i,j, 9) = EBo(bz_,i,j)
          w_w(i,j,10) = EBo(er_,i,j)
          w_w(i,j,11) = EBo(et_,i,j)
          w_w(i,j,12) = EBo(ez_,i,j)
          w_w(i,j,13) = Jcr(jr_,i,j,1)
          w_w(i,j,14) = Jcr(jt_,i,j,1)
          w_w(i,j,15) = Jcr(jz_,i,j,1)
          w_w(i,j,16) = Jcr(jr_,i,j,2)
          w_w(i,j,17) = Jcr(jt_,i,j,2)
          w_w(i,j,18) = Jcr(jz_,i,j,2)
          w_w(i,j,19) = w_p(i,j,1)
          w_w(i,j,20) = w_p(i,j,2)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! --------------------------------- !
    ! --- [3]  Write Field / Fluid  --- !
    ! --------------------------------- !
    FileName = trim(binDirh) // 'Field2D' // cStep // '_' // cRank // '.bin'
    open (50+myRank, file  =trim(FileName), form   ='unformatted', &
         &           status='replace'     , convert='LITTLE_ENDIAN')
    write(50+myRank) w_w
    close(50+myRank)
    
    ! --------------------------------- !
    ! --- [4]  Write  PEinfo.dat    --- !
    ! --------------------------------- !
    if ( myRank.eq.0 ) then
       FileName = trim(binDirh) // 'ijDomain' // cStep //  '.dat'
       open( 50,file=trim(FileName),form='formatted',status='replace' )
       write(50,'(10(a12,1x))') '# PE', 'From', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
       do iPE=0, PEtot-1
          write(50,'(10(i12,1x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
               &                   ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
               &                   ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
       enddo
       close(50)
    endif
    
    
    ! --------------------------------- !
    ! --- [4]  Write NumSP.bin      --- !
    ! --------------------------------- !
    if ( mod( kstep, k_NumSP ).eq.0 ) then
       call nCountUp_SP( NumSP )
       FileName = trim(binDirh) // 'NumSP' // cStep // '_' // cRank // '.bin'
       open( 50, file  =trim(FileName), form   ='unformatted',  &
            &    status='replace'     , convert='LITTLE_ENDIAN'  )
       write(50) NumSP
       close(50)
    endif

    deallocate( w_w )
    
    return
  end subroutine writeField
  
  
  subroutine writePrtcl
    use constants, only : binDir, myRank, cRank
    use variables, only : kstep , pxv
    implicit none
    character(8)       :: cStep
    character(100)     :: FileName, binDirh

    ! --- [1] preparation --- !
    write(cStep,'(i8.8)') kstep
    binDirh = trim(binDir) // 'kstep' // cStep // '/'
    if ( myRank.eq.0 ) call system( 'mkdir -p ' // trim(binDirh) )
    FileName = trim(binDirh) // 'Prtcl' // cStep // '_' // cRank // '.bin'

    ! --- [2] write out --- !
    open( 50, file  =trim(FileName), form   ='unformatted',   &
         &    status='replace'     , convert='LITTLE_ENDIAN' )
    write(50) pxv
    close(50)
    
    return
  end subroutine writePrtcl


  subroutine WriteConst
    use constants
    use variables, only : NspRepres, volRepres, x2Leng, x1Leng
    implicit none
    character(100)     :: FileName

    if ( myRank.eq.0 ) then
       ! ---------------------------------- !
       ! ---  [1]  Write constants.dat  --- !
       ! ---------------------------------- !
       FileName = trim(datDir) // 'constants.dat'
       open (40, file=trim(FileName), form='formatted')
       write(40,'(2(a10,2x),i15)'  ) 'LI'       , 'long'  , LI
       write(40,'(2(a10,2x),i15)'  ) 'LJ'       , 'long'  , LJ
       write(40,'(2(a10,2x),i15)'  ) 'ns'       , 'long'  , ns
       write(40,'(2(a10,2x),i15)'  ) 'nptTot'   , 'long'  , nptTot
       write(40,'(2(a10,2x),i15)'  ) 'IterMax'  , 'long'  , itermax
       
       write(40,'(2(a10,2x),e15.8)') 'dt'       , 'double', dt
       write(40,'(2(a10,2x),e15.8)') 'dr'       , 'double', dr
       write(40,'(2(a10,2x),e15.8)') 'dz'       , 'double', dz

       write(40,'(2(a10,2x),e15.8)') 'npe'      , 'double', dble(npc(1))
       write(40,'(2(a10,2x),e15.8)') 'npeL'     , 'double', dble(np (1))
       write(40,'(2(a10,2x),e15.8)') 'qme'      , 'double', qm(1)
       write(40,'(2(a10,2x),e15.8)') 'wpe'      , 'double', wp(1)
       write(40,'(2(a10,2x),e15.8)') 'rme'      , 'double', rm(1)

       write(40,'(2(a10,2x),e15.8)') 'npi'      , 'double', dble(npc(2))
       write(40,'(2(a10,2x),e15.8)') 'npiL'     , 'double', dble(np (2))
       write(40,'(2(a10,2x),e15.8)') 'qmi'      , 'double', qm(2)
       write(40,'(2(a10,2x),e15.8)') 'wpi'      , 'double', wp(2)
       write(40,'(2(a10,2x),e15.8)') 'rmi'      , 'double', rm(2)

       write(40,'(2(a10,2x),e15.8)') 'mr'       , 'double', mr
       write(40,'(2(a10,2x),e15.8)') 'wpewce'   , 'double', wpewce
       write(40,'(2(a10,2x),e15.8)') 'isoThr'   , 'double', isoThr
       write(40,'(2(a10,2x),e15.8)') 'TiTe'     , 'double', TiTe
       write(40,'(2(a10,2x),e15.8)') 'vthcv'    , 'double', vthcv
       write(40,'(2(a10,2x),e15.8)') 'dr_Debye' , 'double', dr_Debye
       write(40,'(2(a10,2x),e15.8)') 'dz_Debye' , 'double', dz_Debye

       write(40,'(2(a10,2x),e15.8)') 'nuwce'    , 'double', nuwce
       write(40,'(2(a10,2x),e15.8)') 'x2Min'    , 'double', rMin
       write(40,'(2(a10,2x),e15.8)') 'SmCoef'   , 'double', SmCoef
       write(40,'(2(a10,2x),e15.8)') 'NspRepres', 'double', NspRepres
       write(40,'(2(a10,2x),e15.8)') 'volRepres', 'double', volRepres

       write(40,'(2(a10,2x),e15.8)') 'x2Leng'   , 'double', x2Leng
       write(40,'(2(a10,2x),e15.8)') 'x1Leng'   , 'double', x1Leng
       write(40,'(2(a10,2x),e15.8)') 'lDebye'   , 'double', lDebye
       write(40,'(2(a10,2x),e15.8)') 'vthe'     , 'double', vthe
       write(40,'(2(a10,2x),e15.8)') 'vthi'     , 'double', vthi
       close(40)

       ! ---------------------------------- !
       ! --- [2] Write k_Intervals.dat  --- !
       ! ---------------------------------- !
       FileName = trim(datDir) // 'k_Intervals.dat'
       open (40, file=trim(FileName), form='formatted')
       write(40,'(2(a10,2x),i15)'  ) 'k_Enrgy'  , 'long'  , k_Enrgy
       write(40,'(2(a10,2x),i15)'  ) 'k_Field'  , 'long'  , k_Field
       write(40,'(2(a10,2x),i15)'  ) 'k_FlxIp'  , 'long'  , k_FlxIp
       write(40,'(2(a10,2x),i15)'  ) 'k_tShow'  , 'long'  , k_tShow
       write(40,'(2(a10,2x),i15)'  ) 'k_Ecrct'  , 'long'  , k_Ecrct
       write(40,'(2(a10,2x),i15)'  ) 'k_NumSP'  , 'long'  , k_NumSP
       write(40,'(2(a10,2x),i15)'  ) 'k_prtcl'  , 'long'  , k_prtcl
       write(40,'(2(a10,2x),i15)'  ) 'k_SaveP'  , 'long'  , k_SaveP
       write(40,'(2(a10,2x),i15)'  ) 'k_ppcCt'  , 'long'  , k_ppcCt
       close(40)
    endif
    
    return
  end subroutine WriteConst


end module writingMod



  ! subroutine writeField
  !   use constants , only : ns, LIs, LJs, q, binDir, k_NumSP, myRank, cRank, PEtot
  !   use constants , only : er_, et_, ez_, br_, bt_, bz_, jr_, jt_, jz_
  !   use constants, only  : rnk_, frm_, to_, LIs_, LJs_, iFr_, iTo_, inn_, isl_, iel_
  !   use variables , only : kstep, ijDomain
  !   use variables , only : EBf, Jcr, w_p, neMax, neMin
  !   use momentsMod, only : nCountUp_CC, nCountUp_SP
  !   implicit none
  !   include 'mpif.h'
  !   integer             :: i, j, iPE, NumSP(LIs,LJs,ns), ierr
  !   character(8)        :: cStep
  !   character(100)      :: FileName, binDirh
  !   double precision    :: neMin_s, neMax_s

  !   ! --------------------------------- !
  !   ! --- [0]    Preparation        --- !
  !   ! --------------------------------- !
  !   allocate( w_w(LIs,LJs,nWrite) )
  !   write(cStep,'(i8.8)') kstep
  !   binDirh = trim(binDir) // 'kstep' // cStep // '/'
  !   ! call MPI_Barrier( MPI_COMM_WORLD, ierr )
  !   ! if ( myRank.eq.0 ) call system( 'mkdir -p ' // trim(binDirh) )
  !   call system( 'mkdir -p ' // trim(binDirh) )
  !   call MPI_Barrier( MPI_COMM_WORLD, ierr )

  !   ! --------------------------------- !
  !   ! --- [1] Obtain Density / Flow --- !
  !   ! --------------------------------- !
  !   !  -- [1-1] Obtain Density      --  !
  !   call nCountUp_CC( w_p )
  !   !  -- [1-2] Density Max/Min     --  !
  !   neMax_s       = 0.d0
  !   neMin_s       = 1.d8
  !   !$omp parallel default(none) &
  !   !$omp shared(neMax_s,neMin_s,w_p,LIs,LJs) private(i,j)
  !   !$omp do reduction(max:neMax_s) reduction(min:neMin_s)
  !   do j=3, LJs-2
  !      do i=3, LIs-2
  !         neMax_s = max( neMax_s, w_p(i,j,1) )
  !         neMin_s = min( neMin_s, w_p(i,j,1) )
  !      enddo
  !   enddo
  !   !$omp end do
  !   !$omp end parallel
  !   call MPI_AllReduce( neMin_s, neMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
  !   call MPI_AllReduce( neMax_s, neMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    
  !   ! --------------------------------- !
  !   ! --- [2]  Substituition        --- !
  !   ! --------------------------------- !
  !   !$omp parallel default(none) &
  !   !$omp shared(w_w,EBf,Jcr,w_p,q,LIs,LJs) private(i,j)
  !   !$omp do
  !   do j=1, LJs
  !      do i=1, LIs
  !         w_w(i,j, 1) = EBf(br_,i,j)
  !         w_w(i,j, 2) = EBf(bt_,i,j)
  !         w_w(i,j, 3) = EBf(bz_,i,j)
  !         w_w(i,j, 4) = EBf(er_,i,j)
  !         w_w(i,j, 5) = EBf(et_,i,j)
  !         w_w(i,j, 6) = EBf(ez_,i,j)
  !         w_w(i,j, 7) = Jcr(jr_,i,j,0)
  !         w_w(i,j, 8) = Jcr(jt_,i,j,0)
  !         w_w(i,j, 9) = Jcr(jz_,i,j,0)
  !         w_w(i,j,10) = Jcr(jr_,i,j,1)
  !         w_w(i,j,11) = Jcr(jt_,i,j,1)
  !         w_w(i,j,12) = Jcr(jz_,i,j,1)
  !         w_w(i,j,13) = Jcr(jr_,i,j,2)
  !         w_w(i,j,14) = Jcr(jt_,i,j,2)
  !         w_w(i,j,15) = Jcr(jz_,i,j,2)
  !         w_w(i,j,16) = w_p(i,j,1)
  !         w_w(i,j,17) = w_p(i,j,2)
  !         w_w(i,j,18) = q(1)*w_p(i,j,1) + q(2)*w_p(i,j,2)
  !      enddo
  !   enddo
  !   !$omp end do
  !   !$omp end parallel
    
  !   ! --------------------------------- !
  !   ! --- [3]  Write Field / Fluid  --- !
  !   ! --------------------------------- !
  !   FileName = trim(binDirh) // 'Field2D' // cStep // '_' // cRank // '.bin'
  !   open (50+myRank, file  =trim(FileName), form   ='unformatted', &
  !        &           status='replace'     , convert='LITTLE_ENDIAN')
  !   write(50+myRank) w_w
  !   close(50+myRank)
    
  !   ! --------------------------------- !
  !   ! --- [4]  Write  PEinfo.dat    --- !
  !   ! --------------------------------- !
  !   if ( myRank.eq.0 ) then
  !      FileName = trim(binDirh) // 'ijDomain' // cStep //  '.dat'
  !      open( 50,file=trim(FileName),form='formatted',status='replace' )
  !      write(50,'(10(a12,1x))') '# PE', 'From', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
  !      do iPE=0, PEtot-1
  !         write(50,'(10(i12,1x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
  !              &                   ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
  !              &                   ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
  !      enddo
  !      close(50)
  !   endif
    
    
  !   ! --------------------------------- !
  !   ! --- [4]  Write NumSP.bin      --- !
  !   ! --------------------------------- !
  !   if ( mod( kstep, k_NumSP ).eq.0 ) then
  !      call nCountUp_SP( NumSP )
  !      FileName = trim(binDirh) // 'NumSP' // cStep // '_' // cRank // '.bin'
  !      open( 50, file  =trim(FileName), form   ='unformatted',  &
  !           &    status='replace'     , convert='LITTLE_ENDIAN'  )
  !      write(50) NumSP
  !      close(50)
  !   endif

  !   deallocate( w_w )
    
  !   return
  ! end subroutine writeField
