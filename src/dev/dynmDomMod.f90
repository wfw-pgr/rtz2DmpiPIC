module dynmDomMod
  use constants, only : LI, ns, LIs, PEtot
  implicit none
  integer, parameter :: nCmp  = 6
  integer, parameter :: nFrTo = 10
  integer            :: pBuffSize
  integer            :: nDst, nSrc, dst1, src1
  integer            :: FromTo_(0:PEtot-1,nFrTo), exTable(LI,4)
  integer, parameter :: rnk_=1, frm_=2, to_ =3, LIs_=4, LJs_=5
  integer, parameter :: iFr_=6, iTo_=7, inn_=8, isl_=9, iel_=10
  integer, parameter :: PE1_=1, il1_=2, PE2_=3, il2_=4
  logical, parameter :: Flag__Confirmation = .true.
  logical            :: Flag__initilizeDDD = .false.
  integer            :: export_Items(0:PEtot), export_Index(0:PEtot), export_Enemy(0:PEtot)
  integer            :: import_Items(0:PEtot), import_Index(0:PEtot), import_Enemy(0:PEtot)
  integer            :: npLg(LI), npLk(LI,ns)
  integer, allocatable :: npLs(:,:)
  integer, allocatable :: request_Send(:)  , request_Recv(:)
  integer, allocatable ::  status_Send(:,:),  status_Recv(:,:)
  double precision, allocatable :: pSndBuff(:,:,:)
  double precision, allocatable :: pRcvBuff(:,:,:)
contains
  
  subroutine DynamicDomainDecomposition
    
    implicit none
    integer :: ierr
    call npcRedistribute
    call DefineExchangeTable
    call ParticleExchange
    ! call FieldExchange
    ! call RedefineFromTo
    
    call MPI_Finalize( ierr )
    stop
    return
  end subroutine DynamicDomainDecomposition

  
  subroutine npcRedistribute
    use constants , only : LI, LJ, LIs, ns, np, dzInv, OMPNumThreads, myRank, PEtot
    use constants , only : zp_, wp_
    use variables , only : pxv, FromTo
    use ptOrderMod, only : ptReOrdering
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: npLt(LIs) 
    integer             :: npLW(LIs+2,OMPNumThreads)
    integer             :: recvCnt(0:PEtot-1), displs(0:PEtot-1)
    integer             :: i, m, k, ip, ith, iPE, sendCnt, ierr
    integer             :: isl, iel, cnt1, cnt2, is, nptSum
    double precision    :: nptIdeal, nptLevel(0:PEtot-1)

    ! ----------------------------------------------- !
    ! --- [1] Sort Particle according to  i-Grid  --- !
    ! ----------------------------------------------- !
    call ptReOrdering( 'iGrid' )
    
    ! ----------------------------------------------- !
    ! --- [2] Count up #.of pt. on ith-Grid Line  --- !
    ! ----------------------------------------------- !
    !  -- [2-1] Count up k-th species :: npLs     --  !
    allocate( npLs(LIs,ns) )
    npLs = 0
    do k=1, ns
       if ( OMPNumThreads.eq.1 ) ith  = 1
       ! ith  = 1
       npLW = 0
       !$omp parallel default(none) &
       !$omp shared(k,np,pxv,dzInv,npLW) private(m,ip,ith)
       !$ ith = omp_get_thread_num() + 1
       !$omp do
       do m=1, np(k)
          ip           = ceiling( pxv(zp_,m,k)*dzInv ) + 1
          npLW(ip,ith) = npLW(   ip,ith) + 1
       enddo
       !$omp end do
       !$omp end parallel
       do ith=1, OMPNumThreads
          do i=1, LIs
             npLs(i,k) = npLs(i,k) + npLW(i,ith)
          enddo
       enddo
    enddo
    !  -- [2-2] Reduction in species direction    --  !
    npLt = 0
    do k=1, ns
       do i=1, LIs
          npLt(i) = npLt(i) + npLs(i,k)
       enddo
    enddo
    
    ! ----------------------------------------------- !
    ! --- [3] Communication of npLg  ( All PE )   --- !
    ! ----------------------------------------------- !
    !  -- [3-1]   Preparation                     --  !
    isl     = FromTo(myRank,isl_)
    iel     = FromTo(myRank,iel_)
    sendCnt = FromTo(myRank,inn_)
    do iPE=0, PEtot-1
       recvCnt(iPE) = FromTo(iPE,inn_)
       displs (iPE) = FromTo(iPE,iFr_) - 1
    enddo
    !  -- [3-2]   Exchange npLt                   --  !
    npLg   = 0
    call MPI_AllGatherV( npLt(isl:iel)  , sendCnt,         MPI_INTEGER, &
         &               npLg           , recvCnt, displs, MPI_INTEGER, &
         &               MPI_COMM_WORLD , ierr                          )
    do k=1, ns
       call MPI_AllGatherV( npLs(isl:iel,k), sendCnt,         MPI_INTEGER, &
            &               npLk(      :,k), recvCnt, displs, MPI_INTEGER, &
            &               MPI_COMM_WORLD , ierr                          )
    enddo
    
    ! ----------------------------------------------- !
    ! --- [4] Redistribution of LIs by npt        --- !
    ! ----------------------------------------------- !
    !  -- [4-1] Define Criterion Level            --  !
    nptSum   = sum ( npLg   )
    nptIdeal = dble( nptSum ) / dble( PEtot )
    do iPE=0, PEtot-1
       nptLevel(iPE)  = dble(iPE+1) * nptIdeal
    enddo
    nptLevel(PEtot-1) = nptSum
    !  -- [4-2] Determine iFr & iTo               --  !
    cnt1 = 0
    cnt2 = 0
    is   = 1
    do iPE=0, PEtot-1
       do ip=is, LI
          cnt1 = cnt2
          cnt2 = cnt2 + npLg(ip)
          if ( cnt2.ge.nptLevel(iPE) ) exit
       enddo
       ! -- primitive -- !
       ! FromTo_(iPE,iTo_) = ip
       ! --   closer  -- !
       if ( ( nptLevel(iPE)-cnt1 ).lt.( nptLevel(iPE)-cnt2 ) ) then
          FromTo_(iPE,iTo_ ) = ip - 1
          cnt2 = cnt1
       else
          FromTo_(iPE,iTo_ ) = ip
       endif
       is = ip + 1
    enddo
    !  -- [4-3] Determine FromTo                  --  !
    FromTo_(      0,iFr_) = 1
    FromTo_(PEtot-1,iTo_) = LI
    do iPE=1, PEtot-1
       FromTo_(iPE,iFr_) = FromTo_(iPE-1,iTo_ ) + 1
    enddo
    do iPE=0, PEtot-1
       FromTo_(iPE,rnk_) = iPE
       FromTo_(iPE,inn_) = FromTo_(iPE,iTo_) - FromTo_(iPE,iFr_) + 1
       FromTo_(iPE,frm_) = FromTo_(iPE,iFr_) - 1
       FromTo_(iPE,to_ ) = FromTo_(iPE,iTo_) + 1
       FromTo_(iPE,LIs_) = FromTo_(iPE,inn_) + 2
       FromTo_(iPE,LJs_) = LJ
       FromTo_(iPE,isl_) = 2
       FromTo_(iPE,iel_) = FromTo_(iPE,LIs_) - 1  
    enddo
    FromTo_(      0,frm_) = FromTo_(      0,iFr_)
    FromTo_(      0,LIs_) = FromTo_(      0,inn_) + 1
    FromTo_(      0,isl_) = 1
    FromTo_(      0,iel_) = FromTo_(      0,LIs_) - 1
    FromTo_(PEtot-1,to_ ) = FromTo_(PEtot-1,iTo_)
    FromTo_(PEtot-1,LIs_) = FromTo_(PEtot-1,inn_) + 1
    FromTo_(PEtot-1,iel_) = FromTo_(PEtot-1,LIs_)
    
    return
  end subroutine npcRedistribute


  subroutine DefineExchangeTable
    use constants, only : LI, PEtot, myRank
    use variables, only : FromTo
    implicit none
    include 'mpif.h'
    integer            :: i, idx, iPE
    integer            :: iFr, iTo

    ! ---------------------------------- !
    ! --- [1] Definition of exTable  --- !
    ! ---------------------------------- !
    do iPE=0, PEtot-1
       do i=1, FromTo(iPE,inn_)
          idx = FromTo(iPE,iFr_) + (i-1)
          exTable(idx,PE1_) = iPE
          exTable(idx,il1_) = FromTo(iPE,isl_)  + (i-1)
       enddo
    enddo
    do iPE=0, PEtot-1
       do i=1, FromTo_(iPE,inn_)
          idx = FromTo_(iPE,iFr_) + (i-1)
          exTable(idx,PE2_) = iPE
          exTable(idx,il2_) = FromTo_(iPE,isl_) + (i-1)
       enddo
    enddo
    ! ---------------------------------- !
    ! --- [2] Generalized Comm.Table --- !
    ! ---------------------------------- !
    !  -- [2-1] Export Table         --  !
    export_Items = 0
    export_Index = 1
    export_Enemy = 0
    iFr     = FromTo(myRank,iFr_)
    iTo     = FromTo(myRank,iTo_)
    dst1    = exTable(iFr,PE2_)
    nDst    = exTable(iTo,PE2_) - exTable(iFr,PE2_) + 1
    do i=iFr, iTo
       idx  = exTable(i,PE2_) - dst1 + 1
       export_Items(idx) = export_Items(idx) + 1
    enddo
    do i=1, nDst
       export_Enemy(i) = dst1 + (i-1)
       export_Index(i) = export_Index(i-1) + export_Items(i-1)
    enddo
    !  -- [2-2] Import Table         --  !
    import_Items = 0
    import_Index = 1
    import_Enemy = 0
    iFr     = FromTo_(myRank,iFr_)
    iTo     = FromTo_(myRank,iTo_)
    src1    = exTable(iFr,PE1_)
    nSrc    = exTable(iTo,PE1_) - exTable(iFr,PE1_) + 1
    do i=iFr, iTo
       idx  = exTable(i,PE1_) - src1 + 1
       import_Items(idx) = import_Items(idx) + 1
    enddo
    do i=1, nSrc
       import_Enemy(i) = src1 + (i-1)
       import_Index(i) = import_Index(i-1) + import_Items(i-1)
    enddo
    !  -- [2-3] Request / Status     --  !
    allocate( request_Send(nDst), status_Send(MPI_STATUS_SIZE,nDst) )
    allocate( request_Recv(nSrc), status_Recv(MPI_STATUS_SIZE,nSrc) )

    !  -- [2] Confirmation  -- !
    if ( Flag__Confirmation ) then
       if ( myRank.eq.0 ) then
          write(6,*)
          write(6,'(5x,a)') '-- FromTo  ( Old ) --'
          write(6,'(10(a8,1x))') 'rank', 'Frm', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
          do iPE=0, PEtot-1
             write(6,'(10(i8,1x))') FromTo (iPE,rnk_), FromTo (iPE,frm_), FromTo (iPE,to_ ), &
                  &                 FromTo (iPE,LIs_), FromTo (iPE,LJs_), FromTo (iPE,iFr_), FromTo (iPE,iTo_), &
                  &                 FromTo (iPE,inn_), FromTo (iPE,isl_), FromTo (iPE,iel_)
          enddo
          write(6,*)
          write(6,'(5x,a)') '-- FromTo_ ( New ) --'
          write(6,'(10(a8,1x))') 'rank', 'Frm', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
          do iPE=0, PEtot-1
             write(6,'(10(i8,1x))') FromTo_(iPE,rnk_), FromTo_(iPE,frm_), FromTo_(iPE,to_ ), &
                  &                 FromTo_(iPE,LIs_), FromTo_(iPE,LJs_), FromTo_(iPE,iFr_), FromTo_(iPE,iTo_), &
                  &                 FromTo_(iPE,inn_), FromTo_(iPE,isl_), FromTo_(iPE,iel_)
          enddo
          write(6,*)
          write(6,'(5x,a)') '-- ExchangeTable -- '
          write(6,'(5(a10,1x))') 'i(global)', 'PE1', 'il1', 'PE2', 'il2'
          do i=1, LI
             write(6,'(5(i10,1x))') i, exTable(i,PE1_), exTable(i,il1_), exTable(i,PE2_), exTable(i,il2_)
          enddo
       endif
    endif
    
    return
  end subroutine DefineExchangeTable

  
  subroutine FieldExchange
    use constants, only : LIs, LJs, myRank, PEtot
    use variables, only : FromTo, EBf
    implicit none
    include 'mpif.h'
    integer :: ierr
    integer :: i, j, cmp, idx, isl
    integer :: iDst, iSrc, is, ir, len_s, len_r
    double precision, allocatable :: sendBuff(:,:,:), recvBuff(:,:,:)
    
    ! ------------------------- !
    ! --- [1]  Packing      --- !
    ! ------------------------- !
    !  -- [1-1] Allocation  --  !
    allocate( sendBuff(nCmp,LJs,FromTo (myRank,inn_)) )
    allocate( recvBuff(nCmp,LJs,FromTo_(myRank,inn_)) )
    !  -- [1-2] sendBuff    --  !
    isl = FromTo(myRank,isl_)
    do i=1, FromTo(myRank,inn_)
       do j=1, LJs
          do cmp=1, nCmp
             sendBuff(cmp,j,i) = EBf(cmp,isl+(i-1),j)
          enddo
       enddo
    enddo

    ! ------------------------- !
    ! --- [2]  Exchange     --- !
    ! ------------------------- !
    !  -- [2-1] Send        --  !
    do iDst=1, nDst
       is    = export_Index(iDst)
       len_s = export_Items(iDst)*nCmp*LJs
       call MPI_Isend( sendBuff(1,1,is)  , len_s, MPI_DOUBLE_PRECISION, &
            &          export_Enemy(iDst), 0    , MPI_COMM_WORLD      , &
            &          request_Send(iDst), ierr )
    enddo
    !  -- [2-2] Recv        --  !
    do iSrc=1, nSrc
       ir    = import_Index(iSrc)
       len_r = import_Items(iSrc)*nCmp*LJs
       call MPI_Irecv( recvBuff(1,1,ir)  , len_r, MPI_DOUBLE_PRECISION, &
            &          import_Enemy(iSrc), 0    , MPI_COMM_WORLD      , &
            &          request_recv(iSrc), ierr )
    enddo
    
    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll     --  !
    call MPI_WaitAll( nSrc, request_Recv, status_Recv, ierr )
    call MPI_WaitAll( nDst, request_Send, status_Send, ierr )
    !  -- [3-2] Re-Allocate --  !
    LIs = FromTo_(myRank,LIs_)
    LJs = FromTo_(myRank,LJs_)
    deallocate( EBf )
    allocate( EBf(nCmp,LIs,LJs) )
    do j=1, LJs
       do i=1, LIs
          do cmp=1, nCmp
             EBf(cmp,i,j) = 0.d0
          enddo
       enddo
    enddo
    !  -- [3-3] Update EBf  --  !
    ! do i=FromTo_(myRank,iFr_), FromTo_(myRank,iTo_)
    do i=1, FromTo_(myRank,inn_)
       do j=1, LJs
          do cmp=1, nCmp
             idx = FromTo_(myRank,isl_) + (i-1)
             EBf(cmp,idx,j) = recvBuff(cmp,j,i)
          enddo
       enddo
    enddo

    return
  end subroutine FieldExchange


  subroutine ParticleExchange
    use constants, only : ns, nptMax, myRank, PEtot
    use variables, only : FromTo, pxv
    use pMPIComMod, only : ResizeParticleMemory
    use constants, only : ptExpandRatio
    implicit none
    include 'mpif.h'
    integer            :: i, m, iFr, iTo, k, idx, iDst, iSrc, ierr, iPE, cmp, is, ir, len_s, len_r, iter, Niter(ns)
    integer            :: export_ptItm(0:PEtot,ns), export_ptIdx(0:PEtot,ns)
    integer            :: import_ptItm(0:PEtot,ns), import_ptIdx(0:PEtot,ns)
    integer            :: nptIm(ns), nptEx(ns)
    integer, parameter :: ptCmp = 8
    
    !  -- [2-1] Export Table         --  !
    iFr     = FromTo(myRank,iFr_)
    iTo     = FromTo(myRank,iTo_)
    export_ptItm = 0
    export_ptIdx = 1
    nptEx        = 0
    do i=iFr, iTo
       idx = exTable(i,PE2_) - dst1 + 1
       do k=1, ns
          export_ptItm(idx,k)  = export_ptItm(idx,k) + npLk(i,k)
       enddo
    enddo
    do k=1, ns
       do iDst=1, nDst
          export_ptIdx(iDst,k) = export_ptIdx(iDst-1,k) + export_ptItm(iDst-1,k)
          nptEx(k)             = nptEx(k) + export_ptItm(iDst,k)
       enddo
    enddo
    
    !  -- [2-2] Import Table         --  !
    iFr    = FromTo_(myRank,iFr_)
    iTo    = FromTo_(myRank,iTo_)
    import_ptItm = 0
    import_ptIdx = 1
    nptIm        = 0
    src1    = exTable(iFr,PE1_)
    do i=iFr, iTo
       idx = exTable(i,PE1_) - src1 + 1
       do k=1, ns
          import_ptItm(idx,k)  = import_ptItm(idx,k) + npLk(i,k)
       enddo
    enddo
    do k=1, ns
       do iSrc=1, nSrc
          import_ptIdx(iSrc,k) = import_ptIdx(iSrc-1,k) + import_ptItm(iSrc-1,k)
          nptIm(k)             = nptIm(k)               + import_ptItm(iSrc  ,k)
       enddo
    enddo
    pBuffSize = max( nptIm(1), nptIm(2) )
    allocate( pSndBuff(ptCmp,pBuffSize,ns) )
    allocate( pRcvBuff(ptCmp,pBuffSize,ns) )
    if ( nptMax.lt.pBuffSize ) then
       pBuffSize = nint( ptExpandRatio * pBuffSize )
       write(6,'(a,i6,2(i12,1x))') "[ResizeParticleMemory] ", myRank, &
            &     nptMax, pBuffSize
       call ResizeParticleMemory( nptMax, pBuffSize )
    endif
    do k=1, ns
       do iter=1, Niter(k)
          ! -- Send -- !
          do iDst=1, nDst
             is    = export_ptIdx(iDst,k)
             len_s = export_ptItm(iDst,k) * ptCmp
             call MPI_iSend( pSndBuff(1,is,k)  , len_s, MPI_DOUBLE_PRECISION, &
                  &          export_enemy(iDst), 0    , MPI_COMM_WORLD      , &
                  &          request_send(iDst), ierr )
          enddo
          ! -- Recv -- !
          do iSrc=1, nSrc
             ir    = import_ptIdx(iSrc,k)
             len_r = import_ptItm(iSrc,k) * ptCmp
             call MPI_iRecv( pRcvBuff(1,ir,k)  , len_r, MPI_DOUBLE_PRECISION, &
                  &          import_enemy(iSrc), 0    , MPI_COMM_WORLD      , &
                  &          request_recv(iSrc), ierr )
          enddo
          ! -- update -- !
          call MPI_WaitAll( nSrc, request_Recv, status_Recv, ierr )
          call MPI_WaitAll( nDst, request_Send, status_Send, ierr )
          do m=1, nptIm(k)
             do cmp=1, ptCmp
                pxv(cmp,m,k) = pRcvBuff(cmp,m,k)
             enddo
          enddo
          ! -- iter end -- !
       enddo
    enddo

    if ( Flag__Confirmation ) then
       open (40,file='out.dat',form='formatted',status='replace')
       close(40)
       do iPE=0, PEtot-1
          if ( iPE.eq.myRank ) then
             open(40,file='out.dat',status='old',position='append')
             do idx=1, nDst
                write(40,'(7(i8,1x))') myRank, idx, export_ptItm(idx,1), export_ptItm(idx,2), export_ptIdx(idx,1), export_ptIdx(idx,2), pBuffSize
             enddo
             close(40)
          endif
          call MPI_Barrier( MPI_COMM_WORLD, ierr )
       enddo
    endif
    if ( Flag__Confirmation ) then
       open (40,file='recv.dat',form='formatted',status='replace')
       close(40)
       do iPE=0, PEtot-1
          if ( iPE.eq.myRank ) then
             open(40,file='recv.dat',status='old',position='append')
             do idx=1, nSrc
                write(40,'(7(i8,1x))') myRank, idx, import_ptItm(idx,1), import_ptItm(idx,2), import_ptIdx(idx,1), import_ptIdx(idx,2), pBuffSize
             enddo
             close(40)
          endif
          call MPI_Barrier( MPI_COMM_WORLD, ierr )
       enddo
    endif


    return
  end subroutine ParticleExchange


  subroutine RedefineFromTo
    use constants, only : datDir, myRank, PEtot
    use variables, only : FromTo, kstep
    implicit none
    integer            :: cmp, iPE
    character(  8)     :: cStep
    character(100)     :: FileName

    do cmp=1, nFrTo
       do iPE=0, PEtot-1
          FromTo(iPE,cmp) = FromTo_(iPE,cmp)
       enddo
    enddo
    if ( myRank.eq.0 ) then
       write(cStep,'(i8.8)') kstep
       FileName = trim(datDir) // 'PEinfo' // cStep //  '.dat'
       open( 50,file=trim(FileName),form='formatted',status='replace' )
       write(50,'(5(a12,2x))') '# PE', 'From', 'To', 'LIloc', 'LJloc'
       do iPE=0, PEtot-1
          write(50,'(5(i12,2x))') iPE, FromTo(iPE,frm_), FromTo(iPE,to_), FromTo(iPE,LIs_), FromTo(iPE,LJs_)
       enddo
       close(50)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(2x,a10,a45,2x,a6)') "* SAVE :: ", trim(FileName), '[ OK ]'
       write(6,*)
    endif

    
    return
  end subroutine RedefineFromTo
  
end module dynmDomMod
