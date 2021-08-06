module pMPIComMod
  implicit none
  include 'mpif.h'
  ! -- Constants        -- !
  integer         , parameter   :: NeibPEtot_Max = 2
  double precision, parameter   :: buffMargin    = 3.d0
  double precision              :: x1End1, x1End2, x1Mid
  double precision              :: x1slide(neibPEtot_Max)
  ! -- Global Variables -- !
  integer                       :: myRank , PEtot
  integer                       :: CommLen, NeibPEtot   , pBuffSize
  integer                       :: NeibPE(NeibPEtot_Max), dst(2)
  logical                       :: Flag__initialized = .false.
  ! -- Allocatable Var. -- !
  integer         , allocatable :: npComm(:,:)
  integer         , allocatable :: export_Index (:,:)  , import_Index (:,:)
  integer         , allocatable :: export_Items (:,:)  , import_Items (:,:)
  integer         , allocatable :: export_IndexW(:,:,:), export_ItemsW(:,:,:)
  integer         , allocatable :: export_AddrsW(:,:,:), import_Addrs (:,:)
  integer         , allocatable :: import_AddrsW(:,:,:)
  integer         , allocatable :: import_pSeat (:)    , import_iFrom (:,:)
  integer         , allocatable :: import_pSeatW(:,:)
  integer         , allocatable :: export_iBase (:,:,:), import_iBase (:,:)
  integer         , allocatable ::  request_Send(:)    ,  request_Recv(:)
  integer         , allocatable ::   status_Send(:,:)  ,   status_Recv(:,:)
  double precision, allocatable ::      sendBuff(:,:,:),      recvBuff(:,:,:)
contains

  
  subroutine DefineExchangeTable_Particle
    use constants, only          : ns, nptMax, LIs, OMPNumThreads
    use constants, only          : Boundary1__pt
    implicit none
    include 'mpif.h'
    integer                   :: ierr
    ! -------------------------------------------- !
    ! --- [1] Define Grid Number and Partition --- !
    ! -------------------------------------------- !
    !  -- [1-1]  PE infomation                 --  !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )

    ! -------------------------------------------- !
    ! --- [2] Allocate Variables for MPI Comm. --- !
    ! -------------------------------------------- !
    pBuffSize    = nint( buffMargin * nptMax * ( 2.0d0 / dble(LIs) ) )
    allocate( sendBuff(8,pBuffSize,NeibPEtot_Max) )
    allocate( recvBuff(8,pBuffSize,NeibPEtot_Max) )
    allocate( request_Send(NeibPEtot_Max), status_Send(MPI_STATUS_SIZE,NeibPEtot_Max) )
    allocate( request_Recv(NeibPEtot_Max), status_Recv(MPI_STATUS_SIZE,NeibPEtot_Max) )
    allocate( export_Index ( 0:ns,0:NeibPEtot_Max), export_Items(ns,NeibPEtot_Max) )
    allocate( import_Index ( 0:ns, NeibPEtot_Max), import_Items(ns,NeibPEtot_Max) )
    allocate( export_IndexW( 0:ns, NeibPEtot_Max, OMPNumThreads ) )
    allocate( export_ItemsW(   ns, NeibPEtot_Max, OMPNumThreads ) )
    allocate( export_iBase (   ns, NeibPEtot_Max, OMPNumThreads ) )
    allocate( import_pSeatW(   ns, OMPNumThreads ) )
    allocate( import_iBase (   ns, OMPNumThreads ) )
    allocate( export_AddrsW(pBuffSize,NeibPEtot_Max,OMPNumThreads) )
    allocate( import_AddrsW(pBuffSize,ns           ,OMPNumThreads) )
    allocate( import_Addrs (pBuffSize,ns                         ) )
    allocate( import_iFrom(ns,0:NeibPEtot_Max), import_pSeat(ns) )
    import_iFrom  = 0
    import_pSeat  = 0
    import_Index  = 0
    export_Index  = 0
    import_Addrs  = 0
    export_AddrsW = 0
     status_Send  = 0
     status_Recv  = 0
    request_Send  = 0
    request_Recv  = 0
    sendBuff      = 0.d0
    recvBuff      = 0.d0
    export_IndexW = 0
    export_ItemsW = 0
    
    ! -------------------------------------------- !
    ! --- [3]  MPI  Communication  Settings    --- !
    ! -------------------------------------------- !
    !  -- [3-1] Reflect Wall Particle Condition -- !
    if ( trim(Boundary1__pt).eq.'reflect' ) then
       if ( myRank.eq.0       ) then
          NeibPEtot  = 1
          NeibPE(1)  = 1
          dst(1)     = 2
          dst(2)     = 1
       endif
       if ( myRank.eq.PEtot-1 ) then
          NeibPEtot  = 1
          NeibPE(1)  = PEtot-2
          dst(1)     = 1
          dst(2)     = 2
       endif
       if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
          NeibPEtot  = 2
          NeibPE(1)  = myRank-1
          NeibPE(2)  = myRank+1
          dst(1)     = 1
          dst(2)     = 2
       endif
    endif
    !  -- [3-2] Periodic Particle Condition     -- !
    if ( trim(Boundary1__pt).eq.'periodic' ) then
       NeibPEtot  = 2
       NeibPE(1)  = myRank-1
       NeibPE(2)  = myRank+1
       dst(1)     = 1
       dst(2)     = 2
       if ( myRank.eq.0       ) NeibPE(1)  = PEtot-1
       if ( myRank.eq.PEtot-1 ) NeibPE(2)  = 0
    endif
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '[DefineExchangeTable_Particle ]          * Particle  :: ( Initialized )'
    endif
    
    return
  end subroutine DefineExchangeTable_Particle


  subroutine setExchangeGeometry
    use constants, only : Boundary1__pt
    use variables, only : x1Lengloc
    implicit none
    x1End1     = 0.0d0
    x1End2     = x1Lengloc
    x1Mid      = 0.5d0*x1Lengloc
    x1slide(1) = + x1Lengloc
    x1slide(2) = - x1Lengloc
    if ( trim(Boundary1__pt).eq.'reflect' ) then
       if ( myRank.eq.0       ) x1slide(1) = - x1Lengloc
       if ( myRank.eq.PEtot-1 ) x1slide(1) = + x1Lengloc
    endif
    return
  end subroutine setExchangeGeometry

  
  subroutine ParticleExchange
    use constants, only : ns , npt, np, npc, OMPNumThreads
    use constants, only : zp_, wp_
    use variables, only : pxv
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer            :: m, k , ith, idx, gidx, cnt, neib , len_s, len_r, ierr
    integer            :: idx1(OMPNumThreads), idx2(OMPNumThreads)
    integer            :: cnt1(OMPNumThreads), cnt2(OMPNumThreads), cnt3(OMPNumThreads)
    integer            :: iaddr, ibuff, npAdd

    if ( .not.( Flag__initialized ) ) call setExchangeGeometry
    ! ------------------------------ !
    ! --- [1]  Inside / OutSide  --- !
    ! ------------------------------ !
    !  -- [1-1] Initialize variables                       -- !
    idx1(1:OMPNumThreads) = 0
    idx2(1:OMPNumThreads) = 0
    export_IndexW(0:ns,1:NeibPEtot_Max,1:OMPNumThreads)  = 0
    export_ItemsW(1:ns,1:NeibPEtot_Max,1:OMPNumThreads)  = 0
    import_pSeatW         = 0
    !  -- [1-2] Exceeding Particle Search ( OMP ver. (W) ) -- !
    !$omp parallel default(none) &
    !$omp shared(cnt1,cnt2,cnt3,idx1,idx2,pxv,x1End1,x1End2,np,dst,myRank) &
    !$omp shared(export_AddrsW,import_AddrsW,export_ItemsW,export_IndexW,import_pSeatW) &
    !$omp private(m,k,ith)
    !$ ith = omp_get_thread_num() + 1
    !$omp do
    do k=1, ns
       cnt1(ith) = 0
       cnt2(ith) = 0
       cnt3(ith) = 0
       do m=1, np(k)
          if ( pxv(wp_,m,k).gt.0.d0   ) then
             ! -- Particle Exist -- !
             !  -  Left   - !
             if ( pxv(zp_,m,k).le.x1End1 ) then
                write(6,*) ith, cnt1(ith), cnt3(ith), idx1(ith)
                cnt1(ith) = cnt1(ith) + 1
                cnt3(ith) = cnt3(ith) + 1
                idx1(ith) = idx1(ith) + 1
                export_AddrsW(idx1(ith),dst(1),ith) = m
                import_AddrsW(cnt3(ith),k     ,ith) = m
             endif
             !  -  Right  - !
             if ( pxv(zp_,m,k).gt.x1End2 ) then
                write(6,*) ith, cnt3(ith), cnt3(ith), idx2(ith)
                cnt2(ith) = cnt2(ith) + 1
                cnt3(ith) = cnt3(ith) + 1
                idx2(ith) = idx2(ith) + 1
                export_AddrsW(idx2(ith),dst(2),ith) = m
                import_AddrsW(cnt3(ith),k     ,ith) = m
             endif
          else
             ! --  No Particle   -- !
             cnt3(ith)   = cnt3(ith) + 1
             import_AddrsW(cnt3(ith),k,ith) = m
          endif
       enddo
       export_ItemsW(k,dst(1),ith) = cnt1(ith)
       export_ItemsW(k,dst(2),ith) = cnt2(ith)
       export_IndexW(k,dst(1),ith) = idx1(ith)
       export_IndexW(k,dst(2),ith) = idx2(ith)
       import_pSeatW(k       ,ith) = cnt3(ith)
    enddo
    !$omp end do
    !$omp end parallel
    
    ! if ( myRank.eq.0 ) then
    !    do k=1, ns
    !       do neib=1, 2
    !          do ith=1, OMPNumThreads
    !             write(6,*) k, neib, ith, export_ItemsW(k,neib,ith), export_IndexW(k,neib,ith)
    !          enddo
    !       enddo
    !    enddo
    ! endif
    
    !  -- [1-3] Obtain export_Items -- !
    export_iBase(:,:,:) = 0
    export_Items(:,:)   = 0
    export_Index(:,:)   = 0
    import_pSeat(:)     = 0
    import_iBase(:,:)   = 0
    do neib=1, NeibPEtot
       cnt = 0
       do k=1, ns
          do ith=1, OMPNumThreads
             export_Items(k,neib)     = export_Items(k,neib  ) + export_ItemsW(k,neib,ith)
             export_iBase(k,neib,ith) = cnt
             cnt                      = cnt                    + export_ItemsW(k,neib,ith)
          enddo
       enddo
    enddo
    !  -- [1-4] #.of Seat for Import -- !
    do k=1, ns
       cnt = 0
       do ith=1, OMPNumThreads
          import_pSeat(k)             = import_pSeat (k      ) + import_pSeatW(k,ith)
          import_iBase(k,ith)         = cnt
          cnt                         = cnt                    + import_pSeatW(k,ith)
       enddo
    enddo
    !  -- [1-5] export_index for ???  -- !
    do k=1, ns
       do neib=1, NeibPEtot
          export_Index(k,neib)        = export_Index(k,neib-1) + export_Items(k,neib)
       enddo
    enddo
    !$omp parallel default(none) &
    !$omp shared(import_pSeatW,import_iBase,import_Addrs,import_AddrsW) &
    !$omp private(k,idx,gidx,ith)
    !$ ith = omp_get_thread_num() + 1
    !$omp do
    do k=1, ns
       do idx=1,import_pSeatW(k,ith)
          gidx                 = import_iBase(k,ith) + idx
          import_Addrs(gidx,k) = import_AddrsW(idx,k,ith)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! ------------------------------ !
    ! --- [2]  Particle Packing  --- !
    ! ------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,export_IndexW,sendBuff,export_AddrsW,export_iBase,pxv,x1slide,x1Mid) &
    !$omp private(neib,k,idx,gidx,ith)
    !$ ith = omp_get_thread_num() + 1
    !$omp do
    do neib=1, NeibPEtot
       do k=1, ns
          do idx = export_IndexW(k-1,neib,ith)+1, export_IndexW(k,neib,ith)
             gidx = export_iBase(k,neib,ith) + idx
             sendBuff(1,gidx,neib) = pxv( 1, export_AddrsW(idx,neib,ith), k )
             sendBuff(2,gidx,neib) = pxv( 2, export_AddrsW(idx,neib,ith), k ) + x1slide(neib)
             sendBuff(3,gidx,neib) = pxv( 3, export_AddrsW(idx,neib,ith), k )
             sendBuff(4,gidx,neib) = pxv( 4, export_AddrsW(idx,neib,ith), k )
             sendBuff(5,gidx,neib) = pxv( 5, export_AddrsW(idx,neib,ith), k )
             sendBuff(6,gidx,neib) = pxv( 6, export_AddrsW(idx,neib,ith), k )
             sendBuff(7,gidx,neib) = pxv( 7, export_AddrsW(idx,neib,ith), k ) + x1slide(neib)
             sendBuff(8,gidx,neib) = pxv( 8, export_AddrsW(idx,neib,ith), k )
             pxv( 8, export_AddrsW(idx,neib,ith), k ) = 0.d0
             pxv( 2, export_AddrsW(idx,neib,ith), k ) = x1Mid
             pxv( 7, export_AddrsW(idx,neib,ith), k ) = x1Mid
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    ! call MPI_Finalize( ierr )
    ! stop

    ! ------------------------------ !
    ! --- [3]  Exchange pInfo    --- !
    ! ------------------------------ !
    !  -- [3-2] Send    Info. --  !
    do neib=1, NeibPEtot
       call MPI_Isend( export_Items(:,neib), ns, MPI_INTEGER   , NeibPE(neib), 0, &
            &          MPI_COMM_WORLD      , request_Send(neib), ierr )
    enddo
    !  -- [3-3] Receive Info. --  !
    do neib=1, NeibPEtot
       call MPI_Irecv( import_Items(:,neib), ns, MPI_INTEGER   , NeibPE(neib), 0, &
            &          MPI_COMM_WORLD      , request_Recv(neib), ierr )
    enddo
    !  -- [3-4] WaitAll / Update --  !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )
    import_Index = 0
    import_iFrom = 0
    do neib=1, NeibPEtot
       do k=1, ns
          import_Index(k,neib) = import_Index(k-1,neib  ) + import_Items(k,neib)
       enddo
       do k=1, ns
          import_iFrom(k,neib) = import_iFrom(k  ,neib-1) + import_Items(k,neib)
       enddo
    enddo
    
    ! ------------------------- !
    ! --- [4] ISend / IRecv --- !
    ! ------------------------- !
    !  -- [4-1] ISend       --  !
    do neib=1, NeibPEtot
       len_s = export_Index(ns,neib) * 8
       call MPI_Isend( sendBuff(1,1,neib), len_s, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib)      , 0    , MPI_COMM_WORLD,       &
            &          request_Send(neib), ierr                         )
    enddo
    !  -- [4-2] IRecv       --  !
    do neib=1, NeibPEtot
       len_r = import_Index(ns,neib) * 8
       call MPI_Irecv( recvBuff(1,1,neib), len_r, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib)      , 0    , MPI_COMM_WORLD,       &
            &          request_Recv(neib), ierr                         )
    enddo

    ! ------------------------- !
    ! --- [5]  Update       --- !
    ! ------------------------- !
    !  -- [5-1] WaitAll Receive --  !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )
    !  -- [5-2] Update pxv --  !
    do k=1, ns
       npAdd = import_iFrom(k,NeibPEtot) - import_pSeat(k)
       if ( npAdd.ge. 1 ) then
          do idx=1, npAdd
             import_Addrs( import_pSeat(k)+idx,k) = np(k) + idx
          enddo
          np(k) = np(k) + npAdd
       endif
       if ( npAdd.le.-1 ) then
          do idx=import_iFrom(k,NeibPEtot)+1, import_pSeat(k)
             pxv(2,import_Addrs(idx,k),k) = x1Mid
             pxv(7,import_Addrs(idx,k),k) = x1Mid
             pxv(8,import_Addrs(idx,k),k) = 0.d0
          enddo
       endif
       npc(k) = npc(k) + npAdd
    end do
    npt = npc(1) + npc(2)
    do neib=1, NeibPEtot
       do k=1, ns
          do idx=1, Import_Items(k,neib)
             iaddr = import_Addrs(idx+import_iFrom(k  ,neib-1),k)
             ibuff =              idx+import_Index(k-1,neib  )
             pxv(1:8,iaddr,k) = recvBuff(1:8,ibuff,neib)
          enddo
       enddo
    enddo
    !  -- [5-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )
    
    return
  end subroutine ParticleExchange
  
end module pMPIComMod

