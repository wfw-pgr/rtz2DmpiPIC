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
  integer         , allocatable :: export_Index(:,:), import_Index(:,:)
  integer         , allocatable :: export_Items(:,:), import_Items(:,:)
  integer         , allocatable :: export_Addrs(:,:), import_Addrs(:,:)
  integer         , allocatable :: import_pSeat(:)  , import_iFrom(:,:)
  integer         , allocatable :: request_Send(:)  , request_Recv(:)
  integer         , allocatable ::  status_Send(:,:),  status_Recv(:,:)
  double precision, allocatable :: sendBuff(:,:,:)  , recvBuff(:,:,:)
  
contains

  
  subroutine DefineExchangeTable_Particle
    use constants, only        : ns, nptMax, LIs
    use constants, only        : Boundary1__pt
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
    allocate( export_Index(0:ns,NeibPEtot_Max), export_Items(ns,NeibPEtot_Max) )
    allocate( import_Index(0:ns,NeibPEtot_Max), import_Items(ns,NeibPEtot_Max) )
    allocate( export_Addrs(pBuffSize,NeibPEtot_Max) )
    allocate( import_Addrs(pBuffSize,ns)            )
    allocate( import_iFrom(ns,0:NeibPEtot_Max), import_pSeat(ns) )
    import_iFrom = 0
    import_pSeat = 0
    import_Index = 0
    export_Index = 0
    import_Addrs = 0
    export_Addrs = 0
     status_Send = 0
     status_Recv = 0
    request_Send = 0
    request_Recv = 0
    sendBuff     = 0.d0
    recvBuff     = 0.d0

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
    use constants, only : ns , npt, np, npc, zp_, wp_
    use variables, only : pxv
    implicit none
    include 'mpif.h'
    integer            :: m, k , neib , len_s, len_r, ierr
    integer            :: idx  , idx1 , idx2 , cnt1 , cnt2, cnt3
    integer            :: iaddr, ibuff, npAdd

    if ( .not.( Flag__initialized ) ) call setExchangeGeometry
    ! ------------------------------ !
    ! --- [1]  Inside / OutSide  --- !
    ! ------------------------------ !
    !  -- [1-1] Initialize variables      -- !
    idx1       = 0
    idx2       = 0
    export_Index(0:ns,1:NeibPEtot_Max)  = 0
    export_Items(1:ns,1:NeibPEtot_Max)  = 0
    !  -- [1-2] Exceeding Particle Search -- !
    do k=1, ns
       cnt1 = 0 ; cnt2 = 0 ; cnt3 = 0 ;
       do m=1, np(k)
          if ( pxv(wp_,m,k).gt.0.d0   ) then
             ! -- Particle Exist -- !
             !  -  Left   - !
             if ( pxv(zp_,m,k).le.x1End1 ) then
                cnt1 = cnt1 + 1
                cnt3 = cnt3 + 1
                idx1 = idx1 + 1
                export_Addrs(idx1,dst(1)) = m
                import_Addrs(cnt3,k     ) = m
             endif
             !  -  Right  - !
             if ( pxv(zp_,m,k).gt.x1End2 ) then
                cnt2 = cnt2 + 1
                cnt3 = cnt3 + 1
                idx2 = idx2 + 1
                export_Addrs(idx2,dst(2)) = m
                import_Addrs(cnt3,k     ) = m
             endif
          else
             ! --  No Particle   -- !
             cnt3 = cnt3 + 1
             import_Addrs(cnt3,k) = m
          endif
       enddo
       export_Items(k,dst(1)) = cnt1
       export_Items(k,dst(2)) = cnt2
       export_Index(k,dst(1)) = idx1
       export_Index(k,dst(2)) = idx2
       import_pSeat(k)        = cnt3
    enddo

    ! ------------------------------ !
    ! --- [2]  Particle Packing  --- !
    ! ------------------------------ !
    do neib=1, neibPEtot
       do k=1, ns
          do idx= export_Index(k-1,neib)+1, export_Index(k,neib)
             sendBuff(1,idx,neib) = pxv( 1, export_Addrs(idx,neib), k )
             sendBuff(2,idx,neib) = pxv( 2, export_Addrs(idx,neib), k ) + x1slide(neib)
             sendBuff(3,idx,neib) = pxv( 3, export_Addrs(idx,neib), k )
             sendBuff(4,idx,neib) = pxv( 4, export_Addrs(idx,neib), k )
             sendBuff(5,idx,neib) = pxv( 5, export_Addrs(idx,neib), k )
             sendBuff(6,idx,neib) = pxv( 6, export_Addrs(idx,neib), k )
             sendBuff(7,idx,neib) = pxv( 7, export_Addrs(idx,neib), k ) + x1slide(neib)
             sendBuff(8,idx,neib) = pxv( 8, export_Addrs(idx,neib), k )
             pxv( 8, export_Addrs(idx,neib), k ) = 0.d0
             pxv( 2, export_Addrs(idx,neib), k ) = x1Mid
             pxv( 7, export_Addrs(idx,neib), k ) = x1Mid
          enddo
       enddo
    enddo

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

