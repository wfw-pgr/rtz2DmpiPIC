module jMPIComMod
  implicit none
  ! -- Constants        -- !
  integer         , parameter   :: dim           = 2
  integer         , parameter   :: NeibPEtot_Max = 2
  integer         , parameter   :: nLineL        = 2
  integer         , parameter   :: nLineR        = 3
  integer         , parameter   :: nLine         = nLineL+nLineR
  integer                       :: LIloc, LJloc, nsloc
  integer         , parameter   :: id_ = 1,  jd_ = 2
  ! -- Global Variables -- !
  integer                       :: myRank, PEtot, CommLen
  integer                       :: NeibPE(NeibPEtot_Max)
  integer                       :: NeibPEtot
  ! -- Allocatable Var. -- !
  integer         , allocatable :: export_Items(:)  , import_Items(:)
  integer         , allocatable :: export_Index(:)  , import_Index(:)
  integer         , allocatable :: export_Addrs(:,:), import_Addrs(:,:)
  integer         , allocatable :: request_Send(:)  , request_Recv(:)
  integer         , allocatable ::  status_Send(:,:),  status_Recv(:,:)
contains


  subroutine allocExchangeTable_Current( LJs )
    implicit none
    include 'mpif.h'
    integer, intent(in) :: LJs

    ! --- [1]  Size CommLen    --- !
    CommLen = LJs * nLine
    ! --- [2]  Allocation      --- !
    allocate( export_Index(0:NeibPEtot_Max), import_Index(0:NeibPEtot_Max) )
    allocate( export_Items(0:NeibPEtot_Max), import_Items(0:NeibPEtot_Max) )
    allocate( export_Addrs(    CommLen,dim), import_Addrs(    CommLen,dim) )
    allocate( request_Send(NeibPEtot_Max),   request_Recv(NeibPEtot_Max)   )
    allocate(  status_Send(MPI_STATUS_SIZE,NeibPEtot_Max), &
         &     status_Recv(MPI_STATUS_SIZE,NeibPEtot_Max)  )
    export_Index = 0
    import_Index = 0
    export_Items = 0
    import_Items = 0
    export_Addrs = 0
    import_Addrs = 0
    
    return
  end subroutine allocExchangeTable_Current

  
  subroutine setupExchangeTable_Current( LIs, LJs, bctype )
    implicit none
    include 'mpif.h'
    integer      , intent(in) :: LIs, LJs
    character(10), intent(in) :: bctype
    integer                   :: i, j, cnt, idx, ierr
    integer                   :: ex_iidx(nLine), im_iidx(nLine)

    ! -------------------------------------------- !
    ! --- [1] Define Grid Number and Partition --- !
    ! -------------------------------------------- !
    !  -- [1-1]  PE infomation       --- !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )
    LIloc   = LIs
    LJloc   = LJs
    CommLen = LJloc  * nLine
    
    ! -------------------------------------------- !
    ! --- [2] Generate MPI Communication Table --- !
    ! -------------------------------------------- !
    !  -- [2-0] Exception for PEtot=1 -- !
    if ( PEtot.eq.1 ) then
       NeibPEtot = 0
       return
    endif
    ! --------------------------- !
    ! -- [2-1] Current Overlay -- !
    ! --------------------------- !
    if ( trim(bctype).eq.'cWall' ) then
       if ( myRank.eq.0       ) then
          ! -- Neibour Settings -- !
          NeibPEtot            = 1
          NeibPE(1)            = myRank+1
          export_Items(1)      = LJloc * nLineR
          import_Items(1)      = LJloc * nLineL
          do cnt=1, NeibPEtot
             export_Index(cnt) = export_Index(cnt-1) + export_Items(cnt)
             import_Index(cnt) = import_Index(cnt-1) + import_Items(cnt)
          enddo
          ! -- Export -- !
          cnt = 0
          do i=LIloc, LIloc+2
             do j=1, LJloc
                cnt                    = cnt + 1
                export_Addrs(cnt,id_)  = i
                export_Addrs(cnt,jd_)  = j
             enddo
          enddo
          ! -- Import -- !
          cnt = 0
          do i=LIloc-2, LIloc-1
             do j=1, LJloc
                cnt                    = cnt + 1
                import_Addrs(cnt,id_)  = i
                import_Addrs(cnt,jd_)  = j
             enddo
          enddo
       endif
       if ( myRank.eq.PEtot-1 ) then
          ! -- Neibour Settings -- !
          NeibPEtot            = 1
          NeibPE(1)            = myRank-1
          export_Items(1)      = LJloc * nLineL
          import_Items(1)      = LJloc * nLineR
          do cnt=1, NeibPEtot
             export_Index(cnt) = export_Index(cnt-1) + export_Items(cnt)
             import_Index(cnt) = import_Index(cnt-1) + import_Items(cnt)
          enddo
          ! -- Export -- !
          cnt = 0
          do i=0, 1
             do j=1, LJloc
                cnt                    = cnt + 1
                export_Addrs(cnt,id_)  = i
                export_Addrs(cnt,jd_)  = j
             enddo
          enddo
          ! -- Import -- !
          cnt = 0
          do i=2, 4
             do j=1, LJloc
                cnt                    = cnt + 1
                import_Addrs(cnt,id_)  = i
                import_Addrs(cnt,jd_)  = j
             enddo
          enddo
       endif
       if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
          ! -- Neibour Settings -- !
          NeibPEtot            = 2
          NeibPE(1)            = myRank-1
          NeibPE(2)            = myRank+1
          export_Items(1)      = LJloc * nLineL
          export_Items(2)      = LJloc * nLineR
          import_Items(1)      = LJloc * nLineR
          import_Items(2)      = LJloc * nLineL
          do cnt=1, NeibPEtot
             export_Index(cnt) = export_Index(cnt-1) + export_Items(cnt)
             import_Index(cnt) = import_Index(cnt-1) + import_Items(cnt)
          enddo
          ! -- Export -- !
          cnt = 0
          do i=0, 1
             do j=1, LJloc
                cnt                    = cnt + 1
                export_Addrs(cnt,id_)  = i
                export_Addrs(cnt,jd_)  = j
             enddo
          enddo
          do i=LIloc, LIloc+2
             do j=1, LJloc
                cnt                    = cnt + 1
                export_Addrs(cnt,id_)  = i
                export_Addrs(cnt,jd_)  = j
             enddo
          enddo
          ! -- Import -- !
          cnt = 0
          do i=2, 4
             do j=1, LJloc
                cnt                    = cnt + 1
                import_Addrs(cnt,id_)  = i
                import_Addrs(cnt,jd_)  = j
             enddo
          enddo
          do i=LIloc-2, LIloc-1
             do j=1, LJloc
                cnt                    = cnt + 1
                import_Addrs(cnt,id_)  = i
                import_Addrs(cnt,jd_)  = j
             enddo
          enddo
       endif
    endif
    if ( trim(bctype).eq.'periodic' ) then
       ! --  Neibour / Item / Index  Settings  -- !
       NeibPEtot                    = 2
       NeibPE(1)                    = myRank-1
       NeibPE(2)                    = myRank+1
       export_Items(1)              = LJloc * nLineL
       export_Items(2)              = LJloc * nLineR
       import_Items(1)              = LJloc * nLineR
       import_Items(2)              = LJloc * nLineL
       do cnt=1, NeibPEtot
          export_Index(cnt) = export_Index(cnt-1) + export_Items(cnt)
          import_Index(cnt) = import_Index(cnt-1) + import_Items(cnt)
       enddo
       ! --  export Address index i  -- !
       ex_iidx(1) = 0
       ex_iidx(2) = 1
       ex_iidx(3) = LIloc
       ex_iidx(4) = LIloc+1
       ex_iidx(5) = LIloc+2
       ! --  import Address index i  -- !
       im_iidx(1) = 2
       im_iidx(2) = 3
       im_iidx(3) = 4
       im_iidx(4) = LIloc-2
       im_iidx(5) = LIloc-1
       ! --  substituition           -- !
       do i=1, nLine
          do j=1, LJloc
             idx                    = j + LJloc*(i-1)
             export_Addrs(idx,id_)  = ex_iidx(i)
             export_Addrs(idx,jd_)  = j
             import_Addrs(idx,id_)  = im_iidx(i)
             import_Addrs(idx,jd_)  = j
          enddo
       enddo
       if ( myRank.eq.0       ) NeibPE(1) = PEtot-1
       if ( myRank.eq.PEtot-1 ) NeibPE(2) = 0
    endif
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '[DefineExchangeTable_Current  ]          * Current   :: ( Initialized )'
    endif
    
    return
  end subroutine setupExchangeTable_Current


  subroutine CurrentOverlay( xvc )
    implicit none
    include 'mpif.h'
    double precision, intent(inout) :: xvc(3,0:LIloc+2,0:LJloc+2)
    integer                         :: ip, jp, idx, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(3,CommLen)
    double precision                :: recvBuff(3,CommLen)
    
    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,export_Index,sendBuff,export_Addrs,xvc) &
    !$omp private(neib,idx)
    !$omp do
    do neib=1, NeibPEtot
       do idx=export_Index(neib-1)+1, export_Index(neib)
          sendBuff(1,idx) = xvc( 1, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(2,idx) = xvc( 2, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(3,idx) = xvc( 3, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
       enddo
    end do
    !$omp end do
    !$omp end parallel
    
    ! ------------------------- !
    ! --- [2] ISend / IRecv --- !
    ! ------------------------- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot
       is    = export_Index(neib-1) + 1
       len_s = export_Items(neib  ) * 3
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib)      , 0    , MPI_COMM_WORLD,       &
            &          request_Send(neib), ierr                         )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot
       ir    = import_Index(neib-1) + 1
       len_r = import_Items(neib  ) * 3
       call MPI_Irecv( recvBuff(1,ir)    , len_r, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib)      , 0    , MPI_COMM_WORLD,       &
            &          request_Recv(neib), ierr                         )
    enddo

    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll Receive --  !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )
    !  -- [3-2] Update  xVector --  !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,import_Index,import_Addrs,xvc,recvBuff) &
    !$omp private(neib,idx,ip,jp)
    !$omp do
    do neib=1, NeibPEtot
       do idx= import_Index(neib-1)+1, import_Index(neib)
          ip = import_Addrs(idx,id_)
          jp = import_Addrs(idx,jd_)
          xvc(1,ip,jp) = xvc(1,ip,jp) + recvBuff(1,idx)
          xvc(2,ip,jp) = xvc(2,ip,jp) + recvBuff(2,idx)
          xvc(3,ip,jp) = xvc(3,ip,jp) + recvBuff(3,idx)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )

    return
  end subroutine CurrentOverlay


  subroutine DepositOverlay( xvc, nsp )
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: nsp
    double precision, intent(inout) :: xvc(0:LIloc+2,0:LJloc+2,nsp)
    integer                         :: ip, jp, icmp, idx, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(nsp,CommLen)
    double precision                :: recvBuff(nsp,CommLen)
    
    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,export_Index,sendBuff,export_Addrs,xvc,nsp) &
    !$omp private(neib,idx,icmp)
    !$omp do
    do neib=1, NeibPEtot
       do idx=export_Index(neib-1)+1, export_Index(neib)
          do icmp=1, nsp
             sendBuff(icmp,idx) = xvc( export_Addrs(idx,id_), export_Addrs(idx,jd_), icmp )
          enddo
       end do
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ------------------------- !
    ! --- [2] ISend / IRecv --- !
    ! ------------------------- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot
       is    = export_Index(neib-1) + 1
       len_s = export_Items(neib  ) * nsp
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib)      , 0    , MPI_COMM_WORLD,       &
            &          request_Send(neib), ierr                         )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot
       ir    = import_Index(neib-1) + 1
       len_r = import_Items(neib  ) * nsp
       call MPI_Irecv( recvBuff(1,ir)    , len_r, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib)      , 0    , MPI_COMM_WORLD,       &
            &          request_Recv(neib), ierr                         )
    enddo
    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll Receive --  !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )
    !  -- [3-2] Update  xVector --  !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,import_Index,import_Addrs,xvc,recvBuff,nsp) &
    !$omp private(neib,idx,ip,jp,icmp)
    !$omp do
    do neib=1, NeibPEtot
       do idx= import_Index(neib-1)+1, import_Index(neib)
          ip = import_Addrs(idx,id_)
          jp = import_Addrs(idx,jd_)
          do icmp=1, nsp
             xvc(ip,jp,icmp) = xvc(ip,jp,icmp) + recvBuff(icmp,idx)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )

    return
  end subroutine DepositOverlay


end module jMPIComMod
