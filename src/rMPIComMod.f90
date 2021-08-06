module rMPIComMod
  implicit none
  ! -- Constants        -- !
  integer         , parameter   :: dim            = 2
  integer         , parameter   :: NeibPEtot_Max  = 2
  integer         , parameter   :: nLineL         = 2
  integer         , parameter   :: nLineR         = 1
  integer         , parameter   :: nLine          = nLineL+nLineR
  integer                       :: LIloc, LJloc
  integer         , parameter   :: id_  = 1, jd_  = 2
  ! -- Global Variables -- !
  integer                       :: myRank, PEtot, CommLen
  integer                       :: NeibPE(NeibPEtot_Max)
  integer                       :: NeibPEtot
  ! -- Allocatable Var. -- !
  integer         , allocatable :: export_Index(:)  , import_Index(:)
  integer         , allocatable :: export_Items(:)  , import_Items(:)
  integer         , allocatable :: export_Addrs(:,:), import_Addrs(:,:)
  integer         , allocatable :: request_Send(:)  , request_Recv(:)
  integer         , allocatable ::  status_Send(:,:),  status_Recv(:,:)
contains


  subroutine allocExchangeTable_Relocated( LJs )
    implicit none
    include 'mpif.h'
    integer, intent(in) :: LJs

    ! --- [1]  Size CommLen    --- !
    CommLen  = (LJs+1) * 3    
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
  end subroutine allocExchangeTable_Relocated

  
  subroutine setupExchangeTable_Relocated( LIs, LJs, bctype )
    implicit none
    include 'mpif.h'
    integer      , intent(in)   :: LIs, LJs
    character(10), intent(in)   :: bctype 
    integer                     :: i, j, cnt, cnt1, cnt2, cnt3, idx, ierr
    integer                     :: ex_iidx(nLine), im_iidx(nLine)
    
    ! -------------------------------------------- !
    ! --- [1] Define Grid Number and Partition --- !
    ! -------------------------------------------- !
    !  -- [1-1]  PE infomation       --- !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )
    LIloc    =  LIs
    LJloc    =  LJs
    CommLen  = (LJs+1) * 3    
    
    ! -------------------------------------------- !
    ! --- [2] Generate MPI Communication Table --- !
    ! -------------------------------------------- !
    !  -- [2-0] Exception for PEtot=1 -- !
    if ( PEtot.eq.1        ) then
       NeibPEtot = 0
       return
    endif
    ! ------------------------- !
    ! -- [2-1] EBr_Comm      -- !
    ! ------------------------- !
    if ( trim(bctype).eq.'cWall' ) then
       if ( myRank.eq.0       ) then
          NeibPEtot                 = 1
          NeibPE(1)                 = 1
          export_Items(1)           = (LJloc+1)
          import_Items(1)           = (LJloc+1)*2
          export_Index(1)           = (LJloc+1)
          import_Index(1)           = (LJloc+1)*2
          do cnt=1, (LJloc+1)
             cnt1                   = cnt
             cnt2                   = cnt + (LJloc+1)
             export_Addrs(cnt1,id_) = LIloc-1
             export_Addrs(cnt1,jd_) = cnt
             import_Addrs(cnt1,id_) = LIloc
             import_Addrs(cnt1,jd_) = cnt
             import_Addrs(cnt2,id_) = LIloc+1
             import_Addrs(cnt2,jd_) = cnt
          enddo
       endif
       if ( myRank.eq.PEtot-1 ) then
          NeibPEtot                 = 1
          NeibPE(1)                 = PEtot-2
          export_Items(1)           = (LJloc+1)*2
          import_Items(1)           = (LJloc+1)
          export_Index(1)           = (LJloc+1)*2
          import_Index(1)           = (LJloc+1)
          do cnt=1, (LJloc+1)
             cnt1                   = cnt
             cnt2                   = cnt + (LJloc+1)
             export_Addrs(cnt1,id_) = 2
             export_Addrs(cnt1,jd_) = cnt
             import_Addrs(cnt1,id_) = 1
             import_Addrs(cnt1,jd_) = cnt
             export_Addrs(cnt2,id_) = 3
             export_Addrs(cnt2,jd_) = cnt
          enddo
       endif
       if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
          NeibPEtot                 = 2
          NeibPE(1)                 = myRank-1
          NeibPE(2)                 = myRank+1
          export_Items(1)           = (LJloc+1)*2
          export_Items(2)           = (LJloc+1)
          import_Items(1)           = (LJloc+1)
          import_Items(2)           = (LJloc+1)*2
          do cnt=1, NeibPEtot 
             export_Index(cnt)      = export_Index(cnt-1) + export_Items(cnt)
             import_Index(cnt)      = import_Index(cnt-1) + import_Items(cnt)
          enddo
          do cnt=1, (LJloc+1)
             cnt1                   = cnt 
             cnt2                   = cnt + (LJloc+1)
             cnt3                   = cnt + (LJloc+1)*2 
             export_Addrs(cnt1,id_) = 2
             export_Addrs(cnt1,jd_) = cnt
             export_Addrs(cnt2,id_) = 3
             export_Addrs(cnt2,jd_) = cnt
             export_Addrs(cnt3,id_) = LIloc-1
             export_Addrs(cnt3,jd_) = cnt
             import_Addrs(cnt1,id_) = 1
             import_Addrs(cnt1,jd_) = cnt
             import_Addrs(cnt2,id_) = LIloc
             import_Addrs(cnt2,jd_) = cnt
             import_Addrs(cnt3,id_) = LIloc+1
             import_Addrs(cnt3,jd_) = cnt
          enddo
       endif
    endif
    if ( trim(bctype).eq.'periodic' ) then
       NeibPEtot                    = 2
       NeibPE(1)                    = myRank-1
       NeibPE(2)                    = myRank+1
       export_Items(1)              = (LJloc+1)*nLineL
       export_Items(2)              = (LJloc+1)*nLineR
       import_Items(1)              = (LJloc+1)*nLineR
       import_Items(2)              = (LJloc+1)*nLineL
       do cnt=1, NeibPEtot 
          export_Index(cnt)         = export_Index(cnt-1) + export_Items(cnt)
          import_Index(cnt)         = import_Index(cnt-1) + import_Items(cnt)
       enddo
       ! --  export Address index i  -- !
       ex_iidx(1) = 2
       ex_iidx(2) = 3
       ex_iidx(3) = LIloc-1
       ! --  import Address index i  -- !
       im_iidx(1) = 1
       im_iidx(2) = LIloc
       im_iidx(3) = LIloc+1
       ! --  substituition           -- !
       do i=1, nLine
          do j=1, LJloc+1
             idx                    = j + (LJloc+1)*(i-1)
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
       write(6,'(2x,a)') '[DefineExchangeTable_Relocated]          * Relocated :: ( Initialized )'
       write(6,*)
    endif
    
    return
  end subroutine setupExchangeTable_Relocated


  subroutine RelocatedExchange( xvc )
    ! == Communication Routine for xvc(6,LIs,LJs) == !
    implicit none
    include 'mpif.h'
    double precision, intent(inout) :: xvc(6,LIloc+1,LJloc+1)
    integer                         :: idx, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(6,CommLen)
    double precision                :: recvBuff(6,CommLen)
    
    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    !  -- [1-1] Packing -- !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,export_Index,sendBuff,xvc,export_Addrs) &
    !$omp private(neib,idx)
    !$omp do
    do neib=1, NeibPEtot 
       do idx=export_Index(neib-1)+1, export_Index(neib)
          sendBuff(1,idx) = xvc( 1, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(2,idx) = xvc( 2, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(3,idx) = xvc( 3, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(4,idx) = xvc( 4, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(5,idx) = xvc( 5, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(6,idx) = xvc( 6, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! ------------------------- !
    ! --- [2] ISend / IRecv --- !
    ! ------------------------- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot 
       is    = export_Index(neib-1) + 1
       len_s = export_Items(neib  ) * 6
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot 
       ir    = import_Index(neib-1) + 1
       len_r = import_Items(neib  ) * 6
       call MPI_Irecv( recvBuff(1,ir)    , len_r, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Recv(neib), ierr )
    enddo

    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll Receive --  !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )

    !  -- [3-2] Update  xVector --  !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,import_Index,xvc,import_Addrs,recvBuff) &
    !$omp private(neib,idx)
    !$omp do
    do neib=1, NeibPEtot 
       do idx=import_Index(neib-1)+1, import_Index(neib)
          xvc( 1, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(1,idx)
          xvc( 2, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(2,idx)
          xvc( 3, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(3,idx)
          xvc( 4, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(4,idx)
          xvc( 5, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(5,idx)
          xvc( 6, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(6,idx)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )

    return
  end subroutine RelocatedExchange

end module rMPIComMod
