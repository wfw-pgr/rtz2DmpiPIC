module fMPIComMod
  implicit none
  ! -- Constants        -- !
  integer         , parameter   :: dim            = 2
  integer         , parameter   :: NeibPEtot_Max  = 2
  integer         , parameter   :: nLineL         = 1
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


  subroutine allocExchangeTable_Field( LJs )
    implicit none
    include 'mpif.h'
    integer, intent(in) :: LJs
    
    ! --- [1]  Size CommLen    --- !
    CommLen = LJs * NeibPEtot_Max
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
  end subroutine allocExchangeTable_Field

  
  subroutine setupExchangeTable_Field( LIs, LJs, bctype )
    implicit none
    include 'mpif.h'
    integer      , intent(in) :: LIs, LJs
    character(10), intent(in) :: bctype
    integer                   :: i, j, cnt, cnt1, cnt2, idx, ierr
    integer                   :: ex_iidx(nLine), im_iidx(nLine)
    ! -------------------------------------------- !
    ! --- [1] Define Grid Number and Partition --- !
    ! -------------------------------------------- !
    !  -- [1-1]  PE infomation       --- !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )
    LIloc   = LIs
    LJloc   = LJs
    CommLen = LJs * NeibPEtot_Max
    
    ! -------------------------------------------- !
    ! --- [2] Generate MPI Communication Table --- !
    ! -------------------------------------------- !
    !  -- [2-0] Exception for PEtot=1 -- !
    if ( PEtot.eq.1        ) then
       NeibPEtot = 0
       return
    endif
    !  -------------------------  !
    !  -- [2-1] Boundary_Comm --  !
    !  -------------------------  !
    if ( trim(bctype).eq.'cWall' ) then
       if ( myRank.eq.0       ) then
          NeibPEtot                 = 1
          NeibPE(1)                 = 1
          export_Items(1)           = LJloc
          import_Items(1)           = LJloc
          export_Index(1)           = LJloc
          import_Index(1)           = LJloc
          do cnt=1, LJloc
             export_Addrs(cnt,id_)  = LIloc-1
             export_Addrs(cnt,jd_)  = cnt
             import_Addrs(cnt,id_)  = LIloc
             import_Addrs(cnt,jd_)  = cnt
          enddo
       endif
       if ( myRank.eq.PEtot-1 ) then
          NeibPEtot                 = 1
          NeibPE(1)                 = PEtot-2
          export_Items(1)           = LJloc
          import_Items(1)           = LJloc
          export_Index(1)           = LJloc
          import_Index(1)           = LJloc
          do cnt=1, LJloc
             export_Addrs(cnt,id_)  = 2
             export_Addrs(cnt,jd_)  = cnt
             import_Addrs(cnt,id_)  = 1
             import_Addrs(cnt,jd_)  = cnt
          enddo
       endif
       if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
          NeibPEtot                 = 2
          NeibPE(1)                 = myRank-1
          NeibPE(2)                 = myRank+1
          do cnt=1, NeibPEtot
             export_Items(cnt)      = LJloc
             import_Items(cnt)      = LJloc
             export_Index(cnt)      = export_Index(cnt-1) + export_Items(cnt)
             import_Index(cnt)      = import_Index(cnt-1) + import_Items(cnt)
          enddo
          do cnt=1, LJloc
             cnt1                   = cnt
             cnt2                   = cnt + LJloc
             export_Addrs(cnt1,id_) = 2
             export_Addrs(cnt1,jd_) = cnt
             import_Addrs(cnt1,id_) = 1
             import_Addrs(cnt1,jd_) = cnt
             export_Addrs(cnt2,id_) = LIloc-1
             export_Addrs(cnt2,jd_) = cnt
             import_Addrs(cnt2,id_) = LIloc
             import_Addrs(cnt2,jd_) = cnt
          enddo
       endif
    endif
    ! -------------------------- !
    ! -- [2-2]  Periodic_Comm -- !
    ! -------------------------- !
    if ( trim(bctype).eq.'periodic' ) then
       NeibPEtot            = 2
       NeibPE(1)            = myRank-1
       NeibPE(2)            = myRank+1
       export_Items(1:2)    = LJloc
       import_Items(1:2)    = LJloc
       do cnt=1, NeibPEtot
          export_Index(cnt) = export_Index(cnt-1) + export_Items(cnt)
          import_Index(cnt) = import_Index(cnt-1) + import_Items(cnt)
       enddo
       ex_iidx(1) = 2
       ex_iidx(2) = LIloc-1
       im_iidx(1) = 1
       im_iidx(2) = LIloc
       do i=1, nLine
          do j=1, LJloc
             idx  = j + LJloc*(i-1)
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
       write(6,*)
       write(6,'(2x,a)') '---------------         [     MPI-TABLE     ]          ---------------'
       write(6,'(2x,a)') '[DefineExchangeTable_Field    ]          * Field     :: ( Initialized )'
    endif
    
    return
  end subroutine setupExchangeTable_Field


  subroutine FieldCommunicate( xvc, EorB )
    ! == Communication Routine for xvc(6,LIs,LJs) == !
    implicit none
    include 'mpif.h'
    character(1)    , intent(in)    :: EorB
    double precision, intent(inout) :: xvc(6,LIloc,LJloc)
    integer                         :: idx, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(3,CommLen)
    double precision                :: recvBuff(3,CommLen)
    integer                         :: d1_, d2_, d3_

    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    !  -- [1-1] Select Mode -- !
    if ( EorB.eq.'E' ) then
       d1_=1 ; d2_=2 ; d3_=3
    endif
    if ( EorB.eq.'B' ) then
       d1_=4 ; d2_=5 ; d3_=6
    endif
    !  -- [1-2] Packing -- !
    !$omp parallel default(none) &
    !$omp shared (NeibPEtot,export_Index,export_Addrs,sendBuff,xvc,d1_,d2_,d3_) &
    !$omp private(neib,idx)
    !$omp do
    do neib=1, NeibPEtot
       do idx=export_Index(neib-1)+1, export_Index(neib)
          sendBuff(1,idx) = xvc( d1_, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(2,idx) = xvc( d2_, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          sendBuff(3,idx) = xvc( d3_, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
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
       len_s = export_Items(neib  ) * 3
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib), 0   , MPI_COMM_WORLD, &
            &          request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot
       ir    = import_Index(neib-1) + 1
       len_r = import_Items(neib  ) * 3
       call MPI_Irecv( recvBuff(1,ir)    , len_r, MPI_DOUBLE_PRECISION, &
            &          NeibPE(neib), 0   , MPI_COMM_WORLD, &
            &          request_Recv(neib), ierr )
    enddo

    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll Receive / Send --  !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )
    !  -- [3-2] Update  xVector        --  !
    do neib=1, NeibPEtot
       !$omp parallel default(none) &
       !$omp shared (neib,NeibPEtot,import_Index,import_Addrs,recvBuff,xvc,d1_,d2_,d3_) &
       !$omp private(idx)
       !$omp do
       do idx=import_Index(neib-1)+1, import_Index(neib)
          xvc( d1_, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(1,idx)
          xvc( d2_, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(2,idx)
          xvc( d3_, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(3,idx)
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    

    return
  end subroutine FieldCommunicate


  subroutine FluidCommunicate( xvc, nsp )
    ! == Communication Routine for xvc(LIs,LJs,ns) == !
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: nsp
    double precision, intent(inout) :: xvc(0:LIloc+2,0:LJloc+2,nsp)
    integer                         :: idx, isp, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(nsp,CommLen)
    double precision                :: recvBuff(nsp,CommLen)

    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,export_Index,export_Addrs,nsp,sendBuff,xvc) &
    !$omp private(isp,idx,neib)
    !$omp do
    do neib=1, NeibPEtot
       do idx=export_Index(neib-1)+1, export_Index(neib)
          do isp=1, nsp
             sendBuff(isp,idx) = xvc( export_Addrs(idx,id_), export_Addrs(idx,jd_), isp )
          enddo
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
       len_s = export_Items(neib  ) * nsp
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot
       ir    = import_Index(neib-1) + 1
       len_r = import_Items(neib  ) * nsp
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
    !$omp shared(NeibPEtot,import_Index,nsp,xvc,import_Addrs,recvBuff) &
    !$omp private(neib,idx,isp)
    !$omp do
    do neib=1, NeibPEtot
       do idx=import_Index(neib-1)+1, import_Index(neib)
          do isp=1, nsp
             xvc( import_Addrs(idx,id_), import_Addrs(idx,jd_), isp ) = recvBuff(isp,idx)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )

    return
  end subroutine FluidCommunicate


  subroutine CurrentCommunicate( xvc, nCmp )
    ! == Communication Routine for xvc(LIs,LJs,ns) == !
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: nCmp
    double precision, intent(inout) :: xvc(nCmp,0:LIloc+2,0:LJloc+2)
    integer                         :: idx, icmp, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(nCmp,CommLen)
    double precision                :: recvBuff(nCmp,CommLen)

    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    !$omp parallel default(none) &
    !$omp shared(NeibPEtot,export_Index,export_Addrs,nCmp,sendBuff,xvc) &
    !$omp private(icmp,idx,neib)
    !$omp do
    do neib=1, NeibPEtot
       do idx=export_Index(neib-1)+1, export_Index(neib)
          do icmp=1, nCmp
             sendBuff(icmp,idx) = xvc( icmp, export_Addrs(idx,id_), export_Addrs(idx,jd_) )
          enddo
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
       len_s = export_Items(neib  ) * nCmp
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot
       ir    = import_Index(neib-1) + 1
       len_r = import_Items(neib  ) * nCmp
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
    !$omp shared(NeibPEtot,import_Index,nCmp,xvc,import_Addrs,recvBuff) &
    !$omp private(neib,idx,icmp)
    !$omp do
    do neib=1, NeibPEtot
       do idx=import_Index(neib-1)+1, import_Index(neib)
          do icmp=1, nCmp
             xvc( icmp, import_Addrs(idx,id_), import_Addrs(idx,jd_) ) = recvBuff(icmp,idx)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )

    return
  end subroutine CurrentCommunicate

end module fMPIComMod
