module fMPIComMod
  implicit none
  ! -- Constants        -- !
  integer         , parameter   :: dim            = 2
  integer         , parameter   :: NeibPEtot_Max  = 2
  integer         , parameter   :: nType          = 2
  integer                       :: LIloc, LJloc
  integer         , parameter   :: id_  = 1, jd_  = 2
  integer         , parameter   :: bcs_ = 1, prd_ = 2
  ! -- Global Variables -- !
  integer                       :: myRank, PEtot, CommLen
  integer                       :: NeibPE(NeibPEtot_Max,nType)
  integer                       :: NeibPEtot(nType)
  ! -- Allocatable Var. -- !
  integer         , allocatable :: export_Index(:,:)  , import_Index(:,:)
  integer         , allocatable :: export_Items(:,:)  , import_Items(:,:)
  integer         , allocatable :: export_Addrs(:,:,:), import_Addrs(:,:,:)
  integer         , allocatable :: request_Send(:)    , request_Recv(:)
  integer         , allocatable ::  status_Send(:,:)  ,  status_Recv(:,:)
contains
  
  
  subroutine DefineExchangeTable_Field( LIs, LJs )
    implicit none
    include 'mpif.h'
    integer, intent(in) :: LIs, LJs
    integer             :: cnt, cnt1, cnt2, ierr
    ! -------------------------------------------- !
    ! --- [1] Define Grid Number and Partition --- !
    ! -------------------------------------------- !
    !  -- [1-1]  PE infomation       --- !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )
    LIloc   = LIs
    LJloc   = LJs
    CommLen = LJloc * NeibPEtot_Max

    ! -------------------------------------------- !
    ! --- [2] Allocate Variables for MPI Comm. --- !
    ! -------------------------------------------- !
    allocate( export_Index(0:NeibPEtot_Max,nType), import_Index(0:NeibPEtot_Max,nType) )
    allocate( export_Items(0:NeibPEtot_Max,nType), import_Items(0:NeibPEtot_Max,nType) )
    allocate( export_Addrs(    CommLen,dim,nType), import_Addrs(    CommLen,dim,nType) )
    allocate( request_Send(NeibPEtot_Max),   request_Recv(NeibPEtot_Max)   )
    allocate(  status_Send(MPI_STATUS_SIZE,NeibPEtot_Max), &
         &     status_Recv(MPI_STATUS_SIZE,NeibPEtot_Max)  )
    
    ! -------------------------------------------- !
    ! --- [3] Generate MPI Communication Table --- !
    ! -------------------------------------------- !
    !  -- [3-0] Exception for PEtot=1 -- !
    if ( PEtot.eq.1        ) then
       NeibPEtot(1:nType) = 0
       return
    endif
    !  -------------------------  !
    !  -- [3-1] Boundary_Comm --  !
    !  -------------------------  !
    export_Index = 0
    import_Index = 0
    if ( myRank.eq.0       ) then
       NeibPEtot(bcs_)                = 1
       NeibPE(1,bcs_)                 = 1
       export_Items(1,bcs_)           = LJloc
       import_Items(1,bcs_)           = LJloc
       export_Index(1,bcs_)           = LJloc
       import_Index(1,bcs_)           = LJloc
       do cnt=1, LJloc
          export_Addrs(cnt,id_,bcs_)  = LIloc-1
          export_Addrs(cnt,jd_,bcs_)  = cnt
          import_Addrs(cnt,id_,bcs_)  = LIloc
          import_Addrs(cnt,jd_,bcs_)  = cnt
       enddo
    endif
    if ( myRank.eq.PEtot-1 ) then
       NeibPEtot(bcs_)                = 1
       NeibPE(1,bcs_)                 = PEtot-2
       export_Items(1,bcs_)           = LJloc
       import_Items(1,bcs_)           = LJloc
       export_Index(1,bcs_)           = LJloc
       import_Index(1,bcs_)           = LJloc
       do cnt=1, LJloc
          export_Addrs(cnt,id_,bcs_)  = 2
          export_Addrs(cnt,jd_,bcs_)  = cnt
          import_Addrs(cnt,id_,bcs_)  = 1
          import_Addrs(cnt,jd_,bcs_)  = cnt
       enddo
    endif
    if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
       NeibPEtot(bcs_)                = 2
       NeibPE(1,bcs_)                 = myRank-1
       NeibPE(2,bcs_)                 = myRank+1
       export_Items(1,bcs_)           = LJloc
       export_Items(2,bcs_)           = LJloc
       import_Items(1,bcs_)           = LJloc
       import_Items(2,bcs_)           = LJloc
       do cnt=1, NeibPEtot(bcs_)
          export_Index(cnt,bcs_) = export_Index(cnt-1,bcs_) + export_Items(cnt,bcs_)
          import_Index(cnt,bcs_) = import_Index(cnt-1,bcs_) + import_Items(cnt,bcs_)
       enddo
       do cnt=1, LJloc
          cnt1                        = cnt
          cnt2                        = cnt + LJloc
          export_Addrs(cnt1,id_,bcs_) = 2
          export_Addrs(cnt1,jd_,bcs_) = cnt
          import_Addrs(cnt1,id_,bcs_) = 1
          import_Addrs(cnt1,jd_,bcs_) = cnt
          export_Addrs(cnt2,id_,bcs_) = LIloc-1
          export_Addrs(cnt2,jd_,bcs_) = cnt
          import_Addrs(cnt2,id_,bcs_) = LIloc
          import_Addrs(cnt2,jd_,bcs_) = cnt
       enddo
    endif
    ! -------------------------- !
    ! -- [3-2]  Periodic_Comm -- !
    ! -------------------------- !
    NeibPEtot(prd_)                   = 2
    NeibPE(1,prd_)                    = myRank-1
    NeibPE(2,prd_)                    = myRank+1
    export_Items(1,prd_)              = LJloc
    export_Items(2,prd_)              = LJloc
    import_Items(1,prd_)              = LJloc
    import_Items(2,prd_)              = LJloc
    do cnt=1, NeibPEtot(prd_)
       export_Index(cnt,prd_) = export_Index(cnt-1,prd_) + export_Items(cnt,prd_)
       import_Index(cnt,prd_) = import_Index(cnt-1,prd_) + import_Items(cnt,prd_)
    enddo
    do cnt=1, LJloc
       cnt1                           = cnt
       cnt2                           = cnt + LJloc
       export_Addrs(cnt1,id_,prd_)    = 2
       export_Addrs(cnt1,jd_,prd_)    = cnt
       import_Addrs(cnt1,id_,prd_)    = 1
       import_Addrs(cnt1,jd_,prd_)    = cnt
       export_Addrs(cnt2,id_,prd_)    = LIloc-1
       export_Addrs(cnt2,jd_,prd_)    = cnt
       import_Addrs(cnt2,id_,prd_)    = LIloc
       import_Addrs(cnt2,jd_,prd_)    = cnt
    enddo
    if ( myRank.eq.0       ) NeibPE(1,prd_) = PEtot-1
    if ( myRank.eq.PEtot-1 ) NeibPE(2,prd_) = 0

    return
  end subroutine DefineExchangeTable_Field


  subroutine FieldCommunicate( xvc, cType, EorB )
    ! == Communication Routine for xvc(6,LIs,LJs) == !
    implicit none
    include 'mpif.h'
    character(13)   , intent(in)    :: cType
    character(1)    , intent(in)    :: EorB
    double precision, intent(inout) :: xvc(6,LIloc,LJloc)
    integer                         :: idx, cidx, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(3,CommLen)
    double precision                :: recvBuff(3,CommLen)
    integer                         :: d1_, d2_, d3_

    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    !  -- [1-1] Select Mode -- !
    cidx = 0
    if ( trim(cType).eq.'Boundary_Comm' ) cidx = bcs_
    if ( trim(cType).eq.'Periodic_Comm' ) cidx = prd_
    if ( cidx .eq.0              ) stop ' cType ?? -@FieldCommunicate- '
    if ( neibPEtot(cidx).eq.0    ) return
    if ( EorB.eq.'E' ) then
       d1_=1 ; d2_=2 ; d3_=3
    endif
    if ( EorB.eq.'B' ) then
       d1_=4 ; d2_=5 ; d3_=6
    endif
    !  -- [1-2] Packing -- !
    do neib=1, NeibPEtot(cidx)
       do idx=export_Index(neib-1,cidx)+1, export_Index(neib,cidx)
          sendBuff(1,idx) = xvc( d1_, export_Addrs(idx,id_,cidx), export_Addrs(idx,jd_,cidx) )
          sendBuff(2,idx) = xvc( d2_, export_Addrs(idx,id_,cidx), export_Addrs(idx,jd_,cidx) )
          sendBuff(3,idx) = xvc( d3_, export_Addrs(idx,id_,cidx), export_Addrs(idx,jd_,cidx) )
       enddo
    enddo

    ! ------------------------- !
    ! --- [2] ISend / IRecv --- !
    ! ------------------------- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot(cidx)
       is    = export_Index(neib-1,cidx) + 1
       len_s = 3*( export_Index(neib  ,cidx) - export_Index(neib-1,cidx) )
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib,cidx) , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot(cidx)
       ir    = import_Index(neib-1,cidx) + 1
       len_r = 3*( import_Index(neib  ,cidx) - import_Index(neib-1,cidx) )
       call MPI_Irecv( recvBuff(1,ir)    , len_r, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib,cidx) , 0    , MPI_COMM_WORLD &
            &        , request_Recv(neib), ierr )
    enddo

    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll Receive --  !
    call MPI_WaitAll( NeibPEtot(cidx), request_Recv, status_Recv, ierr )
    !  -- [3-2] Update  xVector --  !
    do neib=1, NeibPEtot(cidx)
       do idx=import_Index(neib-1,cidx)+1, import_Index(neib,cidx)
          xvc( d1_, import_Addrs(idx,id_,cidx), import_Addrs(idx,jd_,cidx) ) = recvBuff(1,idx)
          xvc( d2_, import_Addrs(idx,id_,cidx), import_Addrs(idx,jd_,cidx) ) = recvBuff(2,idx)
          xvc( d3_, import_Addrs(idx,id_,cidx), import_Addrs(idx,jd_,cidx) ) = recvBuff(3,idx)
       enddo
    enddo
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot(cidx), request_Send, status_Send, ierr )

    return
  end subroutine FieldCommunicate


  subroutine FluidCommunicate( xvc, nsp, cType )
    ! == Communication Routine for xvc(LIs,LJs,ns) == !
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: nsp
    character(13)   , intent(in)    :: cType
    double precision, intent(inout) :: xvc(0:LIloc+2,0:LJloc+2,nsp)
    integer                         :: idx, isp, cidx, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(nsp,CommLen)
    double precision                :: recvBuff(nsp,CommLen)
    cidx = 0
    if ( cType.eq.'Boundary_Comm' ) cidx = bcs_
    if ( cType.eq.'Periodic_Comm' ) cidx = prd_
    if ( cidx .eq.0               ) write(6,*) ' cType ?? -@FluidCommunicate- ', cType
    if ( neibPEtot(cidx).eq.0     ) return

    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    do neib=1, NeibPEtot(cidx)
       do idx=export_Index(neib-1,cidx)+1, export_Index(neib,cidx)
          do isp=1, nsp
             sendBuff(isp,idx) = xvc( export_Addrs(idx,id_,cidx), export_Addrs(idx,jd_,cidx), isp )
          enddo
       enddo
    end do
    ! ------------------------- !
    ! --- [2] ISend / IRecv --- !
    ! ------------------------- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot(cidx)
       is    =   export_Index(neib-1,cidx) + 1
       len_s = ( export_Index(neib  ,cidx) - export_Index(neib-1,cidx) )*nsp
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib,cidx) , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot(cidx)
       ir    =   import_Index(neib-1,cidx) + 1
       len_r = ( import_Index(neib  ,cidx) - import_Index(neib-1,cidx) )*nsp
       call MPI_Irecv( recvBuff(1,ir)    , len_r, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib,cidx) , 0    , MPI_COMM_WORLD &
            &        , request_Recv(neib), ierr )
    enddo

    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll Receive --  !
    call MPI_WaitAll( NeibPEtot(cidx), request_Recv, status_Recv, ierr )
    !  -- [3-2] Update  xVector --  !
    do neib=1, NeibPEtot(cidx)
       do idx=import_Index(neib-1,cidx)+1, import_Index(neib,cidx)
          do isp=1, nsp
             xvc( import_Addrs(idx,id_,cidx), import_Addrs(idx,jd_,cidx), isp ) = recvBuff(isp,idx)
          enddo
       enddo
    enddo
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot(cidx), request_Send, status_Send, ierr )

    return
  end subroutine FluidCommunicate


  subroutine CurrentCommunicate( xvc, nCmp, cType )
    ! == Communication Routine for xvc(LIs,LJs,ns) == !
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: nCmp
    character(13)   , intent(in)    :: cType
    double precision, intent(inout) :: xvc(nCmp,0:LIloc+2,0:LJloc+2)
    integer                         :: idx, icmp, cidx, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(nCmp,CommLen)
    double precision                :: recvBuff(nCmp,CommLen)
    cidx = 0
    if ( cType.eq.'Boundary_Comm' ) cidx = bcs_
    if ( cType.eq.'Periodic_Comm' ) cidx = prd_
    if ( cidx .eq.0               ) write(6,*) ' cType ?? -@FluidCommunicate- ', cType
    if ( neibPEtot(cidx).eq.0     ) return

    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    do neib=1, NeibPEtot(cidx)
       do idx=export_Index(neib-1,cidx)+1, export_Index(neib,cidx)
          do icmp=1, nCmp
             sendBuff(icmp,idx) = xvc( icmp, export_Addrs(idx,id_,cidx), export_Addrs(idx,jd_,cidx) )
          enddo
       enddo
    end do
    ! ------------------------- !
    ! --- [2] ISend / IRecv --- !
    ! ------------------------- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot(cidx)
       is    = export_Index(neib-1,cidx) + 1
       len_s = export_Items(neib  ,cidx) * nCmp
       call MPI_Isend( sendBuff(1,is)    , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib,cidx) , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot(cidx)
       ir    = import_Index(neib-1,cidx) + 1
       len_r = import_Items(neib  ,cidx) * nCmp
       call MPI_Irecv( recvBuff(1,ir)    , len_r, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib,cidx) , 0    , MPI_COMM_WORLD &
            &        , request_Recv(neib), ierr )
    enddo

    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll Receive --  !
    call MPI_WaitAll( NeibPEtot(cidx), request_Recv, status_Recv, ierr )
    !  -- [3-2] Update  xVector --  !
    do neib=1, NeibPEtot(cidx)
       do idx=import_Index(neib-1,cidx)+1, import_Index(neib,cidx)
          do icmp=1, nCmp
             xvc( icmp, import_Addrs(idx,id_,cidx), import_Addrs(idx,jd_,cidx) ) = recvBuff(icmp,idx)
          enddo
       enddo
    enddo
    !  -- [3-3] WaitAll Send    --  !
    call MPI_WaitAll( NeibPEtot(cidx), request_Send, status_Send, ierr )

    return
  end subroutine CurrentCommunicate

end module fMPIComMod
