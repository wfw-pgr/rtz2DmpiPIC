module nMPIComMod
  implicit none
  ! -- Constants        -- !
  integer         , parameter   :: dim           = 2
  integer         , parameter   :: NeibPEtot_Max = 2
  integer         , parameter   :: nLineL        = 2
  integer         , parameter   :: nLineR        = 3
  integer                       :: nLine
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
  
  
  subroutine DefineExchangeTable_Deposit( LIs, LJs )
    implicit none
    include 'mpif.h'
    integer, intent(in) :: LIs, LJs
    integer             :: i, j, cnt, ierr
    ! -------------------------------------------- !
    ! --- [1] Define Grid Number and Partition --- !
    ! -------------------------------------------- !
    !  -- [1-1]  PE infomation       --- !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )
    LIloc   = LIs
    LJloc   = LJs
    nLine   = nLineL + nLineR
    CommLen = LJloc  * nLine
    
    ! -------------------------------------------- !
    ! --- [2] Allocate Variables for MPI Comm. --- !
    ! -------------------------------------------- !
    allocate( export_Index(0:NeibPEtot_Max), import_Index(0:NeibPEtot_Max) )
    allocate( export_Items(0:NeibPEtot_Max), import_Items(0:NeibPEtot_Max) )
    allocate( export_Addrs(    CommLen,dim), import_Addrs(    CommLen,dim) )
    allocate( request_Send(NeibPEtot_Max)  , request_Recv(NeibPEtot_Max)   )
    allocate(  status_Send(MPI_STATUS_SIZE,NeibPEtot_Max), &
         &     status_Recv(MPI_STATUS_SIZE,NeibPEtot_Max)  )
    
    ! -------------------------------------------- !
    ! --- [3] Generate MPI Communication Table --- !
    ! -------------------------------------------- !
    !  -- [3-0] Exception for PEtot=1 -- !
    if ( PEtot.eq.1        ) then
       NeibPEtot = 0
       return
    endif
    ! --------------------------- !
    ! -- [3-1] Deposit Overlay -- !
    ! --------------------------- !
    export_Items = 0
    import_Items = 0
    export_Index = 0
    import_Index = 0
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
    return
  end subroutine DefineExchangeTable_Deposit



  
end module nMPIComMod
