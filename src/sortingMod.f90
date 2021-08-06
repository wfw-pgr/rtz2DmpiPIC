module sortingMod
  implicit none
  integer, parameter :: NInsertMax = 50
  integer, parameter :: NOMPMax    = 8
  integer, parameter :: nCmp       = 8
contains

  ! ========================== !
  ! === Single Thread Ver. === !
  ! ========================== !
  recursive subroutine sortptcl( nSize, arr, first, last, ptcl )
    implicit none
    integer         , intent(in)    :: first, last, nSize
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(nCmp,nSize)
    integer                         :: ipivot, jpivot
    
    ! ------------------------------- !
    ! --- [1]   Insart  Sorting   --- !
    ! ------------------------------- !
    if ( ( last-first+1 ).lt.NInsertMax ) then
       call InsertSort2( last-first+1, arr(first:last), ptcl(:,first:last) )
       return
    endif
    ! ------------------------------- !
    ! --- [2]   Partition         --- !
    ! ------------------------------- !
    call partition( nSize, arr, first, last, ptcl, ipivot, jpivot )
    ! ------------------------------- !
    ! --- [3]  Sort Each vector   --- !
    ! ------------------------------- !
    if ( first.lt.(ipivot-1) ) call sortptcl( nSize, arr, first   , ipivot-1, ptcl )
    if (  last.gt.(jpivot+1) ) call sortptcl( nSize, arr, jpivot+1, last    , ptcl )
    return
  end subroutine sortptcl
  
  
  ! ========================== !
  ! === Multi. Thread Ver. === !
  ! ========================== !
  subroutine ParallelSort( nSize, arr, first, last, ptcl, Nthreads )
    !$ use omp_lib
    implicit none
    integer         , intent(in)    :: first, last, nSize, Nthreads
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(nCmp,nSize)
    integer                         :: ith
    integer                         :: Nth, step, ompFirst, ompLast
    integer                         :: node(NOMPMax*2)
    ! ---------------------------------------------- !
    ! -- arr(nSize)   :: Order label of sorting   -- !
    ! -- nSize        :: size of a & ptcl         -- !
    ! -- node(m)      :: Start Index of m-th part -- !
    ! -- NOMPMax      :: up to 8                  -- !
    ! ---------------------------------------------- !
    
    ! -------------------------------------- !
    ! --- [1]  Preparation               --- ! 
    ! -------------------------------------- !
    Nth      = min( Nthreads, NOMPMax )
    step     = NOMPMax / Nth * 2
    node( 1) = first
    node(16) = last
    
    ! -------------------------------------- !
    ! --- [2]  pre-Partition ( Manual )  --- !
    ! -------------------------------------- !
    !  -- [2-1] 1/2 partition -- !
    if ( Nth.ge.2 ) then
       call partition( nSize, arr, node(1), node(16), ptcl, node(8), node(9) )
    endif
    !  -- [2-2] 1/4 partition -- !
    if ( Nth.ge.4 ) then
       !$omp parallel num_threads(2) default(none) shared(nSize,arr,ptcl,node) private(ith)
       !$ ith = omp_get_thread_num() + 1
       if ( ith.eq.1 ) call partition( nSize, arr, node( 1)  , node( 8)-1, ptcl, node( 4), node( 5) )
       if ( ith.eq.2 ) call partition( nSize, arr, node( 9)+1, node(16)  , ptcl, node(12), node(13) )
       !$omp end parallel
    endif
    !  -- [2-3] 1/8 partition -- !
    if ( Nth.ge.8 ) then
       !$omp parallel num_threads(4) default(none) shared(nSize,arr,ptcl,node) private(ith)
       !$ ith = omp_get_thread_num() + 1
       if ( ith.eq.1 ) call partition( nSize, arr, node( 1)  , node( 4)-1, ptcl, node( 2), node( 3) )
       if ( ith.eq.2 ) call partition( nSize, arr, node( 5)+1, node( 8)-1, ptcl, node( 6), node( 7) )
       if ( ith.eq.3 ) call partition( nSize, arr, node( 9)+1, node(12)-1, ptcl, node(10), node(11) )
       if ( ith.eq.4 ) call partition( nSize, arr, node(13)+1, node(16)  , ptcl, node(14), node(15) )
       !$omp end parallel
    endif

    ! -------------------------------------- !
    ! --- [3]  Auto (recursive) Sorting  --- !
    ! -------------------------------------- !
    !$omp parallel num_threads(Nth) default(none) &
    !$omp shared(nSize,arr,ptcl,step,node,Nth) private(ith,ompFirst,ompLast)
    !$ ith = omp_get_thread_num() + 1
    ompFirst = node( (ith-1)*step+1 ) + 1
    ompLast  = node( (ith  )*step   ) - 1
    if ( ith.eq.1   ) ompFirst = ompFirst - 1
    if ( ith.eq.Nth ) ompLast  = ompLast  + 1
    call sortptcl( nSize, arr, ompFirst, ompLast, ptcl )
    !$omp end parallel
    
    return
  end subroutine ParallelSort


  subroutine partition( nSize, arr, first, last, ptcl, ipivot, jpivot )
    implicit none
    integer         , intent(in)    :: nSize, first, last
    integer         , intent(inout) :: ipivot, jpivot
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(nCmp,nSize)
    integer                         :: ip, jp, pv, ibuff
    double precision                :: dbuff(nCmp)
    
    ! -- Partition ( [ arr(i)<p ] // [ arr(i)>p ] ) -- !
    pv = arr( (first+last) / 2 )
    ip = first
    jp = last
    do
       do while ( arr(ip).lt.pv )
          ip = ip+1
       enddo
       do while ( pv.lt.arr(jp) )
          jp = jp-1
       enddo
       if ( ip.ge.jp ) exit
       ibuff    =  arr(ip)   ;  arr(ip)   = arr(jp)    ; arr(jp)    = ibuff
       dbuff(:) = ptcl(:,ip) ; ptcl(:,ip) = ptcl(:,jp) ; ptcl(:,jp) = dbuff(:)
       ip       = ip+1
       jp       = jp-1
    enddo
    ipivot = ip
    jpivot = jp
    
    return
  end subroutine partition
  

  ! ----------------------------------------------------- !
  ! -- "Insert Sorting" is faster than "Quick sorting" -- !
  ! --                                for small vector -- !
  ! ----------------------------------------------------- !
  ! === Insert Sorting ( type.1 ) === !
  subroutine InsertSort1( nSize, lbl, ptc )
    implicit none
    integer         , intent(in)    :: nSize
    integer         , intent(inout) :: lbl(nSize)
    double precision, intent(inout) :: ptc(nCmp,nSize)
    integer                         :: ip, jp, lbl_
    double precision                :: ptc_(nCmp)

    ! --- Primitive Insert Sorting --- !
    do ip=2, nSize
       lbl_     = lbl(ip)
       ptc_(:)  = ptc(:,ip)
       jp       = ip-1
       do while( jp.ge.1 )
          if ( lbl(jp).le.lbl_ ) exit
          lbl(  jp+1) = lbl(  jp)
          ptc(:,jp+1) = ptc(:,jp)
          jp          = jp-1
       enddo
       lbl(jp+1)   = lbl_
       ptc(:,jp+1) = ptc_(:)
    enddo
          
    return
  end subroutine InsertSort1


  ! === Insert Sorting ( type.2 : "Index Exchange ver." ) === !
  subroutine InsertSort2( nSize, lbl, ptc )
    implicit none
    integer         , intent(in)    :: nSize
    integer         , intent(inout) :: lbl(nSize)
    double precision, intent(inout) :: ptc(nCmp,nSize)
    integer                         :: ip, jp, lbl_, idx_
    integer                         :: idx(nSize)
    double precision                :: ptc_(nCmp,nSize)
    
    ! --- [1] Make Index & Data Copying --- !
    do ip=1, nSize
       idx (  ip) = ip
       ptc_(:,ip) = ptc(:,ip)
    enddo
    ! --- [2] Insert Sorting ( Main ) --- !
    do ip=2, nSize
       lbl_     = lbl(ip)
       idx_     = idx(ip)
       jp       = ip-1
       do while( jp.ge.1 )
          if ( lbl(jp).le.lbl_ ) exit
          lbl(jp+1) = lbl(jp)
          idx(jp+1) = idx(jp)
          jp        = jp-1
       enddo
       lbl(jp+1)    = lbl_
       idx(jp+1)    = idx_
    enddo
    do ip=1, nSize
       ptc(:,ip) = ptc_(:,idx(ip))
    enddo
    
    return
  end subroutine InsertSort2
  

  ! === Quick Sorting For Index ( integer )      === !
  !  == Sort Ref. :: Data  => arr   ( double  )  ==  !
  !  == target    :: Index => indx  ( integer )  ==  !
  recursive subroutine qsIndex( nSize, arr, first, last, indx )
    implicit none
    integer         , intent(in)    :: first, last, nSize
    integer         , intent(inout) :: indx(nSize)
    double precision, intent(inout) :: arr(nSize)
    integer                         :: ip, jp, ibuff
    double precision                :: dbuff , pivot

    ! --- [1] Preparation    --- !
    pivot = arr( (first+last) / 2 )
    ip    = first
    jp    = last
    ! --- [2] Quick Sorting  --- !
    do
       do while ( arr(ip).lt.pivot  )
          ip = ip+1
       end do
       do while ( arr(jp).gt.pivot )
          jp = jp-1
       end do
       if ( ip.ge.jp ) exit
       dbuff =  arr(ip) ;  arr(ip) =  arr(jp) ;  arr(jp) = dbuff
       ibuff = indx(ip) ; indx(ip) = indx(jp) ; indx(jp) = ibuff
       ip = ip+1
       jp = jp-1
    end do
    ! --- [3] Recursive Call --- !
    if ( first.lt.(ip-1) ) call qsIndex( nSize, arr, first, ip-1, indx )
    if (  last.gt.(jp+1) ) call qsIndex( nSize, arr, jp+1 , last, indx )
    
    return
  end subroutine qsIndex

  
end module sortingMod

