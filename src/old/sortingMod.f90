module sortingMod
  implicit none
  integer, parameter :: NInsertMax = 90
  integer, parameter :: NOMPMax    = 8
contains

  ! ========================== !
  ! === Single Thread Ver. === !
  ! ========================== !
  recursive subroutine sortptcl( nSize, arr, first, last, ptcl )
    implicit none
    integer         , intent(in)    :: first, last, nSize
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(8,nSize)
    integer                         :: pivot
    
    ! ------------------------------- !
    ! --- [1]   Insart  Sorting   --- !
    ! ------------------------------- !
    ! if ( ( last-first+1 ).lt.NInsertMax ) then
    !    call InsertSort1( last-first+1, arr(first:last), ptcl(:,first:last) )
    !    return
    ! endif
    ! ------------------------------- !
    ! --- [2]   Partition         --- !
    ! ------------------------------- !
    call partition( nSize, arr, first, last, ptcl, pivot )
    ! ------------------------------- !
    ! --- [3]  Sort Each vector   --- !
    ! ------------------------------- !
    if ( first.lt.pivot-1 ) call sortptcl( nSize, arr, first, pivot-1, ptcl )
    if ( pivot.lt.last    ) call sortptcl( nSize, arr, pivot, last   , ptcl )
    return
  end subroutine sortptcl
  
  
  ! ========================== !
  ! === Multi. Thread Ver. === !
  ! ========================== !
  subroutine ParallelSort( nSize, arr, first, last, ptcl )
    !$ use omp_lib
    implicit none
    integer         , intent(in)    :: first, last, nSize
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(8,nSize)
    integer                         :: ith, node(NOMPMax+1), Nthreads
    integer                         :: ithF, ithL, stride
    ! ---------------------------------------------- !
    ! -- arr(nSize)   :: Order label of sorting   -- !
    ! -- nSize        :: size of a & ptcl         -- !
    ! -- node(m)      :: Start Index of m-th part -- !
    ! -- NOMPMax      :: up to 8                  -- !
    ! ---------------------------------------------- !

    ! --- [1] Preparation --- ! 
    !$omp parallel default(none) shared(Nthreads)
    !$omp single
    Nthreads = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    Nthreads = min( Nthreads, NOMPMax )
    node(:)  = 0
    node(1)  = first
    node(9)  = last + 1
    stride   = NOMPMax / Nthreads

    ! --- [2] pre-Partition ( Manual ) --- !
    if ( Nthreads.ge.2 ) then
       call partition( nSize, arr, node(1), node(9)-1, ptcl, node(5) )
    endif
    if ( Nthreads.ge.4 ) then
       !$omp parallel num_threads(2) default(none) shared(nSize,arr,ptcl,node) private(ith)
       !$ ith = omp_get_thread_num() + 1
       if ( ith.eq.1 ) call partition( nSize, arr, node(1), node(5)-1, ptcl, node(3) )
       if ( ith.eq.2 ) call partition( nSize, arr, node(5), node(9)-1, ptcl, node(7) )
       !$omp end parallel
    endif
    if ( Nthreads.ge.8 ) then
       !$omp parallel num_threads(4) default(none) shared(nSize,arr,ptcl,node) private(ith)
       !$ ith = omp_get_thread_num() + 1
       if ( ith.eq.1 ) call partition( nSize, arr, node(1), node(3)-1, ptcl, node(2) )
       if ( ith.eq.2 ) call partition( nSize, arr, node(3), node(5)-1, ptcl, node(4) )
       if ( ith.eq.3 ) call partition( nSize, arr, node(5), node(7)-1, ptcl, node(6) )
       if ( ith.eq.4 ) call partition( nSize, arr, node(7), node(9)-1, ptcl, node(8) )
       !$omp end parallel
    endif
    
    ! --- [3] Main Sorting ( Auto: Recursive ) --- !
    !$omp parallel num_threads(Nthreads) default(none) shared(nSize,arr,ptcl,node,stride) private(ith,ithF,ithL)
    !$ ith = omp_get_thread_num() + 1
    ithF   = stride*( ith-1 ) + 1
    ithL   = stride*( ith   ) + 1
    call sortptcl( nSize, arr, node(ithF), node(ithL)-1, ptcl )
    !$omp end parallel
    
    return
  end subroutine ParallelSort
  
  
  subroutine partition( nSize, arr, first, last, ptcl, pivot )
    implicit none
    integer         , intent(in)    :: nSize, first, last
    integer         , intent(inout) :: pivot
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(8,nSize)
    integer                         :: i, j, p, t
    double precision                :: d(8)
    
    ! -- Partition ( [ arr(i)<p ] // [ arr(i)>p ] ) -- !
    p = arr( (first+last) / 2 )
    i = first
    j = last
    do
       do while ( arr(i).lt.p )
          i = i+1
       enddo
       do while ( p.lt.arr(j) )
          j = j-1
       enddo
       if ( i.ge.j ) exit
       t    = arr(i)   ; arr(i)    = arr(j)   ; arr(j)    = t
       d(:) = ptcl(:,i); ptcl(:,i) = ptcl(:,j); ptcl(:,j) = d(:)
       i = i+1
       j = j-1
    enddo
    pivot = i
    
    return
  end subroutine partition
  

  ! ----------------------------------------------------- !
  ! -- "Insert Sorting" is faster than "Quick sorting" -- !
  ! --                                for small vector -- !
  ! ----------------------------------------------------- !
  ! === Insert Sorting ( type.1 ) === !
  subroutine InsertSort1( nSize, arr, ptcl )
    implicit none
    integer         , intent(in)    :: nSize
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(8,nSize)
    integer                         :: i, j, tmp
    double precision                :: dtmp(8)

    do i=2, nSize
       tmp     = arr(i)
       dtmp(:) = ptcl(:,i)
       if ( arr(i-1) .gt. tmp ) then
          j    = i
          do 
             arr(j)    = arr(j-1)
             ptcl(:,j) = ptcl(:,j-1)
             j         = j-1
             if ( ( j.le.1 ).or.( arr(j-1) .le. tmp ) ) exit
          enddo
          arr(j)       = tmp
          ptcl(:,j)    = dtmp(:)
       endif
    enddo
          
    return
  end subroutine InsertSort1


  ! === Insert Sorting ( type.2 ) === !
  subroutine InsertSort2( nSize, arr, ptcl )
    implicit none
    integer         , intent(in)    :: nSize
    integer         , intent(inout) :: arr(nSize)
    double precision, intent(inout) :: ptcl(8,nSize)
    integer                         :: i, j, k, tmp, tid
    integer                         :: indx( NInsertMax )
    double precision                :: ptclcp( 8, NInsertMax )

    do i=1, nSize
       indx(i) = i
       do k = 1, 8
          ptclcp(k,i) = ptcl(k,i)
       enddo
    enddo
    do i=2, nSize
       tmp = arr(i)
       tid = indx(i)
       if ( arr(i-1) .gt. tmp ) then
          j = i
          do 
             arr(j)    = arr(j-1)
             indx(j) = indx(j-1)
             j       = j-1
             if ( ( j.le.1 ).or.( arr(j-1) .le. tmp ) ) exit
          enddo
          arr(j)    = tmp
          indx(j) = tid
       endif
    enddo
    do i=1, nSize
       do k=1, 8
          ptcl(k,i) = ptclcp(k,indx(i))
       enddo
    enddo
    
    return
  end subroutine InsertSort2

  

  ! === Quick Sorting For Index ( integer ) === !
  recursive subroutine qsIndex( nSize, arr, first, last, indx )
    implicit none
    integer         , intent(in)    :: first, last, nSize
    integer         , intent(inout) :: indx(nSize)
    double precision, intent(inout) :: arr(nSize)
    integer                         :: i, j, k
    double precision                :: t, p

    p = arr( (first+last) / 2 )
    i = first
    j = last
    do
       do while ( arr(i) .lt. p )
          i=i+1
       end do
       do while ( p .lt. arr(j) )
          j=j-1
       end do
       if ( i .ge. j ) exit
       t = arr(i)  ; arr(i)  = arr(j)  ; arr(j)  = t
       k = indx(i) ; indx(i) = indx(j) ; indx(j) = k
       i=i+1
       j=j-1
    end do

    if (first < i - 1) call qsIndex(nSize, arr, first, i - 1, indx)
    if (j + 1 < last ) call qsIndex(nSize, arr, j + 1, last , indx)

    return
  end subroutine qsIndex

  
end module sortingMod

