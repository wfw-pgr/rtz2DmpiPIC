module allocatMod
contains
  
  subroutine InitVariables
    use constants , only : myRank  , OMPNumThreads, PEtot
    use constants , only : LI , LJ , LIs, LJs, ns , nptMax
    use constants , only : er_, et_, ez_, br_, bt_, bz_
    use variables , only : pxv, EBf, EBr, EBo, Jcr, JcrW, x1s
    use variables , only : goban, pCellIdx, perCell, cellTop, w_p
    use momentsMod, only : nMoments, fMoments, perCellW
    use momentsMod, only : sumupE, sumupW, sumupP
    implicit none
    include 'mpif.h'
    integer            :: i, j, m, k, n, cmp
    integer, parameter :: OutIdx = (LI-2)*(LJ-2)+1

    ! ============================== !
    ! === [1]     Allocate       === !
    ! ============================== !
    !  -- [1-1]    Particle     -- !
    allocate( pxv(8,nptMax,ns)     )
    allocate( goban(nptMax)        )
    allocate( pCellIdx(nptMax,ns)  )
    !  -- [1-2]    Field        -- !
    allocate( EBf(6,LIs  ,LJs  )   )
    allocate( EBo(6,LIs  ,LJs  )   )
    allocate( EBr(6,LIs+1,LJs+1)   )
    !  -- [1-3]    Current      -- !
    allocate( Jcr( 3,LIs,LJs,0:ns) )
    allocate( JcrW  (3,0:LIs+2,0:LJs+2,ns,OMPNumThreads) )
    !  -- [1-4]    Moment       -- !
    allocate( fMoments(LIs,LJs,ns,nMoments) )
    allocate( sumupW(  0:LIs+2,0:LJs+2,ns,OMPNumThreads) )
    allocate( sumupP(6,0:LIs+2,0:LJs+2,ns,OMPNumThreads) )
    allocate( sumupE(  0:LIs+2,0:LJs+2,ns ) )
    !  -- [1-5]    PerCell      -- !
    allocate( cellTop( LIs,LJs,ns) )
    allocate( perCell( LIs,LJs,ns) )
    allocate( perCellW(LIs,LJs,ns,OMPNumThreads) )
    !  -- [1-6]    Work Space   -- !
    allocate( x1s(LIs) )
    allocate( w_p(LIs,LJs,ns) )

    ! ============================== !
    ! === [2]    Initialize      === !
    ! ============================== !
    !  --------------------------- !
    !  -- [2-1]    Particle     -- !
    !  --------------------------- !
    !$omp parallel do private(m,k)
    do k=1, ns
       do m=1, nptMax
          do cmp=1, 8
             pxv(cmp,m,k) = 0.d0
          enddo
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(m,k)
    do m=1, nptMax
       goban(m)         = OutIdx
       do k=1, ns
          pCellIdx(m,k) = OutIdx
       enddo
    enddo
    !$omp end parallel do

    !  --------------------------- !
    !  -- [2-2]    Field        -- !
    !  --------------------------- !
    !$omp parallel do private(i,j,cmp)
    do j=1, LJs
       do i=1, LIs
          do cmp=1, 6
             EBf(cmp,i,j) = 0.d0
             EBo(cmp,i,j) = 0.d0
          enddo
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(i,j,cmp)
    do j=1, LJs+1
       do i=1, LIs+1
          do cmp=1, 6
             EBr(cmp,i,j) = 0.d0
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    !  --------------------------- !
    !  -- [2-3]    Current      -- !
    !  --------------------------- !
    !$omp parallel do private(i,j,k,cmp)
    do k=0, ns
       do j=1, LJs
          do i=1, LIs
             do cmp=1, 3
                Jcr (cmp,i,j,k) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(i,j,k,n,cmp)
    do n=1, OMPNumThreads
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                do cmp=1, 3
                   JcrW(cmp,i,j,k,n) = 0.d0
                enddo
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    !  --------------------------- !
    !  -- [2-4]    Moment       -- !
    !  --------------------------- !
    !$omp parallel do private(i,j,k,m)
    do m=1, OMPNumThreads
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                sumupW(  i,j,k,m) = 0.d0
                sumupP(1,i,j,k,m) = 0.d0
                sumupP(2,i,j,k,m) = 0.d0
                sumupP(3,i,j,k,m) = 0.d0
                sumupP(4,i,j,k,m) = 0.d0
                sumupP(5,i,j,k,m) = 0.d0
                sumupP(6,i,j,k,m) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(i,j,k)
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             sumupE(i,j,k) = 0.d0
          enddo
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(i,j,k,cmp)
    do cmp=1, nMoments
       do k=1, ns
          do j=1, LJs
             do i=1, LIs
                fMoments(i,j,k,cmp) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do

    !  --------------------------- !
    !  -- [2-5]    PerCell      -- !
    !  --------------------------- !
    !$omp parallel do private(i,j,k)
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             perCell(i,j,k) = 0
             cellTop(i,j,k) = 0
          enddo
       enddo
    enddo
    !$omp end parallel do
    !$omp parallel do private(i,j,k,m)
    do m=1, OMPNumThreads
       do k=1, ns
          do j=1, LJs
             do i=1, LIs
                perCellW(i,j,k,m) = 0
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
    
    !  --------------------------- !
    !  -- [2-6]    Work Space   -- !
    !  --------------------------- !
    !$omp parallel do private(i,j,k)
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             w_p(i,j,k) = 0.d0
          enddo
       enddo
    enddo
    !$omp end parallel do

    
    ! ============================== !
    ! === [3]    Notification    === !
    ! ============================== !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '---------------         [     VARIABLES     ]          ---------------'
       write(6,'(2x,a)') '[InitVariables                ] Particle  variables  :: ( Initialized )'
       write(6,'(2x,a)') '[InitVariables                ] Field     variables  :: ( Initialized )'
    endif
    return
  end subroutine InitVariables
  
  
end module allocatMod
