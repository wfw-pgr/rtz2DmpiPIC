module allocatMod
contains

  ! =================================================================== !
  ! ===  InitVariables  :: Initialize variables                     === !
  ! =================================================================== !
  subroutine InitVariables
    use constants, only : myRank
    implicit none
    
    ! --- [1] Call Allocate Routine --- !
    call allocateParticle
    call allocateField
    call allocateCurrent
    call allocateMoments
    call allocateOthers
    
    ! --- [2] Display --- !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '---------------         [     VARIABLES     ]          ---------------'
       write(6,'(2x,a)') '[InitVariables                ] Particle  variables  :: ( Initialized )'
       write(6,'(2x,a)') '[InitVariables                ] Field     variables  :: ( Initialized )'
    endif
    return
  end subroutine InitVariables


  subroutine allocateParticle
    use constants, only : LI, LJ, ns, nptMax
    use variables, only : pxv, goban, pCellIdx
    implicit none
    integer            :: k, m, cmp
    integer, parameter :: OutIdx = (LI-2)*(LJ-2)+1

    ! --- [1]  Allocation   --- !
    allocate( pxv(8,nptMax,ns)     )
    allocate( goban(nptMax)        )
    allocate( pCellIdx(nptMax,ns)  )
    
    ! --- [2]  Initialize   --- !
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
    return
  end subroutine allocateParticle


  subroutine allocateField
    use constants, only : LIs, LJs
    use variables, only : EBf, EBr, EBo
    implicit none
    integer            :: i, j, cmp
    
    ! --- [1]  Allocation   --- !
    allocate( EBf(6,LIs  ,LJs  )   )
    allocate( EBo(6,LIs  ,LJs  )   )
    allocate( EBr(6,LIs+1,LJs+1)   )
    
    ! --- [2]  Initialize   --- !
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
    return
  end subroutine allocateField


  subroutine allocateCurrent
    use constants, only : LIs, LJs, ns, OMPNumThreads
    use variables, only : Jcr, JcrW
    implicit none
    integer            :: i, j, k, cmp, m
    
    ! --- [1]  Allocation   --- !
    allocate( Jcr (3,LIs,LJs,0:ns)                     )
    allocate( JcrW(3,0:LIs+2,0:LJs+2,ns,OMPNumThreads) )

    ! --- [2]  Initialize   --- !
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
    !$omp parallel do private(i,j,k,m,cmp)
    do m=1, OMPNumThreads
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                do cmp=1, 3
                   JcrW(cmp,i,j,k,m) = 0.d0
                enddo
             enddo
          enddo
       enddo
    enddo
    !$omp end parallel do
    return
  end subroutine allocateCurrent


  subroutine allocateMoments
    use constants , only : LIs, LJs, ns, OMPNumThreads
    use variables , only : sumupE, sumupW, sumupP
    use momentsMod, only : nMoments, fMoments, nPrsQth
    implicit none
    integer             :: i, j, k, m, cmp
    
    ! --- [1]  Allocation   --- !
    allocate( fMoments(LIs,LJs,ns,nMoments) )
    allocate( sumupW(        0:LIs+2,0:LJs+2,ns,OMPNumThreads) )
    allocate( sumupP(nPrsQth,0:LIs+2,0:LJs+2,ns,OMPNumThreads) )
    allocate( sumupE(        0:LIs+2,0:LJs+2,ns ) )
    
    ! --- [2]  Initialize   --- !
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
                sumupP(7,i,j,k,m) = 0.d0
                sumupP(8,i,j,k,m) = 0.d0
                sumupP(9,i,j,k,m) = 0.d0
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
    return
  end subroutine allocateMoments


  subroutine allocateOthers
    use constants , only : LIs, LJs, ns, OMPNumThreads
    use variables , only : perCell, perCellW, cellTop
    use variables , only : w_p
    implicit none
    integer             :: i, j, k, m
    
    ! --- [1]  Allocation   --- !
    allocate( cellTop( LIs,LJs,ns)               )
    allocate( perCell( LIs,LJs,ns)               )
    allocate( perCellW(LIs,LJs,ns,OMPNumThreads) )
    allocate( w_p(LIs,LJs,ns)                    )
    ! --- [2]  Initialize   --- !
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
    !$omp parallel do private(i,j,k)
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             w_p(i,j,k) = 0.d0
          enddo
       enddo
    enddo
    !$omp end parallel do

    return
  end subroutine allocateOthers


  subroutine allocateEBrEBo
    use constants, only : LIs, LJs
    use variables, only : EBr, EBo
    implicit none
    integer            :: i, j, cmp

    ! --- [1]  Allocation   --- !
    allocate( EBo(6,LIs  ,LJs  )   )
    allocate( EBr(6,LIs+1,LJs+1)   )

    ! --- [2]  Initialize   --- !
    !$omp parallel do private(i,j,cmp)
    do j=1, LJs
       do i=1, LIs
          do cmp=1, 6
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
    return
  end subroutine allocateEBrEBo

  
end module allocatMod
