module depBdrcMod
contains

  subroutine depFoldingBoundary( dep )
    use constants,              only : LIs, LJs, ns, myRank, PEtot
    use constants,              only : Boundary1__pt
    use fMPIComMod,             only : FluidCommunicate
    use jMPIComMod,             only : DepositOverlay
    implicit none
    integer                         :: i, j, k
    double precision, intent(inout) :: dep(0:LIs+2,0:LJs+2,ns)

    ! ------------------------------ !
    ! --- [1] Folding Boundary_2 --- !
    ! ------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(dep,LIs,LJs) private(i,k)
    !$omp do
    do k=1, ns
       do i=0, LIs+2
          dep(i,    2,k) = dep(i,    2,k) + dep(i,    1,k)
          dep(i,    3,k) = dep(i,    3,k) + dep(i,    0,k)
          dep(i,LJs-1,k) = dep(i,LJs-1,k) + dep(i,LJs  ,k)
          dep(i,LJs-2,k) = dep(i,LJs-2,k) + dep(i,LJs+1,k)
          dep(i,LJs-3,k) = dep(i,LJs-3,k) + dep(i,LJs+2,k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ------------------------------ !
    ! --- [2]  Inter-PE Overlay  --- !
    ! ------------------------------ !
    call DepositOverlay( dep, ns )

    ! ------------------------------ !
    ! --- [3] Folding Boundary_1 --- !
    ! ------------------------------ !
    !  -- [3-1]  Left  --  !
    if ( trim(Boundary1__pt).eq.'reflect' ) then
       if ( myRank.eq.0       ) then
          !$omp parallel default(none) &
          !$omp shared(dep,LIs,LJs) private(j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                dep(    2,j,k) = dep(    2,j,k) + dep(    1,j,k)
                dep(    3,j,k) = dep(    3,j,k) + dep(    0,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       !  -- [2-2]  Left  --  !
       if ( myRank.eq.PEtot-1 ) then
          !$omp parallel default(none) &
          !$omp shared(dep,LIs,LJs) private(j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                dep(LIs-1,j,k) = dep(LIs-1,j,k) + dep(LIs  ,j,k)
                dep(LIs-2,j,k) = dep(LIs-2,j,k) + dep(LIs+1,j,k)
                dep(LIs-3,j,k) = dep(LIs-3,j,k) + dep(LIs+2,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
    endif

    ! ------------------------------- !
    ! --- [3] Boundary_1 copying  --- !
    ! ------------------------------- !
    !  -- [3-1] Set Neibour Data --   !
    call FluidCommunicate( dep, ns )
    !  -- [3-2] Boundary Copy    --   !
    if ( trim(Boundary1__pt).eq.'reflect' ) then
       if ( myRank.eq.0       ) then
          !$omp parallel default(none) &
          !$omp shared(dep,LIs,LJs) private(j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                dep(  1,j,k) = dep(    2,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       if ( myRank.eq.PEtot-1 ) then
          !$omp parallel default(none) &
          !$omp shared(dep,LIs,LJs) private(j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                dep(LIs,j,k) = dep(LIs-1,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
    endif
    
    ! ------------------------------ !
    ! --- [4] Copying Boundary_2 --- !
    ! ------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(dep,LIs,LJs) private(i,k)
    !$omp do
    do k=1, ns
       do i=1, LIs
          dep(i,  1,k)    = dep(i,    2,k)
          dep(i,LJs,k)    = dep(i,LJs-1,k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine depFoldingBoundary


  subroutine FoldingCopy_CC( FCC )
    use constants, only              : LIs, LJs, myRank, PEtot
    implicit none
    integer                         :: i, j
    double precision, intent(inout) :: FCC(LIs,LJs)

    ! ------------------------------ !
    ! --- [1] Copying Boundary_1 --- !
    ! ------------------------------ !
    !  -- [1-1]  Left Edge        -- !
    if ( myRank.eq.0       ) then
       !$omp parallel default(none) &
       !$omp shared(FCC,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs
          FCC(  1,j) = FCC(    2,j)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    !  -- [1-2] Right Edge        -- !
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(FCC,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs
          FCC(LIs,j) = FCC(LIs-1,j)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    
    ! ------------------------------ !
    ! --- [2] Copying Boundary_2 --- !
    ! ------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(FCC,LIs,LJs) private(i,j)
    !$omp do
    do i=1, LIs
       FCC(i,  1) = FCC(i,    2)
       FCC(i,LJs) = FCC(i,LJs-1)
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine FoldingCopy_CC

  
end module depBdrcMod
