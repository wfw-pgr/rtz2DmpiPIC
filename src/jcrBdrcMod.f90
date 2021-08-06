module jcrBdrcMod
contains
  
  subroutine jcrFoldingBC( bctype, JcrB )
    use constants , only              : LIs, LJs, ns, myRank, PEtot
    use constants , only              : jr_, jt_, jz_
    use jMPIComMod, only              : CurrentOverlay
    use fMPIComMod, only              : CurrentCommunicate
    implicit none
    integer                          :: i, j, k
    character(len=10), intent(in)    :: bctype
    double precision , intent(inout) :: JcrB(3,0:LIs+2,0:LJs+2,ns)
    logical                          :: Flag__Boundary1
    logical                          :: Flag__Boundary2

    ! ---------------------------- !
    ! --- [1] Flag Preparation --- !
    ! ---------------------------- !
    select case( bctype )
    case('Boundary_1')
       Flag__Boundary1 = .true.
       Flag__Boundary2 = .false.
    case('Boundary_2')
       Flag__Boundary1 = .false.
       Flag__Boundary2 = .true.
    case('Boundary_0')
       Flag__Boundary1 = .true.
       Flag__Boundary2 = .true.
    end select
    
    ! ---------------------------- !
    ! --- [2]   Boundary-(2)   --- !
    ! ---------------------------- !
    if ( Flag__Boundary2 ) then
       !$omp parallel default(none) &
       !$omp shared(JcrB,LIs,LJs) private(i,k)
       !$omp do
       do k=1, ns
          do i=0, LIs+2
             ! -- jr_ :: i : Half-Grid, j : Full-Grid -- !
             JcrB(jr_,i,3,k) =   JcrB(jr_,i,3,k) - JcrB(jr_,i,1,k)
             JcrB(jr_,i,4,k) =   JcrB(jr_,i,4,k) - JcrB(jr_,i,0,k)
             ! -- jt_ :: i : Half-Grid, j : Half-Grid -- !
             JcrB(jt_,i,2,k) =   JcrB(jt_,i,2,k) + JcrB(jt_,i,1,k)
             JcrB(jt_,i,3,k) =   JcrB(jt_,i,3,k) + JcrB(jt_,i,0,k)
             ! -- jz_ :: i : Full-Grid, j : Half-Grid -- !
             JcrB(jz_,i,2,k) =   JcrB(jz_,i,2,k) + JcrB(jz_,i,1,k)
             JcrB(jz_,i,3,k) =   JcrB(jz_,i,3,k) + JcrB(jz_,i,0,k)
             ! -- jr  :: jr(r=0) = 0, jr=anti-Symmetric -- !
             JcrB(jr_,i,2,k) =   0.d0
             JcrB(jr_,i,1,k) = - JcrB(jr_,i,3,k)
             JcrB(jr_,i,0,k) = - JcrB(jr_,i,4,k)
             ! -- jt  :: jt=Symmetric -- !
             JcrB(jt_,i,1,k) = + JcrB(jt_,i,2,k)
             JcrB(jt_,i,0,k) = + JcrB(jt_,i,3,k)
             ! -- jz  :: jz=Symmetric -- !
             JcrB(jz_,i,1,k) =   JcrB(jz_,i,2,k)
             JcrB(jz_,i,0,k) =   JcrB(jz_,i,3,k)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif

    ! ------------------------------------------ !
    ! --- [3]   Boundary-(1) :: (i) Overlay  --- !
    ! ------------------------------------------ !
    if ( Flag__Boundary1 ) then
       !  -- [3-1] Overlay Commu. -- ! 
       do k=1, ns
          call CurrentOverlay( JcrB(1:3,0:LIs+2,0:LJs+2,k) )
       enddo
       !  -- [3-2] Left Side Edge -- ! 
       if ( myRank.eq.0 ) then
          !$omp parallel default(none) &
          !$omp shared(JcrB,LIs,LJs) private(j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                ! -- jr_ :: i : Half-Grid, j : Full-Grid -- !
                JcrB(jr_,    2,j,k) = JcrB(jr_,    2,j,k) + JcrB(jr_, 1   ,j,k)
                JcrB(jr_,    3,j,k) = JcrB(jr_,    3,j,k) + JcrB(jr_, 0   ,j,k)
                ! -- jt_ :: i : Half-Grid, j : Half-Grid -- !
                JcrB(jt_,    2,j,k) = JcrB(jt_,    2,j,k) + JcrB(jt_, 1   ,j,k)
                JcrB(jt_,    3,j,k) = JcrB(jt_,    3,j,k) + JcrB(jt_, 0   ,j,k)
                ! -- jz_ :: i : Full-Grid, j : Half-Grid -- !
                JcrB(jz_,    2,j,k) = 0.d0
                JcrB(jz_,    3,j,k) = JcrB(jz_,    3,j,k) - JcrB(jz_, 1   ,j,k)
                JcrB(jz_,    4,j,k) = JcrB(jz_,    4,j,k) - JcrB(jz_, 0   ,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       !  -- [3-3] Right Side Edge -- ! 
       if ( myRank.eq.PEtot-1 ) then
          !$omp parallel default(none) &
          !$omp shared(JcrB,LIs,LJs) private(j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                ! -- jr_ :: i : Half-Grid, j : Full-Grid -- !
                JcrB(jr_,LIs-1,j,k) = JcrB(jr_,LIs-1,j,k) + JcrB(jr_,LIs  ,j,k)
                JcrB(jr_,LIs-2,j,k) = JcrB(jr_,LIs-2,j,k) + JcrB(jr_,LIs+1,j,k)
                ! -- jt_ :: i : Half-Grid, j : Half-Grid -- !
                JcrB(jt_,LIs-1,j,k) = JcrB(jt_,LIs-1,j,k) + JcrB(jt_,LIs  ,j,k)
                JcrB(jt_,LIs-2,j,k) = JcrB(jt_,LIs-2,j,k) + JcrB(jt_,LIs+1,j,k)
                ! -- jz_ :: i : Full-Grid, j : Half-Grid -- !
                JcrB(jz_,LIs  ,j,k) = 0.d0
                JcrB(jz_,LIs-1,j,k) = JcrB(jz_,LIs-1,j,k) - JcrB(jz_,LIs+1,j,k)
                JcrB(jz_,LIs-2,j,k) = JcrB(jz_,LIs-2,j,k) - JcrB(jz_,LIs+2,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif

       ! ------------------------------------------ !
       ! --- [4]   Boundary-(1) :: (ii) Copying --- !
       ! ------------------------------------------ !
       !  -- [4-1] Exchange Boundary Info. --  !
       do k=1, ns
          call CurrentCommunicate( JcrB(1:3,0:LIs+2,0:LJs+2,k), 3 )
       enddo
       !  -- [4-2] Edge Copying            --  !
       !   - (1) Left Edge -   !
       if ( myRank.eq.0       ) then
          !$omp parallel default(none) &
          !$omp shared(JcrB,LIs,LJs) private(i,j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                JcrB(jr_,  1,j,k) =   JcrB(jr_,    2,j,k)
                JcrB(jt_,  1,j,k) =   JcrB(jt_,    2,j,k)
                JcrB(jz_,  1,j,k) = - JcrB(jz_,    3,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       !   - (2) Right Edge -  !
       if ( myRank.eq.PEtot-1 ) then
          !$omp parallel default(none) &
          !$omp shared(JcrB,LIs,LJs) private(i,j,k)
          !$omp do
          do k=1, ns
             do j=1, LJs
                JcrB(jr_,LIs,j,k) =   JcrB(jr_,LIs-1,j,k)
                JcrB(jt_,LIs,j,k) =   JcrB(jt_,LIs-1,j,k)
                ! JcrB(jz_,LIs+1,j,k) = - JcrB(jz_,LIs-1,j,k)
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
    endif
    return
  end subroutine jcrFoldingBC


  subroutine jcrPeriodicBC( bctype, JcrB )
    use constants , only              : LIs, LJs, ns
    use constants , only              : jr_, jt_, jz_
    use jMPIComMod, only              : CurrentOverlay
    use fMPIComMod, only              : CurrentCommunicate
    implicit none
    integer                          :: i, k
    character(len=10), intent(in)    :: bctype
    double precision , intent(inout) :: JcrB(3,0:LIs+2,0:LJs+2,ns)
    logical                          :: Flag__Boundary1, Flag__Boundary2

    ! ---------------------------- !
    ! --- [1] Flag Preparation --- !
    ! ---------------------------- !
    select case( bctype )
    case('Boundary_1')
       Flag__Boundary1 = .true.
       Flag__Boundary2 = .false.
    case('Boundary_2')
       Flag__Boundary1 = .false.
       Flag__Boundary2 = .true.
    case('Boundary_0')
       Flag__Boundary1 = .true.
       Flag__Boundary2 = .true.
    end select
    
    ! ---------------------------- !
    ! --- [2]   Boundary-(2)   --- !
    ! ---------------------------- !
    if ( Flag__Boundary2 ) then
       !$omp parallel default(none) &
       !$omp shared(JcrB,LIs,LJs) private(i,k)
       !$omp do
       do k=1, ns
          do i=0, LIs+2
             ! -- jr_ :: i : Half-Grid, j : Full-Grid -- !
             JcrB(jr_,i,3,k) =   JcrB(jr_,i,3,k) - JcrB(jr_,i,1,k)
             JcrB(jr_,i,4,k) =   JcrB(jr_,i,4,k) - JcrB(jr_,i,0,k)
             ! -- jt_ :: i : Half-Grid, j : Half-Grid -- !
             JcrB(jt_,i,2,k) =   JcrB(jt_,i,2,k) + JcrB(jt_,i,1,k)
             JcrB(jt_,i,3,k) =   JcrB(jt_,i,3,k) + JcrB(jt_,i,0,k)
             ! -- jz_ :: i : Full-Grid, j : Half-Grid -- !
             JcrB(jz_,i,2,k) =   JcrB(jz_,i,2,k) + JcrB(jz_,i,1,k)
             JcrB(jz_,i,3,k) =   JcrB(jz_,i,3,k) + JcrB(jz_,i,0,k)
             ! -- jr  :: jr(r=0) = 0, jr=anti-Symmetric -- !
             JcrB(jr_,i,2,k) =   0.d0
             JcrB(jr_,i,1,k) = - JcrB(jr_,i,3,k)
             JcrB(jr_,i,0,k) = - JcrB(jr_,i,4,k)
             ! -- jt  :: jt=Symmetric -- !
             JcrB(jt_,i,1,k) = + JcrB(jt_,i,2,k)
             JcrB(jt_,i,0,k) = + JcrB(jt_,i,3,k)
             ! -- jz  :: jz=Symmetric -- !
             JcrB(jz_,i,1,k) =   JcrB(jz_,i,2,k)
             JcrB(jz_,i,0,k) =   JcrB(jz_,i,3,k)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif

    ! ---------------------------- !
    ! --- [3]   Boundary-(1)   --- !
    ! ---------------------------- !
    if ( Flag__Boundary1 ) then
       !  -- [3-1] Overlay Commu. -- ! 
       do k=1, ns
          call CurrentOverlay( JcrB(1:3,0:LIs+2,0:LJs+2,k) )
       enddo
       !  -- [3-2] Exchange Boundary Info. --  !
       do k=1, ns
          call CurrentCommunicate( JcrB(1:3,0:LIs+2,0:LJs+2,k), 3 )
       enddo
    endif
    
    return
  end subroutine jcrPeriodicBC
  
end module jcrBdrcMod
