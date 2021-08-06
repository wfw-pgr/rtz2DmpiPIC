module mxwBdrcMod
contains

  ! =================================== !
  ! === Boundary Condition Handlers === !
  ! =================================== !
  subroutine BoundaryB
    use constants, only : Boundary1__EB, Boundary2__EB
    implicit none
    ! -------------------------------- !
    ! --- [1] Boundary - 2         --- !
    ! -------------------------------- !
    if ( trim(Boundary2__EB).eq.'cWall'    ) call cWallBoundary2_B

    ! -------------------------------- !
    ! --- [2] Boundary - 1         --- !
    ! -------------------------------- !
    if ( trim(Boundary1__EB).eq.'cWall'    ) call cWallBoundary1_B
    if ( trim(Boundary1__EB).eq.'periodic' ) call periodicBoundary1_EB( 'B' )
    return
  end subroutine BoundaryB
  

  subroutine BoundaryE
    use constants, only : Boundary1__EB, Boundary2__EB
    implicit none
    ! -------------------------------- !
    ! --- [1] Boundary - 2         --- !
    ! -------------------------------- !
    if ( trim(Boundary2__EB).eq.'cWall'    ) call cWallBoundary2_E
    
    ! -------------------------------- !
    ! --- [2] Boundary - 1         --- !
    ! -------------------------------- !
    if ( trim(Boundary1__EB).eq.'cWall'    ) call cWallBoundary1_E
    if ( trim(Boundary1__EB).eq.'periodic' ) call periodicBoundary1_EB( 'E' )
    return
  end subroutine BoundaryE


  
  ! =================================== !
  ! ===    Boundary Conditions      === !
  ! =================================== !
  subroutine cWallBoundary1_B
    use variables , only : EBf
    use fMPIComMod, only : FieldCommunicate
    implicit none
    
    ! -------------------------------------------- !
    ! --- [1] Boundary_1 (LR) :: ( Bt:Fixed )  --- !
    ! -------------------------------------------- !
    !  -- [1-1] Communication with Neibour   -- !
    call FieldCommunicate( EBf, 'B' )
    ! !  -- [1-2] Fixed Bt ( d/dt = 0 )        -- !
    ! if ( myRank.eq.0       ) then
    !    !$omp parallel default(none) &
    !    !$omp shared(EBf,EBo,LIs,LJs) private(i,j)
    !    !$omp do
    !    do j=1, LJs
    !       EBf(bt_,  2, j) = EBo(bt_,  2, j)
    !    enddo
    !    !$omp end do
    !    !$omp end parallel
    ! endif
    ! if ( myRank.eq.PEtot-1 ) then
    !    !$omp parallel default(none) &
    !    !$omp shared(EBf,EBo,LIs,LJs) private(i,j)
    !    !$omp do
    !    do j=1, LJs
    !       EBf(bt_,LIs, j) = EBo(bt_,LIs, j)
    !    enddo
    !    !$omp end do
    !    !$omp end parallel
    ! endif
    return
  end subroutine cWallBoundary1_B


  subroutine cWallBoundary2_B
    use constants , only : LIs, LJs
    use constants , only : br_, bt_, bz_
    use variables , only : EBf, EBo
    use fMPIComMod, only : FieldCommunicate
    implicit none
    integer             :: i
    
    ! -------------------------------------------- !
    ! --- [1] Boundary_2 (TB) :: ( Bt:Fixed )  --- !
    ! -------------------------------------------- !
    !  -- [1-1] Bt @wall, r=0 :: Fix           --  !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBo,LIs,LJs) private(i)
    !$omp do
    do i=1, LIs
       EBf(bt_,i,  2) = EBo(bt_,i,  2)
       EBf(bt_,i,LJs) = EBo(bt_,i,LJs)
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine cWallBoundary2_B
  

  subroutine cWallBoundary1_E
    use constants , only : LIs, LJs, myRank, PEtot
    use constants , only : er_, et_, ez_
    use variables , only : EBf
    use fMPIComMod, only : FieldCommunicate
    implicit none
    integer             :: j

    ! ----------------------------------------- !
    ! --- [1] Boundary_1 (LR) :: Et=0 @Wall --- !
    ! ----------------------------------------- !
    ! -- [1-1] Communication with Neibour -- !
    call FieldCommunicate( EBf, 'E' )
    ! -- [1-2]  Left Edge Boundary        -- !
    if ( myRank.eq.0       ) then
       !$omp parallel default(none) &
       !$omp shared(EBf,LJs) private(j)
       !$omp do
       do j=1, LJs
          EBf(er_,  1,j) = - EBf(er_,    2,j)
          EBf(et_,  1,j) = - EBf(et_,    2,j)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    ! -- [1-3] Right Edge Boundary        -- !
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(EBf,LIs,LJs) private(j)
       !$omp do
       do j=1, LJs
          EBf(er_,LIs,j) = - EBf(er_,LIs-1,j)
          EBf(et_,LIs,j) = - EBf(et_,LIs-1,j)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    return
  end subroutine cWallBoundary1_E

  
  subroutine cWallBoundary2_E
    use constants , only : LIs, LJs
    use constants , only : er_, et_, ez_
    use variables , only : rh, rhInv
    use variables , only : EBf
    implicit none
    integer             :: i
    double precision    :: factor1, factor2

    ! ----------------------------------------- !
    ! --- [1] Boundary_2 (T) :: Et=0 @Wall  --- !
    ! ----------------------------------------- !
    ! -- [1-1] r=R0 Conductor  ( Et=0 )   -- !
    factor1 = rh(LJs) * rhInv(LJs-1)
    factor2 = 1.d0
    !$omp parallel default(none) &
    !$omp shared(EBf,factor1,factor2,LIs,LJs) private(i)
    !$omp do
    do i=1, LIs
       EBf(er_,i,LJs) =   0.0d0
       EBf(et_,i,LJs) = - factor1 * EBf(et_,i,LJs-1)
       EBf(ez_,i,LJs) = - factor2 * EBf(ez_,i,LJs-1)
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine cWallBoundary2_E

  
  subroutine periodicBoundary1_EB( EorB )
    use variables , only       : EBf
    use fMPIComMod, only       : FieldCommunicate
    implicit none
    character(1), intent(in)  :: EorB     ! 'E'/'B' !
    
    ! ----------------------------------- !
    ! --- [1] Boundary_1 (LR) ::      --- !
    ! ----------------------------------- !
    !  -- [1-1] Field Exchange      -- !
    call FieldCommunicate( EBf, EorB )
    return
  end subroutine periodicBoundary1_EB

  
end module mxwBdrcMod




  ! subroutine conductingWall_B( boundary )
  !   use constants , only : LIs, LJs
  !   use constants , only : br_, bt_, bz_
  !   use variables , only : EBf, EBo
  !   use fMPIComMod, only : FieldCommunicate
  !   implicit none
  !   integer                       :: i
  !   character(len=10), intent(in) :: boundary
  !   logical                       :: Flag__Boundary1, Flag__Boundary2
    
  !   ! ---------------------------- !
  !   ! --- [1] Flag Preparation --- !
  !   ! ---------------------------- !
  !   select case( boundary )
  !   case('Boundary_1')
  !      Flag__Boundary1 = .true.
  !      Flag__Boundary2 = .false.
  !   case('Boundary_2')
  !      Flag__Boundary1 = .false.
  !      Flag__Boundary2 = .true.
  !   case('Boundary_0')
  !      Flag__Boundary1 = .true.
  !      Flag__Boundary2 = .true.
  !   case default
  !      write(6,*) '[WARNING] No B.C. for E '
  !   end select

  !   ! -------------------------------------------- !
  !   ! --- [2] Boundary_2 (TB) :: ( Bt:Fixed )  --- !
  !   ! -------------------------------------------- !
  !   if ( Flag__Boundary2 ) then
  !      ! -- [2-1] Bt @wall, r=0 :: Fix          -- !
  !      !$omp parallel default(none) &
  !      !$omp shared(EBf,EBo,LIs,LJs) private(i)
  !      !$omp do
  !      do i=1, LIs
  !         EBf(bt_,i,  2) = EBo(bt_,i,  2)
  !         EBf(bt_,i,LJs) = EBo(bt_,i,LJs)
  !      enddo
  !      !$omp end do
  !      !$omp end parallel
  !   endif

  !   ! -------------------------------------------- !
  !   ! --- [3] Boundary_1 (LR) :: ( Bt:Fixed )  --- !
  !   ! -------------------------------------------- !
  !   if ( Flag__Boundary1 ) then
  !      !  -- [3-1] Communication with Neibour   -- !
  !      call FieldCommunicate( EBf, 'B' )
  !      ! !  -- [2-2] Fixed Bt ( d/dt = 0 )        -- !
  !      ! if ( myRank.eq.0       ) then
  !      !    !$omp parallel default(none) &
  !      !    !$omp shared(EBf,EBo,LIs,LJs) private(i,j)
  !      !    !$omp do
  !      !    do j=1, LJs
  !      !       EBf(bt_,  2, j) = EBo(bt_,  2, j)
  !      !    enddo
  !      !    !$omp end do
  !      !    !$omp end parallel
  !      ! endif
  !      ! if ( myRank.eq.PEtot-1 ) then
  !      !    !$omp parallel default(none) &
  !      !    !$omp shared(EBf,EBo,LIs,LJs) private(i,j)
  !      !    !$omp do
  !      !    do j=1, LJs
  !      !       EBf(bt_,LIs, j) = EBo(bt_,LIs, j)
  !      !    enddo
  !      !    !$omp end do
  !      !    !$omp end parallel
  !      ! endif
  !   endif

  !   return
  ! end subroutine ConductingWall_B


  ! subroutine conductingWall_E( boundary )
  !   use constants , only : LIs, LJs, myRank, PEtot
  !   use constants , only : er_, et_, ez_
  !   use variables , only : rh, rhInv
  !   use variables , only : EBf
  !   use fMPIComMod, only : FieldCommunicate
  !   implicit none
  !   integer                       :: i, j
  !   double precision              :: factor1, factor2
  !   character(len=10), intent(in) :: boundary  ! - 'Boundary_1' or 'Boundary_2' - !
  !   logical                       :: Flag__Boundary1, Flag__Boundary2

  !   ! ---------------------------- !
  !   ! --- [1] Flag Preparation --- !
  !   ! ---------------------------- !
  !   select case( boundary )
  !   case('Boundary_1')
  !      Flag__Boundary1 = .true.
  !      Flag__Boundary2 = .false.
  !   case('Boundary_2')
  !      Flag__Boundary1 = .false.
  !      Flag__Boundary2 = .true.
  !   case('Boundary_0')
  !      Flag__Boundary1 = .true.
  !      Flag__Boundary2 = .true.
  !   case default
  !      write(6,*) '[WARNING] No B.C. for E '
  !   end select

  !   ! ----------------------------------------- !
  !   ! --- [2] Boundary_2 (T) :: Et=0 @Wall  --- !
  !   ! ----------------------------------------- !
  !   if ( Flag__Boundary2 ) then
  !      ! -- [2-1] r=R0 Conductor  ( Et=0 )   -- !
  !      factor1 = rh(LJs) * rhInv(LJs-1)
  !      factor2 = 1.d0
  !      !$omp parallel default(none) &
  !      !$omp shared(EBf,factor1,factor2,LIs,LJs) private(i)
  !      !$omp do
  !      do i=1, LIs
  !         EBf(er_,i,LJs) =   0.0d0
  !         EBf(et_,i,LJs) = - factor1 * EBf(et_,i,LJs-1)
  !         EBf(ez_,i,LJs) = - factor2 * EBf(ez_,i,LJs-1)
  !      enddo
  !      !$omp end do
  !      !$omp end parallel
  !   endif

  !   ! ----------------------------------------- !
  !   ! --- [3] Boundary_1 (LR) :: Et=0 @Wall --- !
  !   ! ----------------------------------------- !
  !   if ( Flag__Boundary1 ) then
  !      ! -- [3-1] Communication with Neibour -- !
  !      call FieldCommunicate( EBf, 'E' )
  !      ! -- [3-2]  Left Edge Boundary        -- !
  !      if ( myRank.eq.0       ) then
  !         !$omp parallel default(none) &
  !         !$omp shared(EBf,LIs,LJs) private(j)
  !         !$omp do
  !         do j=1, LJs
  !            EBf(er_,  1,j) = - EBf(er_,    2,j)
  !            EBf(et_,  1,j) = - EBf(et_,    2,j)
  !         enddo
  !         !$omp end do
  !         !$omp end parallel
  !      endif
  !      ! -- [3-3] Right Edge Boundary        -- !
  !      if ( myRank.eq.PEtot-1 ) then
  !         !$omp parallel default(none) &
  !         !$omp shared(EBf,LIs,LJs) private(j)
  !         !$omp do
  !         do j=1, LJs
  !            EBf(er_,LIs,j) = - EBf(er_,LIs-1,j)
  !            EBf(et_,LIs,j) = - EBf(et_,LIs-1,j)
  !         enddo
  !         !$omp end do
  !         !$omp end parallel
  !      endif
  !   endif
    
  !   return
  ! end subroutine conductingWall_E
