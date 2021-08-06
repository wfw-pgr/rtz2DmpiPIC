module fBoundaryMod
contains

  
  subroutine periodic( boundary, EorB )
    use constants , only       : LIs, LJs
    use constants , only       : er_, et_, ez_, br_, bt_, bz_
    use variables , only       : EBf
    use fMPIComMod, only       : FieldCommunicate
    implicit none
    integer                   :: i, d1_, d2_, d3_
    character(10), intent(in) :: boundary ! - 'Boundary_1' , 'Boundary_2' - !
    character( 1), intent(in) :: EorB     ! - 'E', 'B'      - !
    
    ! --- [1] For What Component --- !
    if ( EorB.eq.'E' ) then
       d1_ = er_; d2_ = et_; d3_ = ez_
    endif
    if ( EorB.eq.'B' ) then
       d1_ = br_; d2_ = bt_; d3_ = bz_
    endif
    
    ! --- [2] Left & Right  --- !
    call FieldCommunicate( EBf, 'Periodic_Comm', EorB )
    
    if ( boundary.eq.'Boundary_1' ) then
    !    !$omp parallel default(none) &
    !    !$omp shared(EBf,d1_,d2_,d3_,LIs,LJs) private(j)
    !    !$omp do
    !    do j=1, LJs
    !       EBf(d1_,  1, j) = EBf(d1_,LIs-1,j)
    !       EBf(d2_,  1, j) = EBf(d2_,LIs-1,j)
    !       EBf(d3_,  1, j) = EBf(d3_,LIs-1,j)
    !       EBf(d1_,LIs, j) = EBf(d1_,    2,j)
    !       EBf(d2_,LIs, j) = EBf(d2_,    2,j)
    !       EBf(d3_,LIs, j) = EBf(d3_,    2,j)
    !    enddo
    !    !$omp end do
    !    !$omp end parallel
    endif
    ! --- [3] Top & Bottom --- !
    if ( boundary.eq.'Boundary_2' ) then
       !$omp parallel default(none) &
       !$omp shared(EBf,d1_,d2_,d3_,LIs,LJs) private(i)
       !$omp do
       do i=1, LIs
          EBf(d1_, i,  1) = EBf(d1_, i,LJs-1)
          EBf(d2_, i,  1) = EBf(d2_, i,LJs-1)
          EBf(d3_, i,  1) = EBf(d3_, i,LJs-1)
          EBf(d1_, i,LJs) = EBf(d1_, i,    2)
          EBf(d2_, i,LJs) = EBf(d2_, i,    2)
          EBf(d3_, i,LJs) = EBf(d3_, i,    2)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    return
  end subroutine periodic


  subroutine conductingWall_B( boundary )
    use constants , only : LIs, LJs
    use constants , only : br_, bt_, bz_, et_
    use variables , only : EBf, EBo
    use fMPIComMod, only : FieldCommunicate
    implicit none
    integer                       :: i
    character(len=10), intent(in) :: boundary  ! - 'Boundary_1' or 'Boundary_2' - !

    ! --- [1] Communication with Neibour        --- !
    call FieldCommunicate( EBf, 'Boundary_Comm', 'B' )
    
    ! --- [2] Boundary_1 (LR) :: ( Bt:Fixed )   --- !
    if ( boundary.eq.'Boundary_1' ) then
       ! !$omp parallel default(none) &
       ! !$omp shared(EBf,EBo) private(i,j)
       ! !$omp do
       ! do j=1, LJs
       !    EBf(bt_,  2, j) = EBo(bt_,  2, j)
       !    EBf(bt_,LIs, j) = EBo(bt_,LIs, j)
       ! enddo
       ! !$omp end do
       ! !$omp end parallel
    endif

    ! --- [3] Boundary_2 (TB) :: Bt:Fixed  Bz:Ampere --- !
    if ( boundary.eq.'Boundary_2' ) then
       ! -- [2-1] Bt @wall, r=0 :: Fix -- !
       !$omp parallel default(none) &
       !$omp shared(EBf,EBo,LIs,LJs) private(i)
       !$omp do
       do i=1, LIs
          EBf(bt_,i,  2) = EBo(bt_,i,  2)
          EBf(bt_,i,LJs) = EBo(bt_,i,LJs)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    return
  end subroutine ConductingWall_B


  subroutine conductingWall_E( boundary )
    use constants , only : LIs, LJs, myRank, PEtot
    use constants , only : er_, et_, ez_
    use variables , only : rh, rhInv
    use variables , only : EBf
    use fMPIComMod, only : FieldCommunicate
    implicit none
    integer                       :: i, j
    double precision              :: factor1, factor2
    character(len=10), intent(in) :: boundary  ! - 'Boundary_1' or 'Boundary_2' - !
    
    ! --- [1] Communication with Neibour        --- !
    call FieldCommunicate( EBf, 'Boundary_Comm', 'E' )
    

    ! --- [2] Boundary_1 (LR) :: Et=0 @Wall --- !
    if ( boundary.eq.'Boundary_1' ) then
       if ( myRank.eq.0       ) then
          !$omp parallel default(none) &
          !$omp shared(EBf,LIs,LJs) private(j)
          !$omp do
          do j=1, LJs
             EBf(et_,  1,j) = - EBf(et_,    2,j)
             EBf(er_,  1,j) = - EBf(er_,    2,j)
          enddo
          !$omp end do
          !$omp end parallel
       endif
       if ( myRank.eq.PEtot-1 ) then
          !$omp parallel default(none) &
          !$omp shared(EBf,LIs,LJs) private(j)
          !$omp do
          do j=1, LJs
             EBf(et_,LIs,j) = - EBf(et_,LIs-1,j)
             EBf(er_,LIs,j) = - EBf(er_,LIs-1,j)
          enddo
          !$omp end do
          !$omp end parallel
       endif
    endif
    ! --- [3] Boundary_2 (T) :: Et=0 @Wall --- !
    if ( boundary.eq.'Boundary_2' ) then
       factor1 = rh(LJs) * rhInv(LJs-1)
       factor2 = 1.d0
       !$omp parallel default(none) &
       !$omp shared(EBf,factor1,factor2,LIs,LJs) private(i)
       !$omp do
       do i=1, LIs
          EBf(er_,i,LJs) = 0.0d0
          EBf(et_,i,LJs) = - factor1 * EBf(et_,i,LJs-1)
          EBf(ez_,i,LJs) = - factor2 * EBf(ez_,i,LJs-1)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    return
  end subroutine conductingWall_E
  
  
  subroutine conductingWall_divB( EBs )
    use constants, only              : LIs, LJs, dr, dz, myRank, PEtot
    use constants, only              : br_, bt_, bz_
    use variables, only              : EBo, rf , rh, rfInv, rhInv
    implicit none
    integer                         :: i, j
    double precision                :: factor1, factor2, factor3
    double precision, intent(inout) :: EBs(6,LIs,LJs)
    
    ! --- [1] divB = 0 Boundary Conditions      --- !
    !  -- [1-1] z=-Z0,+Z0  Bz :: divB=0         --  !
    !  -- [1-2] r= 0       Br :: anti-Symmetric --  !
    !  -- [1-3] r=R0       Br :: divB=0         --  !
    ! --------------------------------------------- !
    factor1 = dz / dr
    factor2 = rh(LJs-1)*rhInv(LJs)
    factor3 = rf(LJs  )*rhInv(LJs) * dr / dz
    !  -- [1-1] z=-Z0,+Z0  Bz :: divB=0         --  !
    if ( myRank.eq.0       ) then
       !$omp parallel default(none) &
       !$omp shared(EBs,rh,rfInv,factor1,LIs,LJs) private(j)
       !$omp do
       do j=3, LJs-1 ! -- [1-1] ( Left & Right ) : Except for corners :: (1,LJs), (LIs,LJs) -- !
          EBs(bz_,1,j) = EBs(bz_,2,j) + factor1*rfInv(j)*( rh(j)*EBs(br_,2,j) - rh(j-1)*EBs(br_,2,j-1) )
       enddo
       !$omp end do
       !$omp end parallel
    endif
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(EBs,rh,rfInv,factor1,LIs,LJs) private(j)
       !$omp do
       do j=3, LJs-1 ! -- [1-1] ( Left & Right ) : Except for corners :: (1,LJs), (LIs,LJs) -- !
          EBs(bz_,LIs,j) = EBs(bz_,LIs-1,j) - factor1*rfInv(j)*( rh(j)*EBs(br_,LIs,j) - rh(j-1)*EBs(br_,LIs,j-1) )
       enddo
       !$omp end do
       !$omp end parallel
    endif
    !  -- [1-2] r= 0       Br :: anti-Symmetric --  !
    !  -- [1-3] r=R0       Br :: divB=0         --  !
    !$omp parallel default(none) &
    !$omp shared(EBs,EBo,factor2,factor3,LIs,LJs) private(i)
    !$omp do
    do i=2, LIs   ! -- [1-2] ( r=0  ) -- !  !!! sign:[+] because of ( rh(1)=-dr, rh(2)=+dr ) => rBr=[asymm] !!!
       EBs(br_,i,1  ) = + EBs(br_,i,2)
       EBo(br_,i,1  ) = + EBo(br_,i,2)
       EBs(br_,i,LJs) = factor2*EBs(br_,i,LJs-1) - factor3*( EBs(bz_,i,LJs) - EBs(bz_,i-1,LJs) )
    enddo
    !$omp end do
    !$omp end parallel
    
    ! --- [2] Corners --- !
    if ( myRank.eq.0       ) then
       EBs(bz_,  1,  2) = EBs(bz_,2,2) + 4.d0*factor1*EBs(br_,2,2)
       EBs(bz_,  1,LJs) =                        EBs(bz_, 2,LJs  )
       EBs(br_,  2,LJs) = rh(LJs-1)*rhInv(LJs) * EBs(br_, 2,LJs-1)
    endif
    if ( myRank.eq.PEtot-1 ) then
       EBs(bz_,LIs,  2) = EBs(bz_,LIs-1,2) - 4.d0*factor1*EBs(br_,LIs,2)
       EBs(br_,LIs,LJs) = rh(LJs-1)*rhInv(LJs) * EBs(br_,LIs  ,LJs-1)
       EBs(bz_,LIs,LJs) =                        EBs(bz_,LIs-1,LJs  )
    endif
    return
  end subroutine conductingWall_divB


  subroutine conductingWall_E_parallel
    use constants, only : LIs, LJs, myRank, PEtot
    use constants, only : er_, et_, ez_
    use variables, only : rh , rhInv, EBf
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2, tangent
    
    ! --- [1] Boundary Condition for EBf Before Relocate --- !
    factor1 = - rh(    2) * rhInv(  1)
    factor2 = - rh(LJs-1) * rhInv(LJs)
    tangent = - 1.d0
    if ( myRank.eq.0       ) then
       !$omp parallel default(none) &
       !$omp shared(EBf,factor1,factor2,tangent,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs
          EBf(er_,  1,j) = tangent * EBf(er_,    2,j)
          EBf(et_,  1,j) = tangent * EBf(et_,    2,j)
       enddo ! -- z-Boundary -- !
       !$omp end do
       !$omp end parallel
    endif
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(EBf,factor1,factor2,tangent,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs
          EBf(er_,LIs,j) = tangent * EBf(er_,LIs-1,j)
          EBf(et_,LIs,j) = tangent * EBf(et_,LIs-1,j)
       enddo ! -- z-Boundary -- !
       !$omp end do
       !$omp end parallel
    endif
    !$omp parallel default(none) &
    !$omp shared(EBf,factor1,factor2,tangent,LIs,LJs) private(i,j)
    !$omp do
    do i=1, LIs
       EBf(et_,i,  1) = factor1 * EBf(et_,i,    2)
       EBf(ez_,i,  1) = factor1 * EBf(ez_,i,    2)
       EBf(et_,i,LJs) = factor2 * EBf(et_,i,LJs-1)
       EBf(ez_,i,LJs) = tangent * EBf(ez_,i,LJs-1)
    enddo ! -- r-Boundary -- !
    !$omp end do
    !$omp end parallel
    return
  end subroutine conductingWall_E_parallel


  subroutine conductingWall_EBr
    use constants   , only       : LIs, LJs, rMin, myRank, PEtot
    use constants   , only       : br_, bt_, bz_, er_, et_, ez_
    use variables   , only       : EBr, EBf
    implicit none
    integer                     :: i, j
    double precision            :: factor
    double precision, parameter :: normal  = -1.d0
    double precision, parameter :: tangent = +1.d0
    
    ! --- [1] Boundary Condition for EBr --- !
    !  -- [1-1] Bottom j=1 ( r=0 )  --  !
    if ( rMin.eq.0.d0 ) then
       !$omp parallel default(none) &
       !$omp shared(EBr,EBf,LIs,LJs) private(i)
       !$omp do
       do i=2, LIs
          ! EBr(et_,i,2) = EBf(et_,i,2)
          EBr(ez_,i,2) = EBf(ez_,i,2)
       enddo
       !$omp end do
       !$omp do
       do i=1, LIs
          EBr(er_,i,1) = normal  * EBr(er_,i,3)
          EBr(et_,i,1) = tangent * EBr(et_,i,3)
          EBr(ez_,i,1) = tangent * EBr(ez_,i,3)
          EBr(br_,i,1) = normal  * EBr(br_,i,3)
          EBr(bt_,i,1) = tangent * EBr(bt_,i,3)
          EBr(bz_,i,1) = tangent * EBr(bz_,i,3)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    !  -- [1-2] i=1, LIs+1, LJs+1 --  !
    !   - j=LJs+1 <=> j=LJs-1 ( Folding or anti-Symmetric )
    if ( myRank.eq.0 ) then
       !$omp parallel default(none) &
       !$omp shared(EBr,factor,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs
          EBr(er_,   1,j) = EBr(er_,   3,j)
          EBr(et_,   1,j) = EBr(et_,   3,j)
          EBr(ez_,   1,j) = EBr(ez_,   3,j)
          EBr(br_,   1,j) = EBr(br_,   3,j)
          EBr(bt_,   1,j) = EBr(bt_,   3,j)
          EBr(bz_,   1,j) = EBr(bz_,   3,j)
       enddo   !   --- i=1    B.C. --- !
       !$omp end do
       !$omp end parallel
    endif
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(EBr,factor,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs
          EBr(er_,LIs+1,j) = EBr(er_,LIs-1,j)
          EBr(et_,LIs+1,j) = EBr(et_,LIs-1,j)
          EBr(ez_,LIs+1,j) = EBr(ez_,LIs-1,j)
          EBr(br_,LIs+1,j) = EBr(br_,LIs-1,j)
          EBr(bt_,LIs+1,j) = EBr(bt_,LIs-1,j)
          EBr(bz_,LIs+1,j) = EBr(bz_,LIs-1,j)
       enddo   !   --- i=LIs+1 B.C. --- !
       !$omp end do
       !$omp end parallel
    endif
    !$omp parallel default(none) &
    !$omp shared(EBr,factor,LIs,LJs) private(i,j)
    !$omp do
    do i=1, LIs
       EBr(er_,i,LJs+1) = EBr(er_,i,LJs-1)
       EBr(et_,i,LJs+1) = EBr(et_,i,LJs-1)
       EBr(ez_,i,LJs+1) = EBr(ez_,i,LJs-1)
       EBr(br_,i,LJs+1) = EBr(br_,i,LJs-1)
       EBr(bt_,i,LJs+1) = EBr(bt_,i,LJs-1)
       EBr(bz_,i,LJs+1) = EBr(bz_,i,LJs-1)
    enddo  !   --- j=LJs+1 B.C. --- !
    !$omp end do
    !$omp end parallel

    return
  end subroutine conductingWall_EBr


  subroutine ptFoldingBoundary( sumup )
    use constants,             only  : LIs, LJs, ns, myRank, PEtot
    use fMPIComMod,            only  : FluidCommunicate
    implicit none
    integer                         :: i, j, k
    double precision, intent(inout) :: sumup(LIs,LJs,ns)
    
    !  --- [1] Boundary_1 Folding ---  !
    call FluidCommunicate( sumup(1:LIs,1:LJs,1:ns), ns, 'EdgeExchange1' )
    !$omp parallel default(none) &
    !$omp shared(sumup,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          sumup(    2,j,k)  = sumup(    2,j,k) + sumup(  1,j,k)
          sumup(LIs-1,j,k)  = sumup(LIs-1,j,k) + sumup(LIs,j,k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    !  --- [2] Boundary_2 Folding ---  !
    !$omp parallel default(none) &
    !$omp shared(sumup,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do i=2, LIs-1
          sumup(i,    2,k)  = sumup(i,    2,k) + sumup(i,  1,k)
          sumup(i,LJs-1,k)  = sumup(i,LJs-1,k) + sumup(i,LJs,k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    !  --- [3] Boundary_1 copying ---  !
    !   -- [3-1] Set Neibour Data --   !
    call FluidCommunicate( sumup(1:LIs,1:LJs,1:ns), ns, 'EdgeExchange1' )
    !   -- [3-2] Boundary Copy    --   !
    if ( myRank.eq.0       ) then
       !$omp parallel default(none) &
       !$omp shared(sumup,LIs,LJs) private(i,j,k)
       !$omp do
       do k=1, ns
          do j=1, LJs
             sumup(  1,j,k) = sumup(    2,j,k)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(sumup,LIs,LJs) private(i,j,k)
       !$omp do
       do k=1, ns
          do j=1, LJs
             sumup(LIs,j,k) = sumup(LIs-1,j,k)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif
    
    !  --- [4] Boundary_2 copying ---  !
    !$omp parallel default(none) &
    !$omp shared(sumup,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do i=1, LIs
          sumup(i,  1,k)    = sumup(i,    2,k)
          sumup(i,LJs,k)    = sumup(i,LJs-1,k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine ptFoldingBoundary

  
  subroutine FoldingCopy_CC( FCC )
    use constants, only              : LIs, LJs
    implicit none
    integer                         :: i, j
    double precision, intent(inout) :: FCC(LIs,LJs)

    !$omp parallel default(none) &
    !$omp shared(FCC,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       FCC(  1,j) = FCC(    2,j)
       FCC(LIs,j) = FCC(LIs-1,j)
    enddo
    !$omp end do
    !$omp do
    do i=1, LIs
       FCC(i,  1) = FCC(i,    2)
       FCC(i,LJs) = FCC(i,LJs-1)
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine FoldingCopy_CC


end module fBoundaryMod
