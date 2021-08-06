module pcwBdrcMod
contains

  subroutine conductingWall_divB( EBs )
    use constants, only              : LIs, LJs, dr, dz, myRank, PEtot
    use constants, only              : br_, bt_, bz_
    use variables, only              : EBo, rf , rh, rfInv, rhInv
    implicit none
    integer                         :: i, j
    double precision                :: factor1, factor2, factor3
    double precision, intent(inout) :: EBs(6,LIs,LJs)

    ! --------------------------------------------- !
    ! --- [1] divB = 0 Boundary1   (LR)         --- !
    ! --------------------------------------------- !
    factor1 = dz / dr
    factor2 = rh(LJs-1)*rhInv(LJs)
    factor3 = rf(LJs  )*rhInv(LJs) * dr / dz
    !  -- [1-1] z=-Z0  Bz :: divB=0             --  !
    if ( myRank.eq.0       ) then
       !$omp parallel default(none) &
       !$omp shared(EBs,rh,rfInv,factor1,LIs,LJs) private(j)
       !$omp do
       do j=3, LJs-1 ! - : Except for corners :: (1,LJs), (LIs,LJs) - !
          EBs(bz_,  1,j) = EBs(bz_,    2,j) + factor1*rfInv(j)*( rh(j)*EBs(br_,  2,j) - rh(j-1)*EBs(br_,  2,j-1) )
       enddo
       !$omp end do
       !$omp end parallel
    endif
    !  -- [1-2] z=+Z0  Bz :: divB=0             --  !
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(EBs,rh,rfInv,factor1,LIs,LJs) private(j)
       !$omp do
       do j=3, LJs-1 ! - : Except for corners :: (1,LJs), (LIs,LJs) - !
          EBs(bz_,LIs,j) = EBs(bz_,LIs-1,j) - factor1*rfInv(j)*( rh(j)*EBs(br_,LIs,j) - rh(j-1)*EBs(br_,LIs,j-1) )
       enddo
       !$omp end do
       !$omp end parallel
    endif
    
    ! --------------------------------------------- !
    ! --- [2] divB = 0 Boundary2   (BT)         --- !
    ! --------------------------------------------- !
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
    EBs(br_,1,  1) = EBs(br_,1,    2)
    EBo(br_,1,  1) = EBo(br_,1,    2)
    EBs(br_,1,LJs) = EBs(br_,1,LJs-1)
    
    ! --------------------------------------------- !
    ! --- [3] Corners                           --- !
    ! --------------------------------------------- !
    if ( myRank.eq.0       ) then
       EBs(bz_,  1,  2) = EBs(bz_,    2,2) + 4.d0*factor1*EBs(br_,    2,    2)
       EBs(bz_,  1,LJs) =                                 EBs(bz_,    2,LJs  )
       EBs(br_,  2,LJs) = rh(LJs-1)*rhInv(LJs)           *EBs(br_,    2,LJs-1)
    endif
    if ( myRank.eq.PEtot-1 ) then
       EBs(bz_,LIs,  2) = EBs(bz_,LIs-1,2) - 4.d0*factor1*EBs(br_,LIs  ,    2)
       EBs(br_,LIs,LJs) = rh(LJs-1)*rhInv(LJs)           *EBs(br_,LIs  ,LJs-1)
       EBs(bz_,LIs,LJs) =                                 EBs(bz_,LIs-1,LJs  )
    endif
    return
  end subroutine conductingWall_divB


  subroutine conductingWall_Epara
    use constants, only : LIs, LJs, myRank, PEtot
    use constants, only : er_, et_, ez_
    use variables, only : rh , rhInv, EBf
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2, tangent

    ! ---------------------------------------------- !
    ! --- [0]  Preparation ( factor )            --- !
    ! ---------------------------------------------- !
    factor1 = - rh(    2) * rhInv(  1)
    factor2 = - rh(LJs-1) * rhInv(LJs)
    tangent = - 1.d0

    ! ---------------------------------------------- !
    ! --- [1] Parallel E :: Ep = 0 ( Boundary1 ) --- !
    ! ---------------------------------------------- !
    !  -- [1-1]  Left Edge -- !
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
    !  -- [1-2] Right Edge -- !
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

    ! ---------------------------------------------- !
    ! --- [2] Parallel E :: Ep = 0 ( Boundary2 ) --- !
    ! ---------------------------------------------- !
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
  end subroutine conductingWall_Epara


  subroutine conductingWall_EBr
    use constants , only         : LIs, LJs, rMin, myRank, PEtot
    use constants , only         : br_, bt_, bz_, er_, et_, ez_
    use variables , only         : EBr, EBf
    use rMPIComMod, only         : RelocatedExchange
    implicit none
    integer                     :: i, j
    double precision, parameter :: normal  = -1.d0
    double precision, parameter :: tangent = +1.d0

    ! ---------------------------------------- !
    ! --- [1]  Bottom & Top ( Boundary2 )  --- !
    ! ---------------------------------------- !
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
       do i=2, LIs
          EBr(er_,i,    1) =  normal * EBr(er_,i,    3)
          EBr(et_,i,    1) = tangent * EBr(et_,i,    3)
          EBr(ez_,i,    1) = tangent * EBr(ez_,i,    3)
          EBr(br_,i,    1) =  normal * EBr(br_,i,    3)
          EBr(bt_,i,    1) = tangent * EBr(bt_,i,    3)
          EBr(bz_,i,    1) = tangent * EBr(bz_,i,    3)
          EBr(er_,i,LJs+1) =           EBr(er_,i,LJs-1)
          EBr(et_,i,LJs+1) =           EBr(et_,i,LJs-1)
          EBr(ez_,i,LJs+1) =           EBr(ez_,i,LJs-1)
          EBr(br_,i,LJs+1) =           EBr(br_,i,LJs-1)
          EBr(bt_,i,LJs+1) =           EBr(bt_,i,LJs-1)
          EBr(bz_,i,LJs+1) =           EBr(bz_,i,LJs-1)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    
    ! ---------------------------------------- !
    ! --- [2]  Left & Right ( Boundary1 )  --- !
    ! ---------------------------------------- !
    !  -- [2-1]  Left Edge -- !
    if ( myRank.eq.0 ) then
       !$omp parallel default(none) &
       !$omp shared(EBr,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs+1
          EBr(er_,    1,j) = EBr(er_,    3,j)
          EBr(et_,    1,j) = EBr(et_,    3,j)
          EBr(ez_,    1,j) = EBr(ez_,    3,j)
          EBr(br_,    1,j) = EBr(br_,    3,j)
          EBr(bt_,    1,j) = EBr(bt_,    3,j)
          EBr(bz_,    1,j) = EBr(bz_,    3,j)
       enddo   !   --- i=1    B.C. --- !
       !$omp end do
       !$omp end parallel
    endif
    !  -- [2-2] Right Edge -- !
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(EBr,LIs,LJs) private(i,j)
       !$omp do
       do j=1, LJs+1
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
    !  -- [2-3] Exchange -- !
    call RelocatedExchange( EBr )

    return
  end subroutine conductingWall_EBr

end module pcwBdrcMod
