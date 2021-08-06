module maxwellMod
contains
  
  subroutine Bfield
    use constants , only : LIs, LJs, dtdr, dtdz, rMin
    use constants , only : br_, bt_, bz_, er_, et_, ez_
    use variables , only : EBf, EBo, rh, rfinv
    use mxwBdrcMod, only : BoundaryB
    implicit none
    integer             :: i, j
    double precision    :: Factor_B1, Factor_B2, Factor_B3

    ! ------------------------------ !
    ! --- [1] Store OLD B-Field  --- !
    ! ------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(EBo,EBf,LIs,LJs) private(i,j)
    !$omp do 
    do j=1, LJs
       do i=1, LIs
          EBo(br_,i,j) = EBf(br_,i,j)
          EBo(bt_,i,j) = EBf(bt_,i,j)
          EBo(bz_,i,j) = EBf(bz_,i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ------------------------------ !
    ! --- [2]  FDTD ( Faraday )  --- !
    ! ------------------------------ !
    !  -- [2-1] Parameters        -- !
    Factor_B1 =          dtdr
    Factor_B2 =          dtdz
    Factor_B3 = - 4.d0 * dtdr
    !  -- [2-2] Main Region       -- !
    !$omp parallel default(none) &
    !$omp shared(rh,rfinv,EBf,Factor_B1,Factor_B2,LIs,LJs) private(i,j)
    !$omp do
    do j=2, LJs
       do i=2, LIs
          EBf(br_,i,j)  = EBf(br_,i,j) + Factor_B2*(       EBf(et_,i,j) -         EBf(et_,i-1,j  ) )
          EBf(bt_,i,j)  = EBf(bt_,i,j) + Factor_B1*(       EBf(ez_,i,j) -         EBf(ez_,i  ,j-1) ) &
               &                       - Factor_B2*(       EBf(er_,i,j) -         EBf(er_,i-1,j  ) )
          EBf(bz_,i,j)  = EBf(bz_,i,j) - Factor_B1*( rh(j)*EBf(et_,i,j) - rh(j-1)*EBf(et_,i  ,j-1) )*rfinv(j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [2-3] Bz on Axis ( r=0: j=2 )  -- !
    if ( rMin.eq.0.d0 ) then
       !$omp parallel default(none) &
       !$omp shared(EBf,EBo,Factor_B3,LIs,LJs) private(i)
       !$omp do
       do i=2, LIs
          EBf(bz_,i,2) = EBo(bz_,i,2) + Factor_B3 * EBf(et_,i,2)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    ! ------- E-Field Information is NOT sufficient, But OK  ------- !
    ! - OverWrite :: EBf(bt_,:,2), EBf(bt_,  :,LJs)                  !
    !             :: EBf(bt_,2,:), EBf(bt_,LIs,  :)                  !
    !             :: EBf(bz_,:,2)                                    !
    ! -------------------------------------------------------------- !
    
    ! ------------------------------ !
    ! --- [3] Boundary Condition --- !
    ! ------------------------------ !
    call BoundaryB

    return
  end subroutine Bfield

  
  subroutine Efield
    use constants , only : LIs, LJs, myRank, PEtot
    use constants , only : dt, dtdr, dtdz, vAlfe
    use constants , only : br_, bt_, bz_, er_, et_, ez_
    use variables , only : EBf, Jcr, EBo, rf, rhinv
    use mxwBdrcMod, only : BoundaryE
    implicit none
    integer             :: i, j
    double precision    :: Factor_E1, Factor_E2, Factor_E3
    
    ! ------------------------------ !
    ! --- [1] Store OLD E-Field  --- !
    ! ------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(EBo,EBf,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBo(er_,i,j) = EBf(er_,i,j)
          EBo(et_,i,j) = EBf(et_,i,j)
          EBo(ez_,i,j) = EBf(ez_,i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ------------------------------ !
    ! --- [2]  FDTD ( Faraday )  --- !
    ! ------------------------------ !
    !  -- [2-1] Parameter Setting -- !
    Factor_E1 = dtdr
    Factor_E2 = dtdz
    Factor_E3 = dt / vAlfe**2
    ! Factor_E3 = 0.d0
    !  -- [2-2] Main Region       -- !
    !$omp parallel default(none) &
    !$omp shared(rf,rhinv,EBf,Jcr,Factor_E1,Factor_E2,Factor_E3,LIs,LJs) private(i,j)
    !$omp do
    do j=2, LJs-1
       do i=2, LIs-1
          EBf(er_,i,j) = EBf(er_,i,j) - Factor_E2*(         EBf(bt_,i+1,j) -       EBf(bt_,i,j) ) &
               &                      - Factor_E3*Jcr(1,i,j,0)
          EBf(et_,i,j) = EBf(et_,i,j) + Factor_E2*(         EBf(br_,i+1,j) -       EBf(br_,i,j) ) &
               &                      - Factor_E1*(         EBf(bz_,i,j+1) -       EBf(bz_,i,j) ) &
               &                      - Factor_E3*Jcr(2,i,j,0)
          EBf(ez_,i,j) = EBf(ez_,i,j) + Factor_E1*( rf(j+1)*EBf(bt_,i,j+1) - rf(j)*EBf(bt_,i,j) )*rhinv(j) &
               &                      - Factor_E3*Jcr(3,i,j,0)
       enddo
    enddo
    !$omp end do
    !  -- [2-3] Er on Axis ( r=R0: j=LJs )  -- !
    !$omp do
    do i=2, LIs-1
       EBf(er_,i,LJs) = EBf(er_,i,LJs)  - Factor_E2*( EBf(bt_,i+1,LJs) - EBf(bt_,i,LJs) ) &
            &                           - Factor_E3*Jcr(1,i,LJs,0)
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [2-4] Ez( z=Z0:i=LIs )  -- !
    if ( myRank.eq.PEtot-1 ) then
       !$omp parallel default(none) &
       !$omp shared(EBf,Factor_E1,Factor_E3,rf,rhInv,Jcr,LIs,LJs) private(j)
       !$omp do
       do j=2, LJs-1
          EBf(ez_,LIs,j) = EBf(ez_,LIs,j)  + Factor_E1*rhInv(j)*(   rf(j+1)*EBf(bt_,LIs,j+1)   &
               &                                                  - rf(j  )*EBf(bt_,LIs,j  ) ) &
               &                           - Factor_E3*Jcr(3,LIs,j,0)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    
    ! -------------------------------- !
    ! --- [3]  Boundary Condition  --- !
    ! -------------------------------- !
    call BoundaryE

    return
  end subroutine Efield


end module maxwellMod
