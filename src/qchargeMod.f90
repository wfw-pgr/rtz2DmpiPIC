module qchargeMod
contains
  
  subroutine charge
    use constants,    only  : LI, LJ, LIs, LJs, ns, jobDir, myRank, PEtot, k_Ecrct
    use constants,    only  : q, drInv, dzInv, vAlfe
    use constants,    only  : br_, bt_, bz_, er_, et_, ez_
    use variables,    only  : EBf, kstep, rh, rf, rhInv, rfInv, w_p
    use momentsMod,   only  : nCountUp_CC
    use pcwBdrcMod,   only  : conductingWall_divB
    use fMPIComMod,   only  : FieldCommunicate
    implicit none
    include 'mpif.h'
    integer                :: i, j, ierr
    double precision       :: divE, divB, valfe2Inv
    double precision       :: cResave, cResMSE, cRessgm, sum_cRes, sum_cRes2
    double precision       :: divEave, divEMSE, divEsgm, sum_divE, sum_divE2
    double precision       :: divBave, divBMSE, divBsgm, sum_divB, sum_divB2
    double precision       :: rhoC(LIs,LJs), cRes(LIs,LJs), EBs(6,LIs,LJs)
    character(100)         :: FileName
    double precision, save :: cResMSE_0        =  0.d0
    double precision, save :: divBMSE_0        =  0.d0
    double precision, save :: divEMSE_0        =  0.d0
    logical, parameter     :: DensityFiltering = .false.

    ! ------------------------- !
    ! --- [1] calculate w_p --- !
    ! ------------------------- !
    !  -- [1-1] count up          -- !
    call nCountUp_CC( w_p )
    !  -- [1-2] Digital Filtering -- !
    if ( DensityFiltering ) call nTimesFiltering
    !  -- [1-3] rho = e(ni-ne)    -- !
    !$omp parallel default(none) &
    !$omp shared(q,rhoC,w_p,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          w_p(i,j,1) =       q(1) * w_p(i,j,1)
          w_p(i,j,2) =       q(2) * w_p(i,j,2)
          rhoC(i,j)  = w_p(i,j,1) + w_p(i,j,2)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [1-4] Copy :: EBs       -- !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBs(er_,i,j) = EBf(er_,i,j)
          EBs(et_,i,j) = EBf(et_,i,j)
          EBs(ez_,i,j) = EBf(ez_,i,j)
          EBs(br_,i,j) = EBf(br_,i,j)
          EBs(bt_,i,j) = EBf(bt_,i,j)
          EBs(bz_,i,j) = EBf(bz_,i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! ---------------------------------------- !
    ! --- [2] Residual of ( div E - RhoC ) --- !
    ! ---------------------------------------- !
    !  -- [2-1] Sum up divE / cRes -- !
    sum_divE  = 0.d0
    sum_divE2 = 0.d0
    sum_cRes  = 0.d0
    sum_cRes2 = 0.d0
    valfe2Inv = 1.d0 / ( vAlfe**2 )
    !$omp parallel default(none) &
    !$omp shared(cRes,EBs,drInv,dzInv,valfe2Inv,rhoC,rf,rhInv,LIs,LJs) &
    !$omp shared(sum_divE,sum_divE2,sum_cRes,sum_cRes2) &
    !$omp private(i,j,divE)
    !$omp do reduction(+:sum_divE,sum_divE2,sum_cRes,sum_cRes2)
    do j=4, LJs-3
       do i=2, LIs-1
          divE      =    ( rf(j+1)*EBs(er_,i,j+1) - rf(j)*EBs(er_,i,j) )*drInv*rhInv(j) &
               &       + (         EBs(ez_,i+1,j) -       EBs(ez_,i,j) )*dzInv
          cRes(i,j) = divE    - valfe2Inv * rhoC(i,j)
          sum_divE  = sum_divE  + divE
          sum_divE2 = sum_divE2 + divE**2
          sum_cRes  = sum_cRes  + cRes(i,j)
          sum_cRes2 = sum_cRes2 + cRes(i,j)**2
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [2-2] MPI Reduction      -- !
    call MPI_AllReduce( sum_divE , divEave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( sum_divE2, divEMSE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( sum_cRes , cResave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( sum_cRes2, cResMSE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    !  -- [2-3] Average & MSE      -- !
    divEave   = divEave / dble( (LI-2)*(LJ-4) )
    divEMSE   = divEMSE / dble( (LI-2)*(LJ-4) )
    cResave   = cResave / dble( (LI-2)*(LJ-4) )
    cResMSE   = cResMSE / dble( (LI-2)*(LJ-4) )
    !  -- [2-4] Standard Deviation -- !
    divEsgm   = sqrt( divEMSE - divEave**2 )
    cRessgm   = sqrt( cResMSE - cResave**2 )
    if ( kstep.eq.k_Ecrct ) then
       cResMSE_0 = cResMSE
       divEMSE_0 = divEMSE
    endif
    
    ! ------------------- !
    ! --- [3]  div B  --- !
    ! ------------------- !
    !  -- [3-1] Sum up divB        -- !
    call conductingWall_divB( EBs )
    call FieldCommunicate   ( EBs, 'B' )
    sum_divB  = 0.d0
    sum_divB2 = 0.d0
    !$omp parallel default(none) &
    !$omp shared (LIs,LJs,rh,rfInv,drInv,dzInv,EBs) &
    !$omp shared (sum_divB,sum_divB2) &
    !$omp private(i,j,divB)
    !$omp do reduction(+:sum_divB,sum_divB2)
    do j=3, LJs-1
       do i=2, LIs-1
          divB      =   rfInv(j) * drInv * ( rh(j)*EBs(br_,i,j) - rh(j-1)*EBs(br_,i  ,j-1) ) &
               &      +            dzInv * (       EBs(bz_,i,j) -         EBs(bz_,i-1,j  ) )
          sum_divB  = sum_divB  + divB
          sum_divB2 = sum_divB2 + divB**2
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [3-2] MPI Reduction      -- !
    call MPI_AllReduce( sum_divB , divBave, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( sum_divB2, divBMSE, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    !  -- [3-3] Average & MSE      -- !
    divBave   = divBave / dble( (LI-3)*(LJ-3) )
    divBMSE   = divBMSE / dble( (LI-3)*(LJ-3) )
    !  -- [3-4] Standard Deviation -- !
    divBsgm   = sqrt( divBMSE - divBave**2 )
    if ( kstep.eq.k_Ecrct ) divBMSE_0 = divBMSE


    ! -------------------------------- !
    ! --- [4] Summary of Diagnosis --- !
    ! -------------------------------- !
    if ( myRank.eq.0 ) then
       ! --     Display    -- !
       write(6,*) ' ----------------  charge  ------------------ '
       write(6,*) ' Ecorrect :: charge residual ===:   ', cResave
       write(6,*) '             MSE             ===:   ', cResMSE
       write(6,*) '             MSE  (Deviated) ===:   ', cResMSE - cResMSE_0
       write(6,*) '             standard dev.   ===:   ', cRessgm
       write(6,*)
       write(6,*) ' ----------------   divE   ------------------ '
       write(6,*) ' divE     :: divE (Averaged) ===:   ', divEave
       write(6,*) '          :: MSE             ===:   ', divEMSE
       write(6,*) '             MSE  (Deviated) ===:   ', divEMSE - divEMSE_0
       write(6,*) '             standard dev.   ===:   ', divEsgm
       write(6,*)
       write(6,*) ' ----------------   divB   ------------------ '
       write(6,*) ' divB     :: divB (Averaged) ===:   ', divBave
       write(6,*) '          :: MSE             ===:   ', divBMSE
       write(6,*) '             MSE  (Deviated) ===:   ', divBMSE - divBMSE_0
       write(6,*) '             standard dev.   ===:   ', divBsgm
       write(6,*)
       ! -- File Write out -- !
       FileName = trim(jobDir) // 'dat/' // 'charge_residual.dat'
       open( 50, file=trim(FileName), form='formatted', position='append' )
       write(50,'((i10,1x),11(e15.7,1x))') kstep, cResave, cResMSE, cResMSE, cResMSE-cResMSE_0, cRessgm, &
            &                                     divEave, divEMSE, divEsgm, divBave, divBMSE, divBsgm
       close(50)
    endif

    return
  end subroutine charge

  ! ---- NOT in operation --- !
  subroutine Ecorrect
    use constants      , only : LIs, LJs, dr, dz, drInv, dzInv, rmin
    use constants      , only : er_, et_, ez_
    use variables      , only : EBf, rf, rh
    use myBiCGSTAB     , only : BiCGSTABCyl2D
    use lInterpMod     , only : LinearInterpFWD, LinearInterpRET
    implicit none
    integer                  :: i, j
    double precision         :: drLIs, factor1, factor2, rqm, rqp
    double precision         :: rLIs(LJs)
    double precision         :: rhs(LIs,LJs), rhsLIs(LIs,LJs), dphi(LIs,LJs), dphiLIs(LIs,LJs)

    ! ---  [1]  preparation   --- !
    dphi(:,:)   = 0.d0
    dphiLIs(:,:) = 0.d0
    do j=2, LJs-1
       do i=2, LIs-1
          ! rhs(i,j) = cRes(i,j)
       enddo
    enddo
    do i=2, LIs-1
       rhs(i, 1) = 0.d0
       rhs(i,LJs) = 0.d0
    enddo
    do j=1, LJs
       rhs( 1,j) = 0.d0
       rhs(LIs,j) = 0.d0
    enddo

    ! --- [2] rmin = 0.d0 case --- !
    if ( rmin.eq.0.d0 ) then
       ! ---  [2-1] Modify the grid      --- !
       call LinearInterpFWD( rf(2:LJs), rhs(:,2:LJs), rLIs, rhsLIs, LJs-2, LIs )
       drLIs = rLIs(2) - rLIs(1)
       
       ! ---  [2-2] Solve Poisson's eq.  --- !
       call BiCGSTABCyl2D( dphiLIs(2:LIs,2:LJs), rhsLIs(2:LIs,2:LJs), drLIs, dz, LJs-2, LIs-1, rmin )
       
       ! ---  [2-3] Modify the grid      --- !
       call LinearInterpRET( rLIs, dphiLIs, rf(2:LJs), dphi(:,2:LJs), LJs-2, LIs )
    endif

    ! --- [3] rmin > 0.d0 case --- !
    if ( rmin.gt.0.d0 ) then
       call BiCGSTABCyl2D( dphi, rhs, dz, dr, LIs, LJs, rmin )
    endif

    ! --- [4] correction of E-field --- !
    do j=2, LJs
       do i=2, LIs-1
          EBf(er_,i,j) = EBf(er_,i,j) - drInv * ( dphi(i,j) - dphi(i,j-1) )
       enddo
    enddo
    do j=2, LJs-1
       do i=2, LIs
          EBf(ez_,i,j) = EBf(ez_,i,j) - dzInv * ( dphi(i,j) - dphi(i-1,j) )
       enddo
    enddo

    ! --- [5] Boundary Condition  --- !
    rqp     = 0.5d0 * ( rh( 2) + rf(   2) )
    rqm     = 0.5d0 * ( rf( 2) + rf(   1) )
    factor1 =   rqp / rqm
    rqp     = 0.5d0 * ( rf(LJs) + rh(LJs  ) )
    rqm     = 0.5d0 * ( rf(LJs) + rh(LJs-1) )
    factor2 =   rqm / rqp
    do i=1, LIs
       EBf(ez_,i, 1) = - factor1 * EBf(ez_,i,   2)
       EBf(ez_,i,LJs) = - factor2 * EBf(ez_,i,LJs-1)
    enddo
    do j=1, LJs
       EBf(er_, 1,j) = -           EBf(er_,   2,j)
       EBf(er_,LIs,j) = -           EBf(er_,LIs-1,j)
    enddo

    return
  end subroutine Ecorrect


  subroutine nTimesFiltering
    use constants , only : LIs, LJs, ns, SmCnt
    use variables , only : rh, rhInv, w_p
    use smoothfMod, only : smooth_A_linear
    implicit none
    integer            :: i, j, k, iter

    ! --- [1] Cylindrical Correction --- !
    !$omp parallel default(none) &
    !$omp shared(rh,w_p,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             w_p(i,j,k) = w_p(i,j,k) * rh(j)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    ! --- [2] n times Smoothing --- !
    do iter=1, SmCnt
       write(6,'(i4)',advance='no') iter
       call smooth_A_linear( w_p )
    enddo
    write(6,*)
    ! --- [3] Cylindrical Back --- !
    !$omp parallel default(none) &
    !$omp shared(rhInv,w_p,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             w_p(i,j,k) = w_p(i,j,k) * rhInv(j)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine nTimesFiltering
  
  
end module qchargeMod
