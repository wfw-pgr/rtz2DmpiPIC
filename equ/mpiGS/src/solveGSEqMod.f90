module solveGSEqMod
  implicit none
  character(2)    , parameter :: jdg = '1D'
  double precision, parameter :: convergence_factor = 1.e-8
contains

  
  ! ====================================================== !
  ! === solve__GradShafranov                           === !
  ! ====================================================== !
  subroutine solve__GradShafranov
    use variablesMod
    use utilitiesMod, only : print__section, TimeMeasure_MPI
    use PBiCGSTABMod, only : BiCGSTABCyl2D_MPI
    use currentFnMod, only : update__currentDistribution
    use boundaryCMod, only : set__initial_boundaryCondition
    implicit none
    integer        , parameter :: sideLen      = 3
    integer        , parameter :: frameLen     = 70
    integer        , parameter :: Linear_Loop_ = 1
    integer        , parameter :: NonLin_Loop_ = 2
    character(cLen), parameter :: message1 = "Beginning of Main Loop"
    character(cLen), parameter :: message2 = "('Linear Loop = ',(i4),4x,'NonLinear Loop = ',(i4))"
    character(cLen), parameter :: message3 = "Linear Solver Converged...."
    character(cLen)            :: message_
    character(4)               :: citer
    integer                    :: iter1, iter2
    double precision           :: diff1, diff2
    logical                    :: flag__converged1, flag__converged2

    ! ------------------------------------------------------ !
    ! --- [1] display section label                      --- !
    ! ------------------------------------------------------ !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*)
       call print__section( trim(message1), "=", len(trim(message1)), sideLen, frameLen )
       write(6,*)
    endif
    
    ! ------------------------------------------------------ !
    ! --- [2] Non-Linear solver Loop                     --- !
    ! ------------------------------------------------------ !
    iter2            = 0
    flag__converged2 = .false.
    
    do while ( ( .not.( flag__converged2 ) ).and.( iter2.lt.itermax2 ) )
       
       ! ------------------------------------------------------ !
       ! --- [3] Linear solver Loop                         --- !
       ! ------------------------------------------------------ !
       iter1            = 0
       flag__converged1 = .false.
       
       do while ( ( .not.( flag__converged1 ) ).and.( iter1.lt.itermax1 ) )

          ! ------------------------------------------------------ !
          ! --- [3-1] display inner-loop title                 --- !
          ! ------------------------------------------------------ !
          if ( myRank.eq.0 ) then
             write(message_,message2) iter1, iter2
             call print__section( trim(message_), "-", len(trim(message_)), sideLen, frameLen )
          endif
          call TimeMeasure_MPI(3)

          ! ------------------------------------------------------ !
          ! --- [3-2] solve linear equation                    --- !
          ! ------------------------------------------------------ !
          !
          ! if ( solverType.eq."SOR"      ) call solv_poisson
          ! if ( solver.eq.'BCG' ) call BiCGSTABCyl2D    ( psi, rhs, dz, dr, Nz, Nr, rmin )
          !
          if ( solverType.eq."BiCGSTAB" ) then
             call BiCGSTABCyl2D_MPI( psi, rhs, dz, dr, LI, LJ, rMin )
          endif
          call TimeMeasure_MPI(2)

          ! ------------------------------------------------------ !
          ! --- [3-3] determine psi potential                  --- !
          ! ------------------------------------------------------ !
          call determine__psiPotential
          if ( ( iter1.gt.0 ).or.( iter2.gt.0 ) ) then
             call update__PicardIteration
          endif
          
          ! ------------------------------------------------------ !
          ! --- [3-4] check convergence                        --- !
          ! ------------------------------------------------------ !
          call check__convergence( Linear_Loop_, iter1, flag__converged1, convergence1, diff1 )
          iter1 = iter1 + 1
          call TimeMeasure_MPI(3)

          ! ------------------------------------------------------ !
          ! --- [3-5] update current Distribution              --- !
          ! ------------------------------------------------------ !
          call update__currentDistribution
          call TimeMeasure_MPI(4)

          ! ------------------------------------------------------ !
          ! --- [3-6] save log                                 --- !
          ! ------------------------------------------------------ !
          call write__logOfLoop( iter1, iter2, diff1, diff2 )
          call TimeMeasure_MPI(3)
          
       enddo

       ! ------------------------------------------------------ !
       ! --- [4] Display End of Linear solver Loop          --- !
       ! ------------------------------------------------------ !
       if ( myRank.eq.0 ) then
          call print__section( trim(message3), "-", len(trim(message3)), sideLen, frameLen )
       endif
       
       ! ------------------------------------------------------ !
       ! --- [5] check Non-Linear Loop convergence          --- !
       ! ------------------------------------------------------ !
       call check__convergence( NonLin_Loop_, iter2, flag__converged2, convergence2, diff2 )
       iter2 = iter2 + 1
       if ( .not.( flag__converged2 ) ) then
          call set__initial_boundaryCondition
       endif
       call TimeMeasure_MPI(3)
       
    enddo
    
    return
  end subroutine solve__GradShafranov
  
  
  ! ====================================================== !
  ! === determine__psiPotential                        === !
  ! ====================================================== !
  subroutine determine__psiPotential
    use variablesMod
    use cubSplineMod, only : interpolate__cubicSpline1D_point, interpolate__cubicSpline2D_point
    use newtonRapMod, only : solve__cSpline_NewtonRaphson1D, solve__cSpline_NewtonRaphson2D
    implicit none
    integer                 :: i, j, flag
    double precision        :: probe, dpsidr, dpsidz
    double precision        :: dpdr(LJ-1), rAxis_h(LJ-1), rSep(2)
    character(cLen)         :: fmt1, fmt2, fmt3, cvar
    integer, parameter      :: p_ = 3
    
    ! ------------------------------------------------------ !
    ! --- [1] half ver.  Find magnetic Axis              --- !
    ! ------------------------------------------------------ !
    if ( ( trim(boundaryType) == "hlf" ) ) then
       psiMin = 1.d8
       do j=1, LJ
          psiMin = min( psiMin, psi(2,j) )
       enddo
       call interpolate__cubicSpline1D_point( rAxis, psi(2,:), limiter_rPos(1), psiLim, dpsidr )
    endif
    
    ! ------------------------------------------------------ !
    ! --- [2] full ver.  Find magnetic Axis              --- !
    ! ------------------------------------------------------ !
    if ( ( trim(boundaryType) == "pcw" ) ) then
       ! ------------------------------------------------------ !
       ! --- [2-1] 2D routines ( unstable ? )               --- !
       ! ------------------------------------------------------ !
       if ( jdg.eq."2D" ) then
          call solve__cSpline_NewtonRaphson2D( zAxis, rAxis, psi, &
               &                               magAxis(z_), magAxis(r_), psiMin,  &
               &                               flag, "max", convergence_factor )
       endif
       ! ------------------------------------------------------ !
       ! --- [2-2] 1D routines ( relatively stable ? )      --- !
       ! ------------------------------------------------------ !
       if ( jdg.eq."1D" ) then
          !$omp parallel default(none) &
          !$omp shared(LI,LJ,rAxis,rAxis_h,psi,dpdr,drInv) private(j)
          !$omp do
          do j=1, LJ-1
             rAxis_h(j) = ( rAxis(j+1)    + rAxis(j)    ) * 0.5d0
             dpdr(j)    = ( psi(LI/2,j+1) - psi(LI/2,j) ) * drInv
          enddo
          !$omp end do
          !$omp end parallel
          call solve__cSpline_NewtonRaphson1D( rAxis_h, dpdr, magAxis(r_), psiMin, &
               &                               flag, "max", convergence_factor )
          call interpolate__cubicSpline1D_point( rAxis, psi(LI/2,:), magAxis(r_), &
               &                                 psiMin, dpsidr )
          magAxis(z_) = 0.d0
          flag        = 1
       endif
       ! ------------------------------------------------------ !
       ! --- [2-3] Find Minimum ( rough )                   --- !
       ! ------------------------------------------------------ !
       if ( jdg.eq.'mn' ) then
          psiMin = 1.d8
          do j=1, LJ
             if ( psi(LI/2,j).lt.psiMin ) then
                magAxis(r_)  = rAxis(j)
                psiMin       = psi(LI/2,j)
             endif
          enddo
          magAxis(z_) = 0.d0
          flag        = 1
       endif

       ! ------------------------------------------------------ !
       ! --- [2-4] psi @ Limiter                            --- !
       ! ------------------------------------------------------ !
       psiLim = - 1.d5
       do i=1, nLimiter
          call interpolate__cubicSpline2D_point( zAxis, rAxis, psi, &
               &                                 limiter_zPos(i), limiter_rPos(i), &
               &                                 probe, dpsidr, dpsidz )
          limiter_psi(z_,i) = limiter_zPos(i)
          limiter_psi(r_,i) = limiter_rPos(i)
          limiter_psi(p_,i) = probe
          if ( probe > psiLim ) then
             psiLim        = probe
             limiter_index = i
          endif
       enddo
       
    endif

    ! ------------------------------------------------------ !
    ! --- [3] Minor Radius from Separatrix search        --- !
    ! ------------------------------------------------------ !
    rSep(:) = 0.d0
    do j=5, LJ-5, +1
       if ( ( psi(LI/2,j)-psiLim )*( psi(LI/2,j-1)-psiLim ).le.0.d0 ) then
          rSep(1) = rAxis(j)
          exit
       endif
    enddo
    do j=LJ-5, 5, -1
       if ( ( psi(LI/2,j)-psiLim )*( psi(LI/2,j-1)-psiLim ).le.0.d0 ) then
          rSep(2) = rAxis(j)
          exit
       endif
    enddo
    aMinor = 0.5d0 * ( rSep(2)-rSep(1) )

    ! ------------------------------------------------------ !
    ! --- [4] Display Result                             --- !
    ! ------------------------------------------------------ !
    if ( myRank.eq.0 ) then
       write(6,'(4x,a)') '[ Determine Psi     ]'
       if ( psiMin.ne.psiMin ) stop '[ERROR] NaN value -- psiMin [ERROR]'
       
       if (  flag.eq.-1 ) then
          write(6,'( 8x,a)') '====================================================='
          write(6,'(10x,a)') '[CAUTION]   psiMin did not converged   [CAUTION] '
          write(6,'( 8x,a)') '====================================================='
       else
          if ( flag.eq.0 ) write(6,'(8x,a)') 'Converged ---  Saddle Point'
          if ( flag.eq.1 ) write(6,'(8x,a)') 'Converged ---  Extremum minimum value'
          if ( flag.eq.2 ) write(6,'(8x,a)') 'Converged ---  Extremum maximum value'
       endif
       write(6,'(8x,a,2(f8.4,a),e12.5)') &
            & 'Mag Axis   = ( ', magAxis(r_), ',', magAxis(z_), ')  -- psiMin = ', psiMin
       write(6,'(8x,a,2(f8.4,a),e12.5)') &
            & 'Limiter    = ( ', limiter_rPos(limiter_index), ',', limiter_zPos(limiter_index), ')  -- psiLim = ', psiLim
       write(6,'(8x,a,7x,a,e12.5)') &
            & 'Boundary   = Conducting Wall', '-- Bv0    = ', Bv0
    endif

    if ( myRank.eq.0 ) then
       open (lun,file="dat/limiter.dat",status="replace")
       fmt1 = '("limiter_psi(",i0,")")'
       fmt2 = '(a30,4x,"fltarr",4x,"[",2(e12.5,","),e12.5,"]")'
       fmt3 = '(a30,4x,"float",4x,e12.5)'
       do i=1, nLimiter
          write(cvar,trim(fmt1)) i
          write(lun ,trim(fmt2)) trim(cvar), limiter_psi(z_,i), limiter_psi(r_,i), &
               &                             limiter_psi(p_,i)
       enddo
       write(lun,trim(fmt2)) "magAxis", magAxis(z_), magAxis(r_), psiMin
       write(lun,trim(fmt3)) "psiMin", psiMin
       write(lun,trim(fmt3)) "psiLim", psiLim
    end if
    
    return
  end subroutine determine__psiPotential


  ! ====================================================== !
  ! === check convergence of the iteration             === !
  ! ====================================================== !
  subroutine check__convergence( loopType, iter, flag__converged, convergence_criterion, diff )
    use variablesMod
    implicit none
    integer                         :: i, j
    integer         , intent(in)    :: loopType    ! 1:linear_loop, 2:nonLinear_loop
    integer         , intent(inout) :: iter
    logical         , intent(out)   :: flag__converged
    double precision, intent(in)    :: convergence_criterion
    double precision, intent(out)   :: diff
    
    ! ------------------------------------------------------ !
    ! --- [1] Calculation of difference from past_value  --- !
    ! ------------------------------------------------------ !
    diff = 0.0d0
    !$omp parallel default(none) &
    !$omp shared(diff,psi,psi_past,loopType,LI,LJ) private(i,j)
    !$omp do reduction(+:diff)
    do j=1, LJ
       do i=1, LI
          diff = diff + ( psi(i,j) - psi_past(i,j,loopType) )**2
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    diff = sqrt( diff / dble( LI*LJ ) )
    
    ! ------------------------------------------------------ !
    ! --- [2] judge whether converged or not             --- !
    ! ------------------------------------------------------ !
    if ( ( iter.gt.0  ).and.( diff.lt.convergence_criterion ) ) then 
       flag__converged = .false.
    endif
    
    ! ------------------------------------------------------ !
    ! --- [3] save psi as psi_past                       --- !
    ! ------------------------------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(psi_past,psi,loopType,LI,LJ) private(i,j)
    !$omp do
    do j=1, LJ
       do i=1, LI
          psi_past(i,j,loopType) = psi(i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! ------------------------------------------------------ !
    ! --- [4] Output Convergence                         --- !
    ! ------------------------------------------------------ !
    if ( myRank.eq.0 ) then       
       write(6,'(4x,a)') '[ Check Convergence ]'
       write(6,'(8x,a,e12.5)') 'Difference = ', diff
    endif

    return
  end subroutine check__convergence


  ! ====================================================== !
  ! === update__PicardIteration                        === !
  ! ====================================================== !
  subroutine update__PicardIteration
    use variablesMod
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2
    integer, parameter :: linear_loop_ = 1

    ! ------------------------------------------------------ !
    ! --- [1] update Psi using PicardIteration           --- !
    ! ------------------------------------------------------ !
    factor1 =        picard_alphaB
    factor2 = 1.d0 - picard_alphaB
    !$omp parallel default(none) &
    !$omp shared(psi,psi_past,factor1,factor2,LI,LJ) private(i,j)
    !$omp do
    do j=1, LJ
       do i=1, LI
          psi(i,j) = factor1*psi(i,j) + factor2*psi_past(i,j,linear_loop_)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine update__PicardIteration


  ! ====================================================== !
  ! === write__logOfLoop                               === !
  ! ====================================================== !
  subroutine write__logOfLoop( iter1, iter2, diff1, diff2 )
    use variablesMod
    implicit none
    integer         , intent(in) :: iter1, iter2
    double precision, intent(in) :: diff1, diff2
    double precision             :: delta_psi
    logical         , save       :: flag__initialize = .true.

    ! ------------------------------------------------------ !
    ! --- [1] renew file                                 --- !
    ! ------------------------------------------------------ !
    if ( flag__initialize ) then
       if ( myRank.eq.0 ) then
          open (lun,file=trim(convergenceFile),form='formatted',status='replace')
          close(lun)
       endif
       flag__initialize = .false.
    end if
    
    ! ------------------------------------------------------ !
    ! --- [2] write log of a loop                        --- !
    ! ------------------------------------------------------ !
    delta_psi = psiLim - psiMin
    if ( myRank.eq.0 ) then
       open (lun,file=trim(convergenceFile),form='formatted',access='append')
       write(lun,'(2(i4,1x),10(e12.5,1x))') iter1, iter2, diff1, diff2, raxis, zaxis &
            &                              , g0, alpha_g, Bv0, psilim, psimin, delta_psi
       close(lun)
    endif
    
    return
  end subroutine write__logOfLoop

  
end module solveGSEqMod
