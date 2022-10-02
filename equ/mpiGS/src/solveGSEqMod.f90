module solveGSEqMod
  implicit none
  character(2), parameter :: jdg = '1D'
contains

  ! ====================================================== !
  ! === determine__psiPotential                        === !
  ! ====================================================== !
  subroutine determine__psiPotential
    use variablesMod
    use NewtonRaphson, only  : FindExtremum2D, cSplineNewtonRaphson1D, cSpline1D, cSpline2D
    implicit none
    integer                 :: i, j, flag
    double precision        :: probe, det, psizz, dpsidr, dpsidz
    double precision        :: dpdr(LJ-1), rAxis_h(LJ-1), rSep(2)
    
    ! ------------------------------------------------------ !
    ! --- [1] half ver.  Find magnetic Axis              --- !
    ! ------------------------------------------------------ !
    if ( ( trim(boundaryType) == "hlf" ) ) then
       psiMin = 1.d8
       do j=1, LJ
          psiMin = min( psiMin, psi(2,j) )
       enddo
       call cSpline1D( rAxis, psi(2,:), LJ, limiter_rPos(1), psiLim, dpsidr )
    endif
    
    ! ------------------------------------------------------ !
    ! --- [2] full ver.  Find magnetic Axis              --- !
    ! ------------------------------------------------------ !
    if ( ( trim(boundaryType) == "pcw" ) ) then
       ! ------------------------------------------------------ !
       ! --- [2-1] 2D routines ( unstable ? )               --- !
       ! ------------------------------------------------------ !
       if ( jdg.eq."2D" ) then
          call FindExtremum2D( zAxis, rAxis, psi, magAxis(z_), magAxis(r_), &
               &               psiMin,  det, psizz, flag )
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
          call cSplineNewtonRaphson1D( rAxis_h, dpdr, magAxis(r_), psiMin, flag  )
          call cSpline1D( rAxis, psi(LI/2,:), LJ, magAxis(r_), psiMin, dpsidr )
          magAxis(z_) = 0.d0
          det         = 1.d0
          psizz       = 1.d0
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
          det         = 1.d0
          psizz       = 1.d0
          flag        = 1
       endif

       ! ------------------------------------------------------ !
       ! --- [2-4] psi @ Limiter                            --- !
       ! ------------------------------------------------------ !
       psiLim = - 1.d5
       do i=1, nLimiter
          write(6,*) minval( zAxis(:) ), maxval( zAxis(:) )
          write(6,*) minval( rAxis(:) ), maxval( rAxis(:) )
          call cSpline2D( zAxis, rAxis, psi, limiter_zPos(i), limiter_rPos(i), &
               &          probe, dpsidr, dpsidz )
          write(6,*) probe, psiLim, limiter_zpos(i), limiter_rPos(i), minval( psi(:,:) ), maxval(psi(:,:) )
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
    ! write(6,*) " rSep(1), rSep(2), aMinor == ", rSep(1), rSep(2), aMinor

    ! ------------------------------------------------------ !
    ! --- [4] Display Result                             --- !
    ! ------------------------------------------------------ !

    if ( myRank.eq.0 ) then
       write(6,'(4x,a)') '[ Determine Psi     ]'
       if ( psiMin.ne.psiMin ) stop '[ERROR] NaN value -- psiMin [ERROR]'

       if (  flag.eq.0 ) then
          write(6,'( 8x,a)') '====================================================='
          write(6,'(10x,a)') '[CAUTION]   psiMin did not converged   [CAUTION] '
          write(6,'( 8x,a)') '====================================================='
       else
          if ( ( det.gt.0.d0 ).and.( psizz.gt.0.0d0 ) ) write(6,'(8x,a)') 'Converged ---  Extremum minimum value'
          if ( ( det.gt.0.d0 ).and.( psizz.le.0.0d0 ) ) write(6,'(8x,a)') 'Converged ---  Extremum maximum value'
          if (   det.lt.0.d0                          ) write(6,'(8x,a)') 'Converged ---  Saddle Point'
       endif
       write(6,'(8x,a,2(f8.4,a),e12.5)') &
            & 'Mag Axis   = ( ', magAxis(r_), ',', magAxis(z_), ')  -- psiMin = ', psiMin
       write(6,'(8x,a,2(f8.4,a),e12.5)') &
            & 'Limiter    = ( ', limiter_rPos(limiter_index), ',', limiter_zPos(limiter_index), ')  -- psiLim = ', psiLim
       write(6,'(8x,a,7x,a,e12.5)') &
            & 'Boundary   = Conducting Wall', '-- Bv0    = ', Bv0
    endif

    return
  end subroutine determine__psiPotential
  
  
end module solveGSEqMod
