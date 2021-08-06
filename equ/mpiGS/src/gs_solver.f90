module gs_solver
contains

  subroutine SolveGS

    use constants  , only : Nr, Nz, itermax1, itermax2, rmin, solver, jobdir, myRank, PEtot
    use constants  , only : convergence1, convergence2
    use variables  , only : psi, r, z, rhs, dr, dz
    use bdrc       , only : updateBC
    use sorsolver  , only : solv_poisson
    use myPBiCGSTAB, only : BiCGSTABCyl2D_MPI
    use myutil     , only : TimeMeasure
    implicit none
    integer              :: iter1, iter2, cnv_flag1, cnv_flag2
    double precision     :: diff1, diff2

    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*) '  -------------------------------------------  '
       write(6,*) '  ---------      Main    Loop     -----------  '
       write(6,*) '  -------------------------------------------  '
    endif

    !  --- Main Loop ( iterate until psi reach certain convergence ) ---  !
    cnv_flag2    = 0
    iter2        = 0
    if ( myRank.eq.0 ) then
       open (50,file=trim(jobdir)//'convergence.log',form='formatted',status='replace')
       close(50)
    endif
    !  -- Outer Loop :: Update Boundary --  !
    do while( ( cnv_flag2.ne.1 ).and.( iter2.lt.itermax2 ) )

       cnv_flag1 = 0
       iter1     = 0
       !  -- Inner Loop :: Update Current --  !
       do while( ( cnv_flag1.ne.1 ).and.( iter1.lt.itermax1 ) )
          if ( myRank.eq.0 ) then
             write(6,*)
             write(6,*)
             write(6,'(3x,2(a,i8))') 'Iteration1     == ', iter1, '       of Iteration2   == ', iter2
             write(6,*) '--------------------------------------------------------------------'
          endif
          call TimeMeasure(3)
          if ( solver.eq.'SOR' ) call solv_poisson
          ! if ( solver.eq.'BCG' ) call BiCGSTABCyl2D( psi, rhs, dz, dr, Nz, Nr, rmin )
          if ( solver.eq.'BCG' ) call BiCGSTABCyl2D_MPI( psi, rhs, dz, dr, Nz, Nr, rmin )
          call TimeMeasure(2)

          call DeterminePsi
          if ( ( iter1.ne.0 ).or.( iter2.ne.0 ) ) call PicardIteration
          call checkConvg( 1, iter1, cnv_flag1, convergence1, diff1 )
          call TimeMeasure(3)
          call UpdateJphi
          call TimeMeasure(4)
          call WriteLog( iter1, iter2, diff1, diff2 )
          call TimeMeasure(3)

       enddo

       if ( myRank.eq.0 ) then
          write(6,*)
          write(6,*) ' ======================================================================= '
          write(6,'(3x,a)') ' Inner Loop Converged....  '
       endif
       call checkConvg( 2, iter2, cnv_flag2, convergence2, diff2 )
       if ( cnv_flag2 .ne. 1 ) call updateBC
       if ( myRank.eq.0 ) then
          write(6,*) ' ======================================================================= '
          write(6,*)
       endif
       call TimeMeasure(3)

    enddo

    return
  end subroutine SolveGS


  subroutine DeterminePsi

    use constants, only      : Nr, Nz, bctype, myRank
    use constants, only      : Nlim, limposr, limposz, limidx, Bv0
    use variables, only      : r, z, drinv
    use variables, only      : psi, psimin, psilim
    use variables, only      : raxis, zaxis, aMinor
    use NewtonRaphson, only  : FindExtremum2D, cSplineNewtonRaphson1D, cSpline1D, cSpline2D
    implicit none
    integer                 :: i, j, flag
    double precision        :: probe, det, psizz, dpsidr, dpsidz
    double precision        :: dpdr(Nr-1), rh(Nr-1), rSep(2)
    character(2), parameter :: jdg = '1D'
    
    !  --- [1] Get psimin ---  !
    !   -- [1-1] half ver. --  !
    if ( ( bctype.eq.'hlf' ) ) then
       psimin = 1.d8
       do j=1, Nr
          psimin = min( psimin, psi(2,j) )
       enddo
       call cSpline1D( r, psi(2,:), Nr, limposr(1), psilim, dpsidr )
    endif
    
    if ( ( bctype.eq.'pcw' ) ) then
       ! -- [1-2] Magnetic Axis -- !
       if ( jdg.eq.'2D' ) then
          call FindExtremum2D( z, r, psi, zaxis, raxis, psimin, det, psizz, flag )
       endif
       if ( jdg.eq.'1D' ) then
          !$omp parallel default(none) &
          !$omp shared(r,rh,psi,dpdr,drinv) private(j)
          !$omp do
          do j=1, Nr-1
             rh(j)   = ( r(j+1) + r(j) ) * 0.5d0
             dpdr(j) = ( psi(Nz/2,j+1) - psi(Nz/2,j) ) * drinv
          enddo
          !$omp end do
          !$omp end parallel
          call cSplineNewtonRaphson1D( rh, dpdr, raxis, psimin, flag  )
          call cSpline1D( r, psi(Nz/2,:), Nr, raxis, psimin, dpsidr )
          zaxis = 0.d0
          det   = 1.d0
          psizz = 1.d0
          flag  = 1
       endif
       if ( jdg.eq.'mn' ) then
          psimin = 1.d8
          do j=1, Nr
             if ( psi(Nz/2,j).lt.psimin ) then
                raxis  = r(j)
                psimin = psi(Nz/2,j)
             endif
          enddo
          zaxis = 0.d0
          det   = 1.d0
          psizz = 1.d0
          flag  = 1
       endif
       ! -- [1-3] psi @ Limiter -- !
       psilim = - 1.d5
       do i=1, Nlim
          call cSpline2D( z, r, psi, limposz(i), limposr(i), probe, dpsidr, dpsidz )
          if ( probe.gt.psilim ) then
             psilim = probe
             limidx = i
          endif
       enddo
       
    endif
    ! -- [1-4] Minor Radius -- !
    rSep(:) = 0.d0
    do j=5, Nr-5, +1
       if ( ( psi(Nz/2,j)-psilim )*( psi(Nz/2,j-1)-psilim ).le.0.d0 ) then
          rSep(1) = r(j)
          exit
       endif
    enddo
    do j=Nr-5, 5, -1
       if ( ( psi(Nz/2,j)-psilim )*( psi(Nz/2,j-1)-psilim ).le.0.d0 ) then
          rSep(2) = r(j)
          exit
       endif
    enddo
    aMinor = 0.5d0 * ( rSep(2)-rSep(1) )
    ! write(6,*) " rSep(1), rSep(2), aMinor == ", rSep(1), rSep(2), aMinor

    ! --- [2] Display --- !
    if ( myRank.eq.0 ) then       
       write(6,'(4x,a)') '[ Determine Psi     ]'
       if ( psimin.ne.psimin ) stop '[ERROR] NaN value -- psimin [ERROR]'

       if (  flag.eq.0 ) then
          write(6,'( 8x,a)') '====================================================='
          write(6,'(10x,a)') '[CAUTION]   psimin did not converged   [CAUTION] '
          write(6,'( 8x,a)') '====================================================='
       else
          if ( ( det.gt.0.d0 ).and.( psizz.gt.0.0d0 ) ) write(6,'(8x,a)') 'Converged ---  Extremum minimum value'
          if ( ( det.gt.0.d0 ).and.( psizz.le.0.0d0 ) ) write(6,'(8x,a)') 'Converged ---  Extremum maximum value'
          if (   det.lt.0.d0                          ) write(6,'(8x,a)') 'Converged ---  Saddle Point'
       endif
       write(6,'(8x,a,2(f8.4,a),e12.5)') &
            & 'Mag Axis   = ( ', raxis,           ',', zaxis,           ')  -- psiMin = ', psimin
       write(6,'(8x,a,2(f8.4,a),e12.5)') &
            & 'Limiter    = ( ', limposr(limidx), ',', limposz(limidx), ')  -- psiLim = ', psilim
       write(6,'(8x,a,7x,a,e12.5)') &
            & 'Boundary   = Conducting Wall', '-- Bv0    = ', Bv0
    endif

    return
  end subroutine DeterminePsi

  
  subroutine checkConvg( tgt, iter, cnv_flag, convergence, diff )

    use variables,            only : psi, psipast
    use constants,            only : Nr, Nz, myRank
    implicit none
    integer                       :: i, j
    integer, intent(in)           :: tgt
    integer, intent(inout)        :: iter, cnv_flag
    double precision, intent(in)  :: convergence
    double precision, intent(out) :: diff
    
    !  --- [1] Calculation of Difference --- !
    diff = 0.0d0
    !$omp parallel default(none) &
    !$omp shared(diff,psi,psipast,tgt) private(i,j)
    !$omp do reduction(+:diff)
    do j=1, Nr
       do i=1, Nz
          diff = diff + ( psi(i,j) - psipast(i,j,tgt) )**2
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    diff = sqrt( diff / dble( Nz*Nr ) )
    !  --- [2] Judgement of Convergence  --- !
    if ( ( iter.ne.0  ).and.( diff.lt.convergence ) ) then 
       cnv_flag = 1
    endif
    !  --- [3] Copy psi to the past      --- !
    !$omp parallel default(none) &
    !$omp shared(psipast,psi,tgt) private(i,j)
    !$omp do
    do j=1, Nr
       do i=1, Nz
          psipast(i,j,tgt) = psi(i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  --- [4] Output Convergence        --- !
    iter = iter + 1
    if ( myRank.eq.0 ) then       
       write(6,'(4x,a)') '[ Check Convergence ]'
       write(6,'(8x,a,e12.5)') 'Difference = ', diff
    endif
       
    return
  end subroutine checkConvg

  
  subroutine UpdateJphi
    
    !! see Note and 'Computational methods in Plasma Physics' by Stephen Jardin
    ! this is basic method to calculate magnetic equilibrium 
    ! remaining total plasma current Ip and peak value of pressure P0
    ! --- Ip      : conserved
    ! --- alpha   : given parameter
    ! --- epsilon : given parameter
    ! Originally written by Kent in 2015 / 11 / 18
    !                   --- Updated. @ 2018/ 2 / 2
    
    ! (i) current and pressure distribution
    ! consider following pressure and toroidal magnetic field function profile      
    !  (1) --    p ( psi ) = p0 * norm_p ( norm_psi )   --
    !!     -- normalized p ( psi ) = normalized psi^epsilon  ::  p(s) = s^e
    !  (2) -- 1/2 g^2 ( psi ) = 1/2 g0^2 ( 1 + alpha_g * norm_g^2 ( norm_psi ) )
    !! ---- normalized g ( psi ) = normalized psi^alpha
    !! -------------- alpha g is a coefficient to be determined below...
    !! no need for coding ( see Lab Note :: only calculation )

    ! (v) calculate J_phi ( toroidal current ) ( 1st time ( no adjustment for Ip=const.) )
    !! no need for coding ( see Lab Note :: only calculation )
    
    ! (vi) calculate p'(psi) and g'(psi) ( which is normilzed :: p,g,psi )
    
    ! (vii) calculate alpha_g ( coefficient for toroidal field function )
    ! integral over [ Nr, Nz ]
    
    use constants, only : Nr, Nz, solver, equtype, myRank, normInGS, valfe
    use variables, only : dr, dz, r, rinv
    use variables, only : psi, psimin, psilim, rhs, raxis
    use variables, only : Jphi, prs, gfunc, BdrcLR, BdrcTB, Itotal, g0, alpha_g
    use constants, only : epsilon, alpha, Betat, TiTe, vthcv
    use psifunc  , only : pghat
    implicit none
    integer            :: i, j
    double precision   :: psis(Nz,Nr), g2(Nz,Nr), p0, InSep(Nz,Nr)
    double precision   :: dpdp(Nz,Nr), dgdp(Nz,Nr), phat(Nz,Nr), ghat(Nz,Nr)
    double precision   :: Jp, Jg, denom, numer, delta_psi, dpsiinv, Isum
    double precision   :: factor1, factor2, factor3, coef
    
    !  --- [0] preparation  ---  !
    !$omp parallel default(none) &
    !$omp shared(Jphi,prs,gfunc,psis,phat,ghat,dpdp,dgdp,InSep) private(i,j)
    !$omp do
    do j=1, Nr
       do i=1, Nz
          Jphi (i,j) = 0.d0
          prs  (i,j) = 0.d0
          gfunc(i,j) = 0.d0
          psis (i,j) = 0.d0
          InSep(i,j) = 0.d0
          phat (i,j) = 0.d0
          ghat (i,j) = 0.d0
          dpdp (i,j) = 0.d0
          dgdp (i,j) = 0.d0
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    Jp         = 0.d0
    Jg         = 0.d0
    delta_psi  = 0.d0
    alpha_g    = 0.d0
    
    !  --- [1]  Normalized psi,  p(psi), g(psi) ---  !
    p0         = 1.d0
    coef       = Betat
    g0         = sqrt( 2.d0 / coef ) * raxis
    delta_psi  = ( psilim - psimin )
    dpsiinv    = 1.d0 / delta_psi
    !$omp parallel default(none) &
    !$omp shared(psis,psilim,psi,dpsiinv) private(i,j)
    !$omp do
    do j=1, Nr
       do i=1, Nz
          psis(i,j) = ( psilim - psi(i,j) ) * dpsiinv
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    call pghat( psis, phat, ghat, dpdp, dgdp )
    !$omp parallel default(none) &
    !$omp shared(psis,InSep) private(i,j)
    !$omp do
    do j=1, Nr
       do i=1, Nz
          if ( psis(i,j).gt.0.d0 ) then
             InSep(i,j) = 1.d0
          endif
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !  --- [2] Update R.H.S., p(psi), g(psi) ---  !
    !   -- Jphi, g(psi), p(psi) --  !
    !   -- [2-A] ST -- !
    if ( equtype.eq.'STk' ) then
       ! -- [2-A-1] Current Density  -- !
       Jg  = 0.d0
       Jp  = 0.d0
       !$omp parallel default(none) &
       !$omp shared(r,rinv,dpdp,dgdp,psis,InSep,Jg,Jp) private(i,j)
       !$omp do reduction(+:Jg,Jp)
       do j=1, Nr
          do i=1, Nz
             Jp = Jp +    r(j) * dpdp(i,j) * InSep(i,j)
             Jg = Jg + rinv(j) * dgdp(i,j) * InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       ! -- [2-A-2] alpha_g  -- !
       numer   = ( delta_psi * Itotal / ( dr*dz ) - Jp )
       denom   = 0.5d0 * g0**2 * Jg
       alpha_g = numer / denom
       if ( alpha_g .lt. -1.d0 ) then
          ! if ( alpha_g .lt. 0.d0 ) then
          write(6,*)
          write(6,*) ' [CAUTION] alpha_g is less than -1.0 [CAUTION] '
          write(6,*) '               alpha_g === ', alpha_g
          write(6,*) ' See Lab. Note :: Too weak Total Current to confine pressure :: Itotal :: '
          write(6,*) '                         ---> Adjust the poloidal Beta. '
          write(6,*)
          stop
       endif
       ! -- [2-A-3] prs, gfunc -- !
       factor1 =                           dpsiinv
       factor2 = 0.5d0 * alpha_g * g0**2 * dpsiinv
       !$omp parallel default(none) &
       !$omp shared(r,rinv,dpdp,dgdp,Jphi,prs,gfunc,phat,ghat,g0,alpha_g,InSep,factor1,factor2) &
       !$omp private(i,j,Jp,Jg)
       !$omp do
       do j=1, Nr
          do i=1, Nz
             ! if ( psis(i,j).ge.0.d0 ) then
             Jp          = factor1 * r(j)    * dpdp(i,j) * InSep(i,j)
             Jg          = factor2 * rinv(j) * dgdp(i,j) * InSep(i,j)
             
             Jphi(i,j)   = Jp  +  Jg
             prs(i,j)    = phat(i,j) * InSep(i,j)
             gfunc(i,j)  = g0 * sqrt( 1.d0 + alpha_g * ghat(i,j) ) * InSep(i,j)
             ! else
             !    Jphi(i,j)   = 0.d0
             !    prs(i,j)    = 0.d0
             !    gfunc(i,j)  = g0
             ! endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif
    !   -- [2-B] Spheromak --  !
    if ( equtype.eq.'Sph' ) then
       ! -- [2-B-1] Current Density  -- !
       Jg  = 0.d0
       Jp  = 0.d0
       !$omp parallel default(none) &
       !$omp shared(r,rinv,InSep,dpdp,dgdp,Jp,Jg) private(i,j)
       !$omp do reduction(+:Jp,Jg)
       do j=1, Nr
          do i=1, Nz
             Jp = Jp +    r(j) * dpdp(i,j) * InSep(i,j)
             Jg = Jg + rinv(j) * dgdp(i,j) * InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       ! -- [2-B-2] alpha_g  -- !
       numer    = ( delta_psi * Itotal / ( dr*dz ) - Jp )
       denom    = 0.5d0 * Jg
       alpha_g  = numer / denom
       if ( alpha_g .lt. -0.d0 ) then
          ! if ( alpha_g .lt. 0.d0 ) then
          write(6,*)
          write(6,*) ' [CAUTION] alpha_g is less than -1.0 [CAUTION] '
          write(6,*) '               alpha_g === ', alpha_g
          write(6,*) ' See Lab. Note :: Too weak Total Current to confine pressure :: Itotal :: '
          write(6,*) '                         ---> Adjust the poloidal Beta. '
          write(6,*)
          stop
       endif
       ! -- [2-B-3] prs, gfunc -- !
       factor1  =                   dpsiinv
       factor2  = 0.5d0 * alpha_g * dpsiinv
       factor3  = sqrt( alpha_g )
       !$omp parallel default(none) &
       !$omp shared (r,rinv,InSep,Jphi,prs,gfunc,phat,ghat,dpdp,dgdp,factor1,factor2,factor3) &
       !$omp private(i,j,Jp,Jg)
       !$omp do
       do j=1, Nr
          do i=1, Nz
             Jp         = factor1 * r(j)    * dpdp(i,j) * InSep(i,j)
             Jg         = factor2 * rinv(j) * dgdp(i,j) * InSep(i,j)
             
             Jphi (i,j) = Jp  +  Jg
             prs  (i,j) =             phat(i,j) * InSep(i,j)
             gfunc(i,j) = factor3   * ghat(i,j) * InSep(i,j)
             ! else
             !    Jphi(i,j)   = 0.d0
             !    prs(i,j)    = 0.d0
             !    gfunc(i,j)  = 0.d0
             ! endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif
    !   -- [2-C] FRC --  !
    if ( equtype.eq.'FRC' ) then
       !   -- [2-C-1] Jp, Jg  -- !
       Jp  = 0.d0
       !$omp parallel default(none) &
       !$omp shared(Jp,r,dpdp,InSep) private(i,j)
       !$omp do reduction(+:Jp)
       do j=1, Nr
          do i=1, Nz
             Jp = Jp + ( -1.d0 )*r(j)*dpdp(i,j)*InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel

       ! - current conservation - !
       p0      = delta_psi * Itotal / ( dr*dz * Jp )
       ! - pressure top fixed   - !
       ! if ( normInGS.eq.'MHD' ) p0 = 1.d0
       ! if ( normInGS.eq.'PIC' ) p0 = ( 1.d0 + TiTe ) * vthcv**2
       !   -- [2-C-2] Jphi, g(psi), p(psi) --   !
       factor1 = p0 / delta_psi
       !$omp parallel default(none) &
       !$omp shared(r,Jphi,prs,gfunc,p0,phat,dpdp,factor1,InSep) private(i,j,Jp,Jg)
       !$omp do
       do j=1, Nr
          do i=1, Nz
             Jp          = - factor1*r(j)*dpdp(i,j)
             Jg          = 0.d0
             Jphi(i,j)   = Jp + Jg
             prs(i,j)    = p0 * phat(i,j)
             gfunc(i,j)  = 0.d0
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif

    !  --- [3] Post Process  --- !
    !   -- [3-1] R.H.S. --   !
    if ( normInGS.eq.'MHD' ) coef = - 1.d0
    if ( normInGS.eq.'PIC' ) coef = - 1.d0 / valfe**2
    ! if ( normInGS.eq.'PIC' ) coef = 1.d0
    !$omp parallel default(none) &
    !$omp shared(rhs,r,Jphi,coef) private(i,j)
    !$omp do
    do j=1, Nr
       do i=1, Nz
          rhs(i,j) =  coef * r(j) * Jphi(i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !   -- [3-2] Boundary --   !
    if ( solver.eq.'BCG' ) then
       !$omp parallel default(none) &
       !$omp shared(rhs,BdrcLR,BdrcTB) private(i,j)
       !$omp do
       do j=1, Nr
          rhs( 1,j) = BdrcLR(j,1)
          rhs(Nz,j) = BdrcLR(j,2)
       enddo
       !$omp end do
       !$omp do
       do i=1, Nz
          rhs(i,Nr) = BdrcTB(i,1)
          rhs(i, 1) = BdrcTB(i,2)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    !   -- [3-3] Itotal check --  !
    Isum = 0.d0
    !$omp parallel default(none) &
    !$omp shared(Jphi,Isum,dr,dz) private(i,j) 
    !$omp do reduction(+:Isum)
    do j=1, Nr
       do i=1, Nz
          Isum = Isum + Jphi(i,j)*dr*dz
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( myRank.eq.0 ) then
       write(6,'(4x,a)') '[ Update Current    ]'
       write(6,'(8x,a,e12.5)') 'p0         = ', p0
       write(6,'(8x,a,e12.5)') 'alpha_g    = ', alpha_g
       write(6,'(8x,a,e12.5)') 'I_total    = ', Isum
    endif
    
    return
  end subroutine UpdateJphi


  subroutine PicardIteration
    use constants, only : Nr, Nz, picard_alphaB
    use variables, only : psi, psipast
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2

    factor1 =        picard_alphaB
    factor2 = 1.d0 - picard_alphaB
    !$omp parallel default(none) &
    !$omp shared(psi,psipast,factor1,factor2) private(i,j)
    !$omp do
    do j=1, Nr
       do i=1, Nz
          psi(i,j) = factor1*psi(i,j) + factor2*psipast(i,j,1)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine PicardIteration


  subroutine WriteLog( iter1, iter2, diff1, diff2 )
    use variables, only           : zaxis, raxis, g0, alpha_g, psilim, psimin
    use constants, only           : Bv0, jobdir, myRank
    implicit none
    integer         , intent(in) :: iter1, iter2
    double precision, intent(in) :: diff1, diff2
    double precision             :: delta_psi
    
    delta_psi = psilim - psimin
    if ( myRank.eq.0 ) then
       open(50,file=trim(jobdir)//'convergence.log',form='formatted',access='append')
       write(50,'(2(i4,1x),10(e12.5,1x))') iter1, iter2, diff1, diff2, raxis, zaxis &
            &                              , g0, alpha_g, Bv0, psilim, psimin, delta_psi
       close(50)
    endif    
    return
  end subroutine WriteLog

  
end module gs_solver
