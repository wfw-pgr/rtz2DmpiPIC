module currentFnMod
contains

  ! ====================================================== !
  ! === update current distribution                    === !
  ! ====================================================== !
  subroutine update__currentDistribution
    !
    !! ------------------------------------------------------------------------------ !!
    !
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
    ! integral over [ LJ, LI ]
    !
    !! ------------------------------------------------------------------------------ !!
    !
    use variablesMod
    use psiGSFuncMod, only : pghat
    implicit none
    integer            :: i, j
    double precision   :: g2(LI,LJ), p0, InSep(LI,LJ)
    double precision   :: dpdp(LI,LJ), dgdp(LI,LJ), phat(LI,LJ), ghat(LI,LJ)
    double precision   :: Jp, Jg, denom, numer, delta_psi, dpsiInv, Isum
    double precision   :: factor1, factor2, factor3, coef

    ! ------------------------------------------------------ !
    ! --- [1] zero clear of variables                    --- !
    ! ------------------------------------------------------ !
    !$omp parallel default(none) &
    !$omp shared(LI,LJ,current,pressure,gfunc,psis,phat,ghat,dpdp,dgdp,InSep) private(i,j)
    !$omp do
    do j=1, LJ
       do i=1, LI
          current (i,j) = 0.d0
          pressure(i,j) = 0.d0
          gfunc   (i,j) = 0.d0
          psis    (i,j) = 0.d0
          InSep   (i,j) = 0.d0
          phat    (i,j) = 0.d0
          ghat    (i,j) = 0.d0
          dpdp    (i,j) = 0.d0
          dgdp    (i,j) = 0.d0
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    Jp         = 0.d0
    Jg         = 0.d0
    delta_psi  = 0.d0
    alpha_g    = 0.d0

    ! ------------------------------------------------------ !
    ! --- [2] calculate normalized psi                   --- !
    ! ------------------------------------------------------ !
    p0         = 1.d0
    g0         = sqrt( 2.d0 / beta_tor ) * magAxis(r_)
    delta_psi  = ( psiLim - psiMin )
    dpsiInv    = 1.d0 / delta_psi
    !$omp parallel default(none) &
    !$omp shared(LI,LJ,psis,psiLim,psi,dpsiInv) private(i,j)
    !$omp do
    do j=1, LJ
       do i=1, LI
          psis(i,j) = ( psiLim - psi(i,j) ) * dpsiInv
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! ------------------------------------------------------ !
    ! --- [3] p(psi), g(psi)                             --- !
    ! ------------------------------------------------------ !
    call pghat( psis, phat, ghat, dpdp, dgdp )
    !$omp parallel default(none) &
    !$omp shared(LI,LJ,psis,InSep) private(i,j)
    !$omp do
    do j=1, LJ
       do i=1, LI
          if ( psis(i,j).gt.0.d0 ) then
             InSep(i,j) = 1.d0
          endif
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! ------------------------------------------------------ !
    ! --- [4] update R.H.S. (ST)                         --- !
    ! ------------------------------------------------------ !
    if ( trim(equiliType) == 'STk' ) then
       
       ! ------------------------------------------------------ !
       ! --- [4-1] current density                          --- !
       ! ------------------------------------------------------ !
       Jg  = 0.d0
       Jp  = 0.d0
       !$omp parallel default(none) &
       !$omp shared(LI,LJ,rAxis,rInv,dpdp,dgdp,psis,InSep,Jg,Jp) private(i,j)
       !$omp do reduction(+:Jg,Jp)
       do j=1, LJ
          do i=1, LI
             Jp = Jp + rAxis(j) * dpdp(i,j) * InSep(i,j)
             Jg = Jg +  rInv(j) * dgdp(i,j) * InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       
       ! ------------------------------------------------------ !
       ! --- [4-2] alpha_g                                  --- !
       ! ------------------------------------------------------ !
       numer   = ( delta_psi * Itot / ( dr*dz ) - Jp )
       denom   = 0.5d0 * g0**2 * Jg
       alpha_g = numer / denom
       if ( alpha_g .lt. -1.d0 ) then
          ! if ( alpha_g .lt. 0.d0 ) then
          write(6,*)
          write(6,*) ' [CAUTION] alpha_g is less than -1.0 [CAUTION] '
          write(6,*) '               alpha_g === ', alpha_g
          write(6,*) ' See Lab. Note :: Too weak Total Current to confine pressure :: Itot :: '
          write(6,*) '                         ---> Adjust the poloidal Beta. '
          write(6,*)
          stop
       endif
       
       ! ------------------------------------------------------ !
       ! --- [4-3] pressure , gfunc                         --- !
       ! ------------------------------------------------------ !
       factor1 =                           dpsiInv
       factor2 = 0.5d0 * alpha_g * g0**2 * dpsiInv
       !$omp parallel default(none) &
       !$omp shared(LI,LJ,rAxis,rInv,dpdp,dgdp,current,pressure,gfunc,phat,ghat) &
       !$omp shared(g0,alpha_g,InSep,factor1,factor2) &
       !$omp private(i,j,Jp,Jg)
       !$omp do
       do j=1, LJ
          do i=1, LI
             Jp            = factor1 * rAxis(j) * dpdp(i,j) * InSep(i,j)
             Jg            = factor2 * rInv(j)  * dgdp(i,j) * InSep(i,j)
             current(i,j)  = Jp  +  Jg
             pressure(i,j) = phat(i,j) * InSep(i,j)
             gfunc(i,j)    = g0 * sqrt( 1.d0 + alpha_g * ghat(i,j) ) * InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif

    
    ! ------------------------------------------------------ !
    ! --- [5] update R.H.S. (Spheromak)                  --- !
    ! ------------------------------------------------------ !
    if ( trim(equiliType) == 'Sph' ) then
       
       ! ------------------------------------------------------ !
       ! --- [5-1] current density                          --- !
       ! ------------------------------------------------------ !
       Jg  = 0.d0
       Jp  = 0.d0
       !$omp parallel default(none) &
       !$omp shared(LI,LJ,rAxis,rInv,InSep,dpdp,dgdp,Jp,Jg) private(i,j)
       !$omp do reduction(+:Jp,Jg)
       do j=1, LJ
          do i=1, LI
             Jp = Jp +    rAxis(j) * dpdp(i,j) * InSep(i,j)
             Jg = Jg + rInv(j) * dgdp(i,j) * InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       
       ! ------------------------------------------------------ !
       ! --- [5-2] alpha_g                                  --- !
       ! ------------------------------------------------------ !
       numer    = ( delta_psi * Itot / ( dr*dz ) - Jp )
       denom    = 0.5d0 * Jg
       alpha_g  = numer / denom
       if ( alpha_g .lt. -0.d0 ) then
          ! if ( alpha_g .lt. 0.d0 ) then
          write(6,*)
          write(6,*) ' [CAUTION] alpha_g is less than -1.0 [CAUTION] '
          write(6,*) '               alpha_g === ', alpha_g
          write(6,*) ' See Lab. Note :: Too weak Total Current to confine pressure :: Itot :: '
          write(6,*) '                         ---> Adjust the poloidal Beta. '
          write(6,*)
          stop
       endif

       ! ------------------------------------------------------ !
       ! --- [5-3] pressure , gfunc                         --- !
       ! ------------------------------------------------------ !
       factor1  =                   dpsiInv
       factor2  = 0.5d0 * alpha_g * dpsiInv
       factor3  = sqrt( alpha_g )
       !$omp parallel default(none) &
       !$omp shared (LI,LJ,rAxis,rInv,InSep,current,pressure,gfunc,phat,ghat) &
       !$omp shared (dpdp,dgdp,factor1,factor2,factor3) &
       !$omp private(i,j,Jp,Jg)
       !$omp do
       do j=1, LJ
          do i=1, LI
             Jp            = factor1 * rAxis(j) * dpdp(i,j) * InSep(i,j)
             Jg            = factor2 * rInv(j)  * dgdp(i,j) * InSep(i,j)
             current (i,j) = Jp  +  Jg
             pressure(i,j) =             phat(i,j) * InSep(i,j)
             gfunc   (i,j) = factor3   * ghat(i,j) * InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif

    
    ! ------------------------------------------------------ !
    ! --- [6] update R.H.S. (FRC)                        --- !
    ! ------------------------------------------------------ !
    if ( trim(equiliType) == 'FRC' ) then
       
       ! ------------------------------------------------------ !
       ! --- [6-1] current density                          --- !
       ! ------------------------------------------------------ !
       Jp  = 0.d0
       !$omp parallel default(none) &
       !$omp shared(LI,LJ,Jp,rAxis,dpdp,InSep) private(i,j)
       !$omp do reduction(+:Jp)
       do j=1, LJ
          do i=1, LI
             Jp = Jp + ( -1.d0 )*rAxis(j)*dpdp(i,j)*InSep(i,j)
          enddo
       enddo
       !$omp end do
       !$omp end parallel

       ! ------------------------------------------------------ !
       ! --- [6-2] current conservation                     --- !
       ! ------------------------------------------------------ !
       ! - current conservation - !
       p0      = delta_psi * Itot / ( dr*dz * Jp )
       ! - pressure top fixed   - !
       ! if ( normInGS == 'MHD' ) p0 = 1.d0
       ! if ( normInGS == 'PIC' ) p0 = ( 1.d0 + TiTe ) * vthcv**2

       ! ------------------------------------------------------ !
       ! --- [6-3] current, pressure, gfunc                 --- !
       ! ------------------------------------------------------ !
       factor1 = p0 / delta_psi
       !$omp parallel default(none) &
       !$omp shared(LI,LJ,rAxis,current,pressure,gfunc,p0,phat,dpdp,factor1,InSep) private(i,j,Jp,Jg)
       !$omp do
       do j=1, LJ
          do i=1, LI
             Jp            = - factor1*rAxis(j)*dpdp(i,j)
             Jg            = 0.d0
             current(i,j)  = Jp + Jg
             pressure(i,j) = p0 * phat(i,j)
             gfunc(i,j)    = 0.d0
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif

    
    ! ------------------------------------------------------ !
    ! --- [7] post Process 1 (r.h.s.)                    --- !
    ! ------------------------------------------------------ !
    if ( trim(normalizationInGS) == "MHD" ) coef = - 1.d0
    if ( trim(normalizationInGS) == "PIC" ) coef = - 1.d0 / valfe**2
    ! if ( trim(normalizationInGS) == "PIC" ) coef =   1.d0
    !$omp parallel default(none) &
    !$omp shared(LI,LJ,rhs,rAxis,current,coef) private(i,j)
    !$omp do
    do j=1, LJ
       do i=1, LI
          rhs(i,j) =  coef * rAxis(j) * current(i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ------------------------------------------------------ !
    ! --- [8] post Process 2 (boundary)                  --- !
    ! ------------------------------------------------------ !
    if ( trim(solverType) == "BCG" ) then
       !$omp parallel default(none) &
       !$omp shared(LI,LJ,rhs,bdrcLR,bdrcTB) private(i,j)
       !$omp do
       do j=1, LJ
          rhs( 1,j) = bdrcLR(j,lft_)
          rhs(LI,j) = bdrcLR(j,rgt_)
       enddo
       !$omp end do
       !$omp do
       do i=1, LI
          rhs(i,LJ) = bdrcTB(i,top_)
          rhs(i, 1) = bdrcTB(i,bot_)
       enddo
       !$omp end do
       !$omp end parallel
    endif
    
    ! ------------------------------------------------------ !
    ! --- [9] post Process 3 ( check Itot )              --- !
    ! ------------------------------------------------------ !
    Isum = 0.d0
    !$omp parallel default(none) &
    !$omp shared(LI,LJ,current,Isum,dr,dz) private(i,j) 
    !$omp do reduction(+:Isum)
    do j=1, LJ
       do i=1, LI
          Isum = Isum + current(i,j)*dr*dz
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( myRank == 0 ) then
       write(6,'(4x,a)') '[ Update Current    ]'
       write(6,'(8x,a,e12.5)') 'p0         = ', p0
       write(6,'(8x,a,e12.5)') 'alpha_g    = ', alpha_g
       write(6,'(8x,a,e12.5)') 'Itotal     = ', Isum
    endif
    
    return
  end subroutine update__currentDistribution

  
  ! ====================================================== !
  ! === set__initial_currentDistribution               === !
  ! ====================================================== !
  subroutine set__initial_currentDistribution( mode )
    use variablesMod
    use saveFieldMod
    use diagnosisMod, only : measure__perimeter
    implicit none
    character(8), intent(in) :: mode        ! -- [ arbitral:1st, pressure:2nd ] -- !
    integer                  :: i, j
    double precision         :: rCent, zCent, aMinorInv, radius
    double precision         :: Ptot , Stot, Jtot, dS, pLength, dpsiInv, coef

    ! ------------------------------------------------------ !
    ! --- [1] (trial) arbitral cylindrical current       --- !
    ! ------------------------------------------------------ !
    if ( trim(mode) == "arbitral" ) then

       ! ------------------------------------------------------ !
       ! --- [1-1] consider to set cylindrical current      --- !
       ! ------------------------------------------------------ !
       aMinor  = min( 0.5d0*(rMax-rMin), 0.5d0*(zMax-zMin) ) * 0.5d0
       rCent   = 0.5d0 * ( rMax + rMin )
       zCent   = 0.5d0 * ( zMax + zMin )
       if ( boundaryType == 'hlf' ) zCent = 0.d0

       ! ------------------------------------------------------ !
       ! --- [1-2] set cylindrical current dist.            --- !
       ! ------------------------------------------------------ !
       aMinorInv = 1.d0 / aMinor
       do j=1, LJ
          do i=1, LI
             radius = sqrt( ( rAxis(j) - rCent )**2 + ( zAxis(i) - zCent )**2 )
             if ( radius < aMinor ) then
                current (i,j) = ( 1.d0 - radius**2/aMinor**2 ) * rInv(j)
                psis    (i,j) = ( aMinor - radius ) * aMinorInv
                pressure(i,j) = pTop * psis(i,j)     ! -- Assumption : p(psi) = p0 psi  -- !
             else
                current (i,j) =   0.d0
                psis    (i,j) =   0.d0
             endif
          enddo
       enddo
    endif
    
    ! ------------------------------------------------------ !
    ! --- [2] normalization of psi ( psis )              --- !
    ! ------------------------------------------------------ !
    if ( trim(mode) == "pressure" ) then
       dpsiInv =  1.d0 / ( psiMin - psiLim )
       do j=1, LJ
          do i=1, LI
             psis(i,j) = ( psi(i,j) - psiLim ) * dpsiInv
          enddo
       enddo
    endif

    ! ------------------------------------------------------ !
    ! --- [3] area integration of total Pressure         --- !
    ! ------------------------------------------------------ !
    Ptot   = 0.d0
    Stot   = 0.d0
    dS     = dr*dz
    do j=1, LJ
       do i=1, LI
          if ( psis(i,j) > 0.d0 ) then
             Ptot  = Ptot + dS*pressure(i,j)
             Stot  = Stot + dS
          endif
       enddo
    enddo
    Ptot  = Ptot / Stot
    
    ! ------------------------------------------------------ !
    ! --- [4] peripheral length                          --- !
    ! ------------------------------------------------------ !
    if      ( trim(mode) == "pressure" ) then
       pLength = measure__perimeter( psis, LI, LJ, dz, dr )
       
    else if ( trim(mode) == "arbitral"   ) then
       pLength = 2.d0 * pi * aMinor
       
    endif

    ! ------------------------------------------------------ !
    ! --- [5] count up all current :: Jtot               --- !
    ! ------------------------------------------------------ !
    Jtot = 0.d0
    do j=1, LJ
       do i=1, LI
          Jtot = Jtot + current(i,j)*dS
       enddo
    enddo
    
    ! ------------------------------------------------------ !
    ! --- [6] Define Itot                                --- !
    ! ------------------------------------------------------ !
    if ( trim(normalizationInGS) == "MHD" ) then
       coef =  - 1.d0
       ! coef =  + 1.d0
       Itot =  pLength * sqrt( 2.d0 * Ptot / beta_pol )
    endif
    if ( trim(normalizationInGS) == "PIC" ) then
       coef =  - 1.d0 / valfe**2
       ! coef =  + 1.d0 / valfe**2
       Itot =  pLength * sqrt( 2.d0 * Ptot / beta_pol ) * valfe
       ! Itot =  pLength * sqrt( 2.d0 * Ptot / beta_pol )
    endif

    ! ------------------------------------------------------ !
    ! --- [7] normalization by Itot                      --- !
    ! ------------------------------------------------------ !
    Jtot   = Itot / Jtot
    do j=1, LJ
       do i=1, LI
          current(i,j) =            Jtot * current(i,j)
          rhs    (i,j) = coef * rAxis(j) * current(i,j)
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [8] display of the initial guess equ.          --- !
    ! ------------------------------------------------------ !
    if ( myRank == 0 ) then
       write(6,*)
       write(6,'(6x,a)'      ) '  --- Initial Guess by Cylindrical Tube Current --- '
       write(6,'(6x,a,a10  )') ' mode      :: ', trim( mode )
       write(6,'(6x,a,f12.7)') ' Beta_pol  :: ', beta_pol
       write(6,'(6x,a,f12.7)') ' pLength   :: ', pLength
       write(6,'(6x,a,f12.7)') ' Iplasma   :: ', Itot
       write(6,'(6x,a,f12.7)') ' Area      :: ', Stot
       write(6,'(6x,a,f12.7)') ' Paverage  :: ', Ptot
       write(6,'(6x,a)'      ) '  ------------------------------------------------- '
       write(6,*)
    endif

    return
  end subroutine set__initial_currentDistribution
  
end module currentFnMod




  ! ! ====================================================== !
  ! ! === set__initial_currentDistribution               === !
  ! ! ====================================================== !
  ! subroutine set__initial_currentDistribution
  !   use variablesMod
  !   implicit none
  !   integer            :: i, j
  !   double precision   :: Ptot, psihat, area, dS, coef
  !   double precision   :: radius, rCent, zCent, Jtot, aMinorInv

  !   ! ------------------------------------------------------ !
  !   ! --- [4] count up all current :: Jtot               --- !
  !   ! ------------------------------------------------------ !
  !   Jtot = 0.d0
  !   do j=1, LJ
  !      do i=1, LI
  !         Jtot = Jtot + current(i,j)*dS
  !      enddo
  !   enddo
    
  !   ! ------------------------------------------------------ !
  !   ! --- [5] Define Itot                                --- !
  !   ! ------------------------------------------------------ !
  !   if ( trim(normalizationInGS) == 'MHD' ) then
  !      coef =  - 1.d0
  !      Itot =    2.d0*pi*aMinor * sqrt( 2.d0 * Ptot / beta_pol )
  !   endif
  !   if ( trim(normalizationInGS) == 'PIC' ) then
  !      coef =  - 1.d0 / valfe**2
  !      Itot = valfe * 2.d0*pi*aMinor * sqrt( 2.d0 * Ptot / beta_pol )
  !      ! Itot = 2.d0*pi*aMinor * sqrt( 2.d0 * Ptot / Betap )
  !   endif

  !   ! ------------------------------------------------------ !
  !   ! --- [6] normalization by Itot                      --- !
  !   ! ------------------------------------------------------ !
  !   Jtot   = Itot / Jtot
  !   do j=1, LJ
  !      do i=1, LI
  !         current(i,j) = Jtot * current(i,j)
  !         rhs (i,j) = coef * rAxis(j) * current(i,j)
  !      enddo
  !   enddo

  !   ! ------------------------------------------------------ !
  !   ! --- [7] display of the initial guess equ.          --- !
  !   ! ------------------------------------------------------ !
  !   if ( myRank == 0 ) then
  !      write(6,*)
  !      write(6,'(6x,a)'      ) '  --- Initial Guess by Cylindrical Tube Current --- '
  !      write(6,'(6x,a,f12.7)') ' Beta_p  :: ', beta_pol
  !      write(6,'(6x,a,f12.7)') ' Iplasma :: ', Itot
  !      write(6,'(6x,a,f12.7)') ' Area    :: ', area
  !      write(6,'(6x,a,f12.7)') ' p_ave   :: ', Ptot
  !      write(6,'(6x,a)'      ) '  ------------------------------------------------- '
  !      write(6,*)
  !   endif

  !   return
  ! end subroutine set__initial_currentDistribution

