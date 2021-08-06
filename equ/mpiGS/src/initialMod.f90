module initialMod
contains

  subroutine InitCond

    use constants  , only : Nz, Nr, rMin
    use variables  , only : psi, rhs, dz, dr
    use gs_solver  , only : DeterminePsi, UpdateJphi
    use myPBiCGSTAB, only : BiCGSTABCyl2D_MPI
    implicit none
    integer :: i,j 
    
    ! --- [1] Initialize Variables           --- !
    call InitVariables
    ! --- [2] Try to Find Approx.GS-Solution --- !
    call SetTrialCurrent
    call SetInitialBC
    call BiCGSTABCyl2D_MPI( psi, rhs, dz, dr, Nz, Nr, rMin )
    call DeterminePsi
    call UpdateJphi
    ! --- [3] Determine Initial Condition    --- !
    call SetInitialCurrent
    call SetInitialBC
    
    return
  end subroutine InitCond

  
  subroutine InitVariables

    use constants, only : Nr, Nz, Nmid, myRank, normType
    use constants, only : vthcv, wpewce, dr_Debye, dz_Debye
    use constants, only : rMax, rMin, zMax, zMin, Bv0, AspectR
    use constants, only : Nlim, limposr, limposz, Bsensor_r1, Bsensor_r2
    use variables, only : dr, dz, drinv, dzinv, r, z, rinv
    use variables, only : raxis, zaxis
    use variables, only : psi, psipast, rhs
    use variables, only : Jphi, prs, gfunc
    use variables, only : psilim, psimin
    implicit none
    integer            :: i, j
    double precision   :: lDebye, rleng, zleng

    ! --- [1] Initialize x1, x2 grid --- !
    !  -- [1-1] Grid Parameters  --  !
    if ( normType.eq.'PIC' ) then
       lDebye     =    vthcv / wpewce
       dr         = dr_Debye * lDebye
       dz         = dz_Debye * lDebye
       drinv      = 1.d0 / dr
       dzinv      = 1.d0 / dz
       rleng      = dble( Nr-1 ) * dr
       zleng      = dble( Nz-1 ) * dz
       if ( AspectR.gt.1.d0 ) then
          rMin    =   rleng * ( AspectR - 1.d0 ) * ( Bsensor_r1 + Bsensor_r2 )*0.5d0
       else
          rMin    = 0.d0
       endif
       zMin       = - dz*dble(Nz-1) * 0.5d0
       zMax       = + dz*dble(Nz-1) * 0.5d0
       rMax       = + dr*dble(Nr-1) + rMin
    endif
    if ( normType.eq.'MHD' ) then
       rleng      = rMax - rMin
       zleng      = zMax - zMin
       dz         = ( zMax - zMin ) / dble( Nz-1 )
       dr         = ( rMax - rMin ) / dble( Nr-1 )
       drinv      = 1.d0 / dr
       dzinv      = 1.d0 / dz
    endif
    do i=1, Nz
       z(i) = dz*dble( i-1) + zMin
    enddo
    do j=1, Nr
       r(j) = dr*dble( j-1) + rMin
    enddo

    
    !  --- [1] Grid Initialize    --- !
    !   -- [1-2] rMin / zMin -- !
    !   -- [1-3]  r/z grid   -- !
    do j=1, Nr
       r(j)    = dr*dble(j-1) + rMin
       if ( r(j).ne.0.d0 ) then
          rinv(j) = 1.d0 / r(j)
       else
          rinv(j) = 0.d0
       endif
    enddo
    do i=1, Nz
       z(i)    = dz*dble(i-1) + zMin
    enddo
    rMax       = r(Nr)
    zMax       = z(Nz)
    raxis      = 0.5d0 * ( rMax + rMin )
    zaxis      = 0.5d0 * ( zMax + zMin )
    !   -- [1-4]  Limiter & Sensor Position --  !
    Bsensor_r1 = Bsensor_r1 * rleng + rMin
    Bsensor_r2 = Bsensor_r2 * rleng + rMin
    do i=1, Nlim
       limposr(i) = limposr(i) * rleng + rMin
       limposz(i) = limposz(i) * zleng + zMin
    enddo

    !  --- [2] Display Settings --- !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*)
       write(6,*) '  -------------------------------------------  '
       write(6,*) '  ---------    Grid   Setting     -----------  '
       write(6,*) '  -------------------------------------------  '
       write(6,*)
       write(6,'(6x,3(a8,1x))'            ) 'Nr', 'Nz', 'Nmid'
       write(6,'(6x,3(i8,1x))'            )  Nr ,  Nz ,  Nmid
       write(6,*)
       write(6,'(6x,a48)'                 ) '  rMin  ---  rMax    ||     zMin   ---   zMax'
       write(6,'(6x,2(2(f10.5,2x),4x))'   )  rMin, rMax, zMin, zMax
       write(6,*)
       write(6,'(6x,2(a10  ,2x))'         ) 'rleng', 'zleng'
       write(6,'(6x,2(f10.5,2x))'         )  rleng ,  zleng
       write(6,*)
       write(6,'(6x,2(1(a10  ,2x),4x))'   ) 'limposr', 'limposz'
       do i=1, Nlim
          write(6,'(6x,2(1(f10.5,2x),4x))')  limposr(i),  limposz(i)
       enddo
       write(6,*)
       write(6,'(6x,2(1(a10  ,2x),4x))'   ) 'Bsensor_r1', 'Bsensor_r2'
       write(6,'(6x,2(1(f10.5,2x),4x))'   )  Bsensor_r1 ,  Bsensor_r2
       write(6,*)
       write(6,*)
    endif
    
    !  --- [3] Initialize variables --- !
    psilim  = 0.d0
    psimin  = 0.d0
    do j=1, Nr
       do i=1, Nz
          psi(i,j)       = 0.d0
          psipast(i,j,1) = 0.d0
          psipast(i,j,2) = 0.d0
          rhs(i,j)       = 0.d0
          Jphi(i,j)      = 0.d0
          prs(i,j)       = 0.d0
          gfunc(i,j)     = 0.d0
       enddo
    enddo
    
    return
  end subroutine InitVariables


  subroutine SetTrialCurrent
    
    use constants, only : Nr, Nz, bctype, myRank, normInGS
    use constants, only : rMax, rMin, zMax, zMin
    use constants, only : Betap, Bv0, valfe, pTop
    use variables, only : r, z, dr, dz, aMinor
    use variables, only : Jphi, rhs, prs, psi, Itotal
    implicit none
    integer            :: i, j
    double precision   :: prs_t, psihat, prs_h, area, dS, coef
    double precision   :: radius, rcent, zcent, Jtot, aMinorinv
    double precision, parameter :: pi = 4.d0 * atan( 1.d0 )
    
    !  --- [1]   Initial guess : circular shape current distribution  --- !
    !   -- [1-1] Geometry -- !
    aMinor  = min( 0.5d0*(rMax-rMin), 0.5d0*(zMax-zMin) ) * 0.5d0
    rcent   = 0.5d0 * ( rMax + rMin )
    zcent   = 0.5d0 * ( zMax + zMin )
    if ( bctype.eq.'hlf' ) zcent = 0.d0
    !  -- [1-2] Set Arbitral Current --  !
    do j=1, Nr
       do i=1, Nz
          radius = sqrt( ( r(j) - rcent )**2 + ( z(i) - zcent )**2 )
          if ( radius.lt.aMinor ) then
             Jphi(i,j) = ( 1.d0 - radius**2/aMinor**2 ) / r(j)
          else
             Jphi(i,j) = 0.d0
          endif
       enddo
    enddo
    !  -- [1-3] Area Integrated Pressure  --  !
    prs_t     = 0.d0
    area      = 0.d0
    aMinorinv = 1.d0 / aMinor
    dS        = dr*dz
    do j=1, Nr
       do i=1, Nz
          radius    = sqrt( ( r(j) - rcent )**2 + ( z(i) - zcent )**2 )
          if ( radius.lt.aMinor ) then
             psihat = ( aMinor - radius ) * aMinorinv
             prs_h  = pTop  * psihat           ! -- Assumption : p(psi) = p0 psi -- !
             prs_t  = prs_t + prs_h*dS
             area   = area  + dS
          endif
       enddo
    enddo
    prs_t     = prs_t / area
    
    !  --- [2] Total Current is given value :: Jtot => Itotal  ---  !
    !   -- [2-1] Sum up Current --  !
    Jtot = 0.d0
    do j=1, Nr
       do i=1, Nz
          Jtot = Jtot + Jphi(i,j)*dS
       enddo
    enddo
    !   -- [3-2] Define Itotal and Adjust to it --  !
    if ( normInGS.eq.'MHD' ) then
       Itotal =         2.d0*pi*aMinor * sqrt( 2.d0 * prs_t / Betap )
       coef   =  - 1.d0
    endif
    if ( normInGS.eq.'PIC' ) then
       Itotal = valfe * 2.d0*pi*aMinor * sqrt( 2.d0 * prs_t / Betap )
       ! Itotal = 2.d0*pi*aMinor * sqrt( 2.d0 * prs_t / Betap )
       coef   =  - 1.d0 / valfe**2
    endif
    
    Jtot   = Itotal / Jtot
    do j=1, Nr
       do i=1, Nz
          Jphi(i,j) = Jtot * Jphi(i,j)
          rhs(i,j)  = coef * r(j) * Jphi(i,j)
       enddo
    enddo

    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'(6x,a)'      ) '  --- Initial Guess by Cylindrical Tube Current --- '
       write(6,'(6x,a,f12.7)') ' Beta_p  :: ', betap
       write(6,'(6x,a,f12.7)') ' Iplasma :: ', Itotal
       write(6,'(6x,a,f12.7)') ' Area    :: ', area
       write(6,'(6x,a,f12.7)') ' prs_ave :: ', prs_t
       write(6,'(6x,a)'      ) '  ------------------------------------------------- '
       write(6,*)
    endif
    
    return
  end subroutine SetTrialCurrent


  subroutine SetInitialCurrent
    
    use constants, only : Nr, Nz, solver, bctype, myRank, normInGS
    use constants, only : rMax, rMin, zMax, zMin, Bv0
    use constants, only : Betap, valfe
    use variables, only : r, z, dr, dz, aMinor
    use variables, only : Jphi, rhs, prs, psi, psilim, psimin, Itotal
    use diagnosis, only : Perimeter
    implicit none
    integer            :: i, j
    double precision   :: prs_t, area, dS, pLength, Jtot, dpsiinv, coef
    double precision   :: psis(Nz,Nr)
    
    !  --- [1] Initial Current : From Trial GSS-Estimation  --- !
    !   -- [1-1] Normalized psi --  !
    dpsiinv =  1.d0 / ( psilim - psimin )
    do j=1, Nr
       do i=1, Nz
          psis(i,j) = ( psilim - psi(i,j) ) * dpsiinv
       enddo
    enddo
    !   -- [1-2] Total Pressure --  !
    prs_t  = 0.d0
    area   = 0.d0
    dS     = dr*dz
    do j=1, Nr
       do i=1, Nz
          if ( psis(i,j).gt.0.d0 ) then
             prs_t = prs_t + prs(i,j)*dS
             area  = area  + dS
          endif
       enddo
    enddo
    prs_t  = prs_t / area
    ! - [1-2] Peripheral Length - !
    pLength = Perimeter( psis, Nz, Nr, dz, dr )
    
    !  --- [2] Total Current is given value :: Jtot => Itotal  ---  !
    !   -- [2-1] Sum up Current --  !
    Jtot = 0.d0
    do j=1, Nr
       do i=1, Nz
          Jtot = Jtot + Jphi(i,j)*dS
       enddo
    enddo
    !   -- [2-2] Define Itotal and Adjust to it --  !
    if ( normInGS.eq.'MHD' ) then
       Itotal =       pLength * sqrt( 2.d0 * prs_t / Betap )
       coef   =  - 1.d0
    endif
    if ( normInGS.eq.'PIC' ) then
       Itotal = valfe*pLength * sqrt( 2.d0 * prs_t / Betap )
       coef   =  - 1.d0 / valfe**2
    endif
    Jtot   = Itotal / Jtot
    do j=1, Nr
       do i=1, Nz
          Jphi(i,j) = Jtot * Jphi(i,j)
          rhs(i,j)  = coef * r(j) * Jphi(i,j)
       enddo
    enddo

    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'(6x,a,f12.7)') ' Beta_p  :: ', betap
       write(6,'(6x,a,f12.7)') ' pLength :: ', pLength
       write(6,'(6x,a,f12.7)') ' Iplasma :: ', Itotal
       write(6,'(6x,a,f12.7)') ' Area    :: ', area
       write(6,'(6x,a,f12.7)') ' prs_ave :: ', prs_t
       write(6,*)
    endif
    return
  end subroutine SetInitialCurrent

  
  subroutine SetInitialBC
    use constants, only : Nr, Nz, Bv0,   solver, betap, bctype, myRank, normInGS, valfe
    use variables, only : r, rhs, BdrcLR, BdrcTB,  Itotal, raxis, aMinor
    implicit none
    integer                     :: i, j
    double precision            :: coef
    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    
    !  --- [1] Initial Guess of Bv0  ---  !
    if ( normInGS.eq.'MHD' ) &
         & Bv0 =  - Itotal / (            4.d0 * pi * raxis ) * ( log( 8.d0*raxis / aMinor ) - 1.5d0 + betap )
    if ( normInGS.eq.'PIC' ) &
         & Bv0 =  - Itotal / ( valfe**2 * 4.d0 * pi * raxis ) * ( log( 8.d0*raxis / aMinor ) - 1.5d0 + betap )
    ! if ( normInGS.eq.'PIC' ) &
    !      & Bv0 =  Itotal / ( 4.d0 * pi * raxis ) * ( log( 8.d0*raxis / aMinor ) - 1.5d0 + betap )
    if ( myRank.eq.0 ) then
       write(6,'(6x,a,f12.7)') ' Itotal  :: ', Itotal
       write(6,'(6x,a,f12.7)') ' raxis   :: ', raxis
       write(6,'(6x,a,f12.7)') ' aMinor  :: ', aMinor
       write(6,'(6x,a,f12.7)') ' betap   :: ', betap
       write(6,'(6x,a,f12.7)') ' Bv0     :: ', Bv0
    endif
    
    !  --- [2] Store Boundary Condition --- !
    coef = +0.5d0
    if ( bctype.eq.'pcw' ) then
       do j=1, Nr
          BdrcLR(j,1) = coef * r( j)**2 * Bv0
          BdrcLR(j,2) = coef * r( j)**2 * Bv0
       enddo
       do i=1, Nz
          BdrcTB(i,1) = coef * r(Nr)**2 * Bv0
          BdrcTB(i,2) = coef * r( 1)**2 * Bv0
       enddo
    endif
    if ( bctype.eq.'hlf' ) then
       do j=1, Nr
          BdrcLR(j,1) = 0.0d0 
          BdrcLR(j,2) = coef * r( j)**2 * Bv0
       enddo
       do i=1, Nz
          BdrcTB(i,1) = coef * r(Nr)**2 * Bv0
          BdrcTB(i,2) = coef * r( 1)**2 * Bv0
       enddo
    endif
    if ( solver.eq.'BCG' ) then
       do j=1, Nr
          rhs( 1,j) = BdrcLR(j,1)
          rhs(Nz,j) = BdrcLR(j,2)
       enddo
       do i=1, Nz
          rhs(i,Nr) = BdrcTB(i,1)
          rhs(i, 1) = BdrcTB(i,2)
       enddo
    endif
    
    return
  end subroutine SetInitialBC
  
end module initialMod
