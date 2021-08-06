module Field
  use constants, only : Nr, Nz
  implicit none
  integer, parameter :: LIg = Nz+1
  integer, parameter :: LJg = Nr+1
  double precision   :: Br(0:Nz,0:Nr), Bz(0:Nz,0:Nr), Bt(0:Nz,0:Nr)
  double precision   :: Jr(Nz,Nr), Jz(Nz,Nr), Jt(Nz,Nr)
  double precision   :: rf(Nr), rh(0:Nr), rfinv(Nr), rhinv(0:Nr)
  double precision   :: BgB(3,LIg,LJg),  JgB(3,LIg,LJg), fgB(LIg,LJg,3), ugB(LIg,LJg,3)
contains
  
  subroutine calc_Bfield_Egrid

    use constants, only : Nr, Nz, rmin, rmax, Bmax, valfe, vthcv, TiTe, rhofloor, normSW
    use variables, only : r, z, dr, dz
    use variables, only : psi, prs, gfunc, rho
    implicit none
    integer            :: i, j
    double precision   :: dzx2, drx2, factor1, factor2
    double precision   :: drinv, dzinv, gfcnt(0:Nz,0:Nr)
    double precision   :: absBmax, prscoef, valfe2, Ttot, Ttotinv, rhof

    !  --- [0] preparation  ---  !
    !   -- Initialization   --   !
    Br(:,:) = 0.d0
    Bt(:,:) = 0.d0
    Bz(:,:) = 0.d0
    !   -- Coordinate       --   !
    drinv   = 1.d0 / dr
    dzinv   = 1.d0 / dz
    do j=1, Nr
       rf(j)    = dr * dble(j-1) + rmin
       rfinv(j) = 1.d0 / rf(j)
    enddo
    do j=1, Nr-1
       rh(j)    = 0.5d0 * ( rf(j) + rf(j+1) )
       rhinv(j) = 1.0d0 / rh(j)
    enddo
    if ( rmin.eq.0.d0 ) then
       rf(1)    = 0.d0
       rfinv(1) = 0.d0
    endif
    rh( 0)      = rmin  - 0.5d0 * dr
    rhinv( 0)   = 1.0d0 / rh( 0)
    rh(Nr)      = rf(Nr)+ 0.5d0 * dr
    rhinv(Nr)   = 1.0d0 / rh(Nr)

    !  --- [1] Interpolation --- !
    do j=1, Nr-1
       do i=1, Nz-1
          gfcnt(i,j) = 0.25d0 * ( gfunc(i,j) + gfunc(i+1,j) + gfunc(i,j+1) + gfunc(i+1,j+1) )
       enddo
    enddo
    !  --- [2] calculate B  ---  !
    do j=1, Nr-1
       do i=1, Nz-1
          Br(i,j) =   ( psi(i+1,j) - psi( i,j) ) * rfinv( j) * dzinv
          Bt(i,j) =     gfcnt(i,j) * rhinv(j)
          Bz(i,j) = - ( psi(i,j+1) - psi( i,j) ) * rhinv( j) * drinv
       enddo
    enddo
    !  --- [3] calculate B on boundary ---  !
    do j=1, Nr-1
       Bz(Nz,j) = - ( psi(Nz,j+1) - psi(Nz,j) ) * rhinv( j) * drinv
    enddo
    do i=1, Nz-1
       Br(i,Nr) = + ( psi(i+1,Nr) - psi(i,Nr) ) * rfinv(Nr) * dzinv
    enddo
    !   -- [3-1] Boundary Condition -- !
    Bz(:, 0) = Bz(:,   1)
    Bz(:,Nr) = Bz(:,Nr-1)
    Br( 0,:) = Br(   1,:)
    Br(Nz,:) = Br(Nz-1,:)
    factor1  = ( rmin + 0.5d0*dr ) / ( rmin - 0.5d0*dr )
    factor2  = ( rmax - 0.5d0*dr ) / ( rmax + 0.5d0*dr )
    Bt( 0,:) = Bt(   1,:)
    Bt(Nz,:) = Bt(Nz-1,:)
    Bt(:, 0) = Bt(:,   1) * factor1
    Bt(:,Nr) = Bt(:,Nr-1) * factor2

    !  --- [4] Adjusting to Bmax  ---  !
    if ( normSW ) then
       absBmax = 0.d0
       do j=1, Nr
          do i=1, Nz
             absBmax = max( absBmax, Br(i,j)**2 + Bt(i,j)**2 + Bz(i,j)**2 )
          enddo
       enddo
       absBmax = Bmax / sqrt( absBmax )
       prscoef = absBmax**2
       do j=0, Nr
          do i=0, Nz
             Br(i,j)  = absBmax *  Br(i,j)
             Bt(i,j)  = absBmax *  Bt(i,j)
             Bz(i,j)  = absBmax *  Bz(i,j)
          enddo
       enddo
       do j=1, Nr
          do i=1, Nz
             prs(i,j) = prscoef * prs(i,j)
          enddo
       enddo
    endif

    !  --- [5] Set prs, rho  --- !
    ! prscoef = valfe**2 * absBmax**2
    ! do j=1, Nr
    !    do i=1, Nz
    !       prs(i,j) = prscoef * prs(i,j)
    !    enddo
    ! enddo
    Ttot    = ( 1.d0 + TiTe ) * vthcv**2
    Ttotinv =   1.d0 / Ttot
    rhof    = maxval( prs ) * Ttotinv * rhofloor
    do j=1, Nr
       do i=1, Nz
          rho(i,j) = prs(i,j) * Ttotinv + rhof
       enddo
    enddo

    !  --- [5] Adjusting to Bmax  ---  !
    Jr(:,:) = 0.d0
    Jt(:,:) = 0.d0
    Jz(:,:) = 0.d0
    valfe2  = valfe**2    
    do j=1, Nr
       do i=1, Nz
          Jr(i,j) = - valfe2 *   ( Bt(i,j) - Bt(i-1,j) ) * dzinv
          Jt(i,j) =   valfe2 * ( ( Br(i,j) - Br(i-1,j) ) * dzinv &
               &               - ( Bz(i,j) - Bz(i,j-1) ) * drinv )
          Jz(i,j) =   valfe2 * ( rh(j)*Bt(i,j) - rh(j-1)*Bt(i,j-1) ) * drinv * rfinv(j)
       enddo
    enddo

    return
  end subroutine calc_Bfield_Egrid

  
  subroutine calc_Bfield_Bgrid

    use constants, only : Bmax, valfe, TiTe, vthcv, rhofloor, rmin, normSW
    use variables, only : dr, dz, psi, gfunc, prs
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2, valfe2, drinv, dzinv
    double precision   :: absBmax, prscoef, Ttot, Ttotinv, rhof
    double precision   :: psiext(LIg+1,LJg+1), gfcext(LIg,LJg)
    double precision   :: rf(LJg), rfinv(LJg), rh(LJg), rhinv(LJg)

    !  --- [1] Set Grid  ---   !
    drinv       = 1.d0 / dr
    dzinv       = 1.d0 / dz
    do j=1, LJg
       rf(j)    = dr * dble(j-2) + rmin
       rfinv(j) = 1.d0 / rf(j)
    enddo
    if ( rmin .eq. 0.d0 ) then
       rf(2)    = 0.d0
       rfinv(2) = 0.d0
    endif
    do j=1, LJg-1
       rh(j)    = 0.5d0 * ( rf(j) + rf(j+1) )
       rhinv(j) = 1.d0 / rh(j)
    enddo
    rh(LJg)     = rf(LJg) + 0.5d0*dr
    rhinv(LJg)  = 1.d0 / rh(LJg)
        
    !  --- [2]  Copy and reflect Psi  ---  !
    !  --- [2-1] Main Region  ---  !
    do j=1, Nr          ! j=2, LJg   for j+1
       do i=1, Nz       ! i=2, LIg   for i+1
          psiext(i+1,j+1) =   psi(i,j)
          gfcext(i+1,j+1) = gfunc(i,j)          
       enddo
    enddo
    
    !  --- [2-2] j=1, LJ+1 boundary copy ---  !
    factor1     =   rh(  1) * rhinv(    2)
    factor2     =   rh(LJg) * rhinv(LJg-1)
    do i=2, LIg
       psiext(i,    1) = ( 1.d0+factor1 )*psiext(i,  2) - factor1*psiext(i,    3)
       psiext(i,LJg+1) = ( 1.d0+factor2 )*psiext(i,LJg) - factor2*psiext(i,LJg-1)
    enddo
    !  --- [2-3] i=1, LI+1 boundary copy ---  !
    do j=1, LJg+1
       psiext(    1,j) = 2.d0 * psiext(  2,j) - psiext(    3,j)
       psiext(LIg+1,j) = 2.d0 * psiext(LIg,j) - psiext(LIg-1,j)
    enddo
    
    !  --- [3]   Interpolation of psi    ---  !
    do j=1, LJg
       do i=1, LIg
          fgB(i,j,1) = 0.25d0 * ( + psiext(i,j  ) + psiext(i+1,j  ) &
               &                  + psiext(i,j+1) + psiext(i+1,j+1) )
       enddo
    enddo
    
    !  --- [4]    Calculate Br & Bz      ---  !
    BgB(:,:,:)  = 0.d0
    do j=1, LJg
       do i=2, LIg
          BgB(1,i,j)   = + ( fgB(i,j,1) - fgB(i-1,j,1) ) * rhinv(j) * dzinv
       enddo
    enddo ! ~~ [4-1] Br ~~ !
    do j=2, LJg
       do i=2, LIg
          BgB(2,i,j)   = + rfinv(j) * gfcext(i,j)
       enddo
    enddo ! ~~ [4-2] Bt ~~ !
    do j=2, LJg
       do i=1, LIg
          BgB(3,i,j)   = - ( fgB(i,j,1) - fgB(i,j-1,1) ) * rfinv(j) * drinv
       enddo
    enddo ! ~~ [4-3] Bz ~~ !
    !   -- [4-4] Boundary Condition B -- !
    BgB(1,1,:)  = BgB(1,2,:)
    BgB(2,1,:)  = BgB(2,2,:)
    if ( rmin.eq.0.d0 ) then
       BgB(1,:,1)      = - BgB(1,:,2)
       BgB(2,:,1)      = - BgB(2,:,3)
       BgB(2,:,2)      = 0.d0
       BgB(3,:,1)      = + BgB(3,:,3)
       BgB(3,:,2)      = + BgB(3,:,3)
    else
       BgB(2,:,1)      = rf(2) / rf(1) * BgB(2,:,2)
       BgB(3,:,1)      = BgB(3,:,2)
    endif

    !  --- [5] Adjust to Bmax  ---  !
    if ( normSW ) then
       absBmax = 0.d0
       do j=1, LJg
          do i=1, LIg
             absBmax   = max( absBmax, BgB(1,i,j)**2 + BgB(2,i,j)**2 + BgB(3,i,j)**2 )
          enddo
       enddo
       absBmax = Bmax / sqrt( absBmax )
       prscoef = absBmax**2
       do j=1, LJg
          do i=1, LIg
             BgB(1,i,j)   = absBmax *   BgB(1,i,j)
             BgB(2,i,j)   = absBmax *   BgB(2,i,j)
             BgB(3,i,j)   = absBmax *   BgB(3,i,j)
             fgB(i,j,3) = prscoef * fgB(i,j,3)
          enddo
       enddo
    endif
    
    !  --- [6] pressure & rho  ---  !
    !   -- [6-1] Copy FRP -- !
    do j=1, Nr
       do i=1, Nz
          fgB(i+1,j+1,3) = prs(i,j)
       enddo
    enddo
    fgB(1,:,3) = fgB(2,:,3)
    fgB(:,1,3) = fgB(:,2,3)
    !   -- [6-2] Limit rho to rhoFloor -- !
    Ttot    = ( 1.d0 + TiTe ) * vthcv**2
    Ttotinv =   1.d0 / Ttot
    rhof    = maxval( prs ) * Ttotinv * rhofloor
    do j=1, LJg
       do i=1, LIg
          fgB(i,j,2) = fgB(i,j,3) * Ttotinv + rhof
          fgB(i,j,3) = fgB(i,j,2) * Ttot
       enddo
    enddo
    
    !  --- [7]   Calculate Current Density ---  !
    !   -- [7-1] Current = rot B -- !
    valfe2  = valfe**2
    do j=1, LJg-1
       do i=1, LIg-1
          JgB(1,i,j)   = - valfe2 *   ( BgB(2,i+1,j) - BgB(2,i,j) ) * dzinv
          JgB(2,i,j)   =   valfe2 * ( ( BgB(1,i+1,j) - BgB(1,i,j) ) * dzinv &
               &         - ( BgB(3,i,j+1) - BgB(3,i,j) ) * drinv )
          JgB(3,i,j)   =   valfe2 *   ( rf(j+1)*BgB(2,i,j+1) - rf(j)*BgB(2,i,j) ) * drinv * rhinv(j)
       enddo
    enddo
    !   -- [7-2] Boundary Condition -- !
    JgB(1,:,LJg)       = 0.d0
    JgB(2,:,LJg)       = 0.d0
    JgB(3,:,LJg)       = 0.d0
    JgB(1,LIg,:)       = 0.d0
    JgB(2,LIg,:)       = 0.d0
    JgB(3,LIg,:)       = 0.d0
    JgB(3,:,  1)       = 0.d0
    
    !  --- [8] Get Electron Flow Velocity ---  !
    do j=1, LJg
       do i=1, LIg
          ugB(i,j,1)   = - JgB(1,i,j) / fgB(i,j,2)
          ugB(i,j,2)   = - JgB(2,i,j) / fgB(i,j,2)
          ugB(i,j,3)   = - JgB(3,i,j) / fgB(i,j,2)
       enddo
    enddo

    return
  end subroutine calc_Bfield_Bgrid




end module Field
