module BgridXZMod
contains

  subroutine CalculateField( coordinate, onGrid )
    implicit none
    character(1), intent(in) :: onGrid
    character(2), intent(in) :: coordinate
    
    if   ( coordinate .eq. 'XZ' )                           call calcBJ_XZ
    if ( ( coordinate .eq. 'RZ' ).and.( onGrid .eq. 'E' ) ) call calcBJ_RZ_EGrid
    if ( ( coordinate .eq. 'RZ' ).and.( onGrid .eq. 'B' ) ) call calcBJ_RZ_BGrid
    call setPrsrRho
    
    return
  end subroutine CalculateField

  
  subroutine calcBJ_XZ
    
    use constants, only : N1 , N2, Bmax, valfe
    use variables, only : Avp, x1, x2
    use variables, only : Bfd, Jcr
    implicit none
    integer            :: i, j
    double precision   :: dx1inv, dx2inv, absB2Max, valfe2

    !  --- [1] Magnetic Field  :: B --- !
    do j=1, N2-1
       dx2inv = 1.0d0 / ( x2(j+1) - x2(j) )
       do i=1, N1-1
          dx1inv     = 1.0d0 / ( x1(i+1) - x1(i) )
          Bfd(1,i,j) = + ( Avp(i,j+1) - Avp(i,j) ) * dx2inv
          Bfd(2,i,j) = - ( Avp(i+1,j) - Avp(i,j) ) * dx1inv
       enddo
    enddo
    call PCW_bc

    !  --- [2]  Normalize B-field  ---- !
    absB2Max = 0.d0
    do j=1, N2
       do i=1, N1
          absB2Max = max( absB2Max, Bfd(1,i,j)**2+Bfd(2,i,j)**2+Bfd(3,i,j)**2 )
       enddo
    enddo
    absB2Max = Bmax / sqrt( absB2Max )
    do j=1, N2
       do i=1, N1
          Bfd(1,i,j) = absB2Max * Bfd(1,i,j)
          Bfd(2,i,j) = absB2Max * Bfd(2,i,j)
          Bfd(3,i,j) = absB2Max * Bfd(3,i,j)
       enddo
    enddo
    
    !  --- [3] Current Density :: J --- !
    valfe2 = valfe**2
    do j=1, N2-1
       dx2inv = 1.0d0 / ( x2(j+1) - x2(j) )
       do i=1, N1-1
          dx1inv     = 1.0d0 / ( x1(i+1) - x1(i) )

          Jcr(1,i,j) = valfe2 * ( + ( Bfd(3,i,j+1) - Bfd(3,i,j-1) ) * dx2inv )
          Jcr(2,i,j) = valfe2 * ( - ( Bfd(3,i+1,j) - Bfd(3,i-1,j) ) * dx1inv )
          Jcr(3,i,j) = valfe2 * ( + ( Bfd(2,i+1,j) - Bfd(2,i-1,j) ) * dx1inv &
               &       - ( Bfd(1,i,j+1) - Bfd(1,i,j-1) ) * dx2inv )
       enddo
    enddo
    ! -- B.C. -- !
    do j=1, N2
       Jcr(1, 1,j) = 0.d0
       Jcr(1,N1,j) = 0.d0
       Jcr(2, 1,j) = 0.d0
       Jcr(2, 2,j) = 0.d0
       Jcr(2,N1,j) = 0.d0
       Jcr(3, 1,j) = 0.d0
       Jcr(3, 2,j) = 0.d0
       Jcr(3,N1,j) = 0.d0
    enddo
    do i=1, N1
       Jcr(1,i, 1) = 0.d0
       Jcr(1,i, 2) = 0.d0
       Jcr(1,i,N2) = 0.d0
       Jcr(2,i, 1) = 0.d0
       Jcr(2,i,N2) = 0.d0
       Jcr(3,i, 1) = 0.d0
       Jcr(3,i, 2) = 0.d0
       Jcr(3,i,N2) = 0.d0
    enddo

    return 
  end subroutine calcBJ_XZ

  
  subroutine calcBJ_RZ_EGrid
    
    use constants, only : N1 , N2, Bmax, valfe, x2min
    use variables, only : Avp, psi, pcf, x1, x2
    use variables, only : Bfd, Jcr
    implicit none
    integer            :: i, j
    double precision   :: dx1inv, x2inv, dx2inv, absB2Max, valfe2

    !  --- [1] Magnetic Fluc Function :: psi --- !
    do j=1, N2
       do i=1, N1
          psi(i,j) = x2(j) * Avp(i,j)
       enddo
    enddo
    
    !  --- [2] Magnetic Field  :: B --- !
    !     --- ( Bt is initially stored in Bfd(3,:,:), thus relocate to the center of a grid. ) ---
    do j=1, N2-1
       do i=1, N1-1
          Bfd(2,i,j) = + 0.25d0 * ( Bfd(3,i,j) + Bfd(3,i+1,j) + Bfd(3,i,j+1) + Bfd(3,i+1,j+1) )
       enddo
    enddo !   --  Bt  --  !
    
    do j=1, N2-1
       x2inv  =  2.d0 / ( x2(j+1) + x2(j) )
       dx2inv = x2inv / ( x2(j+1) - x2(j) )
       do i=1, N1-1
          dx1inv     = 1.d0 / ( ( x1(i+1) - x1(i) )*x2(j) )
          
          Bfd(1,i,j) = - ( psi(i+1,j) - psi(i,j) ) * dx1inv
          Bfd(3,i,j) = + ( psi(i,j+1) - psi(i,j) ) * dx2inv
       enddo
    enddo !   -- Br & Bz  -- !
    
    if ( x2min .eq. 0.d0 ) then
       do i=1, N1-1
          Bfd(1,i,1) = 0.d0
       enddo
    endif !   -- r=0 ver. --   !
    
    do j=1, N2-1
       Bfd(1,N1,j) = Bfd(1,N1-1,j)
       Bfd(2,N1,j) = Bfd(2,N1-1,j)
       Bfd(3,N1,j) = + ( psi(i,j+1) - psi(i,j) ) / ( 0.5d0*( x2(j+1)+x2(j) ) * ( x2(j+1)-x2(j) ) )
    enddo
    do i=1, N1
       Bfd(1,i,N2) = 0.d0
       Bfd(2,i,N2) = x2(N2-1) / x2(N2) * Bfd(2,i,N2-1)
       Bfd(3,i,N2) = Bfd(3,i,N2-1)
    enddo !   --    b.c.  --   !
    
    !  --- [3]  Normalize B-field  ---- !
    absB2Max = 0.d0
    do j=1, N2
       do i=1, N1
          absB2Max = max( absB2Max, Bfd(1,i,j)**2+Bfd(2,i,j)**2+Bfd(3,i,j)**2 )
       enddo
    enddo
    absB2Max = Bmax / sqrt( absB2Max )
    do j=1, N2
       do i=1, N1
          Bfd(1,i,j) = absB2Max * Bfd(1,i,j)
          Bfd(2,i,j) = absB2Max * Bfd(2,i,j)
          Bfd(3,i,j) = absB2Max * Bfd(3,i,j)
       enddo
    enddo
    
    !  --- [4] Poloidal Current Function :: F --- !
    do j=1, N2
       do i=1, N1
          pcf(i,j) = 0.5d0 * ( x2(j+1) + x2(j) ) * Bfd(2,i,j)
       enddo
    enddo
    
    !  --- [5] Current Density :: J --- !
    valfe2 = valfe**2
    do j=1, N2-1
       x2inv  = 1.d0 / x2(j)
       dx2inv = 1.d0 / ( x2(j+1) - x2(j) )
       do i=1, N1-1
          dx1inv     = 1.0d0 / ( x1(i+1) - x1(i) )
          
          Jcr(1,i,j) = valfe2 * ( - ( Bfd(2,i+1,j) - Bfd(2,i,j) ) * dx1inv )
          Jcr(2,i,j) = valfe2 * ( + ( Bfd(1,i+1,j) - Bfd(1,i,j) ) * dx1inv &
               &                  - ( Bfd(3,i,j+1) - Bfd(3,i,j) ) * dx2inv )
          Jcr(3,i,j) = valfe2 * ( + ( pcf(i,j+1) - pcf(i,j) ) * dx2inv * x2inv )
       enddo
    enddo
    ! -- B.C. -- !
    do j=1, N2
       Jcr(1,N1,j) = 0.d0
       Jcr(2,N1,j) = 0.d0
       Jcr(3,N1,j) = 0.d0
    enddo
    do i=1, N1
       Jcr(1,i,N2) = 0.d0
       Jcr(2,i,N2) = 0.d0
       Jcr(3,i,N2) = 0.d0
       Jcr(3,i, 1) = 0.d0
    enddo
    
    return 
  end subroutine calcBJ_RZ_EGrid


  subroutine calcBJ_RZ_BGrid
    
    use constants, only : N1 , N2, Bmax, valfe, x2min, BtPattern, polarity
    use variables, only : Avp, psi, pcf, x2, dx1, dx2
    use variables, only : Bfd, Jcr
    implicit none
    integer            :: i, j
    double precision   :: rh(0:N2), rf(0:N2+1), rhinv(0:N2), rfinv(0:N2+1), drinv, dzinv
    double precision   :: absB2Max, valfe2, BtSign1, BtSign2
    double precision   :: Aext(0:N1+1,0:N2+1), Aintp(0:N1,0:N2), psih(0:N1,0:N2)
    double precision   :: factor1, factor2, factor3, factor4
    double precision   :: Reverse = -1.d0

    !  --- [1] Grid Making  ---   !
    dzinv = 1.d0 / dx1
    drinv = 1.d0 / dx2
    
    do j=1, N2
       rf(j) =  x2(j)
    enddo
    rf(0)    = x2(1)  - dx2
    rf(N2+1) = x2(N2) + dx2
    do j=0, N2+1
       rfinv(j) = 1.d0 / rf(j)
    enddo
    if ( x2min.eq.0.d0 ) rfinv(1) = 0.d0
    
    do j=0, N2
       rh(j)    = ( rf(j) + rf(j+1) ) * 0.5d0
       rhinv(j) = 1.d0 / rh(j)
    enddo
    
    
    !  --- [2] extend Avp  ---   !
    do j=1, N2
       do i=1, N1
          Aext(i,j) = Avp(i,j)
       enddo
    enddo
    ! factor1 = + rf(1   ) * rfinv(   0) * ( 1.d0 + rh( 0)*rhinv(   1) )
    ! factor2 = - rf(2   ) * rfinv(   0) * rh( 0) * rhinv(1)
    factor1 =   0.d0
    factor2 = - 1.d0
    factor3 = + rf(N2  ) * rfinv(N2+1) * ( 1.d0 + rh(N2)*rhinv(N2-1) )
    factor4 = - rf(N2-1) * rfinv(N2+1) * rh(N2) * rhinv(N2-1)
    do i=1, N1
       Aext(i,   0) = factor1 * Aext(i, 1) + factor2 * Aext(i,   2)
       Aext(i,N2+1) = factor3 * Aext(i,N2) + factor4 * Aext(i,N2-1)
    enddo
    do j=0, N2+1
       Aext(   0,j) = 2.d0 * Aext( 1,j) - Aext(   2,j)
       Aext(N1+1,j) = 2.d0 * Aext(N1,j) - Aext(N1-1,j)
    enddo
          
    !  --- [3] Interpolate At  --- !
    do j=0, N2
       do i=0, N1
          Aintp(i,j) = 0.25d0 * rhinv(j) * ( rf(j)   * ( Aext(i,j  ) + Aext(i+1,j  ) ) &
               &                           + rf(j+1) * ( Aext(i,j+1) + Aext(i+1,j+1) ) )
          psih(i,j)  =  rh(j) * Aintp(i,j)
       enddo
    enddo
    
    !  --- [4] Magnetic Flux Function :: psi --- !
    do j=0, N2
       do i=0, N1
          psi(i,j) = psih(i,j)
       enddo
    enddo
    
    !  --- [5] Magnetic Field  :: B --- !
    !     --- ( Bt is initially stored in Bfd(3,:,:), thus relocate to the center of a grid. ) ---
    Bfd(2,:,:) = 0.d0
    if ( BtPattern.eq.'Co-' ) then       
       BtSign1 = + 1.d0
       BtSign2 = + 1.d0
    endif
    if ( BtPattern.eq.'Ctr' ) then
       if ( polarity.eq.'CaseI' ) then
          BtSign1 = + 1.d0
          Btsign2 = - 1.d0
       endif
       if ( polarity.eq.'CaseO' ) then
          BtSign1 = - 1.d0
          Btsign2 = + 1.d0
       endif
    endif
    BtSign1 = BtSign1*Reverse
    BtSign2 = BtSign2*Reverse
    do j=1, N2
       do i=1, N1/2
          Bfd(2,i,j) = BtSign1 * Bfd(3,i,j)
       enddo
    enddo
    do j=1, N2
       do i=N1/2+1, N1
          Bfd(2,i,j) = BtSign2 * Bfd(3,i,j)
       enddo       
    enddo !   --  Bt  --  !

    Bfd(1,:,:) = 0.d0
    do j=0, N2
       do i=1, N1
          Bfd(1,i,j) = - rhinv(j) * dzinv * ( psih(i,j) - psih(i-1,j) )
       enddo
    enddo 
    Bfd(3,:,:) = 0.d0
    do j=1, N2
       do i=0, N1
          Bfd(3,i,j) = + rfinv(j) * drinv * ( psih(i,j) - psih(i,j-1) )
       enddo
    enddo!   --  Br, Bz  --  !

    
    !  --- [6]  Normalize B-field  ---- !
    absB2Max = 0.d0
    do j=0, N2
       do i=0, N1
          absB2Max = max( absB2Max, Bfd(1,i,j)**2+Bfd(2,i,j)**2+Bfd(3,i,j)**2 )
       enddo
    enddo
    absB2Max = Bmax / sqrt( absB2Max )
    do j=0, N2
       do i=0, N1
          Bfd(1,i,j) = absB2Max * Bfd(1,i,j)
          Bfd(2,i,j) = absB2Max * Bfd(2,i,j)
          Bfd(3,i,j) = absB2Max * Bfd(3,i,j)
       enddo
    enddo

    !  --- [7] Boundary Condition --- !
    do j=0, N2
       Bfd(1,0,j) = Bfd(1,1,j)
       Bfd(2,0,j) = Bfd(2,1,j)
    enddo
    ! factor1 = rf(1) * rfinv(0)
    factor1 = - 1.d0
    do i=0, N1
       Bfd(2,i,0) = - Bfd(2,i,2)
       Bfd(3,i,1) = + Bfd(3,i,2)
       Bfd(3,i,0) = + Bfd(3,i,2)
    enddo
    
    !  --- [8] Poloidal Current Function :: F --- !
    do j=1, N2
       do i=1, N1
          pcf(i,j) = rf(j) * Bfd(2,i,j)
       enddo
    enddo
    
    !  --- [9] Current Density :: J --- !
    valfe2 = valfe**2
    do j=1, N2-1
       do i=1, N1-1
          Jcr(1,i,j) = valfe2 * ( - dzinv * ( Bfd(2,i+1,j) - Bfd(2,i,j) ) )
          Jcr(2,i,j) = valfe2 * ( + dzinv * ( Bfd(1,i+1,j) - Bfd(1,i,j) ) &
               &                  - drinv * ( Bfd(3,i,j+1) - Bfd(3,i,j) ) )
          Jcr(3,i,j) = valfe2 * ( + drinv * rhinv(j) * ( pcf(i,j+1) - pcf(i,j) ) )
       enddo
    enddo
    ! -- B.C. -- !
    factor1 = rh(   2) / rh( 1)
    factor2 = rh(N2-1) / rh(N2)
    do i=1, N1-1
       Jcr(1,i, 0) = 0.d0
       Jcr(1,i, 1) = 0.d0
       Jcr(1,i,N2) = 0.d0
       Jcr(2,i, 0) = - Jcr(2,i,   1)
       Jcr(2,i,N2) = - Jcr(2,i,N2-1) * factor2
       Jcr(3,i, 0) = + Jcr(3,i,   1)
       Jcr(3,i,N2) = + Jcr(3,i,N2-1)
    enddo
    do j=0, N2
       Jcr(1, 0,j) = - Jcr(1,   1,j)
       Jcr(1,N1,j) = - Jcr(1,N1-1,j)
       Jcr(2, 0,j) = - Jcr(2,   1,j)
       Jcr(2,N1,j) = - Jcr(2,N1-1,j)
       Jcr(3, 0,j) = 0.d0
       Jcr(3, 1,j) = 0.d0
       Jcr(3,N1,j) = 0.d0
    enddo
    
    return 
  end subroutine calcBJ_RZ_BGrid

  
  subroutine PCW_bc
    use constants, only : N1, N2
    use variables, only : Bfd
    implicit none
    integer            :: i, j

    !  -- Left & Right -- !
    do i=1, N1
       Bfd(1,i, 1) = Bfd(1,i,   2)
       Bfd(1,i,N2) = Bfd(1,i,N2-1)
       Bfd(2,i, 1) = 0.d0
       Bfd(2,i,N2) = 0.d0
       Bfd(3,i, 1) = Bfd(3,i,   2)
       Bfd(3,i,N2) = Bfd(3,i,N2-1)
    enddo
    !  -- Top & Bottom -- !
    do j=1, N2
       Bfd(1, 1,j) = 0.d0
       Bfd(1,N1,j) = 0.d0
       Bfd(2, 1,j) = Bfd(2,   2,j)
       Bfd(2,N1,j) = Bfd(2,N1-1,j)
       Bfd(3, 1,j) = Bfd(3,   2,j)
       Bfd(3,N1,j) = Bfd(3,N1-1,j)
    enddo

    return
  end subroutine PCW_bc

  
  subroutine setPrsrRho
    use constants, only : N1, N2
    use constants, only : vthcv, TiTe, mr, desiredBetapol, valfe
    use variables, only : Bfd, Jcr, prs, rho, uvc, frp, psi
    implicit none
    integer            :: i, j, k
    double precision   :: maxBpol2, coef, uvcMax, uvcMin
    double precision   :: q(2)
    
    ! --- [1] Find Max Bpol --- !
    maxBpol2 = 0.d0
    do j=1, N2
       do i=1, N1
          maxBpol2 = max( maxBpol2, Bfd(1,i,j)**2 + Bfd(2,i,j)**2 )
       enddo
    enddo
    
    ! --- [2] set Pressure  --- !
    do j=0, N2
       do i=0, N1
          prs(i,j) = desiredBetapol * 0.5d0 * valfe**2 * maxBpol2
       enddo
    enddo
    write(6,'(a,f12.6)') ' Set Pressure as constant :: P0 ===  ', prs(1,1)

    ! --- [3] set rho --- !
    coef = 1.d0 / ( vthcv**2 * ( 1.d0 + TiTe ) )
    do j=0, N2
       do i=0, N1
          rho(i,j) = coef * prs(i,j)
       enddo
    enddo
    write(6,'(a,f12.6)') ' Set density  as constant :: n0 ===  ', rho(1,1)

    ! --- [4] set velocity --- !
    q(1) = - 1.d0
    q(2) = + 1.d0
    uvcMax = -1.d0
    uvcMin = +1.d8
    do k=1, 2
       do j=0, N2
          do i=0, N1
             uvc(1,i,j,k) = Jcr(1,i,j) / ( q(k)*rho(i,j) )
             uvc(2,i,j,k) = Jcr(2,i,j) / ( q(k)*rho(i,j) )
             uvc(3,i,j,k) = Jcr(3,i,j) / ( q(k)*rho(i,j) )
             uvcMax       = max( uvcMax, uvc(1,i,j,k)**2 + uvc(2,i,j,k)**2 + uvc(3,i,j,k)**2 )
             uvcMin       = min( uvcMin, uvc(1,i,j,k)**2 + uvc(2,i,j,k)**2 + uvc(3,i,j,k)**2 )
          enddo
       enddo
    enddo       
    write(6,'(a,f12.6)') ' Minimum velocity  :: vMin ===  ', sqrt( uvcMin )
    write(6,'(a,f12.6)') ' Maximum velocity  :: vMax ===  ', sqrt( uvcMax )

    ! --- [5] set frp --- !
    do j=0, N2
       do i=0, N1
          frp(i,j,1) = psi(i,j)
          frp(i,j,2) = rho(i,j)
          frp(i,j,3) = prs(i,j)
       enddo
    enddo
    
    return
  end subroutine setPrsrRho

end module BgridXZMod






  ! subroutine calcBJ_RZ_BGrid_org
    
  !   use constants, only : N1 , N2, Bmax, valfe, x2min
  !   use variables, only : Avp, psi, pcf, x1, x2, dx1, dx2
  !   use variables, only : Bfd, Jcr, frp
  !   implicit none
  !   integer            :: i, j
  !   double precision   :: rh(0:N2), rf(0:N2+1), rhinv(0:N2), rfinv(0:N2+1), drinv, dzinv
  !   double precision   :: dx1inv, x2inv, dx2inv, absB2Max, valfe2
  !   double precision   :: Aext(0:N1+1,0:N2+1), Aintp(0:N1,0:N2), psih(0:N1,0:N2)
  !   double precision   :: factor1, factor2, factor3, factor4

  !   !  --- [1] Grid Making  ---   !
  !   if ( x2min.ne.0.d0 ) write(6,*) ' [ERROR] Not Spheromak [ERROR] '
  !   dzinv = 1.d0 / dx1
  !   drinv = 1.d0 / dx2
  !   !  -- rf / rfinv -- !
  !   do j=1, N2
  !      rf(j) =  x2(j)
  !   enddo
  !   rf(0)    = x2(1)  - dx2
  !   rf(N2+1) = x2(N2) + dx2
  !   do j=0, N2+1
  !      rfinv(j) = 1.d0 / rf(j)
  !   enddo
  !   !  -- rh / rhinv -- !
  !   do j=0, N2
  !      rh(j)    = ( rf(j) + rf(j+1) ) * 0.5d0
  !      rhinv(j) = 1.d0 / rh(j)
  !   enddo
    
  !   !  --- [2] extend Avp  ---   !
  !   !   --    Bulk    --  !
  !   do j=1, N2
  !      do i=1, N1
  !         Aext(i,j) = Avp(i,j)
  !      enddo
  !   enddo
  !   !   -- Peripheral --  !
  !   factor1 = + rf(1   ) * rfinv(   0) * ( 1.d0 + rh( 0)*rhinv(   1) )
  !   factor2 = - rf(2   ) * rfinv(   0) * rh( 0) * rhinv(1)
  !   factor3 = + rf(N2  ) * rfinv(N2+1) * ( 1.d0 + rh(N2)*rhinv(N2-1) )
  !   factor4 = - rf(N2-1) * rfinv(N2+1) * rh(N2) * rhinv(N2-1)
  !   do i=1, N1
  !      Aext(i,   0) = factor1 * Aext(i, 1) + factor2 * Aext(i,   2)
  !      Aext(i,N2+1) = factor3 * Aext(i,N2) + factor4 * Aext(i,N2-1)
  !   enddo
  !   do j=0, N2+1
  !      Aext(   0,j) = 2.d0 * Aext( 1,j) - Aext(   2,j)
  !      Aext(N1+1,j) = 2.d0 * Aext(N1,j) - Aext(N1-1,j)
  !   enddo

  !   !  --- [3] Interpolation At  --- !
  !   do j=0, N2
  !      do i=0, N1
  !         Aintp(i,j) = 0.25d0 * rhinv(j) * ( rf(j)   * ( Aext(i,j  ) + Aext(i+1,j  ) ) &
  !              &                           + rf(j+1) * ( Aext(i,j+1) + Aext(i+1,j+1) ) )
  !         psih(i,j)  =  rh(j) * Aintp(i,j)
  !      enddo
  !   enddo
    
  !   !  --- [4] Magnetic Flux Function :: psi --- !
  !   do j=0, N2
  !      do i=0, N1
  !         psi(i,j) = psih(i,j)
  !      enddo
  !   enddo
    
  !   !  --- [5] Magnetic Field  :: B --- !
  !   !     --- ( Bt is initially stored in Bfd(3,:,:), thus relocate to the center of a grid. ) ---
  !   Bfd(2,:,:) = 0.d0
  !   do j=1, N2
  !      do i=1, N1/2
  !         Bfd(2,i,j) = Bfd(3,i,j)
  !      enddo
  !   enddo
  !   do j=1, N2
  !      do i=N1/2+1, N1
  !         Bfd(2,i,j) = - Bfd(3,i,j)
  !      enddo       
  !   enddo !   --  Bt  --  !

  !   Bfd(1,:,:) = 0.d0
  !   do j=0, N2
  !      do i=1, N1
  !         Bfd(1,i,j) = - rhinv(j) * dzinv * ( psih(i,j) - psih(i-1,j) )
  !      enddo
  !   enddo !   --  Br, Bz  --  !

  !   Bfd(3,:,:) = 0.d0
  !   do j=1, N2
  !      do i=0, N1
  !         Bfd(3,i,j) = + rfinv(j) * drinv * ( psih(i,j) - psih(i,j-1) )
  !      enddo
  !   enddo
    
  !   !  --- [6]  Normalize B-field  ---- !
  !   absB2Max = 0.d0
  !   do j=0, N2
  !      do i=0, N1
  !         absB2Max = max( absB2Max, Bfd(1,i,j)**2+Bfd(2,i,j)**2+Bfd(3,i,j)**2 )
  !      enddo
  !   enddo
  !   absB2Max = Bmax / sqrt( absB2Max )
  !   do j=0, N2
  !      do i=0, N1
  !         Bfd(1,i,j) = absB2Max * Bfd(1,i,j)
  !         Bfd(2,i,j) = absB2Max * Bfd(2,i,j)
  !         Bfd(3,i,j) = absB2Max * Bfd(3,i,j)
  !      enddo
  !   enddo

  !   !  --- [7] Boundary Condition --- !
  !   do j=0, N2
  !      Bfd(1,0,j) = Bfd(1,1,j)
  !      Bfd(2,0,j) = Bfd(2,1,j)
  !   enddo
  !   factor1 = rf(1) * rfinv(0)
  !   do i=0, N1
  !      Bfd(2,i,0) = factor1 * Bfd(2,i,2)
  !      Bfd(3,i,0) = factor1 * Bfd(3,i,2)
  !   enddo
  !   !  --- [8] Poloidal Current Function :: F --- !
  !   do j=1, N2
  !      do i=1, N1
  !         pcf(i,j) = rf(j) * Bfd(2,i,j)
  !      enddo
  !   enddo
    
  !   !  --- [9] Current Density :: J --- !
  !   valfe2 = valfe**2
  !   do j=1, N2-1
  !      do i=1, N1-1
  !         Jcr(1,i,j) = valfe2 * ( - dzinv * ( Bfd(2,i+1,j) - Bfd(2,i,j) ) )
  !         Jcr(2,i,j) = valfe2 * ( + dzinv * ( Bfd(1,i+1,j) - Bfd(1,i,j) ) &
  !              &                  - drinv * ( Bfd(3,i,j+1) - Bfd(3,i,j) ) )
  !         Jcr(3,i,j) = valfe2 * ( + drinv * rhinv(j) * ( pcf(i,j+1) - pcf(i,j) ) )
  !      enddo
  !   enddo
  !   !   -- B.C. -- !
  !   factor1 = rh(   2) / rh( 1)
  !   factor2 = rh(N2-1) / rh(N2)
  !   do i=1, N1-1
  !      Jcr(1,i, 0) = 0.d0
  !      Jcr(1,i, 1) = 0.d0
  !      Jcr(1,i,N2) = 0.d0
  !      Jcr(2,i, 0) = - Jcr(2,i,   1) * factor1
  !      Jcr(2,i,N2) = - Jcr(2,i,N2-1) * factor2
  !      Jcr(3,i, 0) = - Jcr(3,i,   1)
  !      Jcr(3,i,N2) = - Jcr(3,i,N2-1)
  !   enddo
  !   do j=0, N2
  !      Jcr(1, 0,j) = - Jcr(1,   1,j)
  !      Jcr(1,N1,j) = - Jcr(1,N1-1,j)
  !      Jcr(2, 0,j) = - Jcr(2,   1,j)
  !      Jcr(2,N1,j) = - Jcr(2,N1-1,j)
  !      Jcr(3, 0,j) = 0.d0
  !      Jcr(3, 1,j) = 0.d0
  !      Jcr(3,N1,j) = 0.d0
  !   enddo
    
  !   return 
  ! end subroutine calcBJ_RZ_BGrid_org

