module EgridEBMod
  use constants, only : N1, N2
  implicit none
  integer, parameter :: LIg = N1+1
  integer, parameter :: LJg = N2+1
  double precision   :: rf(LJg+1), rh(LJg+1), rfinv(LJg+1), rhinv(LJg+1)
  double precision   :: BgB(3,LIg,LJg), JgB(3,LIg,LJg), fgB(LIg,LJg,3), ugB(LIg,LJg,3)
contains

  subroutine EgridField( coordinate, onGrid )
    implicit none
    character(1), intent(in) :: onGrid
    character(2), intent(in) :: coordinate
    
    if ( ( coordinate.eq.'XZ' ).and.( onGrid.eq.'E' ) ) call Cartesian__EgridField
    if ( ( coordinate.eq.'RZ' ).and.( onGrid.eq.'E' ) ) call calcBJ_RZ_EGrid
    call setPrsrRho
    
    return
  end subroutine EgridField

  
  subroutine Cartesian__EgridField
    use constants, only : Bmax, valfe, normType, MergingType, normSW, symmBreak, Flag__ReverseBt
    use variables, only : dx1, dx2, Avp, Bfd
    implicit none
    integer            :: i, j
    double precision   :: dx1inv, dx2inv, coef, absB2Max, BtSign1, BtSign2
    double precision   :: AvpExt(LIg+1,LJg+1)

    !  --- [1] Extend Avp Field     --- !
    !   -- [1-1] Change Polarity of Bt -- !
    select case( MergingType )
    case( 'CoHx' )
       BtSign1 = + ( 1.d0 + symmBreak )
       BtSign2 = + ( 1.d0 - symmBreak )
    case( 'CtrI' )
       BtSign1 = + ( 1.d0 + symmBreak )
       Btsign2 = - ( 1.d0 - symmBreak )
    case( 'CtrO' )
       BtSign1 = - ( 1.d0 + symmBreak )
       Btsign2 = + ( 1.d0 - symmBreak )
    case default
       write(6,*) '[Cartesian_EgridField-@EgridRZMod.f90] MergingType == ??', MergingType
       stop
    end select
    if ( Flag__ReverseBt ) then
       BtSign1 = -1.d0 * BtSign1
       BtSign2 = -1.d0 * BtSign2
    endif
    do j=1, N2
       do i=1, N1/2
          Bfd(3,i,j) = BtSign1 * Bfd(3,i,j)
       enddo
    enddo
    do j=1, N2
       do i=N1/2+1, N1
          Bfd(3,i,j) = BtSign2 * Bfd(3,i,j)
       enddo
    enddo !   --  Bt  --  !
    !   -- [1-2] Get AvpExt, BgB       -- !
    do j=2, LJg
       do i=2, LIg
          AvpExt(i,j) = Avp(  i-1,j-1)
          BgB (3,i,j) = 0.25d0*( Bfd(3,i-1,j-1) + Bfd(3,i  ,j-1) &
               &               + Bfd(3,i-1,j  ) + Bfd(3,i  ,j  ) )
       enddo
    enddo
    AvpExt(    1,:) = AvpExt(  2,:) + ( AvpExt(  2,:) - AvpExt(    3,:) )
    AvpExt(LIg+1,:) = AvpExt(LIg,:) + ( AvpExt(LIg,:) - AvpExt(LIg-1,:) )
    AvpExt(:,    1) = AvpExt(:,  2) + ( AvpExt(:,  2) - AvpExt(:,    3) )
    AvpExt(:,LJg+1) = AvpExt(:,LJg) + ( AvpExt(:,LJg) - AvpExt(:,LJg-1) )
    !   -- [1-3] fgB :: Avp -- !
    fgB(:,:,1)      = AvpExt(1:LIg,1:LJg)
    
    !  --- [2] Magnetic Field  :: B --- !
    !   -- [2-1] Inner Region       --  !
    dx1inv = 1.0d0 / dx1
    dx2inv = 1.0d0 / dx2
    do j=2, LJg
       do i=2, LIg
          BgB(1,i,j) = + ( AvpExt(i,j+1) - AvpExt(i,j) ) * dx2inv
          BgB(2,i,j) = - ( AvpExt(i+1,j) - AvpExt(i,j) ) * dx1inv
       enddo
    enddo
    !   -- [2-2] Boundary Region    --  !
    !    - Left & Right - !
    do j=1, LJg
       BgB(1,  1,j) = BgB(1,    2,j)
       BgB(2,  1,j) = BgB(2,    2,j)
       BgB(2,LIg,j) = BgB(2,LIg-1,j)
       BgB(3,  1,j) = BgB(3,    2,j)
       BgB(3,LIg,j) = BgB(3,LIg-1,j)
    enddo
    !    - Top & Bottom - !
    do i=1, LIg
       BgB(1,i,  1) = BgB(1,i,    2)
       BgB(1,i,LJg) = BgB(1,i,LJg-1)
       BgB(2,i,  1) = BgB(2,i,    2)
       BgB(3,i,  1) = BgB(3,i,    2)
    enddo
    !   -- [2-3] Normalization    --  !
    if ( normSW ) then
       absB2Max = 0.d0
       do j=2, LJg
          do i=2, LIg
             absB2Max = max( absB2Max, BgB(1,i,j)**2+BgB(2,i,j)**2+BgB(3,i,j)**2 )
          enddo
       enddo
       absB2Max = Bmax / sqrt( absB2Max )
       do j=1, LJg
          do i=1, LIg
             BgB(1,i,j) = absB2Max * BgB(1,i,j)
             BgB(2,i,j) = absB2Max * BgB(2,i,j)
             BgB(3,i,j) = absB2Max * BgB(3,i,j)
             fgB(i,j,1) = absB2Max * fgB(i,j,1)
          enddo
       enddo
    endif

    !  --- [3] Current Density :: J --- !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    do j=3, LJg-1
       do i=3, LIg-1
          JgB(1,i,j) = coef * ( + ( BgB(3,i  ,j  ) - BgB(3,  i,j-1) ) * dx2inv )
          JgB(2,i,j) = coef * ( - ( BgB(3,i  ,j  ) - BgB(3,i-1,j  ) ) * dx1inv )
          JgB(3,i,j) = coef * ( + ( BgB(2,i  ,j  ) - BgB(2,i-1,j  ) ) * dx1inv &
               &                - ( BgB(1,i  ,j  ) - BgB(1,i  ,j-1) ) * dx2inv )
       enddo
    enddo
    ! -- B.C. -- !
    JgB(:,  1,  :) = 0.d0
    JgB(:,  :,  1) = 0.d0
    JgB(:,  2,  :) = 0.d0
    JgB(:,  :,  2) = 0.d0
    JgB(:,LIg,  :) = 0.d0
    JgB(:,  :,LJg) = 0.d0

    !  --- [4] Rho and Pressure :: rho, p --- !
    call setPrsrRho

    return 
  end subroutine Cartesian__EgridField

  
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


  subroutine setPrsrRho
    use constants, only : vthcv, TiTe, mr, desiredBetapol, valfe, normType
    implicit none
    integer            :: i, j
    double precision   :: maxBpol2, pMax, Ttot, Ttotinv
    
    ! --- [1] Find Max Bpol --- !
    maxBpol2 = 0.d0
    do j=2, LJg
       do i=2, LIg
          maxBpol2 = max( maxBpol2, BgB(1,i,j)**2 + BgB(2,i,j)**2 )
       enddo
    enddo
    
    ! --- [2] set Pressure  --- !
    do j=1, LJg
       do i=1, LIg
          fgB(i,j,3) = desiredBetapol * 0.5d0 * maxBpol2
       enddo
    enddo
    write(6,'(a,f12.6)') ' Set Pressure as constant :: P0 ===  ', fgB(2,2,3)
    
    ! --- [3] set rho --- !
    pMax    = maxval( fgB(:,:,3) )
    if ( normType.eq.'PIC' ) Ttot = ( 1.d0 + TiTe ) * vthcv**2
    if ( normType.eq.'MHD' ) Ttot = pMax
    Ttotinv =  1.d0 / Ttot
    do j=1, LJg
       do i=1, LIg
          fgB(i,j,2) = Ttotinv * fgB(i,j,3)
       enddo
    enddo
    write(6,'(a,f12.6)') ' Set density  as constant :: n0 ===  ', fgB(2,2,2)
    
    ! --- [4] Get Electron Flow Velocity ---  !
    do j=1, LJg
       do i=1, LIg
          ugB(i,j,1)   = - JgB(1,i,j) / fgB(i,j,2)
          ugB(i,j,2)   = - JgB(2,i,j) / fgB(i,j,2)
          ugB(i,j,3)   = - JgB(3,i,j) / fgB(i,j,2)
       enddo
    enddo
    
    return
  end subroutine setPrsrRho

end module EgridEBMod
