module BgridRZMod
  use constants, only : N1, N2
  implicit none
  integer, parameter :: LIg = N1+1
  integer, parameter :: LJg = N2+1
  integer            :: bx_, by_, bz_, jx_, jy_, jz_
  integer            :: br_, bt_, jr_, jt_
  integer            :: at_ = 1, rh_ = 2, pr_ = 3
  double precision   :: rf(LJg+1), rfInv(LJg+1)
  double precision   :: rh(LJg+1), rhInv(LJg+1)
  double precision   :: BgB(3,LIg,LJg), JgB(3,LIg,LJg)
  double precision   :: fgB(3,LIg,LJg), ugB(3,LIg,LJg)
contains

  
  ! =================================================================== !
  ! ===  BgridField  ::  Bgrid B-Field Calculator                   === !
  ! =================================================================== !
  subroutine BgridField( coordinate, onGrid )
    implicit none
    character(1), intent(in) :: onGrid
    character(2), intent(in) :: coordinate

    if ( coordinate.eq.'XZ' ) then
       bx_=1; by_=2; bz_=3
       jx_=1; jy_=2; jz_=3
    endif
    if ( coordinate.eq.'RZ' ) then
       br_=1; bt_=2; bz_=3
       jr_=1; jt_=2; jz_=3
    endif
    if ( ( coordinate.eq.'XZ' ).and.( onGrid.eq.'B' ) ) call Cartesian__BgridField
    if ( ( coordinate.eq.'RZ' ).and.( onGrid.eq.'B' ) ) call Cylindrical__BgridField
    
    return
  end subroutine BgridField


  ! =================================================================== !
  ! ===  Cartesian__BgridField :: Cartesian Bgrid Field Calculator  === !
  ! =================================================================== !
  subroutine Cartesian__BgridField
    use constants, only : N1, N2, Bmax, valfe, normType, MergingType, normSW, symmBreak
    use variables, only : dx1, dx2, Avp, Bfd
    implicit none
    integer            :: i, j
    double precision   :: dx1Inv, dx2Inv, coef, absB2Max, BtSign1, BtSign2
    double precision   :: AvpExt(LIg,LJg)
    logical, parameter :: Flag__ReverseBt = .true.

    ! ------------------------------------- !
    ! --- [1] Extend Avp Field          --- !
    ! ------------------------------------- !
    !  -- [1-1] Interpolate AvpExt      --  !
    do j=2, LJg-1
       do i=2, LIg-1
          AvpExt (i,j) = 0.25d0 * ( Avp(i-1,j  ) + Avp(i  ,j  ) &
               &                  + Avp(i-1,j-1) + Avp(i  ,j-1) )
       enddo
    enddo
    !  -- [1-2] Substitute Bfd          --  !
    do j=2, LJg
       do i=2, LIg
          BgB(by_,i,j) = Bfd(3,i-1,j-1)
       enddo
    enddo
    !  -- [1-3] Boundary Treatment      --  !
    AvpExt(  1,  :) = AvpExt(    2,    :) + ( AvpExt(    2,    :) - AvpExt(    3,    :) )
    AvpExt(LIg,  :) = AvpExt(LIg-1,    :) + ( AvpExt(LIg-1,    :) - AvpExt(LIg-2,    :) )
    AvpExt(  :,  1) = AvpExt(    :,    2) + ( AvpExt(    :,    2) - AvpExt(    :,    3) )
    AvpExt(  :,LJg) = AvpExt(    :,LJg-1) + ( AvpExt(    :,LJg-1) - AvpExt(    :,LJg-2) )
    !  -- [1-4] Substitute Avp->fgB     --  !
    fgB(at_,:,:)    = AvpExt(1:LIg,1:LJg)
    
    ! ------------------------------------- !
    ! --- [2] Magnetic Field  :: B      --- !
    ! ------------------------------------- !
    !  -- [2-1] Determine Merging Mode  --  !
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
    !  -- [2-2] Instantly Reverse Bt    --  !
    if ( Flag__ReverseBt ) then
       BtSign1 = -1.d0 * BtSign1
       BtSign2 = -1.d0 * BtSign2
    endif
    !  -- [2-3] Apply BtSign            --  !
    do j=1, N2
       do i=1, N1/2
          Bfd(3,i,j) = BtSign1 * Bfd(3,i,j)
       enddo
    enddo
    do j=1, N2
       do i=N1/2+1, N1
          Bfd(3,i,j) = BtSign2 * Bfd(3,i,j)
       enddo
    enddo
    !  -- [2-4] Inner Region            --  !
    dx1Inv = 1.0d0 / dx1
    dx2Inv = 1.0d0 / dx2
    do j=2, LJg
       do i=2, LIg
          BgB(bx_,i,j) = - ( AvpExt(i,j) - AvpExt(i-1,j  ) ) * dx1Inv
          BgB(bz_,i,j) = + ( AvpExt(i,j) - AvpExt(i  ,j-1) ) * dx2Inv
       enddo
    enddo
    !  -- [2-5] Boundary Region         --  !
    !   - Left & Right - !
    BgB(bz_,  1,:) = BgB(bz_,    2,:)
    BgB(bx_,  1,:) = BgB(bx_,    2,:)
    BgB(by_,  1,:) = BgB(by_,    2,:)
    BgB(bz_,LIg,:) = BgB(bz_,LIg-1,:)
    !   - Top & Bottom - !
    BgB(bz_,:,  1) = BgB(bz_,:,    2)
    BgB(bx_,:,  1) = BgB(bx_,:,    2)
    BgB(by_,:,  1) = BgB(by_,:,    2)
    BgB(bx_,:,LJg) = BgB(bx_,:,LJg-1)
    !  -- [2-6] Normalization           --  !
    if ( normSW ) then
       absB2Max = 0.d0
       do j=2, LJg
          do i=2, LIg
             absB2Max = max( absB2Max, BgB(bx_,i,j)**2 + BgB(by_,i,j)**2 + BgB(bz_,i,j)**2 )
          enddo
       enddo
       absB2Max = Bmax / sqrt( absB2Max )
       do j=1, LJg
          do i=1, LIg
             BgB(bx_,i,j) = absB2Max * BgB(bx_,i,j)
             BgB(by_,i,j) = absB2Max * BgB(by_,i,j)
             BgB(bz_,i,j) = absB2Max * BgB(bz_,i,j)
             fgB(at_,i,j) = absB2Max * fgB(at_,i,j)
          enddo
       enddo
    endif

    ! ------------------------------------- !
    ! --- [3] Current Density :: J      --- !
    ! ------------------------------------- !
    !  -- [3-1] Inner Region            --  !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    do j=2, LJg-1
       do i=2, LIg-1
          JgB(jx_,i,j) = coef * ( - ( BgB(by_,i+1,j  ) - BgB(by_,i  ,j  ) ) * dx1Inv )
          JgB(jy_,i,j) = coef * ( + ( BgB(bx_,i+1,j  ) - BgB(bx_,i  ,j  ) ) * dx1Inv &
               &                  - ( BgB(bz_,i  ,j+1) - BgB(bz_,i  ,j  ) ) * dx2Inv )
          JgB(jz_,i,j) = coef * ( + ( BgB(by_,i+1,j  ) - BgB(by_,i  ,j  ) ) * dx2Inv )
       enddo
    enddo
    !  -- [3-2] Boundary Region         --  !
    JgB(  :,  1,  :) = 0.d0
    JgB(  :,  :,  1) = 0.d0
    JgB(jz_,  2,  :) = 0.d0
    JgB(jx_,  :,  2) = 0.d0
    JgB(  :,LIg,  :) = 0.d0
    JgB(  :,  :,LJg) = 0.d0

    ! ------------------------------------- !
    ! --- [4] Pressure & Rho            --- !
    ! ------------------------------------- !
    call setPrsrRho

    return 
  end subroutine Cartesian__BgridField


  ! =================================================================== !
  ! ===  Cylindrical__BgridField  :: Cylindrical Bgrid Field        === !
  ! =================================================================== !
  subroutine Cylindrical__BgridField
    use constants, only : N1, N2, normSW, MergingType, normType, symmBreak, Flag__ReverseBt
    use constants, only : x2Min, Bmax, valfe
    use variables, only : dx1, dx2, x2
    use variables, only : Avp, Bfd
    implicit none
    integer            :: i, j
    double precision   :: dz, dr, dzInv, drInv
    double precision   :: psiExt(1:LIg+1,1:LJg+1), AvpExt(1:LIg,1:LJg)
    double precision   :: BtSign1, BtSign2, absB2Max, coef

    ! ------------------------------------- !
    ! --- [0] Grid Making               --- !
    ! ------------------------------------- !
    dz    = dx1
    dr    = dx2
    dzInv = 1.d0 / dz
    drInv = 1.d0 / dr
    do j=2, LJg
       rf(j) =  x2(j-1)
    enddo
    rf(LJg+1) = rf(LJg) + dr
    if ( x2min.eq.0.d0 ) then
       rfInv(1) = 1.d0 / rf(1)
       rfInv(2) = 0.d0
       do j=3, LJg+1
          rfInv(j) = 1.d0 / rf(j)
       enddo
    else
       do j=1, LJg+1
          rfInv(j) = 1.d0 / rf(j)
       enddo
    endif
    do j=1, LJg+1
       rh(j)    = rf(j) + 0.5d0*dr
       rhInv(j) = 1.d0 / rh(j)
    enddo
    ! ------------------------------------- !
    ! --- [1] Extend Avp Field          --- !
    ! ------------------------------------- !
    !  -- [1-1] Interpolate AvpExt      --  !    
    do j=2, LJg
       do i=2, LIg
          psiExt(i,j) = rf(j)*Avp(i-1,j-1)
       enddo
    enddo
    psiExt(    1,:) = psiExt(  2,:) + ( psiExt(  2,:) - psiExt(    3,:) )
    psiExt(LIg+1,:) = psiExt(LIg,:) + ( psiExt(LIg,:) - psiExt(LIg-1,:) )
    psiExt(:,    1) = psiExt(:,  2) + ( psiExt(:,  2) - psiExt(:,    3) )
    psiExt(:,LJg+1) = psiExt(:,LJg) + ( psiExt(:,LJg) - psiExt(:,LJg-1) )
    do j=1, LJg
       do i=1, LIg
          fgB(at_,i,j) = 0.25d0 * ( psiExt(i  ,j+1) + psiExt(i+1,j+1) &
               &                  + psiExt(i  ,j  ) + psiExt(i+1,j  ) )
          AvpExt(i,j)  = rhInv(j) * fgB(at_,i,j)
       enddo
    enddo
    ! ------------------------------------- !
    ! --- [2] Magnetic Field  :: B      --- !
    ! ------------------------------------- !
    !  -- [2-1] Determine Merging Mode  --  !
    BgB(:,:,:) = 0.d0
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
    !  -- [2-2] Instantly Reverse Bt    --  !
    if ( Flag__ReverseBt ) then
       BtSign1 = -1.d0 * BtSign1
       BtSign2 = -1.d0 * BtSign2
    endif
    !  -- [2-3] Apply BtSign            --  !
    do j=1, N2
       do i=1, N1/2
          BgB(bt_,i+1,j+1) = BtSign1 * Bfd(3,i,j)
       enddo
    enddo
    do j=1, N2
       do i=N1/2+1, N1
          BgB(bt_,i+1,j+1) = BtSign2 * Bfd(3,i,j)
       enddo       
    enddo
    !  -- [2-4] Inner Region            --  !
    do j=2, LJg
       do i=2, LIg
          BgB(br_,i,j) = - rhinv(j) * dzinv * ( fgB(at_,i,j) - fgB(at_,i-1,j) )
       enddo
    enddo ! - Br - !
    do j=2, LJg
       do i=2, LIg
          BgB(bz_,i,j) = + rfinv(j) * drinv * ( fgB(at_,i,j) - fgB(at_,i,j-1) )
       enddo
    enddo ! - Bz - !
    !  -- [2-5] Boundary Region         --  !
    !   -  Inner  -   !
    if ( x2Min.eq.0.d0 ) then ! -( Stokes's Theorem )- !
       ! do i=2, LIg
       !    BgB(bz_,i,2) = 4.d0 * drinv * AvpExt(i,2)
       ! enddo !! --- Gauge invariant ??? --- !!
       do i=2, LIg
          BgB(bz_,i,2) = BgB(bz_,i,3)
       enddo
    endif
    if ( x2Min.eq.0.d0 ) then
       do i=2, LIg
          BgB(br_,i,1) = - BgB(br_,i,2)
          BgB(bt_,i,1) = - BgB(bt_,i,2)
          BgB(bz_,i,1) = + BgB(bz_,i,3)
       enddo
    else
       do i=2, LIg
          BgB(br_,i,1) = - rhinv(1) * dzinv * ( fgB(at_,i,1) - fgB(at_,i-1,1) )
          BgB(bt_,i,1) =   rhinv(1) * rh(2) *   BgB(bt_,i,2)
          BgB(bz_,i,1) =                        BgB(bz_,i,2)
       enddo
    endif
    !   - Left -   !
    do j=1, LJg
       BgB(br_,1,j) = BgB(br_,2,j)
       BgB(bt_,1,j) = BgB(bt_,2,j)
       BgB(bz_,1,j) = BgB(bz_,2,j)
    enddo
    !  -- [2-6] fgB substitution        --  !
    do j=1, LJg
       do i=1, LIg
          fgB(at_,i,j) = rhInv(j) * fgB(at_,i,j)
       enddo
    enddo
    !  -- [2-7] Normalization           --  !
    if ( normSW ) then
       absB2Max = 0.d0
       do j=2, LJg
          do i=2, LIg
             absB2Max = max( absB2Max, BgB(br_,i,j)**2+BgB(bt_,i,j)**2+BgB(bz_,i,j)**2 )
          enddo
       enddo
       absB2Max = Bmax / sqrt( absB2Max )
       do j=1, LJg
          do i=1, LIg
             BgB(br_,i,j) = absB2Max * BgB(br_,i,j)
             BgB(bt_,i,j) = absB2Max * BgB(bt_,i,j)
             BgB(bz_,i,j) = absB2Max * BgB(bz_,i,j)
             fgB(at_,i,j) = absB2Max * fgB(at_,i,j)
          enddo
       enddo
    endif

    ! ------------------------------------- !
    ! --- [3] Current Density :: J      --- !
    ! ------------------------------------- !
    !  -- [3-1] Inner Region            --  !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    do j=2, LJg-1
       do i=2, LIg-1
          JgB(jr_,i,j) = coef * ( - (         BgB(bt_,i+1,j  ) -         BgB(bt_,i  ,j  ) ) * dzInv )
          JgB(jt_,i,j) = coef * ( + (         BgB(br_,i+1,j  ) -         BgB(br_,i  ,j  ) ) * dzInv &
               &                  - (         BgB(bz_,i  ,j+1) -         BgB(bz_,i  ,j  ) ) * drInv )
          JgB(jz_,i,j) = coef * ( + ( rf(j+1)*BgB(bt_,i  ,j+1) - rf(j  )*BgB(bt_,i  ,j  ) ) * drInv * rhInv(j) )
       enddo
    enddo
    !  -- [3-2] Boundary Region         --  !
    JgB(  :,  1,  :) = 0.d0
    JgB(  :,  :,  1) = 0.d0
    JgB(jr_,  :,  2) = 0.d0
    JgB(jt_,  :,  2) = 0.d0
    JgB(  :,LIg,  :) = 0.d0
    JgB(  :,  :,LJg) = 0.d0
    JgB(jz_,  2,  :) = 0.d0
    JgB(jt_,  2,  :) = 0.d0
    if ( x2Min.eq.0.d0 ) then
       JgB(jr_,:, 2) = 0.d0
       JgB(jt_,:, 2) = 0.d0
    endif

    ! ------------------------------------- !
    ! --- [4] Pressure & Rho            --- !
    ! ------------------------------------- !    
    call setPrsrRho

    return 
  end subroutine Cylindrical__BgridField


  ! =================================================================== !
  ! ===  setPrsrRho  ::  set Pressure & Rho                         === !
  ! =================================================================== !
  subroutine setPrsrRho
    use constants, only          : vthcv, TiTe, desiredBeta, normType, valfe, BetaMode, coordinate, myRank
    use variables, only          : dx1, dx2
    implicit none
    integer                     :: i, j
    double precision            :: PMag, Svol, pMax, Ttot, coef, vol(LJg)
    double precision, parameter :: pi       = 4.d0*atan( 1.d0 )
    double precision, parameter :: onethird = 1.d0 / 3.d0

    ! ------------------------------------- !
    ! --- [1] Define Magnetic Pressure  --- !
    ! ------------------------------------- !
    if ( trim(BetaMode).eq."MaxOfabsB" ) then
       !  -- [1-1] Find Max |B|^2 --  !
       PMag = 0.d0
       do j=2, LJg
          do i=2, LIg
             PMag = max( PMag, BgB(1,i,j)**2 + BgB(2,i,j)**2 + BgB(3,i,j)**2 )
          enddo
       enddo
    endif
    if ( trim(BetaMode).eq."EnergyOfB" ) then    
       !  -- [1-2] Calc. Total Energy --  !
       vol(:) = 0.d0
       if ( coordinate.eq."XZ" ) then
          do j=2, LJg
             vol(j) = dx1*dx2
          enddo
       endif
       if ( coordinate.eq."RZ" ) then
          do j=2, LJg
             vol(j) = dx1*dx2*2.d0*pi*rh(j)
          enddo
       endif
       do j=2, LJg
          do i=2, LIg
             if ( BgB(1,i,j).ne.BgB(1,i,j) ) write(6,*) "Nan Detec. :: ( i, j, BgB1 )", i, j, BgB(1,i,j)
             if ( BgB(2,i,j).ne.BgB(2,i,j) ) write(6,*) "Nan Detec. :: ( i, j, BgB2 )", i, j, BgB(2,i,j)
             if ( BgB(3,i,j).ne.BgB(3,i,j) ) write(6,*) "Nan Detec. :: ( i, j, BgB3 )", i, j, BgB(3,i,j)
          enddo
       enddo
       Svol = 0.d0
       PMag = 0.d0
       do j=2, LJg
          do i=2, LIg
             Svol = Svol + vol(j)
             PMag = PMag + ( BgB(1,i,j)**2 + BgB(2,i,j)**2 + BgB(3,i,j)**2 )*vol(j)
          enddo
       enddo
       PMag = PMag / Svol
    endif
    ! ------------------------------------- !
    ! --- [2] set Pressure Profile      --- !
    ! ------------------------------------- !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    pMax  = desiredBeta * onethird * PMag * coef
    do j=1, LJg
       do i=1, LIg
          fgB(pr_,i,j) = pMax
       enddo
    enddo
    ! ------------------------------------- !
    ! --- [3] set Rho Profile           --- !
    ! ------------------------------------- !
    if ( normType.eq.'PIC' ) Ttot = ( 1.d0 + TiTe ) * vthcv**2
    if ( normType.eq.'MHD' ) Ttot = pMax
    coef    = 1.d0 / Ttot
    do j=1, LJg
       do i=1, LIg
          fgB(rh_,i,j) = coef * fgB(pr_,i,j)
       enddo
    enddo
    ! ------------------------------------- !
    ! --- [4] set Flow Velocity (e)     --- !
    ! ------------------------------------- !
    do j=1, LJg
       do i=1, LIg
          ugB(1,i,j)   = - JgB(1,i,j) / fgB(rh_,i,j)
          ugB(2,i,j)   = - JgB(2,i,j) / fgB(rh_,i,j)
          ugB(3,i,j)   = - JgB(3,i,j) / fgB(rh_,i,j)
       enddo
    enddo
    ! ------------------------------------- !
    ! --- [5] Display Summary           --- !
    ! ------------------------------------- !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'      ) '[ setPrsrRho      @ BgridMod ]'
       write(6,'(10x,a,f12.6)') 'Beta               :: beta        ===', desiredBeta
       write(6,'(10x,a,f12.6)') 'Thermal  Pressure  :: P(Thermal ) ===', pMax
       write(6,'(10x,a,f12.6)') 'Magnetic Pressure  :: P(Magnetic) ===', pMag
       write(6,'(10x,a,f12.6)') 'Total Temperature  :: T           ===', Ttot
       write(6,'(10x,a,f12.6)') '                   :: 1 / T       ===', 1.d0 / Ttot
       write(6,'(10x,a,f12.6)') 'Density            :: n0          ===', fgB(rh_,LIg/2,LJg/2)
    endif
    
    return
  end subroutine setPrsrRho
  

end module BgridRZMod
