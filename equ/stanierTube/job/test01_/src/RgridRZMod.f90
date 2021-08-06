module RgridRZMod
  use constants, only : N1, N2
  implicit none
  integer, parameter :: LIr = N1+1
  integer, parameter :: LJr = N2+1
  integer            :: bx_, by_, bz_, jx_, jy_, jz_
  integer            :: br_, bt_, jr_, jt_
  integer            :: at_ = 1, rh_ = 2, pr_ = 3
  double precision   :: rf(LJr+1), rfInv(LJr+1)
  double precision   :: rh(LJr+1), rhInv(LJr+1)
  double precision   :: BgR(3,LIr,LJr), JgR(3,LIr,LJr)
  double precision   :: fgR(3,LIr,LJr), ugR(3,LIr,LJr)
contains


  ! =================================================================== !
  ! ===  RgridField  :: B-Field on Rgrid calculator                 === !
  ! =================================================================== !
  subroutine RgridField( coordinate, onGrid )
    implicit none
    character(1), intent(in) :: onGrid
    character(2), intent(in) :: coordinate

    ! ------------------------------------- !
    ! --- [1] coordinate Settings       --- !
    ! ------------------------------------- !
    if ( coordinate.eq.'XZ' ) then
       bx_=1; by_=2; bz_=3
       jx_=1; jy_=2; jz_=3
    endif
    if ( coordinate.eq.'RZ' ) then
       br_=1; bt_=2; bz_=3
       jr_=1; jt_=2; jz_=3
    endif
    ! ------------------------------------- !
    ! --- [2] call Calculator           --- !
    ! ------------------------------------- !
    if ( ( coordinate.eq.'XZ' ).and.( onGrid.eq.'R' ) ) call Cartesian__RegularGrid
    if ( ( coordinate.eq.'RZ' ).and.( onGrid.eq.'R' ) ) call Cylindric__RegularGrid
    return
  end subroutine RgridField
  

  ! =================================================================== !
  ! ===  Cartesian__RegularGrid  ::  B,J on Rgrid ( cartesian )     === !
  ! =================================================================== !
  subroutine Cartesian__RegularGrid
    use constants, only : Bmax, valfe, normType, normSW, symmBreak
    use variables, only : dx1, dx2, Avp, Bfd
    implicit none
    integer            :: i, j
    double precision   :: dx1Inv, dx2Inv, coef, absB2Max, BtSign1, BtSign2
    double precision   :: AvpExt(LIr+1,LJr+1)
    double precision   :: Bfd_p(LIr+1,LJr+1), Avp_p(LIr+1,LJr+1)
    
    ! ------------------------------------- !
    ! --- [1] Extend Avp Field          --- !
    ! ------------------------------------- !
    !  -- [1-1] Change Polarity (Bt)    --  !
    call changePolarityBt( Avp_p, Bfd_p )
    !  -- [1-2] Get AvpExt              --  !
    do j=2, LJr
       do i=2, LIr
          AvpExt( i,j) = Avp_p(i,j)
          BgR(by_,i,j) = Bfd_p(i,j)
       enddo
    enddo
    AvpExt(    1,:) = AvpExt(  2,:) + ( AvpExt(  2,:) - AvpExt(    3,:) )
    AvpExt(LIr+1,:) = AvpExt(LIr,:) + ( AvpExt(LIr,:) - AvpExt(LIr-1,:) )
    AvpExt(:,    1) = AvpExt(:,  2) + ( AvpExt(:,  2) - AvpExt(:,    3) )
    AvpExt(:,LJr+1) = AvpExt(:,LJr) + ( AvpExt(:,LJr) - AvpExt(:,LJr-1) )
    BgR(by_,1,:)    = BgR(by_,2,:)
    BgR(by_,:,1)    = BgR(by_,:,2)
    !   -- [1-3] fgR :: Flux            --  !
    fgR(at_,:,:)    = AvpExt(1:LIr,1:LJr)

    ! ------------------------------------- !
    ! --- [2] Magnetic Field :: B       --- !
    ! ------------------------------------- !
    !  -- [2-1] Calculate B-field       --  !
    dx1Inv = 1.d0 / ( 2.d0 * dx1 )
    dx2Inv = 1.d0 / ( 2.d0 * dx2 )
    do j=2, LJr
       do i=2, LIr
          BgR(bx_,i,j) = - ( AvpExt(i+1,j) - AvpExt(i-1,j) ) * dx1Inv
          BgR(bz_,i,j) = + ( AvpExt(i,j+1) - AvpExt(i,j-1) ) * dx2Inv
       enddo
    enddo
    BgR(1:2,  2,  2) = 0.d0
    BgR(1:2,  2,LJr) = 0.d0
    BgR(1:2,LIr,  2) = 0.d0
    BgR(1:2,LIr,LJr) = 0.d0
    !  -- [2-2] Outer Boundary          --  !
    BgR(bx_,1,:) = BgR(bx_,2,:)
    BgR(bz_,1,:) = BgR(bz_,2,:)
    BgR(by_,1,:) = BgR(by_,2,:)
    BgR(bx_,:,1) = BgR(bx_,:,2)
    BgR(bz_,:,1) = BgR(bz_,:,2)
    BgR(by_,:,1) = BgR(by_,:,2)
    BgR(by_,:,1) = BgR(by_,:,2)
    !  -- [2-3] Normalize B-field       --  !
    if ( normSW ) call normalizeField

    ! ------------------------------------- !
    ! --- [3] Current Density :: J      --- !
    ! ------------------------------------- !
    !  -- [3-1] Inner Region            --  !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    do j=3, LJr-1
       do i=3, LIr-1
          JgR(jx_,i,j) = coef * ( - ( BgR(by_,i+1,j  ) - BgR(by_,i-1,j  ) ) * dx1Inv )
          JgR(jy_,i,j) = coef * ( + ( BgR(bz_,i  ,j+1) - BgR(bz_,i  ,j-1) ) * dx2Inv &
               &                  - ( BgR(bx_,i+1,j  ) - BgR(bx_,i-1,j  ) ) * dx1Inv )
          JgR(jz_,i,j) = coef * ( + ( BgR(by_,i  ,j+1) - BgR(by_,i  ,j-1) ) * dx2Inv )
       enddo
    enddo
    !  -- [3-2] Boundary Region         --  !
    JgR(:,  1,:) = 0.d0
    JgR(:,  2,:) = 0.d0
    JgR(:,LIr,:) = 0.d0
    JgR(:,:,  1) = 0.d0
    JgR(:,:,  2) = 0.d0
    JgR(:,:,LJr) = 0.d0

    ! ------------------------------------- !
    ! --- [4] Rho and Pressure          --- !
    ! ------------------------------------- !
    call setPrsrRho

    return
  end subroutine Cartesian__RegularGrid


  ! =================================================================== !
  ! ===  Cylindric__RegularGrid  :: J, B on Rgrid ( cylindrical )   === !
  ! =================================================================== !
  subroutine Cylindric__RegularGrid
    use constants, only : Bmax, valfe, x2min, normType, normSW
    use variables, only : dx1, dx2, Avp, Bfd
    implicit none
    integer            :: i, j, ip, jp
    double precision   :: dx1Inv, dx2Inv, coef
    double precision   :: AvpExt(LIr+1,LJr+1)
    double precision   :: Avp_p(LIr+1,LJr+1), Bfd_p(LIr+1,LJr+1)

    ! ------------------------------------- !
    ! --- [0] Grid Making               --- !
    ! ------------------------------------- !
    call gridMaking
    ! ------------------------------------- !
    ! --- [1] Bt Polarity               --- !
    ! ------------------------------------- !
    call changePolarityBt( Avp_p, Bfd_p )
    ! ------------------------------------- !
    ! --- [2] Extend Avp Field          --- !
    ! ------------------------------------- !
    !  -- [2-1] Main Region             --  !
    do j=2, LJr
       do i=2, LIr
          AvpExt (i,j) = 0.25d0*(   rf(j  )*( Avp_p(i+1,j  ) + Avp_p(i,j  ) ) &
               &                  + rf(j+1)*( Avp_p(i+1,j+1) + Avp_p(i,j+1) ) ) * rhInv(j)
          BgR(bt_,i,j) = 0.25d0*(   rf(j  )*( Bfd_p(i+1,j  ) + Bfd_p(i,j  ) ) &
               &                  + rf(j+1)*( Bfd_p(i+1,j+1) + Bfd_p(i,j+1) ) ) * rhInv(j)
       enddo
    enddo
    !  -- [2-2] Boundary Region         --  !
    if ( x2min.eq.0.d0 ) then
       AvpExt(    1,:) =   AvpExt(  2,:) + ( AvpExt(  2,:) - AvpExt(    3,:) )
       AvpExt(    :,1) =   AvpExt(  :,2)
       AvpExt(:,LJr+1) =   ( rf(LJr)*AvpExt(:,LJr) &
            &            + ( rf(LJr)*AvpExt(:,LJr) - rf(LJr-1)*AvpExt(:,LJr-1) ) ) * rfInv(LJr+1)
       AvpExt(LIr+1,:) =   AvpExt(LIr,:) + ( AvpExt(LIr,:) - AvpExt(LIr-1,:) )
       BgR(bt_,  1,:)  =   BgR(bt_,    2,:)
       BgR(bt_,:,  1)  = + BgR(bt_,:,    2)
       ! BgR(bt_,LIr,:)  =   BgR(bt_,LIr-1,:)
       ! BgR(bt_,:,LJr)  = + BgR(bt_,:,LJr-1)
    else
       AvpExt(:,    1) =   ( rf(  2)*AvpExt(:,  2) &
            &            + ( rf(  2)*AvpExt(:,  2) - rf(    3)*AvpExt(:,    3) ) ) * rfInv(    1)
       AvpExt(:,LJr+1) =   ( rf(LJr)*AvpExt(:,LJr) &
            &            + ( rf(LJr)*AvpExt(:,LJr) - rf(LJr-1)*AvpExt(:,LJr-1) ) ) * rfInv(LJr+1)
       AvpExt(    1,:) = AvpExt(  2,:) + ( AvpExt(  2,:) - AvpExt(    3,:) )
       AvpExt(LIr+1,:) = AvpExt(LIr,:) + ( AvpExt(LIr,:) - AvpExt(LIr-1,:) )
       BgR(bt_,1,:)    = BgR(bt_,2,:)
       BgR(bt_,:,1)    = rf(2) * rfInv(1) * BgR(bt_,:,2)
    endif
    !  -- [2-3] Avp substitution        --  !
    do j=1, LJr
       do i=1, LIr
          fgR(at_,i,j) = AvpExt(i,j)
       enddo
    enddo
    fgR(at_,  :,LJr)   = fgR(at_,    :,LJr-1)
    fgR(at_,LIr,  :)   = fgR(at_,LIr-1,    :)
    fgR(at_,  :,  1)   = fgR(at_,    :,    2)
    fgR(at_,  1,  :)   = fgR(at_,    2,    :)
    
    ! ------------------------------------- !
    ! --- [3] Magnetic Field :: B       --- !
    ! ------------------------------------- !
    !  -- [3-1] Calculate B-field       --  !
    dx1Inv = 1.d0 / ( 2.d0 * dx1 )
    dx2Inv = 1.d0 / ( 2.d0 * dx2 )
    do j=2, LJr
       do i=2, LIr
          BgR(br_,i,j) = - (         AvpExt(i+1,j) -         AvpExt(i-1,j) )            * dx1Inv
          BgR(bz_,i,j) = + ( rh(j+1)*AvpExt(i,j+1) - rh(j-1)*AvpExt(i,j-1) ) * rhInv(j) * dx2Inv
       enddo
    enddo
    !  -- [3-2] Boundary Region         --  !
    if ( rf(2).eq.0.d0 ) then
       do i=2, LIr
          BgR(bz_,i,2) = ( rh(3)*AvpExt(i,3) - rh(2)*AvpExt(i,2) ) * rfInv(3) * dx2Inv*2.d0
       enddo
    endif
    BgR(br_,1,:) = BgR(br_,2,:)
    BgR(bt_,1,:) = BgR(bt_,2,:)
    BgR(bz_,1,:) = BgR(bz_,2,:)
    BgR(br_,:,1) = BgR(br_,:,2)
    BgR(bt_,:,1) = BgR(bt_,:,2)
    BgR(bz_,:,1) = BgR(bz_,:,3)
    !  -- [3-3] Normalization           --  !
    if ( normSW ) call normalizeField

    ! ------------------------------------- !
    ! --- [4] Current Density :: J      --- !
    ! ------------------------------------- !
    !  -- [4-1] Main Region             --  !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    do j=3, LJr-1
       do i=3, LIr-1
          JgR(jr_,i,j) = coef * ( - dx1Inv * ( BgR(bt_,i+1,j) - BgR(bt_,i-1,j) ) )
          JgR(jt_,i,j) = coef * ( + dx1Inv * ( BgR(br_,i+1,j) - BgR(br_,i-1,j) ) &
               &                  - dx2Inv * ( BgR(bz_,i,j+1) - BgR(bz_,i,j-1) ) )
          JgR(jz_,i,j) = coef * (   + dx2Inv * rfInv(j) &
               &              * (   rf(j+1)*BgR(bt_,i,j+1) - rf(j-1)*BgR(bt_,i,j-1) ) )
       enddo
    enddo
    !  -- [4-2] Boundary Region         --  !
    JgR(:,  1,:) = 0.d0
    JgR(:,  2,:) = 0.d0
    JgR(:,LIr,:) = 0.d0
    JgR(:,:,  1) = 0.d0
    JgR(:,:,  2) = 0.d0
    JgR(:,:,LJr) = 0.d0

    ! ------------------------------------- !
    ! --- [5] Rho and Pressure          --- !
    ! ------------------------------------- !
    call setPrsrRho
    
    return
  end subroutine Cylindric__RegularGrid
  

  ! =================================================================== !
  ! ===  setPrsrRho  ::  set Pressure & Rho                         === !
  ! =================================================================== !
  subroutine setPrsrRho
    use constants, only          : vthcv, TiTe, desiredBeta, normType, valfe, BetaMode, coordinate, myRank
    use variables, only          : dx1, dx2
    implicit none
    integer                     :: i, j
    double precision            :: PMag, Svol, pMax, Ttot, coef, vol(LJr)
    double precision, parameter :: pi       = 4.d0*atan( 1.d0 )
    double precision, parameter :: onethird = 1.d0 / 3.d0

    ! ------------------------------------- !
    ! --- [1] Define Magnetic Pressure  --- !
    ! ------------------------------------- !
    if ( trim(BetaMode).eq."MaxOfabsB" ) then
       !  -- [1-1] Find Max |B|^2 --  !
       PMag = 0.d0
       do j=2, LJr
          do i=2, LIr
             PMag = max( PMag, BgR(1,i,j)**2 + BgR(2,i,j)**2 + BgR(3,i,j)**2 )
          enddo
       enddo
    endif
    if ( trim(BetaMode).eq."EnergyOfB" ) then    
       !  -- [1-2] Calc. Total Energy --  !
       vol(:) = 0.d0
       if ( coordinate.eq."XZ" ) then
          do j=2, LJr
             vol(j) = dx1*dx2
          enddo
       endif
       if ( coordinate.eq."RZ" ) then
          do j=2, LJr
             vol(j) = dx1*dx2*2.d0*pi*rh(j)
          enddo
       endif
       Svol = 0.d0
       PMag = 0.d0
       do j=2, LJr
          do i=2, LIr
             Svol = Svol + vol(j)
             PMag = PMag + ( BgR(1,i,j)**2 + BgR(2,i,j)**2 + BgR(3,i,j)**2 )*vol(j)
          enddo
       enddo
       PMag = PMag / Svol
    endif
    ! ------------------------------------- !
    ! --- [2] set Pressure Profile      --- !
    ! ------------------------------------- !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    pMax  = desiredBeta * onethird*PMag * coef
    do j=1, LJr
       do i=1, LIr
          fgR(pr_,i,j) = pMax
       enddo
    enddo
    ! ------------------------------------- !
    ! --- [3] set Rho Profile           --- !
    ! ------------------------------------- !
    if ( normType.eq.'PIC' ) Ttot = ( 1.d0 + TiTe ) * vthcv**2
    if ( normType.eq.'MHD' ) Ttot = pMax
    coef =  1.d0 / Ttot
    do j=1, LJr
       do i=1, LIr
          fgR(rh_,i,j) = coef * fgR(pr_,i,j)
       enddo
    enddo
    ! ------------------------------------- !
    ! --- [4] set Flow Velocity (e)     --- !
    ! ------------------------------------- !
    do j=1, LJr
       do i=1, LIr
          ugR(1,i,j)   = - JgR(1,i,j) / fgR(rh_,i,j)
          ugR(2,i,j)   = - JgR(2,i,j) / fgR(rh_,i,j)
          ugR(3,i,j)   = - JgR(3,i,j) / fgR(rh_,i,j)
       enddo
    enddo

    ! ------------------------------------- !
    ! --- [5] Display Summary           --- !
    ! ------------------------------------- !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'      ) '[ setPrsrRho      @ RgridMod ]'
       write(6,'(2x,a)'      ) '[ setPrsrRho      @ RgridMod ]'
       write(6,'(10x,a,f12.6)') 'Beta               :: beta        ===', desiredBeta
       write(6,'(10x,a,f12.6)') 'Thermal  Pressure  :: P(Thermal ) ===', pMax
       write(6,'(10x,a,f12.6)') 'Magnetic Pressure  :: P(Magnetic) ===', pMag
       write(6,'(10x,a,f12.6)') 'Total Temperature  :: T           ===', Ttot
       write(6,'(10x,a,f12.6)') '                   :: 1 / T       ===', 1.d0 / Ttot
       write(6,'(10x,a,f12.6)') 'Pressure (const.)  :: P0          ===', fgR(pr_,LIr/2,LJr/2)
       write(6,'(10x,a,f12.6)') 'Density  (const.)  :: n0          ===', fgR(rh_,LIr/2,LJr/2)
    endif

    return
  end subroutine setPrsrRho


  ! =================================================================== !
  ! ===  gridMaking  ::  make cylindrical Grid Axis                 === !
  ! =================================================================== !
  subroutine gridMaking
    use constants, only : x2min
    use variables, only : dx2
    implicit none
    integer            :: j
    
    do j=1, LJr+1
       rf(j) = dx2*dble(j-2) + x2min
       rh(j) = dx2*dble(j-2) + x2min + 0.5d0*dx2
    enddo
    rfInv(:) = 0.d0
    rhInv(:) = 0.d0
    do j=1, LJr+1
       if ( rf(j).ne.0.d0 ) rfInv(j) = 1.d0 / rf(j)
       if ( rh(j).ne.0.d0 ) rhInv(j) = 1.d0 / rh(j)
    enddo
    return
  end subroutine gridMaking


  ! =================================================================== !
  ! ===  normalizeField  :: normalization of ( BgR/fgR )            === !
  ! =================================================================== !
  subroutine normalizeField
    use constants, only : Bmax
    implicit none
    integer            :: i, j
    double precision   :: absB2Max

    absB2Max = 0.d0
    do j=1, LJr
       do i=1, LIr
          absB2Max = max( absB2Max, BgR(1,i,j)**2 + BgR(2,i,j)**2 + BgR(3,i,j)**2 )
       enddo
    enddo
    absB2Max = Bmax / sqrt( absB2Max )
    do j=1, LJr
       do i=1, LIr
          BgR(  1,i,j) = absB2Max * BgR(  1,i,j)
          BgR(  2,i,j) = absB2Max * BgR(  2,i,j)
          BgR(  3,i,j) = absB2Max * BgR(  3,i,j)
          fgR(at_,i,j) = absB2Max * fgR(at_,i,j)
       enddo
    enddo
    return
  end subroutine normalizeField


  ! =================================================================== !
  ! ===  changePolarityBt  :: arrange Bt :: CaseIO of CHSM / CoHx   === !
  ! =================================================================== !
  subroutine changePolarityBt( Avp_p, Bfd_p )
    use constants, only : N1, N2, symmBreak, MergingType, Flag__ReverseBt
    use variables, only : Avp, Bfd
    implicit none
    integer                       :: i, j, ip, jp
    double precision              :: BtSign1, BtSign2
    double precision, intent(out) :: Avp_p(LIr+1,LJr+1), Bfd_p(LIr+1,LJr+1)
    
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
          ip = i+1
          jp = j+1
          Bfd_p(ip,jp) = BtSign1 * Bfd(3,i,j)
       enddo
    enddo
    do j=1, N2
       do i=N1/2+1, N1
          ip = i+1
          jp = j+1
          Bfd_p(ip,jp) = BtSign2 * Bfd(3,i,j)
       enddo
    enddo
    do j=1, N2
       do i=1, N1
          ip = i+1
          jp = j+1
          Avp_p(ip,jp) = Avp(i,j)
       enddo
    enddo

    return
  end subroutine changePolarityBt

  
end module RgridRZMod
