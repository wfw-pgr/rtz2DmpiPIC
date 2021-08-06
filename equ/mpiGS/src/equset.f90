module equset
  use constants, only : Nr, Nz, Nmid
  implicit none
  integer, parameter :: LI = 2*Nz + Nmid + 1
  integer, parameter :: LJ = Nr+1
  double precision   :: Bfd(3,LI,LJ), Jcr(3,LI,LJ)
  double precision   :: frp(LI,LJ,3), uvc(LI,LJ,3) 
  double precision   :: rf(LJ), rh(LJ), rfinv(LJ), rhinv(LJ)
  ! ------------------------------------ !
  ! ---------     Variables    --------- !
  ! ------------------------------------ !
  !                                      !
  ! --- Bfd -- Bfield             ------ !
  ! --- Jcr -- J current          ------ !
  ! --- frp -- ( Flux, Rho, Pressure )-- !
  ! --- uvc -- u vector           ------ !
  ! ------------------------------------ !
contains

  subroutine Set2GS_Egrid

    use constants, only : Nr, Nz, Nmid, rmin, rhofloor, TiTe, vthcv, valfe
    use variables, only : dr, dz, drinv, dzinv, prs
    use Field    , only : Br, Bz, Bt
    implicit none
    integer            :: i, j, ip2
    double precision   :: dBr, dBt, Nmidinv, valfe2
    double precision   :: Ttot, Ttotinv
    
    !  --- [0] Preparation   ---  !
    ip2         = Nz + Nmid
    Nmidinv     = 1.d0 / dble( Nmid )
    do j=1, LJ
       rf(j)    = dr * dble(j-2) + rmin
       rfinv(j) = 1.d0 / rf(j)
    enddo
    do j=1, LJ-1
       rh(j)    = 0.5d0 * ( rf(j) + rf(j+1) )
       rhinv(j) = 1.d0 / rh(j)
    enddo
    if ( rmin.eq.0.d0 ) then
       rf(2)    = 0.d0
       rfinv(2) = 0.d0
    end if
    rh(LJ)      = rf(LJ) + 0.5d0 * dr
    rhinv(LJ)   = 1.d0 / rh(LJ)
    Bfd(:,:,:)  = 0.d0
    frp(:,:,:)  = 0.d0
    
    !  --- [1] Set Equilibrium plasma  ---  !
    do j=1, Nr+1
       do i=1, Nz
          Bfd(1,    i,j) =  Br(i-1,j-1)
          Bfd(1,ip2+i,j) =  Br(i-1,j-1)
          Bfd(2,    i,j) =  Bt(i-1,j-1)
          Bfd(2,ip2+i,j) =  Bt(i-1,j-1)
          Bfd(3,    i,j) =  Bz(i-1,j-1)
          Bfd(3,ip2+i,j) =  Bz(i-1,j-1)
       enddo          
    enddo
    do j=2, Nr+1
       do i=2, Nz+1
          frp(    i,j,3) = prs(i-1,j-1)
          frp(ip2+i,j,3) = prs(i-1,j-1)
       enddo
    enddo

    !  --- [2] Linear Interpolation :: Br  ---    !
    do j=1, Nr+1
       
       dBr = ( Bfd(1,Nz+Nmid+1,j) - Bfd(1,Nz,j) ) * Nmidinv
       dBt = ( Bfd(2,Nz+Nmid+1,j) - Bfd(2,Nz,j) ) * Nmidinv
       
       do i=Nz+1, Nz+Nmid
          Bfd(1,i,j) = Bfd(1,Nz,j) + dBr * dble(i-(Nz+1))
          Bfd(2,i,j) = Bfd(2,Nz,j) + dBt * dble(i-(Nz+1))
       enddo
    enddo
    
    !  --- [3] psi Calculation from z = zmin ---  !
    do j=2, LJ-1
       frp(2,j+1,1)    = frp(2,j,1) - rh(j) * Bfd(3,2,j) * dr
    enddo
    do j=1, LJ
       do i=2, LI-1
          frp(i+1,j,1) = frp(i,j,1) + rf(j) * Bfd(1,i,j) * dz
       enddo
    enddo
    frp(1,:,1) = frp(2,:,1)
    
    !  --- [4] Calculate Br & Bz  ---             !
    !   -- [4-1] gradPsi x t -- !
    do j=1, LJ-1
       do i=1, LI-1
          Bfd(1,i,j) = + ( frp(i+1,j,1) - frp(i,j,1) ) * rfinv(j) * dzinv
          Bfd(3,i,j) = - ( frp(i,j+1,1) - frp(i,j,1) ) * rhinv(j) * drinv
       enddo
    enddo
    !   -- [4-2] Boundary Conditions -- !
    do i=1, LI-1
       Bfd(1,i,LJ) = + ( frp(i+1,LJ,1) - frp(i,LJ,1) ) * rfinv(LJ) * dzinv
    enddo
    do j=1, LJ-1
       Bfd(3,LI,j) = - ( frp(LI,j+1,1) - frp(LI,j,1) ) * rhinv( j) * drinv
    enddo
    Bfd(1, 1,:) = Bfd(1,   2,:)
    Bfd(1,LI,:) = Bfd(1,LI-1,:)
    Bfd(2, 1,:) = Bfd(2,   2,:)
    Bfd(2,LI,:) = Bfd(2,LI-1,:)
    Bfd(3,:, 1) = Bfd(3,:,   2)
    Bfd(3,:,LJ) = Bfd(3,:,LJ-1)
    if ( rmin.eq.0.d0 ) then
       Bfd(1,i,1) = - Bfd(1,i,3)
       Bfd(1,i,2) = 0.d0
       Bfd(2,i,1) = + Bfd(2,i,2)
       Bfd(3,i,1) = + Bfd(3,i,2)
    endif

    !  --- [5] Set rho, prs  --- !    
    Ttot    = ( 1.d0 + TiTe ) * vthcv**2
    Ttotinv =   1.d0 / Ttot
    do j=1, LJ
       do i=1, LI
          frp(i,j,2) = frp(i,j,3) * Ttotinv + rhofloor
          frp(i,j,3) = frp(i,j,2) * Ttot
       enddo
    enddo

    !  --- [6] Set Current   --- !
    valfe2 = valfe**2
    do j=2, LJ
       do i=2, LI
          Jcr(1,i,j) = - valfe2 *   ( Bfd(2,i,j) - Bfd(2,i-1,j) ) * dzinv
          Jcr(2,i,j) =   valfe2 * ( ( Bfd(1,i,j) - Bfd(1,i-1,j) ) * dzinv &
               &                  - ( Bfd(3,i,j) - Bfd(3,i,j-1) ) * drinv )
          Jcr(3,i,j) =   valfe2 *   ( rh(j)*Bfd(2,i,j) - rh(j-1)*Bfd(2,i,j-1) ) * drinv * rfinv(j)
       enddo
    enddo
    Jcr(1,:,1) = 0.d0
    Jcr(2,:,1) = 0.d0
    Jcr(3,:,1) = 0.d0
    Jcr(1,1,:) = 0.d0
    Jcr(2,1,:) = 0.d0
    Jcr(3,1,:) = 0.d0

    !  --- [7] Get Electron Flow Velocity ---  !
    do j=1, LJ
       do i=1, LI
          uvc(i,j,1) = - Jcr(1,i,j) / frp(i,j,2)
          uvc(i,j,2) = - Jcr(2,i,j) / frp(i,j,2)
          uvc(i,j,3) = - Jcr(3,i,j) / frp(i,j,2)
       enddo
    enddo

    return
  end subroutine Set2GS_Egrid


  subroutine Set2GS_Bgrid

    use constants, only : Nmid, rmin, rhofloor, TiTe, vthcv, valfe, CHSM, CaseIO
    use variables, only : dr, dz, drinv, dzinv
    use Field    , only : LIg, LJg, BgB, fgB
    implicit none
    integer            :: i, j, ip2, ip2m1
    double precision   :: dBr, dBt, dpr, Nmidinv, valfe2
    double precision   :: Ttot, Ttotinv, polarity(2)
    
    !  --- [0] Preparation   ---  !
    !    -- [0-1] grid --   !
    ip2         = LIg + Nmid
    ip2m1       = LIg + Nmid - 1
    Nmidinv     = 1.d0 / dble( Nmid )
    do j=1, LJ
       rf(j)    = dr * dble(j-2) + rmin
       rfinv(j) = 1.d0 / rf(j)
    enddo
    do j=1, LJ-1
       rh(j)    = 0.5d0 * ( rf(j) + rf(j+1) )
       rhinv(j) = 1.d0 / rh(j)
    enddo
    if ( rmin.eq.0.d0 ) then
       rf(2)    = 0.d0
       rfinv(2) = 0.d0
    end if
    rh(LJ)      = rf(LJ) + 0.5d0 * dr
    rhinv(LJ)   = 1.d0 / rh(LJ)

    !    -- [0-2] Initialize --   !
    Bfd(:,:,:)  = 0.d0
    frp(:,:,:)  = 0.d0

    !    -- [0-3] Polarity   --   !
    if (   CHSM.eq.'co-' ) then
       polarity(1) = + 1.d0
       polarity(2) = + 1.d0
    endif
    if ( ( CHSM.eq.'ctr' ).and.( CaseIO.eq.'CaseI' ) ) then
       polarity(1) = + 1.d0
       polarity(2) = - 1.d0
    endif
    if ( ( CHSM.eq.'ctr' ).and.( CaseIO.eq.'CaseO' ) ) then
       polarity(1) = - 1.d0
       polarity(2) = + 1.d0
    endif
    
    !  --- [1] Set Equilibrium plasma  ---  !
    do j=1, LJg
       do i=1, LIg
          Bfd(1,      i,j) =  BgB(1,i,j)
          Bfd(1,ip2m1+i,j) =  BgB(1,i,j)
          Bfd(2,      i,j) =  BgB(2,i,j) * polarity(1)
          Bfd(2,ip2m1+i,j) =  BgB(2,i,j) * polarity(2)
          Bfd(3,      i,j) =  BgB(3,i,j)
          Bfd(3,ip2m1+i,j) =  BgB(3,i,j)
       enddo
    enddo
    do j=1, LJg
       do i=1, LIg
          frp(      i,j,3) = fgB(i,j,3)
          frp(ip2m1+i,j,3) = fgB(i,j,3)
       enddo
    enddo

    !  --- [2] Linear Interpolation :: Br  ---    !
    do j=1, LJg
       dBr = ( Bfd(1,ip2+1,j) - Bfd(1,LIg,j) ) * Nmidinv
       dBt = ( Bfd(2,ip2+1,j) - Bfd(2,LIg,j) ) * Nmidinv
       dpr = ( frp(ip2+1,j,3) - frp(LIg,j,3) ) * Nmidinv       
       do i=LIg+1, ip2
          Bfd(1,i,j) = Bfd(1,LIg,j) + dBr * dble(i-LIg)
          Bfd(2,i,j) = Bfd(2,LIg,j) + dBt * dble(i-LIg)
          frp(i,j,3) = frp(LIg,j,3) + dBt * dble(i-LIg)
       enddo
    enddo
    
    !  --- [3] psi Calculation from z = zmin ---  !
    do j=2, LJ-1
       frp(2,j+1,1)    = frp(2,j,1) - ( rf(j+1)*Bfd(3,2,j+1) ) * dr
    enddo ! --  i = 2     -- !
    do j=2, LJ
       do i=2, LI-1
          frp(i+1,j,1) = frp(i,j,1) + rh(j)*( Bfd(1,i+1,j) )*dz
       enddo
    enddo ! --  i = 3-LI  -- !
    do j=2, LJ
       frp(1,j,1) = frp(2,j,1) - rh(j)*Bfd(1,2,j)*dz
    enddo ! --  i = 1     -- !
    frp(:,1,1) = frp(:,2,1)    
    
    !  --- [4] Calculate Br & Bz  ---  !
    do j=2, LJ-1
       do i=2, LI-1
          Bfd(1,i,j) = + ( frp(i,j,1) - frp(i-1,j,1) ) * rhinv( j) * dzinv
          Bfd(3,i,j) = - ( frp(i,j,1) - frp(i,j-1,1) ) * rfinv( j) * drinv
       enddo
    enddo
    do j=2, LJ
       Bfd(1,LI,j) = + ( frp(LI,j,1) - frp(LI-1,j,1) ) * rhinv( j) * dzinv
    enddo
    do i=2, LI
       Bfd(3,i,LJ) = - ( frp(i,LJ,1) - frp(i,LJ-1,1) ) * rfinv(LJ) * drinv
    enddo
    !  -- b.c. --  !
    Bfd(1, 1,:) = + Bfd(1,   2,:)
    Bfd(2, 1,:) = + Bfd(2,   2,:)
    Bfd(3,:, 1) = + Bfd(3,   :,3)
    Bfd(3,:, 2) = + Bfd(3,   :,3)
    if ( rmin.eq.0.d0 ) then
       Bfd(1,i,1) = - Bfd(1,i,3)
       Bfd(1,i,2) = 0.d0
       Bfd(2,i,1) = + Bfd(2,i,2)
       Bfd(3,i,1) = + Bfd(3,i,2)
    endif
    
    !  --- [5] Set rho, prs  --- !
    Ttot    = ( 1.d0 + TiTe ) * vthcv**2
    Ttotinv =   1.d0 / Ttot
    do j=1, LJ
       do i=1, LI
          frp(i,j,2) = frp(i,j,3) * Ttotinv + rhofloor
          frp(i,j,3) = frp(i,j,2) * Ttot
       enddo
    enddo

    !  --- [6] Set Current   --- !
    valfe2 = valfe**2
    do j=2, LJ-1
       do i=2, LI-1
          Jcr(1,i,j) = - valfe2 *   ( Bfd(2,i+1,j) - Bfd(2,i,j) ) * dzinv
          Jcr(2,i,j) =   valfe2 * ( ( Bfd(1,i+1,j) - Bfd(1,i,j) ) * dzinv &
               &                  - ( Bfd(3,i,j+1) - Bfd(3,i,j) ) * drinv )
          Jcr(3,i,j) =   valfe2 * ( rf(j+1)*Bfd(2,i,j+1) - rf(j)*Bfd(2,i,j) )*drinv*rhinv(j)
       enddo
    enddo
    !  -- b.c. --  !
    Jcr(1, 1, :) = 0.d0
    Jcr(2, 1, :) = 0.d0
    Jcr(3, 1, :) = 0.d0
    Jcr(2, 2, :) = 0.d0
    Jcr(3, 2, :) = 0.d0
    Jcr(1, :, 1) = 0.d0
    Jcr(2, :, 1) = 0.d0
    Jcr(3, :, 1) = 0.d0
    Jcr(1, :, 2) = 0.d0
    Jcr(2, :, 2) = 0.d0
    Jcr(1,LI, :) = 0.d0
    Jcr(2,LI, :) = 0.d0
    Jcr(3,LI, :) = 0.d0
    Jcr(1, :,LJ) = 0.d0
    Jcr(2, :,LJ) = 0.d0
    Jcr(3, :,LJ) = 0.d0    

    !  --- [7] Get Electron Flow Velocity ---  !
    do j=1, LJ
       do i=1, LI
          uvc(i,j,1) = - Jcr(1,i,j) / frp(i,j,2)
          uvc(i,j,2) = - Jcr(2,i,j) / frp(i,j,2)
          uvc(i,j,3) = - Jcr(3,i,j) / frp(i,j,2)
       enddo
    enddo
    
    return
  end subroutine Set2GS_Bgrid


  subroutine BgridConversion

    use constants, only : valfe
    use variables, only : dr, drinv, dzinv
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2, valfe2
    double precision   :: psiext(LI+1,LJ+1)
    double precision   :: Btfext(0:LI,0:LJ)

    !  --- [1]  Copy and reflect Psi  ---  !
    !  --- [1-1] Main Region  ---  !
    do j=1, LJ
       do i=1, LI
          psiext(i,j) = frp(i,j,1)
          Btfext(i,j) = Bfd(2,i,j)
       enddo
    enddo
    
    !  --- [1-2] j=1, LJ+1 boundary copy ---  !
    factor1 =   rh( 1)        * rhinv( 2)
    factor2 = ( rh(LJ) + dr ) * rhinv(LJ)
    ! ! factor3 =   rh( 2)        * rhinv( 1)
    ! factor4 =   rh(LJ) / (  rh(LJ) + dr )
    do i=2, LI
       psiext(i,   1) = ( 1.d0+factor1 )*psiext(i, 2) - factor1*psiext(i,   3)
       psiext(i,LJ+1) = ( 1.d0+factor2 )*psiext(i,LJ) - factor2*psiext(i,LJ-1)
       ! Btfext(i,   1) =        factor3  *Btfext(i, 2)
       ! Btfext(i,LJ+1) =        factor4  *Btfext(i,LJ)
    enddo

    !  --- [1-3] i=1, LI+1 boundary copy ---  !
    do j=1, LJ+1
       psiext(   1,j) = 2.d0 * psiext( 2,j) - psiext(   3,j)
       psiext(LI+1,j) = 2.d0 * psiext(LI,j) - psiext(LI-1,j)
       ! ! Btfext(   1,j) =        Btfext( 2,j)
       ! Btfext(LI+1,j) =        Btfext(LI,j)
    enddo

    !  --- [2]   Interpolation of psi    ---  !
    do j=1, LJ
       do i=1, LI
          frp(i,j,1) = 0.25d0 * ( psiext(i,j  ) + psiext(i+1,j  ) &
               &                + psiext(i,j+1) + psiext(i+1,j+1) )
       enddo
    enddo
    
    !  --- [3]    Calculate Br & Bz      ---  !
    Bfd(:,:,:) = 0.d0
    do j=1, LJ
       do i=2, LI
          Bfd(1,i,j) = + ( frp(i,j,1) - frp(i-1,j,1) ) * rhinv(j) * dzinv
       enddo
    enddo
    do j=2, LJ
       do i=2, LI
          Bfd(2,i,j) = 0.25d0 * ( rh(j)   * ( Btfext(i,j  ) + Btfext(i-1,j  ) ) &
               &                + rh(j-1) * ( Btfext(i,j-1) + Btfext(i-1,j-1) ) ) * rfinv(j)
       enddo
    enddo
    do j=2, LJ
       do i=1, LI
          Bfd(3,i,j) = - ( frp(i,j,1) - frp(i,j-1,1) ) * rfinv(j) * drinv
       enddo
    enddo
    Bfd(1,1,:) = Bfd(1,2,:)
    Bfd(2,1,:) = Bfd(2,2,:)
    Bfd(2,:,1) = Bfd(2,:,2)
    Bfd(3,:,1) = Bfd(3,:,2)
    
    !  --- [4]   Calculate Current Density ---  !
    valfe2 = valfe**2
    do j=1, LJ-1
       do i=1, LI-1
          Jcr(1,i,j) = - valfe2 *   ( Bfd(2,i+1,j) - Bfd(2,i,j) ) * dzinv
          Jcr(2,i,j) =   valfe2 * ( ( Bfd(1,i+1,j) - Bfd(1,i,j) ) * dzinv &
               &       - ( Bfd(3,i,j+1) - Bfd(3,i,j) ) * drinv )
          Jcr(3,i,j) =   valfe2 *   ( rf(j+1)*Bfd(2,i,j+1) - rf(j)*Bfd(2,i,j) ) * drinv * rhinv(j)
       enddo
    enddo
    Jcr(1,:,LJ) = 0.d0
    Jcr(2,:,LJ) = 0.d0
    Jcr(3,:,LJ) = 0.d0
    Jcr(1,LI,:) = 0.d0
    Jcr(2,LI,:) = 0.d0
    Jcr(3,LI,:) = 0.d0
    
    return
  end subroutine BgridConversion

end module equset
