module BgridMod
  use constants, only : Nz, Nr
  implicit none
  integer, parameter :: LIg = Nz+1
  integer, parameter :: LJg = Nr+1
  integer, parameter :: br_ = 1, bt_ = 2, bz_ = 3
  integer, parameter :: jr_ = 1, jt_ = 2, jz_ = 3
  integer, parameter :: at_ = 1, rh_ = 2, pr_ = 3
  double precision   :: rf(LJg+1), rh(LJg+1), rfInv(LJg+1), rhInv(LJg+1)
  double precision   :: BgB(3,LIg,LJg), JgB(3,LIg,LJg), fgB(3,LIg,LJg), ugB(3,LIg,LJg)
contains
  
  subroutine BgridField( coordinate, onGrid )
    implicit none
    character(1), intent(in) :: onGrid
    character(2), intent(in) :: coordinate
    
    if ( ( coordinate.eq.'RZ' ).and.( onGrid.eq.'B' ) ) call Cylindric__BgridField
    return
  end subroutine BgridField
  
  
  subroutine Cylindric__BgridField
    use constants, only : Bmax, valfe, TiTe, vthcv, rhofloor, rMin, normSW, normType
    use variables, only : dr, dz, psi, gfunc, prs
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2, valfe2, drInv, dzInv
    double precision   :: absB, absBmax, prscoef, Bfdcoef, Ttot, TInv, rhof, pMax, coef
    double precision   :: psiExt(LIg+1,LJg+1), gfcExt(LIg+1,LJg+1)
    double precision   :: divB, sum_divB, sum_divB2, divBave, divBMSE

    ! ------------------------------------ !
    ! --- [1]  Preparation             --- !
    ! ------------------------------------ !
    !  -- [1-1] Grid Making  -- !
    drInv    = 1.d0 / dr
    dzInv    = 1.d0 / dz
    do j=1, LJg+1
       rf(j) = dr*dble(j-2) + rMin
       rh(j) = dr*dble(j-2) + rMin + 0.5d0*dr
    enddo
    rfInv(:) = 0.d0
    rhInv(:) = 0.d0
    do j=1, LJg+1
       if ( rf(j).ne.0.d0 ) rfInv(j) = 1.d0 / rf(j)
       if ( rh(j).ne.0.d0 ) rhInv(j) = 1.d0 / rh(j)
    enddo
    
    ! ------------------------------------ !
    !  --- [2]  Copy and reflect Psi  ---  !
    ! ------------------------------------ !
    !  --- [2-1] Main Region  ---  !
    do j=2, LJg          ! j=2, LJg   for j+1
       do i=2, LIg       ! i=2, LIg   for i+1
          psiExt(i,j)  = - psi  (i-1,j-1)
          gfcExt(i,j)  = + gfunc(i-1,j-1)
       enddo
    enddo
    !  --- [2-2] j=1, LJ+1 Boundary copy :: ( psiExt ) ---  !
    factor1 = rh(  2)*rhInv(    3)
    factor2 = rh(LJg)*rhInv(LJg-1)
    do i=2, LIg
       psiExt(i,    1) = psiExt(i,  2) + factor1 * ( psiExt(i,2  ) - psiExt(i,3    ) )
       psiExt(i,LJg+1) = psiExt(i,LJg) + factor2 * ( psiExt(i,LJg) - psiExt(i,LJg-1) )
    enddo
    !  --- [2-3] i=1, LI+1 Boundary copy :: ( psiExt ) ---  !
    do j=1, LJg+1
       psiExt(    1,j) = psiExt(  2,j) + ( psiExt(2  ,j) - psiExt(    3,j) )
       psiExt(LIg+1,j) = psiExt(LIg,j) + ( psiExt(LIg,j) - psiExt(LIg-1,j) )
    enddo
    !  --- [2-4] i=1, j=1  Boundary copy :: ( gfcExt ) ---  !
    gfcExt(1,:)        = gfcExt(  2,:) + ( gfcExt(  2,:) - gfcExt(    3,:) )
    gfcExt(:,1)        = gfcExt(  :,2) + ( gfcExt(  :,2) - gfcExt(    :,3) )
    
    ! ------------------------------------ !
    !  --- [3]  Interpolation of psi  ---  !
    ! ------------------------------------ !
    do j=1, LJg
       do i=1, LIg
          fgB(at_,i,j) = 0.25d0 * ( + psiExt(i,j+1) + psiExt(i+1,j+1) &
               &                    + psiExt(i,j  ) + psiExt(i+1,j  ) )
       enddo
    enddo

    ! ------------------------------------ !
    !  --- [4]   Calculate  Br & Bz   ---  !
    ! ------------------------------------ !
    BgB(:,:,:) = 0.d0
    do j=2, LJg
       do i=2, LIg
          BgB(br_,i,j) = - ( fgB(at_,i,j) - fgB(at_,i-1,j) ) * rhInv(j) * dzInv
          BgB(bt_,i,j) = +   gfcExt(i,j)                     * rfInv(j)
          BgB(bz_,i,j) = + ( fgB(at_,i,j) - fgB(at_,i,j-1) ) * rfInv(j) * drInv
       enddo
    enddo ! ~~ [4-1] Br, Bt, Bz ( 2~LIg, 2~LJg ) ~~ !
    !   -- [4-2] Boundary Condition B -- !
    if ( rMin.eq.0.d0 ) then
       ! - i=1 : dB// / dn = 0 - !
       BgB(br_,  1,  :) = + BgB(br_,    2,    :)
       BgB(bt_,  1,  :) = + BgB(bt_,    2,    :)
       BgB(bz_,  1,  :) = + BgB(bz_,    2,    :)
       ! - j=2 : pole           - !
       BgB(bt_,  :,  2) =   0.d0
       BgB(bz_,  :,  2) = + BgB(bz_,    :,    3)
       ! - j=1 : anti-symmetric - !
       BgB(br_,  :,  1) = - BgB(br_,    :,    2)
       BgB(bt_,  :,  1) = - BgB(bt_,    :,    3)
       BgB(bz_,  :,  1) = + BgB(bz_,    :,    2)
    else
       ! - i=1 : dB// / dn = 0 - !
       BgB(br_,  1,  :) = + BgB(br_,    2,    :)
       BgB(bt_,  1,  :) = + BgB(bt_,    2,    :)
       BgB(bz_,  1,  :) = + BgB(bz_,    2,    :)
       ! - j=1 : Vacuum Field  - !
       BgB(br_,  :,  1) = + BgB(br_,    :,    2) * rf(2) * rfInv(1)
       BgB(bt_,  :,  1) = + BgB(bt_,    :,    2) * rf(2) * rfInv(1)
       BgB(bz_,  :,  1) = + BgB(bz_,    :,    2) * rf(2) * rfInv(1)
    endif

    ! ! -------------- !
    ! sum_divB  = 0.d0
    ! sum_divB2 = 0.d0
    ! do j=3, LJg-1
    !    do i=3, LIg-1
    !       divB      =   rfInv(j) * drInv * ( rh(j)*BgB(br_,i,j) - rh(j-1)*BgB(br_,i  ,j-1) ) &
    !            &      +            dzInv * (       BgB(bz_,i,j) -         BgB(bz_,i-1,j  ) )
    !       sum_divB  = sum_divB  + divB
    !       sum_divB2 = sum_divB2 + divB**2
    !    enddo
    ! enddo
    ! write(6,*) "divB Ave  :: ", sum_divB  / dble( (LIg-3)*(LJg-3) )
    ! write(6,*) "divB MSE  :: ", sum_divB2 / dble( (LIg-3)*(LJg-3) )
    ! stop
    ! ! -------------- !
    

    !   -- [5-1] Copy FRP -- !
    do j=2, LJg
       do i=2, LIg
          fgB(pr_,i,j) = prs(i-1,j-1)
          fgB(at_,i,j) = fgB(at_,i,j)*rhInv(j)
       enddo
    enddo
    fgB(pr_,  1,  :) = fgB(pr_,    2,    :)
    fgB(pr_,  :,  1) = fgB(pr_,    :,    3)
    fgB(at_,  1,  :) = fgB(at_,    2,    :)
    fgB(at_,  :,  1) = fgB(at_,    :,    2)
    
    !  --- [6] Adjust to Bmax  ---  !
    if ( normSW ) then
       absBmax = 0.d0
       do j=2, LJg
          do i=2, LIg
             absBmax   = max( absBmax, BgB(br_,i,j)**2 + BgB(bt_,i,j)**2 + BgB(bz_,i,j)**2 )
          enddo
       enddo
       absBmax =    sqrt( absBmax )
       Bfdcoef =   Bmax / absBmax
       prscoef = ( Bmax / absBmax )**2
       write(6,*) ' absBmax == ', absBmax
       write(6,*) ' Bfdcoef == ', Bfdcoef
       write(6,*) ' prscoef == ', prscoef
       do j=1, LJg
          do i=1, LIg
             BgB(br_,i,j) = Bfdcoef * BgB(br_,i,j)
             BgB(bt_,i,j) = Bfdcoef * BgB(bt_,i,j)
             BgB(bz_,i,j) = Bfdcoef * BgB(bz_,i,j)
             fgB(at_,i,j) = Bfdcoef * fgB(at_,i,j)
             fgB(pr_,i,j) = prscoef * fgB(pr_,i,j)
          enddo
       enddo
    endif
    
    !  --- [6] pressure & rho  ---  !
    !   -- [6-2] Limit rho to rhoFloor -- !
    pMax    = maxval( fgB(pr_,:,:) )
    if ( normType.eq.'PIC' ) then
       Ttot = ( 1.d0 + TiTe ) * vthcv**2 ! * valfe**2
       TInv = 1.d0 / Ttot
    endif
    if ( normType.eq.'MHD' ) then
       Ttot = pMax
       TInv = 1.d0 / Ttot
    endif
    rhof    =  pMax * TInv * rhofloor
    do j=1, LJg
       do i=1, LIg
          fgB(rh_,i,j) = fgB(pr_,i,j) * TInv + rhof
          fgB(pr_,i,j) = fgB(rh_,i,j) * Ttot
       enddo
    enddo

    !  --- [7]   Calculate Current Density ---  !
    !   -- [7-1] Current = rot B -- !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    do j=2, LJg-1
       do i=2, LIg-1
          JgB(jr_,i,j) = coef * ( - (         BgB(bt_,i+1,j  ) -       BgB(bt_,i,j) ) * dzInv )
          JgB(jt_,i,j) = coef * ( + (         BgB(br_,i+1,j  ) -       BgB(br_,i,j) ) * dzInv &
               &                  - (         BgB(bz_,i  ,j+1) -       BgB(bz_,i,j) ) * drInv )
          JgB(jz_,i,j) = coef * ( + ( rf(j+1)*BgB(bt_,i  ,j+1) - rf(j)*BgB(bt_,i,j) ) * drInv * rhInv(j) )
       enddo
    enddo
    ! -- B.C. -- !
    JgB(  :,  1,  :) = JgB(:  ,2,:)
    JgB(  :,  :,  1) = JgB(:  ,:,2)
    JgB(jz_,  2,  :) = 0.d0
    JgB(jt_,  2,  :) = 0.d0
    JgB(  :,LIg,  :) = JgB(:,LIg-1,:)
    JgB(  :,  :,LJg) = JgB(:,:,LJg-1)
    if ( rMin.eq.0.d0 ) then
       JgB(jr_,:, 2) = 0.d0
       JgB(jt_,:, 2) = 0.d0
    endif
    ! JgB(  :,  1,  :) = 0.d0
    ! JgB(  :,  :,  1) = 0.d0
    ! JgB(jz_,  2,  :) = 0.d0
    ! JgB(jt_,  2,  :) = 0.d0
    ! JgB(  :,LIg,  :) = 0.d0
    ! JgB(  :,  :,LJg) = 0.d0
    ! if ( rMin.eq.0.d0 ) then
    !    JgB(jr_,:, 2) = 0.d0
    !    JgB(jt_,:, 2) = 0.d0
    ! endif
   
    !  --- [8] Get Electron Flow Velocity ---  !
    do j=1, LJg
       do i=1, LIg
          ugB(jr_,i,j)   = - JgB(jr_,i,j) / fgB(rh_,i,j)
          ugB(jt_,i,j)   = - JgB(jt_,i,j) / fgB(rh_,i,j)
          ugB(jz_,i,j)   = - JgB(jz_,i,j) / fgB(rh_,i,j)
       enddo
    enddo
    
    return
  end subroutine Cylindric__BgridField
  
end module BgridMod
