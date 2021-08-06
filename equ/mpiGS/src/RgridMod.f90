module RgridMod
  use constants, only : Nz, Nr
  implicit none
  integer, parameter :: LIr = Nz+1
  integer, parameter :: LJr = Nr+1
  integer, parameter :: br_ = 1, bt_ = 2, bz_ = 3
  integer, parameter :: jr_ = 1, jt_ = 2, jz_ = 3
  integer, parameter :: at_ = 1, rh_ = 2, pr_ = 3
  double precision   :: rf(LJr+1), rh(LJr+1), rfinv(LJr+1), rhinv(LJr+1)
  double precision   :: BgR(3,LIr,LJr), JgR(3,LIr,LJr), fgR(LIr,LJr,3), ugR(LIr,LJr,3)
contains


  subroutine RgridField( coordinate, onGrid )
    implicit none
    character(1), intent(in) :: onGrid
    character(2), intent(in) :: coordinate
    
    if ( ( coordinate.eq.'RZ' ).and.( onGrid.eq.'R' ) ) call Cylindric__RegularGrid
    return
  end subroutine RgridField

  subroutine Cylindric__RegularGrid

    use constants, only : Bmax, valfe, TiTe, vthcv, rhofloor, rMin, normSW, normType
    use variables, only : dr, dz, psi, gfunc, prs
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2, valfe2, drinv, dzinv
    double precision   :: absB, absBmax, prscoef, Bfdcoef, Ttot, Tinv, rhof, pMax, coef
    double precision   :: psiExt(LIr+1,LJr+1), gfcExt(LIr+1,LJr+1)

    ! --- [1] Preparation   --- !
    !  -- [1-1] Grid Making  -- !
    drinv       = 1.d0 / ( 2.d0*dr )
    dzinv       = 1.d0 / ( 2.d0*dz )
    do j=1, LJr+1
       rf(j) = dr*dble(j-2) + rMin
       rh(j) = dr*dble(j-2) + rMin + 0.5d0*dr
    enddo
    rfinv(:) = 0.d0
    rhinv(:) = 0.d0
    do j=1, LJr+1
       if ( rf(j).ne.0.d0 ) rfinv(j) = 1.d0 / rf(j)
       if ( rh(j).ne.0.d0 ) rhinv(j) = 1.d0 / rh(j)
    enddo
    
    !  --- [2]  Copy and reflect Psi  ---  !
    !  --- [2-1] Main Region  ---  !
    do j=2, LJr          ! j=2, LJr   for j+1
       do i=2, LIr       ! i=2, LIr   for i+1
          psiExt(i,j) = psi  (i-1,j-1)
          gfcExt(i,j) = gfunc(i-1,j-1)
       enddo
    enddo
    !  --- [2-2] j=1, LJ+1 boundary copy :: ( psiExt ) ---  !
    factor1 = rh(  2) * rhInv(    3)
    factor2 = rh(LJr) * rhInv(LJr-1)
    do i=2, LIr
       psiExt(i,    1) = psiExt(i,  2) + factor1 * ( psiExt(i,2  ) - psiExt(i,3    ) )
       psiExt(i,LJr+1) = psiExt(i,LJr) + factor2 * ( psiExt(i,LJr) - psiExt(i,LJr-1) )
    enddo
    !  --- [2-3] i=1, LI+1 boundary copy :: ( psiExt ) ---  !
    do j=1, LJr+1
       psiExt(    1,j) = psiExt(  2,j) + ( psiExt(2  ,j) - psiExt(    3,j) )
       psiExt(LIr+1,j) = psiExt(LIr,j) + ( psiExt(LIr,j) - psiExt(LIr-1,j) )
    enddo
    !  --- [2-4] i=1, j=1 boundary copy :: ( gfcExt ) ---  !
    gfcExt(1,:)        = gfcExt(  2,:) + ( gfcExt(  2,:) - gfcExt(    3,:) )
    gfcExt(:,1)        = gfcExt(  :,2) + ( gfcExt(  :,2) - gfcExt(    :,3) )
    
    !  --- [3]   Interpolation of psi    ---  !
    do j=1, LJr
       do i=1, LIr
          fgR(i,j,at_) = psiExt(i,j) * rfInv(j)
       enddo
    enddo

    !  --- [4]    Calculate Br & Bz      ---  !
    BgR(:,:,:) = 0.d0
    do j=2, LJr
       do i=2, LIr
          BgR(br_,i,j) = - ( psiExt(i+1,j) - psiExt(i-1,j) ) * rfinv(j) * dzinv
          BgR(bt_,i,j) = +   gfcExt(i  ,j)                   * rfinv(j)
          BgR(bz_,i,j) = + ( psiExt(i,j+1) - psiExt(i,j-1) ) * rfinv(j) * drinv
       enddo
    enddo ! ~~ [4-1] Br, Bt, Bz ( 2~LIr, 2~LJr ) ~~ !
    !   -- [4-2] Boundary Condition B -- !
    if ( rMin.eq.0.d0 ) then
       ! - i=1 : dB// / dn = 0 - !
       BgR(br_,  1,  :) = + BgR(br_,    3,    :)
       BgR(bt_,  1,  :) = + BgR(bt_,    3,    :)
       BgR(bz_,  1,  :) = + BgR(bz_,    2,    :)
       ! - j=1 : anti-symmetric - !
       BgR(br_,  :,  1) = - BgR(br_,    :,    3)
       BgR(bt_,  :,  1) = - BgR(bt_,    :,    3)
       BgR(bz_,  :,  1) = + BgR(bz_,    :,    3)
       ! - j=2 : pole           - !
       BgR(br_,  :,  2) = 0.d0
       BgR(bt_,  :,  2) = 0.d0
       BgR(bz_,  :,  2) = BgR(bz_,:,3)
    else
       ! - i=1 : dB// / dn = 0 - !
       BgR(br_,  1,  :) = + BgR(br_,    3,    :)
       BgR(bt_,  1,  :) = + BgR(bt_,    3,    :)
       BgR(bz_,  1,  :) = + BgR(bz_,    2,    :)
       ! - j=1 : Vacuum Field  - !
       BgR(br_,  :,  1) = + BgR(br_,    :,    2) * rf(2) * rfInv(1)
       BgR(bt_,  :,  1) = + BgR(bt_,    :,    2) * rf(2) * rfInv(1)
       BgR(bz_,  :,  1) = + BgR(bz_,    :,    2) * rf(2) * rfInv(1)
    endif

    !   -- [5-1] Copy FRP -- !
    do j=2, LJr
       do i=2, LIr
          fgR(i,j,pr_) = prs(i-1,j-1)
       enddo
    enddo
    fgR(  1,  :,pr_) = fgR(    3,    :,pr_)
    fgR(  :,  1,pr_) = fgR(    :,    3,pr_)
    
    !  --- [6] Adjust to Bmax  ---  !
    if ( normSW ) then
       absBmax = 0.d0
       do j=2, LJr
          do i=2, LIr
             absB      = sqrt( BgR(br_,i,j)**2 + BgR(bt_,i,j)**2 + BgR(bz_,i,j)**2 )
             absBmax   = max( absBmax, absB )
          enddo
       enddo
       Bfdcoef =   Bmax / absBmax
       prscoef = ( Bmax / absBmax )**2
       write(6,*) ' absBmax == ', absBmax
       write(6,*) ' Bfdcoef == ', Bfdcoef
       write(6,*) ' prscoef == ', prscoef
       do j=1, LJr
          do i=1, LIr
             BgR(br_,i,j) = Bfdcoef * BgR(br_,i,j)
             BgR(bt_,i,j) = Bfdcoef * BgR(bt_,i,j)
             BgR(bz_,i,j) = Bfdcoef * BgR(bz_,i,j)
             fgR(i,j,at_) = Bfdcoef * fgR(i,j,at_)
             fgR(i,j,pr_) = prscoef * fgR(i,j,pr_)
          enddo
       enddo
    endif
    
    !  --- [6] pressure & rho  ---  !
    !   -- [6-2] Limit rho to rhoFloor -- !
    pMax    = maxval( fgR(:,:,pr_) )
    if ( normType.eq.'PIC' ) then
       Ttot = ( 1.d0 + TiTe ) * vthcv**2 * valfe**2
       Tinv = 1.d0 / Ttot
    endif
    if ( normType.eq.'MHD' ) then
       Ttot = pMax
       Tinv = 1.d0 / (            Ttot )
    endif
    rhof    =  pMax * Tinv * rhofloor
    do j=1, LJr
       do i=1, LIr
          fgR(i,j,rh_) = fgR(i,j,pr_) * Tinv + rhof
          fgR(i,j,pr_) = fgR(i,j,rh_) * Ttot
       enddo
    enddo
    
    !  --- [7]   Calculate Current Density ---  !
    !   -- [7-1] Current = rot B -- !
    if ( normType.eq.'PIC' ) coef = valfe**2
    if ( normType.eq.'MHD' ) coef = 1.d0
    do j=2, LJr-1
       do i=2, LIr-1
          JgR(jr_,i,j)   = coef * ( - (         BgR(bt_,i+1,j  ) -         BgR(bt_,i-1,j  ) ) ) * dzinv
          JgR(jt_,i,j)   = coef * ( + (         BgR(br_,i+1,j  ) -         BgR(br_,i-1,j  ) )   * dzinv &
               &                    - (         BgR(bz_,i  ,j+1) -         BgR(bz_,i  ,j-1) )   * drinv )
          JgR(jz_,i,j)   = coef * ( + ( rf(j+1)*BgR(bt_,i  ,j+1) - rf(j-1)*BgR(bt_,i  ,j-1) ) ) * drinv * rfinv(j)
       enddo
    enddo
    !   -- [7-2] Boundary Condition -- !
    JgR(jr_,  :,  1)     = 0.d0
    JgR(jt_,  :,  1)     = 0.d0
    JgR(jz_,  :,  1)     = 0.d0
    JgR(jr_,  1,  :)     = 0.d0
    JgR(jt_,  1,  :)     = 0.d0
    JgR(jz_,  1,  :)     = 0.d0
    JgR(jr_,  :,LJr)     = 0.d0
    JgR(jt_,  :,LJr)     = 0.d0
    JgR(jz_,  :,LJr)     = 0.d0
    JgR(jr_,LIr,  :)     = 0.d0
    JgR(jt_,LIr,  :)     = 0.d0
    JgR(jz_,LIr,  :)     = 0.d0
    if ( rMin.eq.0.d0 ) then
       JgR(jr_,:,2)      = 0.d0
       JgR(jt_,:,2)      = 0.d0
       JgR(jz_,:,2)      = 2.d0 * coef * drinv * BgR(bt_,:,3)
    endif
   
    !  --- [8] Get Electron Flow Velocity ---  !
    do j=1, LJr
       do i=1, LIr
          ugR(i,j,jr_)   = - JgR(jr_,i,j) / fgR(i,j,rh_)
          ugR(i,j,jt_)   = - JgR(jt_,i,j) / fgR(i,j,rh_)
          ugR(i,j,jz_)   = - JgR(jz_,i,j) / fgR(i,j,rh_)
       enddo
    enddo
    
    return
  end subroutine Cylindric__RegularGrid
  
end module RgridMod
