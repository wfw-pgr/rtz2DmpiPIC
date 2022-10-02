module boundaryCMod
contains

  ! ====================================================== !
  ! === set__initial_boundaryCondition                 === !
  ! ====================================================== !
  subroutine set__initial_boundaryCondition
    use variablesMod
    implicit none
    integer                     :: i, j
    double precision, parameter :: coef = 0.5d0

    ! ------------------------------------------------------ !
    ! --- [1]  approximate Vetical Field (EF) value      --- !
    ! ------------------------------------------------------ !
    !  -- Shafranov's Eq. -- 
    if      ( trim( normalizationInGS ) == "MHD" ) then
       Bv0 =  - Itot / (            4.d0 * pi * magAxis(r_) ) &
            &          * ( log( 8.d0*magAxis(r_)/aMinor ) - 1.5d0 + beta_pol )
       
    else if ( trim( normalizationInGS ) == "PIC" ) then
       Bv0 =  - Itot / ( valfe**2 * 4.d0 * pi * magAxis(r_) ) &
            &          * ( log( 8.d0*magAxis(r_)/aMinor ) - 1.5d0 + beta_pol )
       !  Bv0 =    Itot / ( 4.d0 * pi * magAxis(r_) ) * ( log( 8.d0*magAxis(r_) / aMinor ) - 1.5d0 + beta_pol )
    endif

    ! ------------------------------------------------------ !
    ! --- [2] set perfect conducting wall boundary       --- !
    ! ------------------------------------------------------ !
    if ( trim( boundaryType ) == "pcw" ) then
       do j=1, LJ
          bdrcLR(j,lft_) = coef * rAxis( j)**2 * Bv0
          bdrcLR(j,rgt_) = coef * rAxis( j)**2 * Bv0
       enddo
       do i=1, LI
          bdrcTB(i,top_) = coef * rAxis(LJ)**2 * Bv0
          bdrcTB(i,bot_) = coef * rAxis( 1)**2 * Bv0
       enddo
    endif

    ! ------------------------------------------------------ !
    ! --- [3] set symmetry pcw boundary (symm.@left )    --- !
    ! ------------------------------------------------------ !
    if ( trim( boundaryType ) == "hlf" ) then
       do j=1, LJ
          bdrcLR(j,lft_) = 0.d0
          bdrcLR(j,rgt_) = coef * rAxis( j)**2 * Bv0
       enddo
       do i=1, LI
          bdrcTB(i,top_) = coef * rAxis(LJ)**2 * Bv0
          bdrcTB(i,bot_) = coef * rAxis( 1)**2 * Bv0
       enddo
    endif

    ! ------------------------------------------------------ !
    ! --- [4] set r.h.s. for BiCG solver                 --- !
    ! ------------------------------------------------------ !
    if ( trim( solverType ) == "BCG" ) then
       do j=1, LJ
          rhs( 1,j) = bdrcLR(j,lft_)
          rhs(LI,j) = bdrcLR(j,rgt_)
       enddo
       do i=1, LI
          rhs(i,LJ) = bdrcTB(i,top_)
          rhs(i, 1) = bdrcTB(i,bot_)
       enddo
    endif

    return
  end subroutine set__initial_boundaryCondition

  
end module boundaryCMod


