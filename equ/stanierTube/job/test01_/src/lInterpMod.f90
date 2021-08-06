module lInterpMod
contains
  
  ! =================================================================== !
  ! ===  LinearInterpFWD :: Linear Interpolation                    === !
  ! =================================================================== !
  subroutine LinearInterpFWD( y1, fxy1, y2, fxy2, LI, LJ, sftype )
    implicit none
    integer                       :: i, j
    double precision              :: y1Len, dy1, dy1Inv, dy2
    integer                       :: LJp
    integer         , intent(in)  :: LI, LJ
    double precision, intent(in)  :: y1(1:LJ), fxy1(LI,1:LJ)
    double precision, intent(out) :: y2(0:LJ), fxy2(LI,0:LJ)
    character(3)    , intent(in)  :: sftype
    
    ! ------------------------------------- !
    ! --- [1] New Grid Making           --- !
    ! ------------------------------------- !
    if ( y1(1).ne.0.d0 ) stop '[ERROR] x(1) != 0.0'
    !  - y1 Grid Info.   - !
    y1Len  = y1(LJ) - y1(1)
    dy1    = y1( 2) - y1(1)
    dy1Inv =  1.d0 / dy1
    LJp    = LJ + 1
    !  - Define y2 Grid  - !
    dy2    = y1Len / ( ( dble(LJ)-1.d0 ) + 0.5d0 )
    do j=0, LJ
       y2(j) = dy2 * ( dble(j) - 0.5d0 )
    enddo
    ! ------------------------------------- !
    ! --- [2] Interpolate Forward       --- !
    ! ------------------------------------- !
    do i=1, LI
       call SF1DInterpolate( y1(1:LJ), fxy1(i,1:LJ), y2(1:LJ), fxy2(i,1:LJ), LJ, LJ, sftype )
    enddo
    do i=1, LI
       fxy2(i,0) = fxy2(i,1)
    enddo
    
    return
  end subroutine LinearInterpFWD

  
  ! =================================================================== !
  ! ===  LinearInterpRET  :: back transform y2->y1                  === !
  ! =================================================================== !
  subroutine LinearInterpRET( x, fin, y, fout, N, M )
    implicit none
    integer                       :: i, j, jp1, jp2
    double precision              :: dx, dxInv, coef(2), ydx
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: x(0:M), y(1:M), fin( N,0:M)
    double precision, intent(out) :: fout(N,1:M)

    ! ------------------------------------- !
    ! --- [1] Grid Settings             --- !
    ! ------------------------------------- !
    dx    = x(2) - x(1)
    dxInv = 1.d0 / dx
    ! ------------------------------------- !
    ! --- [2] Back Transform            --- !
    ! ------------------------------------- !
    do j=1, M
       ydx             = ( y(M)-y(j) ) * dxInv
       jp1             = min( M - ceiling( ydx ), M-1 )
       jp2             = jp1 + 1

       coef(2)         = dble( M - jp1 ) - ydx
       coef(1)         = 1.d0 - coef(2)
       do i=1, N
          fout(i,j) = coef(1)*fin(i,jp1) + coef(2)*fin(i,jp2)
       enddo
    enddo
    return
  end subroutine LinearInterpRET


  ! =================================================================== !
  ! ===  shapef  :: Shape Function                                  === !
  ! =================================================================== !
  function shapef( iposit, rposit, drinv, sftype )
    implicit none
    integer         , intent(in) :: iposit
    double precision, intent(in) :: rposit
    double precision, intent(in) :: drinv
    character(3)    , intent(in) :: sftype
    double precision             :: shapef(5)
    double precision             :: delta

    select case ( sftype )
    case( "CIC" )
       ! ------------- CIC ------------ !
       delta     = rposit*drinv - dble( iposit )
       shapef(1) = 0.d0
       shapef(2) = 0.d0
       shapef(3) = 1.d0 - abs(delta)
       shapef(4) = 0.d0
       shapef(5) = 0.d0
       shapef( 3 + int( sign( 1.d0, delta ) ) ) = abs(delta)

    case( "2nd" )
       ! ----------- 2nd-b-spline ---------- !
       delta       = rposit*drinv - dble( iposit )
       shapef(1) = 0.d0
       shapef(2) = 0.50d0*( 0.5d0-delta )**2
       shapef(3) = 0.75d0 - delta*delta
       shapef(4) = 0.50d0*( 0.5d0+delta )**2
       shapef(5) = 0.d0

    end select

    return
  end function shapef


  ! =================================================================== !
  ! ===  SF1DInterpolate  ::  Interpolation using Shape Function    === !
  ! =================================================================== !
  subroutine SF1DInterpolate( x1, y1, x2, y2, LJ1, LJ2, sftype )
    implicit none
    integer         , intent(in   ) :: LJ1    , LJ2
    double precision, intent(in   ) :: x1(LJ1), y1(LJ1)
    double precision, intent(inout) :: x2(LJ2), y2(LJ2)
    character(3)    , intent(in)    :: sftype
    integer                         :: j, jg, jp, jpp
    double precision                :: x1Floor, dx1Inv, retg
    double precision                :: sf(-2:2)

    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    x1Floor = x1(1)
    dx1Inv  = 1.d0 / ( x1(2) - x1(1) )
    ! ------------------------------------- !
    ! --- [2] Interpolation             --- !
    ! ------------------------------------- !
    do jg=1, LJ2

       jp   = nint( ( x2(jg) - x1Floor ) * dx1Inv )
       sf   = shapef( jp, x2(jg), dx1Inv, sftype )
       jp   = jp + 1
       retg = 0.d0
       do j=-2, 2
          jpp  = min( LJ1, max( 1, jp+j ) )
          retg = retg + sf(j)*y1(jpp)
       enddo
       y2(jg) = retg

    enddo
    return
  end subroutine SF1DInterpolate

  
end module lInterpMod



! subroutine rct2rct( xp, yp, data, ret, xg, yg, LIp, LJp, LIg, LJg, sftype, ierr )
!   implicit none
!   integer         , intent(in)  :: LIp, LJp, LIg, LJg
!   double precision, intent(in)  :: xp(LIp), yp(LJp)
!   double precision, intent(in)  :: xg(LIg), yg(LJg)
!   double precision, intent(in)  :: data(LIp,LJp)
!   character(3)    , intent(in)  :: sftype
!   double precision, intent(out) ::  ret(LIg,LJg)
!   integer         , intent(out) :: ierr
!   integer                       :: i, j, ig, jg, ip, jp, ipp, jpp
!   double precision              :: dx, dy, dxInv, dyInv, xMin, yMin, xMax, yMax
!   double precision              :: xn(LIg), yn(LJg), sfx(-2:2), sfy(-2:2), sfxy, retg
!   logical                       :: chk

!   ! --- [0] Preparation --- !
!   ierr  = 0
!   dx    = xp(2) - xp(1)
!   dy    = yp(2) - yp(1)
!   dxInv = 1.d0 / dx
!   dyInv = 1.d0 / dy
!   xMin  = xp(1)
!   yMin  = yp(1)
!   ! -- ( xn, yn = xg, yg - xMin, yMin ) -- !
!   do ig=1, LIg
!      xn(ig) = xg(ig) - xMin
!   enddo
!   do jg=1, LJg
!      yn(jg) = yg(jg) - yMin
!   enddo
!   ! -- Boundary Check :: x -- !
!   xMin  = 0.d0
!   xMax  = xp(LIp) - xp(1)
!   yMin  = 0.d0
!   yMax  = yp(LJp) - yp(1)
!   do ig=1, LIg
!      if ( xn(ig).lt.xMin ) then
!         chk    = .true.
!         xn(ig) = xMin
!      endif
!      if ( xn(ig).gt.xMax ) then
!         chk    = .true.
!         xn(ig) = xMax
!      endif
!   enddo
!   ! -- Boundary Check :: y -- !
!   do jg=1, LJg
!      if ( yn(jg).lt.yMin ) then
!         chk   = .true.
!         yn(jg) = yMin
!      endif
!      if ( yn(jg).gt.yMax ) then
!         chk   = .true.
!         yn(jg) = yMax
!      endif
!   enddo
!   if ( chk ) ierr = 1

!   return
  
! contains
! end subroutine rct2rct


    ! !  - x & fin --> y, itrp - !
    ! do j=1, M
    !    y1dy1           = ( y(M)-y(j) ) * dxinv
    !    jp1             = min( max( M - ceiling( ydx ), 1 ), M-1 )
    !    jp2             = jp1 + 1

    !    coef(2)         = dble( M - jp1 ) - ydx
    !    coef(1)         = 1.d0 - coef(2)
    !    do i=1, N
    !       fout(i,j) = coef(1)*fin(i,jp1) + coef(2)*fin(i,jp2)
    !    enddo
    ! enddo

  
  ! subroutine LinearInterpRET( x, fin, y, fout, N, M )
  !   implicit none
  !   integer                       :: i, j, jp1, jp2
  !   double precision              :: dx, dxInv, coef(2), ydx
  !   integer         , intent(in)  :: N, M
  !   double precision, intent(in)  :: x(0:M), y(1:M), fin( N,0:M)
  !   double precision, intent(out) :: fout(N,1:M)

  !   ! -- dx, dxInv -- !
  !   dx    = x(2) - x(1)
  !   dxInv = 1.d0 / dx
  !   ! -- Return -- !
  !   do j=1, M
  !      ydx             = ( y(M)-y(j) ) * dxInv
  !      jp1             = min( M - ceiling( ydx ), M-1 )
  !      jp2             = jp1 + 1

  !      coef(2)         = dble( M - jp1 ) - ydx
  !      coef(1)         = 1.d0 - coef(2)
  !      do i=1, N
  !         fout(i,j) = coef(1)*fin(i,jp1) + coef(2)*fin(i,jp2)
  !      enddo
  !   enddo

  !   return
  ! end subroutine LinearInterpRET

  ! subroutine LinearInterpRET( y2, fxy2, y1, fxy1, LI, LJ, sftype )
  !   implicit none
  !   integer                         :: i
  !   integer         , intent(in)    :: LI, LJ
  !   double precision, intent(inout) :: y2(0:LJ), fxy2(LI,0:LJ)
  !   double precision, intent(inout) :: y1(1:LJ), fxy1(LI,1:LJ)
  !   character(3)    , intent(in)    :: sftype

  !   ! -- Transform Back Grid  ( y2 -> y1 ) -- !
  !   ! do i=1, LI
  !   !    fxy2(i,0) = - fxy2(i,1)
  !   ! enddo
  !   do i=1, LI
  !      call SF1DInterpolate( y2(1:LJ), fxy2(i,1:LJ), y1(1:LJ), fxy1(i,1:LJ), LJ, LJ, sftype )
  !   enddo
  !   do i=1, LI
  !      fxy1(i, 1) = fxy2(i, 1)
  !      fxy1(i,LJ) = fxy2(i,LJ)
  !   enddo

  !   return
  ! end subroutine LinearInterpRET


