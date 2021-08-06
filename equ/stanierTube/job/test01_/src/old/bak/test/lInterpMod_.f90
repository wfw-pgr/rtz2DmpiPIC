module lInterpMod
contains

  subroutine LinearInterpFWD( x, fin, y, fout, N, M )
    implicit none
    integer                       :: i, j, jp1, jp2
    double precision              :: dx, dy, dxinv, coef(2), ydx
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: x(1:M), fin( N,1:M)
    double precision, intent(out) :: y(0:M), fout(N,0:M)

    ! -- New Grid -- !
    if ( y1(1).ne.0.d0 ) stop '[ERROR] x(1) != 0.0'
    !  - y1 Grid Info.   - !
    y1Len  = y1(M) - y1(1)
    dy1    = y1(2) - y1(1)
    dy1Inv =  1.d0 / dy1
    !  - Define y2 Grid  - !
    dy2    = y1Len / ( ( dble(M)-1.d0 ) + 0.5d0 )
    do j=0, M
       y2(j) = dy2 * ( dble(j) - 0.5d0 )
    enddo
    
    
    return
  end subroutine LinearInterpFWD

  
  subroutine LinearInterpRET( x, fin, y, fout, N, M )
    implicit none
    integer                       :: i, j, jp1, jp2
    double precision              :: dx, dxInv, coef(2), ydx
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: x(0:M), y(1:M), fin( N,0:M)
    double precision, intent(out) :: fout(N,1:M)

    ! -- dx, dxInv -- !
    dx    = x(2) - x(1)
    dxInv = 1.d0 / dx
    ! -- Return -- !
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

  
  subroutine SF1DInterpolate( x1, y1, x2, y2, LJ1, LJ2 )
    implicit none
    integer         , intent(in   ) :: LJ1    , LJ2
    double precision, intent(in   ) :: x1(LJ1), y1(LJ1)
    double precision, intent(inout) :: x2(LJ2), y2(LJ2)
    integer                         :: j, jg, jp
    double precision                :: x1Floor, dx1Inv, retg
    double precision                :: sf(-2:2)

    ! -- Preparation   -- !
    x1Floor = x1(1)
    dx1Inv  = x1(2) - x1(1)
    ! -- Interpolation -- !
    do jg=1, LJ1
       
       jp   = nint( ( x2(jg) - x1Floor ) * dx1Inv )
       sf   = shapef( jp, x2(jg), dx1Inv, sftype )
       jp   = jp + 1
       retg = 0.d0
       do j=-2, 2
          jpp  = min( LJ1, max( 1, jp+j ) )
          retg = retg + sf(j)*y1(jpp)
       enddo
       ret(jg) = retg
       
    enddo

    return
  end subroutine SF1DInterpolate

  
end module lInterpMod



subroutine rct2rct( xp, yp, data, ret, xg, yg, LIp, LJp, LIg, LJg, sftype, ierr )
  implicit none
  integer         , intent(in)  :: LIp, LJp, LIg, LJg
  double precision, intent(in)  :: xp(LIp), yp(LJp)
  double precision, intent(in)  :: xg(LIg), yg(LJg)
  double precision, intent(in)  :: data(LIp,LJp)
  character(3)    , intent(in)  :: sftype
  double precision, intent(out) ::  ret(LIg,LJg)
  integer         , intent(out) :: ierr
  integer                       :: i, j, ig, jg, ip, jp, ipp, jpp
  double precision              :: dx, dy, dxInv, dyInv, xMin, yMin, xMax, yMax
  double precision              :: xn(LIg), yn(LJg), sfx(-2:2), sfy(-2:2), sfxy, retg
  logical                       :: chk

  ! --- [0] Preparation --- !
  ierr  = 0
  dx    = xp(2) - xp(1)
  dy    = yp(2) - yp(1)
  dxInv = 1.d0 / dx
  dyInv = 1.d0 / dy
  xMin  = xp(1)
  yMin  = yp(1)
  ! -- ( xn, yn = xg, yg - xMin, yMin ) -- !
  do ig=1, LIg
     xn(ig) = xg(ig) - xMin
  enddo
  do jg=1, LJg
     yn(jg) = yg(jg) - yMin
  enddo
  ! -- Boundary Check :: x -- !
  xMin  = 0.d0
  xMax  = xp(LIp) - xp(1)
  yMin  = 0.d0
  yMax  = yp(LJp) - yp(1)
  do ig=1, LIg
     if ( xn(ig).lt.xMin ) then
        chk    = .true.
        xn(ig) = xMin
     endif
     if ( xn(ig).gt.xMax ) then
        chk    = .true.
        xn(ig) = xMax
     endif
  enddo
  ! -- Boundary Check :: y -- !
  do jg=1, LJg
     if ( yn(jg).lt.yMin ) then
        chk   = .true.
        yn(jg) = yMin
     endif
     if ( yn(jg).gt.yMax ) then
        chk   = .true.
        yn(jg) = yMax
     endif
  enddo
  if ( chk ) ierr = 1

  return
  
contains
end subroutine rct2rct


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

