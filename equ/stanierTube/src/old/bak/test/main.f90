program main

  implicit none
  integer, parameter :: LJ1   = 30
  integer, parameter :: LJ2   = 60
  double precision   :: x1Len = 20.d0
  double precision   :: x2Len = 20.d0
  double precision   :: dx1, dx2
  double precision   :: x1(LJ1), x2(LJ2), y1(LJ1), y2(LJ2)
  integer            :: i, j

  dx1 = x1Len / dble( LJ1-1 )
  dx2 = x2Len / dble( LJ2-1 )
  do i=1, LJ1
     x1(i) = dx1 * dble(i-1)
  enddo
  do i=1, LJ2
     x2(i) = dx2 * dble(i-1)
  enddo
  do i=1, LJ1
     y1(i) = x1(i)**2
  enddo

  call SF1DInterpolate( x1, y1, x2, y2, LJ1, LJ2, '2nd' )
  
  do i=1, LJ1
     write(6,'(i10,1x,2(f10.4,1x),2(f10.4))') i, x1(i), y1(i), x2(i), y2(i)
  enddo
  do i=LJ1+1, LJ2
     write(6,'(i10,1x,2( 10x ,1x),2(f10.4))') i,               x2(i), y2(i)
  enddo

  

contains

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


  subroutine SF1DInterpolate( x1, y1, x2, y2, LJ1, LJ2, sftype )
    implicit none
    integer         , intent(in   ) :: LJ1    , LJ2
    double precision, intent(in   ) :: x1(LJ1), y1(LJ1)
    double precision, intent(inout) :: x2(LJ2), y2(LJ2)
    character(3)    , intent(in)    :: sftype
    integer                         :: j, jg, jp, jpp
    double precision                :: x1Floor, dx1Inv, retg
    double precision                :: sf(-2:2)

    ! -- Preparation   -- !
    x1Floor = x1(1)
    dx1Inv  = 1.d0 / ( x1(2) - x1(1) )
    ! -- Interpolation -- !
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


end program main
