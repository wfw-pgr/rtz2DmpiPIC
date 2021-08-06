module lInterpMod
contains

  subroutine LinearInterpFWD( x, fin, y, fout, N, M )
    implicit none
    integer                       :: i, j, jp1, jp2
    double precision              :: dx, dy, dxinv, coef(2), ydx
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: x(1:N+1), fin( M,1:N+1)
    double precision, intent(out) :: y(0:N+1), fout(M,0:N+1)

    ! -- dx, dy, y-coordinate -- !
    dx    = x(2) - x(1)
    dxinv = 1.d0 / dx
    dy    = 2.d0 * x(N+1) / ( 2*N+1 )
    do j=0, N+1
       y(j) = dy * ( dble(j) - 0.5d0 )
    enddo

    ! -- x & fin --> y, itrp -- !
    do j=1, N+1
       ydx             = ( y(N+1)-y(j) ) * dxinv
       jp1             = min( max( N+1 - ceiling( ydx ), 1 ), N )
       jp2             = jp1 + 1

       coef(2)         = dble( N+1 - jp1 ) - ydx
       coef(1)         = 1.d0 - coef(2)
       do i=1, M
          fout(i,j) = coef(1)*fin(i,jp1) + coef(2)*fin(i,jp2)
       enddo
    enddo

    return
  end subroutine LinearInterpFWD

  
  subroutine LinearInterpRET( x, fin, y, fout, N, M )
    implicit none
    integer                       :: i, j, jp1, jp2
    double precision              :: dx, dxinv, coef(2), ydx
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: x(0:N+1), y(1:N+1), fin( M,0:N+1)
    double precision, intent(out) :: fout(M,1:N+1)

    ! -- dx, dxinv -- !
    dx    = x(2) - x(1)
    dxinv = 1.d0 / dx
    ! -- Return    -- !
    do j=1, N+1
       ydx             = ( y(N+1)-y(j) ) * dxinv
       jp1             = min( N+1 - ceiling( ydx ), N )
       jp2             = jp1 + 1

       coef(2)         = dble( N+1 - jp1 ) - ydx
       coef(1)         = 1.d0 - coef(2)
       do i=1, M
          fout(i,j) = coef(1)*fin(i,jp1) + coef(2)*fin(i,jp2)
       enddo
    enddo

    return
  end subroutine LinearInterpRET

end module lInterpMod
