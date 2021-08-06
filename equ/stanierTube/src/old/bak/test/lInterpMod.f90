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
    if ( x(1).ne.0.d0 ) stop '[ERROR] x(1) != 0.0'
    dx    = x(2) - x(1)
    dxinv = 1.d0 / dx
    dy    = 2.d0 * x(M) / dble( 2*M-1 )
    do j=0, M
       y(j) = dy * ( dble(j) - 0.5d0 )
    enddo
    ! -- Interpolation -- !
    !  - x & fin --> y, itrp - !
    do j=1, M
       ydx             = ( y(M)-y(j) ) * dxinv
       jp1             = min( max( M - ceiling( ydx ), 1 ), M-1 )
       jp2             = jp1 + 1

       coef(2)         = dble( M - jp1 ) - ydx
       coef(1)         = 1.d0 - coef(2)
       do i=1, N
          fout(i,j) = coef(1)*fin(i,jp1) + coef(2)*fin(i,jp2)
       enddo
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

end module lInterpMod

