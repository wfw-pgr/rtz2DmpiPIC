module NewtonRaphson
contains

  subroutine LinearNewtonRaphson1D( fx, rix, flag )
    !  --  Newton-Raphson with Linear Interpolation  --  !
    implicit none
    integer         , intent(inout) :: flag
    double precision, intent(in)    :: fx(:) ! input function
    double precision, intent(inout) :: rix   ! initial / return index of extremum value
    integer         , parameter     :: itermax = 1000
    double precision, parameter     :: convg   = 1.d-2
    double precision, parameter     :: Fract   = 0.05d0
    double precision                :: chMax   = 0.2d0
    integer                         :: k, ix, Nx
    double precision                :: fk, dfdx, alpha, deltax, chng

    !  --- [0] Preparation  --- !
    Nx    = size( fx )
    chMax = chMax * Nx
    flag  = 0

    !  --- [1] Main Loop    --- !
    do k=1, itermax
       ! ---  (1) Calculate Derivative  --- !
       ix     = ceiling( rix ) - 1
       alpha  = rix - dble( ix )
       fk     = alpha * fx(ix+1) + ( 1.d0 - alpha ) * fx(ix)
       dfdx   = fx(ix+1) - fx(ix)

       ! ---  (2) calculate the delta_x --- ! 
       deltax = - fk / dfdx
100    continue
       if ( deltax.ne.deltax ) then
          flag = -1
          write(6,*) ' [CAUTION]  DeltaX == NaN  :: Return Former rix [CAUTION] '
          exit
       endif
       if ( deltax.gt.chMax  ) then
          deltax = deltax * Fract
          goto 100
       endif

       ! ---  (3) step forward          --- !
       rix    = rix + deltax

       ! ---  (4) check the convergence --- !
       chng   = sqrt( deltax**2 )
       if ( chng.le.convg ) then
          flag = 1
          exit
       else if ( ( rix.le.0.d0 ).or.( rix.ge.dble(Nx) ) ) then
          flag = -1
          write(6,*) ' [CAUTION]  :: No solution --- return rix = 0.d0 --- :: [CAUTION] '
          exit
       endif
    enddo

    return
  end subroutine LinearNewtonRaphson1D
  

  subroutine cSplineNewtonRaphson1D( x, fx, xk, fk, flag )
    !  --  Newton-Raphson with cSpline Interpolation  --  !
    implicit none
    integer         , intent(inout) :: flag
    double precision, intent(in)    :: x(:), fx(:)    ! input function
    double precision, intent(inout) :: xk, fk   ! initial / return x, y-value
    integer         , parameter     :: itermax = 1000
    double precision, parameter     :: convg   = 1.d-2
    double precision, parameter     :: Fract   = 0.05d0
    double precision                :: chMax   = 0.2d0
    integer                         :: k, ix, Nx
    double precision                :: dfdx, alpha, deltax, chng

    !  --- [0] Preparation  --- !
    Nx    = size( fx )
    if ( size(x).ne.Nx ) stop ' [ERROR] x & fx are incompatible size @ cSplineNewtonRaphson1D [ERROR] '
    chMax = chMax * Nx * (x(2)-x(1))
    flag  = 0

    !  --- [1] Main Loop    --- !
    do k=1, itermax
       ! ---  (1) Calculate Derivative  --- !
       call cSpline1D( x, fx, Nx, xk, fk, dfdx )

       ! ---  (2) calculate the delta_x --- ! 
       deltax = - fk / dfdx
100    continue
       if ( deltax.ne.deltax ) then
          flag = -1
          write(6,*) ' [CAUTION]  DeltaX == NaN  :: Return Former xk [CAUTION] '
          exit
       endif
       if ( deltax.gt.chMax  ) then
          deltax = deltax * Fract
          goto 100
       endif

       ! ---  (3) step forward          --- !
       xk    = xk + deltax

       ! ---  (4) check the convergence --- !
       chng   = sqrt( deltax**2 )
       if ( chng.le.convg ) then
          flag = 1
          exit
       else if ( ( xk.lt.x(1) ).or.( xk.gt.x(Nx) ) ) then
          flag = -1
          write(6,*) ' [ERROR]  :: No solution --- Stop --- :: [ERROR] '
          stop
       endif
    enddo

    return
  end subroutine CSplineNewtonRaphson1D


  subroutine NewtonRaphson2D( x, y, f, g, xk, yk, dfdx, dfdy, dgdx, dgdy, flag )
    ! -- Find Extremum by means of "Newton-Raphson with 2D spline Interpolation" -- !
    implicit none
    double precision, intent(in)    :: x(:), y(:)
    double precision, intent(in)    :: f(:,:), g(:,:)
    double precision, intent(inout) :: xk, yk
    double precision, intent(out)   :: dfdx, dfdy, dgdx, dgdy
    integer         , intent(out)   :: flag
    integer         , parameter     :: iteratemax = 1000
    double precision, parameter     :: convg      = 5.d-2
    double precision, parameter     :: alpha      = 0.05d0
    double precision                :: chMax      = 0.3d0
    integer          :: k, Nx, Ny
    double precision :: fk, gk
    double precision :: deltax, deltay
    double precision :: det, chng, crit
    
    !  --- [0] preparation --- !
    Nx    = size( x )
    Ny    = size( y )
    chMax = chMax*max( x(Nx)-x(2),y(Ny)-y(2) )
    flag  = 0
    crit  = convg*min( x(2)-x(1), y(2)-y(1) )

    !  --- [1] Main Loop  --- !
    do k=0, iteratemax
       !  --- [1-1] make Jacobian  --- !
       call cSpline2D( x, y, f, xk, yk, fk, dfdx, dfdy )
       call cSpline2D( x, y, g, xk, yk, gk, dgdx, dgdy )
       det    = 1.d0 / ( (dfdx * dgdy) - ( dfdy * dgdx ) )

       !  --- [1-2] Calculate Delta_xy --- !
       deltax = ( (dfdy * gk) - (dgdy * fk) )*det
       deltay = ( (dgdx * fk) - (dfdx * gk) )*det
       chng   = sqrt( deltax**2 + deltay**2 )
       if ( chng.gt.chMax ) then
          deltax = alpha * deltax
          deltay = alpha * deltay
       endif
       !  --- [1-3] Step Forward       --- !
       xk     = xk + deltax*alpha
       yk     = yk + deltay*alpha
       !  --- [1-4] Check Convergence  --- !
       if ( chng.le.crit ) then
          flag= 1        ! converge          
       endif
       if ( ( xk.lt.x(1) ).or.( yk.lt.y(1) ).or.( xk.gt.x(Nx) ).or.( yk.gt.y(Ny) ) ) then
          flag=-1        ! out of boundary
          ! write(6,*) ' [ERROR] 2D Newton-Raphson Solver :: Out of Boundary [ERROR] '
       endif
    enddo
    
    return
  end subroutine NewtonRaphson2D

  
  subroutine FindExtremum2D( x, y, f, xk, yk, fext, det, fxx, flag )
    ! -- Wrapper for "2D-spline-Newton-Raphson" -- !
    implicit none
    double precision, intent(in)    :: x(:), y(:)
    double precision, intent(in)    :: f(:,:)
    double precision, intent(inout) :: xk, yk
    double precision, intent(out)   :: fext, det, fxx
    integer         , intent(out)   :: flag
    integer                         :: i, j, Nx, Ny
    double precision                :: fxy, fyx, fyy
    double precision                :: dfextdx, dfextdy
    double precision, allocatable   :: dfdx(:,:), dfdy(:,:)
    double precision, allocatable   :: dxinv(:), dyinv(:)
    
    !  --- [0] preparation --- !
    !    -- [0-1] Allocation -- !
    Nx = size( x )
    Ny = size( y )
    allocate( dxinv(Nx-1), dyinv(Ny-1) )
    allocate( dfdx(2:Nx-1,2:Ny-1), dfdy(2:Nx-1,2:Ny-1) )
    !    -- [0-2] dx, dy -- !
    do i=1,Nx-1
       dxinv(i) = 1.d0 / ( x(i+1) - x(i) )
    enddo
    do j=1,Ny-1
       dyinv(j) = 1.d0 / ( y(j+1) - y(j) )
    enddo
    !    -- [0-3] Calculate Partial Derivative of f -- !
    do j=2, Ny-1
       do i=2, Nx-1
          dfdx(i,j) = ( f(i+1,j) - f(i-1,j) )*dxinv(i)
       enddo
    enddo
    do j=2, Ny-1
       do i=2, Nx-1
          dfdy(i,j) = ( f(i,j+1) - f(i,j-1) )*dyinv(j)
       enddo
    enddo
    
    !  --- [1] Newton-Raphson 2D --- !
    call NewtonRaphson2D( x(2:Nx-1), y(2:Ny-1), dfdx, dfdy, xk, yk, fxx, fxy, fyx, fyy, flag )
    if ( flag.eq.1 ) then 
       det = fxx*fyy - fxy*fyx
    endif
    
    !  --- [2] Interpolation of Answer --- !
    call cSpline2D( x, y, f, xk, yk, fext, dfextdx, dfextdy )
    
    return
  end subroutine FindExtremum2D

  
  subroutine FindExtremum2D_ZigZag( fxy, rik, rjk, flag )
    ! -- Find Extremum of fxy(i,j) by means of "Iterative Zig-Zag 1D search" -- !
    implicit none
    integer                         :: i, j, ik, jk, iter, Nx, Ny, inflag, flagsum
    double precision                :: rik_old, rjk_old
    double precision, allocatable   :: dfdx(:,:), dfdy(:,:), fx(:), fy(:)
    integer         , intent(out)   :: flag
    double precision, intent(in)    :: fxy(:,:)
    double precision, intent(inout) :: rik, rjk
    double precision, parameter     :: convg   = 1.d-3
    integer         , parameter     :: itermax = 100

    !  --- [1] Derivative --- !
    Nx      = ubound( fxy, 1 )
    Ny      = ubound( fxy, 2 )
    flag    = 0
    inflag  = 0
    flagsum = 0
    allocate( dfdx(Nx-1,Ny-1), dfdy(Nx-1,Ny-1), fx(Nx-1), fy(Ny-1) )
    do j=2, Ny-1
       do i=2, Nx-1
          dfdx(i,j) = fxy(i+1,j) - fxy(i-1,j)
          dfdy(i,j) = fxy(i,j+1) - fxy(i,j-1)
       enddo
    enddo
    rik_old = - rik
    rjk_old = - rjk

    !  --- [2] Main Loop  --- !
    do iter=1, itermax
       jk      = nint( rjk )
       call LinearNewtonRaphson1D( dfdx(2:Nx-1,jk), rik, inflag )
       flagsum = flagsum + max( - inflag, 0 )
       ik      = nint( rik )
       call LinearNewtonRaphson1D( dfdy(ik,2:Ny-1), rjk, inflag )
       flagsum = flagsum + max( - inflag, 0 )

       if ( flagsum.eq.2 ) then
          flag = -1
          write(6,*) ' [CAUTION] Failed to converge in Newton-Raphson1D [CAUTION] '
          return
       endif
       if ( ( ( rik-rik_old )**2 + ( rjk-rjk_old )**2 ).lt.convg ) then
          flag = +1
          return
       end if
       rik_old = rik
       rjk_old = rjk
    enddo

    return
  end subroutine FindExtremum2D_ZigZag
  

  ! ====================================================== !
  ! === cubic spline 1D                                === !
  ! ====================================================== !
  subroutine cSpline1D( x, y, M , xnew, ynew, ynewprime)
    implicit none
    integer         , intent(in)    :: M
    double precision, intent(in)    :: x(M), y(M), xnew ! [caution] x(N+1) !!!!
    double precision, intent(inout) :: ynew, ynewprime
    integer                         :: i, j, sizexnew
    double precision                :: h(M-1), hInv(M-1), v(2:M-1), u(M)
    double precision                :: mat(2:M-1,M)
    double precision                :: a(M-1), b(M-1), c(M-1), d(M-1)

    ! ------------------------------------------------------ !
    ! --- [1] prepare :: h(j) : dx(j)                    --- !
    ! ------------------------------------------------------ !
    do j=1, M-1
       h(j)    = x(j+1) - x(j)
       if ( h(j) == 0.d0 ) then
          write(6,*) "[cSpline1D @ NewtonRaphson.f90] x-vector is not Monotonic... [ERROR] "
          do i=1, M
             write(6,*) i, j, x(i), y(i)
          enddo
          stop
       else
          hInv(j) = 1.d0 / h(j)
       endif
    enddo
    ! ------------------------------------------------------ !
    ! --- [2] prepare :: v(j) : derivative = dy/dx       --- !
    ! ------------------------------------------------------ !
    do j=2, M-1
       v(j) = 6.0d0*( ( y(j+1)-y(j  ) )*hInv(j  ) - &
            &         ( y(j  )-y(j-1) )*hInv(j-1)   )
    enddo
    ! ------------------------------------------------------ !
    ! --- [3] prepare :: u(j) : d2y / dx2 = f''          --- !
    ! ------------------------------------------------------ !
    u(:) = 0.0d0
    !    -- [2-1] B.C. Disignation ( f''= 0 (natural spline) ) -- !
    u(1) = 0.0d0
    u(M) = 0.0d0

    ! ------------------------------------------------------ !
    ! --- [4] Solve equations for u(j)                   --- !
    ! ------------------------------------------------------ !
    !  -- [4-1] Preparation :: Make Matrix               --  !
    mat(:,:) = 0.0d0
    do i=2, M-1
       mat(i,i-1) = h(i-1)
       mat(i,i)   = 2.0d0 * ( h(i-1)+h(i) )
       mat(i,i+1) = h(i)
    enddo
    !  -- [4-2] Solve Equation                           --  !
    call TriDiag( M-1, u, h, mat )
    !  -- [4-3] Obtain Coefficients                      --  !
    call GetCoef( u, x, y, a, b, c, d)
    
    ! ------------------------------------------------------ !
    ! --- [5] Interpolation using coefficients           --- !
    ! ------------------------------------------------------ !
    ynew      = SplInterp_Value( x, y, xnew, a, b, c, d )
    ynewprime = SplInterp_Deriv( x,    xnew, a, b, c, d )

    return
  end subroutine cSpline1D


  ! ====================================================== !
  ! === cubic-spline 2D ( for grided 2D Data:fxy )     === !
  ! ====================================================== !
  subroutine cSpline2D( x, y, fxy, xk, yk, fk, dfdx, dfdy )
    implicit none
    double precision, intent(in)  :: xk, yk
    double precision, intent(in)  :: x(:), y(:), fxy(:,:)
    double precision, intent(out) :: fk, dfdx, dfdy
    integer                       :: i, j, Nx, Ny
    double precision              :: fnew, fnewprime
    double precision              :: xtoy, ytox
    double precision, allocatable :: xinterpd(:), yinterpd(:)

    ! ------------------------------------------------------ !
    ! --- [1] Preparation                                --- !
    ! ------------------------------------------------------ !
    Nx = size( x )
    Ny = size( y )
    allocate( xinterpd(Ny), yinterpd(Nx) )

    ! ------------------------------------------------------ !
    ! --- [2] calculate df/dy at (xk,yk)                 --- !
    ! ------------------------------------------------------ !
    do j=1, Ny
       call cSpline1D( x, fxy(:,j), Nx, xk, fnew, fnewprime )
       xinterpd(j) = fnew
    enddo
    call cSpline1D( y, xinterpd, Ny, yk, fnew, fnewprime )
    xtoy = fnew
    dfdy = fnewprime

    ! ------------------------------------------------------ !
    ! --- [3] calculate df/dx at  (xk,yk)                --- !
    ! ------------------------------------------------------ !
    do i=1, Nx
       call cSpline1D( y, fxy(i,:), Ny, yk, fnew, fnewprime )
       yinterpd(i) = fnew
    enddo
    call cSpline1D( x, yinterpd, Nx, xk, fnew, fnewprime )
    ytox = fnew
    dfdx = fnewprime

    ! ------------------------------------------------------ !
    ! --- [4] choose whether xtoy / ytox                 --- !
    ! ------------------------------------------------------ !
    fk   = xtoy
    ! fk   = ytox

    return
  end subroutine cSpline2D

  
  function SplInterp_Value( x, y, xnew, a, b, c, d )
    
    implicit none
    integer                      :: i, j, count, N
    double precision             :: SplInterp_Value
    double precision             :: ynew, Sprime
    double precision, intent(in) :: x(0:), y(0:), xnew
    double precision, intent(in) :: a(0:), b(0:), c(0:), d(0:)

    !  --- [0] Preparation   --- !
    N = size( x ) - 1
    
    !  --- [1] Interpolation --- !
    if ( xnew.lt.x(0) ) then
       Sprime = 3.d0*a(0)*( xnew-x(0) )**2 &
            & + 2.d0*b(0)*( xnew-x(0) ) + c(0)
       ynew   = Sprime*( xnew-x(0) ) + y(0)
    else if ( ( xnew.ge.x(0) ).and.( xnew.lt.x(N) ) ) then
       do j=0, N
          if ( xnew.ge.x(j).and.( xnew.lt.x(j+1)) ) then
             ynew = a(j)*( xnew-x(j) )**3 &
                  + b(j)*( xnew-x(j) )**2 &
                  + c(j)*( xnew-x(j) ) + d(j)
             exit
          endif
       enddo
    else
       Sprime = 3.d0*a(N)*( xnew-x(N) )**2 &
            & + 2.d0*b(N)*( xnew-x(N) ) + c(N)
       ynew = Sprime*( xnew-x(N) ) + y(N)
    endif
    
    !  --- [2] Return Value  --- !
    SplInterp_Value = ynew
        
    return
  end function SplInterp_Value

  
  function SplInterp_Deriv( x, xnew, a, b, c, d )
    
    implicit none
    integer                      :: i, j, count, N
    double precision             :: SplInterp_Deriv
    double precision             :: ynewprime, Sprime
    double precision, intent(in) :: x(0:), xnew
    double precision, intent(in) :: a(0:), b(0:), c(0:), d(0:)

    !  --- [0] Preparation   --- !
    N = size( x ) - 1
    
    !  --- [1] Interpolation --- !
    if ( xnew .lt. x(0) ) then
       Sprime = 3.d0*a(0)*( xnew-x(0) )**2 &
            & + 2.d0*b(0)*( xnew-x(0) ) + c(0)
       ynewprime = Sprime
       
    else if ( ( xnew.ge.x(0) ).and.( xnew.lt.x(N) ) ) then
       do j=0, N
          if ( xnew.ge.x(j).and.( xnew.lt.x(j+1)) ) then
             ynewprime = 3.d0*a(j)*( xnew-x(j) )**2 &
                  + 2.d0*b(j)*( xnew-x(j) ) &
                  + c(j)
             exit
          endif
       enddo

    else
       Sprime = 3.d0*a(N)*( xnew-x(N) )**2 &
            & + 2.d0*b(N)*( xnew-x(N) ) + c(N)
       ynewprime = Sprime
    endif
    
    !  --- [2] Return Value  --- !
    SplInterp_Deriv = ynewprime
        
    return
  end function SplInterp_Deriv

  ! poisson solver using tridiagonal algorithm
  subroutine TriDiag( N, ux, rhs, mat )
    implicit none
    integer         , intent(in)    :: N
    double precision, intent(in)    :: rhs(1:N-1)
    double precision, intent(in)    :: mat(1:N-1,0:N)
    double precision, intent(inout) :: ux(0:N)
    double precision                :: A(0:N), B(0:N), C(0:N), D(0:N), E(0:N), F(0:N)

    call SetMatrix( N, ux, rhs, mat, A, B, C, D, E, F )
    call calc_EF  ( N, A , B, C, D, E, F )
    call calc_ux  ( N, ux, B, C, D, E, F )

    return
  end subroutine TriDiag

  
  subroutine GetCoef( u, x, y, a, b, c, d )
    
    implicit none
    double precision, intent(in)  :: u(0:), x(0:), y(0:)
    double precision, intent(out) :: a(0:), b(0:), c(0:), d(0:)    
    integer                       :: j, N

    N = size(u)-1
    do j=0, N-1
       b(j) = 0.5d0*u(j)
       a(j) = ( u(j+1)-u(j) )/( 6.0d0 * ( x(j+1)-x(j) ) )
       d(j) = y(j)
       c(j) = ( y(j+1)-y(j) )/( x(j+1)-x(j) ) &
            & - ( x(j+1)-x(j) )*( 2.d0 * u(j) +u(j+1) )/6.d0
    enddo
    
    return
  end subroutine GetCoef

  ! set initial condition
  subroutine SetMatrix( N, ux, rhs, mat, A, B, C, D, E, F )

    implicit none
    integer         , intent(in)    :: N
    double precision, intent(in)    :: ux(0:N), rhs(1:N-1), mat(1:N-1,0:N)
    double precision, intent(inout) :: A(0:N), B(0:N), C(0:N), D(0:N)
    double precision, intent(inout) :: E(0:N), F(0:N)
    integer                         :: i

    !  --- [1] Initialize --- !
    E(:) = 0.0d0
    F(:) = 0.0d0
    
    !  --- [2] Set Coefficients --- !
    !    --  [2-1] For i=0      -- !
    A(0) =  0.0d0
    B(0) = -1.0d0
    C(0) =  0.0d0
    D(0) =  ux(0)
    !    --  [2-2] For i=1->N-1 -- !
    do i=1, N-1
       A(i) =  mat(i,i+1)
       B(i) = -mat(i,i)
       C(i) =  mat(i,i-1)
       D(i) =  rhs(i)
    enddo    
    !    --  [2-3] For i=N      -- !
    A(N) = 0.0d0
    B(N) = -1.0d0
    C(N) = 0.0d0
    D(N) = ux(N)
    
    return
  end subroutine SetMatrix

  ! calculate coefficients of E and F
  subroutine calc_EF( N, A, B, C, D, E, F )

    implicit none
    integer         , intent(in)    :: N
    double precision, intent(in)    :: A(0:N), B(0:N), C(0:N), D(0:N)
    double precision, intent(inout) :: E(0:N), F(0:N)
    integer                         :: i
    double precision                :: inv

    !  --- [1] i = 0 --- !
    E(0) =   A(0)/B(0)
    F(0) = - D(0)/B(0)
    
    !  --- [2] i = [ 1, N ] --- !
    do i = 1, N
       inv  = 1.0d0 / ( B(i) - C(i)*E(i-1) )
       E(i) =   A(i)*inv
       F(i) = ( C(i)*F(i-1)-D(i) ) * inv
    enddo

    return
  end subroutine calc_EF

  ! solve equation by stepping forward recursion equation
  subroutine calc_ux( N, ux, B, C, D, E, F )

    implicit none
    integer         , intent(in)    :: N
    double precision, intent(inout) :: ux(0:N)
    double precision, intent(in)    :: B(0:N), C(0:N), D(0:N), E(0:N), F(0:N)
    integer                         :: i
    
    !  --- [1] For i = N-1 --- !
    ux( N-1 ) = ( B(N)*F(N-1) - D(N)*E(N-1) ) / ( B(N) - E(N-1)*C(N) )
    
    !  --- [2] For i < N-1 --- !
    do i = N-2, 1, -1
       ux(i)  = E(i)*ux(i+1) + F(i)
    enddo

    return
  end subroutine calc_ux
  

end module NewtonRaphson
