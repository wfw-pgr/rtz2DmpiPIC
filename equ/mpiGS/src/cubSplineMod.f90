module cubSplineMod
contains

  ! ====================================================== !
  ! === description                                    === !
  ! ====================================================== !
  ! ---                                                --- !
  ! ---  * call interpolate__cubicSpline1D_point       --- !
  ! ---         for 1-point data ( double )            --- !
  ! ---  * call interpolate__cubicSpline1D for Array   --- !
  ! ---  * call solve__cubicSpline to obtain coef      --- !
  ! ---  * call return__value_cSpline to obtain value  --- !
  ! ---         using coef calculated before.          --- !
  ! ---  * call interpolate__cubicSpline2D for 2D      --- !
  ! ---                                                --- !
  ! ====================================================== !


  ! ====================================================== !
  ! ===  interpolate__cubicSpline1D_point              === !
  ! ====================================================== !
  subroutine interpolate__cubicSpline1D_point( xD, yD, xI, yI, dydxI )
    implicit none
    double precision, intent(in)    :: xD(:), yD(:)
    double precision, intent(in)    :: xI
    double precision, intent(inout) :: yI, dydxI
    double precision, allocatable   :: coef(:,:)
    integer                         :: LD, LI
    double precision                :: xI_(1), yI_(1), dydxI_(1)

    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !
    LD = size( xD )
    if ( .not.( allocated( coef ) ) ) then
       allocate( coef(4,LD), source=0.d0 )
    endif
    ! ------------------------------------------------------ !
    ! --- [2] solve & interpolate                        --- !
    ! ------------------------------------------------------ !
    xI_(1) = xI
    call solve__cubicSpline   ( xD, yD, coef, LD )
    call return__value_cSpline( xD, yD, xI_, yI_, dydxI_, coef, LD, 1 )

    yI    =    yI_(1)
    dydxI = dydxI_(1)
    
    return
  end subroutine interpolate__cubicSpline1D_point

  
  ! ====================================================== !
  ! ===  interpolate__cubicSpline1D                    === !
  ! ====================================================== !
  subroutine interpolate__cubicSpline1D( xD, yD, xI, yI, dydxI )
    implicit none
    double precision, intent(in)    :: xD(:), yD(:)
    double precision, intent(inout) :: xI(:), yI(:), dydxI(:)
    double precision, allocatable   :: coef(:,:)
    integer                         :: LD, LI

    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !
    LD = size( xD )
    LI = size( xI )
    if ( .not.( allocated( coef ) ) ) then
       allocate( coef(4,LD), source=0.d0 )
    endif
    ! ------------------------------------------------------ !
    ! --- [2] solve & interpolate                        --- !
    ! ------------------------------------------------------ !
    call solve__cubicSpline   ( xD, yD, coef, LD )
    call return__value_cSpline( xD, yD, xI, yI, dydxI, coef, LD, LI )
    return
  end subroutine interpolate__cubicSpline1D


  ! ====================================================== !
  ! === interpolate cubicSpline (from grided 2D Data)  === !
  ! ====================================================== !
  subroutine interpolate__cubicSpline2D_point( xD, yD, fxy, xk, yk, fk, dfdx, dfdy )
    implicit none
    double precision, intent(in)  :: xk, yk
    double precision, intent(in)  :: xD(:), yD(:), fxy(:,:)
    double precision, intent(out) :: fk, dfdx, dfdy
    integer                       :: i, j, Nx, Ny
    double precision              :: xI(1), yI(1), dfdxI(1), dfdyI(1), fI(1)
    double precision              :: xtoy, ytox
    double precision, allocatable :: fI_x(:), fI_y(:)
    character(4)    , parameter   :: type = "mean"       ! [ "mean", "xtoy", "ytox" ]
    
    ! ------------------------------------------------------ !
    ! --- [1] Preparation                                --- !
    ! ------------------------------------------------------ !
    Nx = size( xD )
    Ny = size( yD )
    allocate( fI_x(Nx), fI_y(Ny), source=0.d0 )

    ! ------------------------------------------------------ !
    ! --- [2] calculate df/dy at (xk,yk)                 --- !
    ! ------------------------------------------------------ !
    xI(1) = xk
    yI(1) = 0.d0
    do j=1, Ny
       call interpolate__cubicSpline1D( xD, fxy(:,j), xI, yI, dfdxI )
       fI_y(j) = yI(1)
    enddo
    yI(1) = yk
    fI(1) = 0.d0
    call interpolate__cubicSpline1D( yD, fI_y, yI, fI, dfdyI )
    xtoy  = fI(1)
    dfdy  = dfdyI(1)
    
    ! ------------------------------------------------------ !
    ! --- [3] calculate df/dx at (xk,yk)                 --- !
    ! ------------------------------------------------------ !
    yI(1) = yk
    xI(1) = 0.d0
    do i=1, Nx
       call interpolate__cubicSpline1D( yD, fxy(i,:), yI, xI, dfdyI )
       fI_x(i) = xI(1)
    enddo
    xI(1) = xk
    fI(1) = 0.d0
    call interpolate__cubicSpline1D( xD, fI_x, xI, fI, dfdxI )
    ytox  = fI(1)
    dfdx  = dfdxI(1)

    ! ------------------------------------------------------ !
    ! --- [4] choose return value [ xtoy / ytox / mean ] --- !
    ! ------------------------------------------------------ !
    if      ( trim(type) == "mean" ) then
       fk   = 0.5d0 * ( xtoy + ytox )
    else if ( trim(type) == "xtoy" ) then
       fk   = xtoy
    else if ( trim(type) == "ytox" ) then
       fk   = ytox
    endif
    ! return ( fk, dfdx, dfdy )
    return
  end subroutine interpolate__cubicSpline2D_point


  ! ====================================================== !
  ! === interpolate cubicSpline (from grided 2D Data)  === !
  ! ====================================================== !
  subroutine interpolate__cubicSpline2D( xD, yD, fxy, xk, yk, fk, dfdx, dfdy )
    implicit none
    double precision, intent(in)  :: xk(:), yk(:)
    double precision, intent(in)  :: xD(:), yD(:), fxy(:,:)
    double precision, intent(out) :: fk(:), dfdx(:), dfdy(:)
    integer                       :: ik, nData

    ! ------------------------------------------------------ !
    ! --- [1] size check                                 --- !
    ! ------------------------------------------------------ !
    nData = size( xk )
    if (   ( size(  yk).ne.nData ).or.( size(  fk).ne.nData ).or.&
         & ( size(dfdx).ne.nData ).or.( size(dfdy).ne.nData )      ) then
       write(6,*) "[interpolate__cubicSpline2D] incompatible size :: xk, yk, fk, dfdx, dfdy... "
       stop
    endif

    ! ------------------------------------------------------ !
    ! --- [2] Main Loop (iterate to call point routine)  --- !
    ! ------------------------------------------------------ !
    do ik=1, nData
       call interpolate__cubicSpline2D_point( xD, yD, fxy, xk(ik), yk(ik), fk(ik), &
            &                                 dfdx(ik), dfdy(ik) )
    enddo
    return
  end subroutine interpolate__cubicSpline2D
  
  
  ! ! ====================================================== !
  ! ! === cubic-spline 2D ( for grided 2D Data:fxy )     === !
  ! ! ====================================================== !
  ! subroutine cSpline2D( xD, yD, fxy, xk, yk, fk, dfdx, dfdy )
  !   implicit none
  !   double precision, intent(in)  :: xk, yk
  !   double precision, intent(in)  :: xD(:), yD(:), fxy(:,:)
  !   double precision, intent(out) :: fk, dfdx, dfdy
  !   integer                       :: i, j, Nx, Ny
  !   double precision              :: fnew, fnewprime
  !   double precision              :: xtoy, ytox
  !   double precision, allocatable :: xinterpd(:), yinterpd(:), xi_(:), yi_(:), coef(:,:)

  !   ! ------------------------------------------------------ !
  !   ! --- [1] Preparation                                --- !
  !   ! ------------------------------------------------------ !
  !   Nx = size( xD )
  !   Ny = size( yD )
  !   allocate( xI(1), fI_x(Ny), dfdxI(Ny), coef_x(4,Nx), source=0.d0 )
  !   allocate( yI(1), fI_y(Ny), dfdyI(Ny), coef_y(4,Nx), source=0.d0 )


  !   ! ------------------------------------------------------ !
  !   ! --- [2] calculate df/dy at (xk,yk)                 --- !
  !   ! ------------------------------------------------------ !
  !   xI(1) = xk
  !   yI(1) = yk
  !   do j=1, Ny
  !      call cSpline1D( xD, fxy(:,j), xI, fI_x(j), dfdxI(j), coef_x, Nx, 1 )
  !   enddo
  !   call cSpline1D( yD, fI_x, yI, fnew, dfdyI(j), coef_y, Ny, 1 )
  !   xtoy = fnew
  !   dfdy = fnewprime

  !   ! ! ------------------------------------------------------ !
  !   ! ! --- [3] calculate df/dx at  (xk,yk)                --- !
  !   ! ! ------------------------------------------------------ !
  !   ! do i=1, Nx
  !   !    call cSpline1D( yD, fxy(i,:), Ny, yk, fnew, fnewprime )
  !   !    yinterpd(i) = fnew
  !   ! enddo
  !   ! call cSpline1D( xD, yinterpd, Nx, xk, fnew, fnewprime )
  !   ! ytox = fnew
  !   ! dfdx = fnewprime

  !   ! ------------------------------------------------------ !
  !   ! --- [4] choose whether xtoy / ytox                 --- !
  !   ! ------------------------------------------------------ !
  !   fk   = xtoy
  !   ! fk   = ytox

  !   return
  ! end subroutine cSpline2D


  ! ====================================================== !
  ! === solve cubic spline                             === !
  ! ====================================================== !
  subroutine solve__cubicSpline( xD, yD, coef, LD )
    implicit none
    integer         , intent(in)    :: LD
    double precision, intent(in)    :: xD(LD), yD(LD)
    double precision, intent(inout) :: coef(4,LD)
    integer                         :: i, j
    double precision                :: hx(LD), hxInv(LD), vx(LD), ux(LD)
    double precision                :: mat(LD,-1:+1)
    integer         , parameter     :: a_=1, b_=2, c_=3, d_=4
    
    ! ---------------------------------------------- !
    ! --  [ range ]                               -- !
    ! --       * xD :: x   :: LD   => (1,LD  )    -- !
    ! --       * hx :: f'  :: LD-1 => (1,LD-1)    -- !
    ! --       * vx :: f'' :: LD-2 => (2,LD-1)    -- !
    ! --       * ux :: f'' :: LD-2 => (2,LD-1)    -- !
    ! --                                          -- !
    ! --      x1  x2  x3       xN-2 xN-1 xN       -- !
    ! -- xD   *---*---*---*-....-*---*---*        -- !
    ! --        |   |   |          |   |          -- !
    ! -- hx     h1  h2  h3       hN-2 hN-1 hN     -- |
    ! --          *   *   * .... *   *   *        -- !
    ! -- ux       u2  u3  u3   uN-2 uN-1          -- |
    ! -- vx       v2  v3  v3   vN-2 vN-1          -- |
    ! --                                          -- !
    ! ---------------------------------------------- !

    ! ------------------------------------------------------ !
    ! --- [1] prepare :: h(x_j) : Delta ( x_j+1 - x_j )  --- !
    ! ------------------------------------------------------ !
    hx(:) = 0.d0
    do j=1, LD-1
       hx(j)    = xD(j+1) - xD(j)
       if ( hx(j) == 0.d0 ) then
          write(6,*) "[solve__cubicSpline @ cubSplineMod] x-vector is not Monotonic... [ERROR] "
          do i=1, LD
             write(6,*) i, j, xD(i), yD(i)
          enddo
          stop
       else
          hxInv(j) = 1.d0 / hx(j)
       endif
    enddo
    
    ! ------------------------------------------------------ !
    ! --- [2] prepare :: v(j) : derivative = d2y/dx2     --- !
    ! ------------------------------------------------------ !
    vx(:) = 0.d0
    do j=2, LD-1
       vx(j) = 6.0d0*( ( yD(j+1)-yD(j  ) )*hxInv(j  ) - &
            &          ( yD(j  )-yD(j-1) )*hxInv(j-1)   )
    enddo
    
    ! ------------------------------------------------------ !
    ! --- [3] prepare :: u(j) : d2y / dx2 = f''          --- !
    ! ------------------------------------------------------ !
    !  -- [3-1] initialize                               --  !
    ux(:)  = 0.0d0
    !  -- [3-2] Boundary condition                       --  !
    !  --       ( natural spline )                       --  !
    ux( 1) = 0.0d0
    ux(LD) = 0.0d0

    ! ------------------------------------------------------ !
    ! --- [4] Solve equations for u(j)                   --- !
    ! ------------------------------------------------------ !
    !  -- [4-1] Preparation :: Make Matrix               --  !
    do i=2, LD-1
       mat(i,-1) =         hx(i-1)
       mat(i, 0) = 2.0d0*( hx(i-1) + hx(i) )
       mat(i,+1) =                   hx(i)
    enddo
    !  -- [4-2] Solve Equation                           --  !
    call solve__triDiagonalEqs( ux, vx, mat, LD )
    !  -- [4-3] Obtain Coefficients                      --  !
    call get__coefficients( ux, xD, yD, coef, LD )
        
    return
  end subroutine solve__cubicSpline


  ! ====================================================== !
  ! ===  interpolate cSpline 1D                        === !
  ! ====================================================== !
  subroutine return__value_cSpline( xD, yD, xI, yI, dydxI, coef, LD, LI )
    implicit none
    integer                         :: i, j, ik, count
    integer         , intent(in)    :: LD, LI
    double precision, intent(in)    :: xD(LD), yD(LD), coef(4,LD)
    double precision, intent(inout) :: xI(LI), yI(LI), dydxI(LI)
    integer         , parameter     :: a_=1, b_=2, c_=3, d_=4

    ! ------------------------------------------------------ !
    ! --- [1] initialize variables                       --- !
    ! ------------------------------------------------------ !
    yI   (:) = 0.d0
    dydxI(:) = 0.d0

    ! ------------------------------------------------------ !
    ! --- [2] do loop begin                              --- !
    ! ------------------------------------------------------ !
    do ik=1, LI
       
       if ( xI(ik).lt.xD(1) ) then
          ! s_prime =   3.d0*coef(a_,1)*( xD(1)-xD(1) )**2 &
          !      &    + 2.d0*coef(b_,1)*( xD(1)-xD(1) )    &
          !      &    +      coef(c_,1)
          dydxI(ik) = coef(c_,1)
          yI(ik)    = dydxI(ik)*( xI(ik)-xD(1) ) + yD(1)
       else if ( ( xI(ik).ge.xD(1) ).and.( xI(ik).lt.xD(LD) ) ) then
          do j=1, LD-1
             if ( xI(ik).ge.xD(j).and.( xI(ik).lt.xD(j+1)) ) then
                dydxI(ik) =   3.d0*coef(a_,j)*( xI(ik)-xD(j) )**2 &
                     &      + 2.d0*coef(b_,j)*( xI(ik)-xD(j) )    &
                     &      +      coef(c_,j)
                yI(ik)    =        coef(a_,j)*( xI(ik)-xD(j) )**3 &
                     &      +      coef(b_,j)*( xI(ik)-xD(j) )**2 &
                     &      +      coef(c_,j)*( xI(ik)-xD(j) )    &
                     &      +      coef(d_,j)
                exit
             endif
          enddo
       else
          ! -- max index of coef :: LD-1 &   -- !
          ! -- max index of   xD :: LD       -- !
          dydxI(ik) =   3.d0*coef(a_,LD-1)*( xD(LD)-xD(LD-1) )**2 &
               &      + 2.d0*coef(b_,LD-1)*( xD(LD)-xD(LD-1) )    &
               &      +      coef(c_,LD-1)
          yI(ik)    = dydxI(ik)*( xI(ik)-xD(LD) ) + yD(LD)
       endif
       
    end do
    return
  end subroutine return__value_cSpline


  ! ====================================================== !
  ! === get__coefficients :: calculate coef of spline  === !
  ! ====================================================== !
  subroutine get__coefficients( ux, xD, yD, coef, LD )
    implicit none
    integer         , intent(in)  :: LD
    double precision, intent(in)  :: ux(LD), xD(LD), yD(LD)
    double precision, intent(out) :: coef(4,LD)
    integer                       :: j
    double precision              :: dInv
    integer         , parameter   :: a_=1, b_=2, c_=3, d_=4
    double precision, parameter   :: onesixth = 1.d0 / 6.d0

    do j=1, LD-1
       dInv = 1.d0 / ( xD(j+1)-xD(j) )
       coef(b_,j) = 0.5d0*ux(j)
       coef(a_,j) = ( ux(j+1)-ux(j) ) * onesixth * dInv
       coef(d_,j) = yD(j)
       coef(c_,j) = ( yD(j+1)-yD(j) ) * dInv &
            &     - ( xD(j+1)-xD(j) ) * ( 2.d0*ux(j)+ux(j+1) ) * onesixth
    enddo
    return
  end subroutine get__coefficients

  
  ! ====================================================== !
  ! ===  Tridiagonal solover algorithm                 === !
  ! ====================================================== !
  subroutine solve__triDiagonalEqs( ux, vx, mat, LD )
    implicit none
    integer         , intent(in)    :: LD
    double precision, intent(in)    :: vx(LD), mat(LD,-1:+1)
    double precision, intent(inout) :: ux(LD)
    integer                         :: i
    double precision                :: inv
    double precision                :: triParts(6,LD)
    integer         , parameter     :: A_=1, B_=2, C_=3, D_=4, E_=5, F_=6

    ! ------------------------------------------------------ !
    ! --- [1] set matrix coefficients                    --- !
    ! ------------------------------------------------------ !
    call set__cSplineMatrix( ux, vx, mat, triParts, LD )

    ! ------------------------------------------------------ !
    ! --- [2] calculate E and F of the eqs.              --- !
    ! ------------------------------------------------------ !
    triParts(E_,1) =   triParts(A_,1)/triParts(B_,1)
    triParts(F_,1) = - triParts(D_,1)/triParts(B_,1)
    do i=2, LD-1
       inv            = 1.0d0 / ( triParts(B_,i) - triParts(C_,i)*triParts(E_,i-1) )
       triParts(E_,i) =   triParts(A_,i)*inv
       triParts(F_,i) = ( triParts(C_,i)*triParts(F_,i-1)-triParts(D_,i) ) * inv
    enddo

    ! ------------------------------------------------------ !
    ! --- [3] calculate u(x) :: A u = v                  --- !
    ! ------------------------------------------------------ !
    ux( LD-1 ) =   ( triParts(B_,LD)*triParts(F_,LD-1) - triParts(D_,LD)*triParts(E_,LD-1) ) &
         &       / ( triParts(B_,LD) - triParts(E_,LD-1)*triParts(C_,LD) )
    do i = LD-2, 2, -1
       ux(i)  = triParts(E_,i)*ux(i+1) + triParts(F_,i)
    enddo

    return
  end subroutine solve__triDiagonalEqs

  
  ! ====================================================== !
  ! ===  set matrix coefficients                       === !
  ! ====================================================== !
  subroutine set__cSplineMatrix( ux, vx, mat, triParts, LD )
    ! --------------------------------------------------------- !
    ! --- see Yamamoto Masashi's web page or others....     --- !
    ! --- http://www.yamamo10.jp/yamamoto/lecture/2006/5E/  -->
    ! ---   --> interpolation/interpolation_html/node3.html --- !
    ! --------------------------------------------------------- !
    implicit none
    integer         , intent(in)    :: LD
    double precision, intent(in)    :: ux(LD), vx(LD), mat(LD,-1:+1)
    double precision, intent(inout) :: triParts(6,LD)
    integer                         :: i
    integer         , parameter     :: A_=1, B_=2, C_=3, D_=4, E_=5, F_=6

    ! ------------------------------------------------------ !
    ! --- [1] initialize                                 --- !
    ! ------------------------------------------------------ !
    triParts(E_,:)  = 0.0d0
    triParts(F_,:)  = 0.0d0
    ! ------------------------------------------------------ !
    ! --- [2] Set Coefficients                           --- !
    ! ------------------------------------------------------ !
    !  -- [2-1] for i=0, LD                              --  !
    triParts(A_,1)  =  0.0d0
    triParts(A_,LD) =  0.0d0
    triParts(B_,1)  = -1.0d0
    triParts(B_,LD) = -1.0d0
    triParts(C_,1)  =  0.0d0
    triParts(C_,LD) =  0.0d0
    triParts(D_,1)  =  ux( 1)
    triParts(D_,LD) =  ux(LD)
    ! -- [2-2] for i=1 >> LD-1                           --  !
    do i=2, LD-1
       triParts(A_,i) =  mat(i,+1)
       triParts(B_,i) = -mat(i, 0)
       triParts(C_,i) =  mat(i,-1)
       triParts(D_,i) =   vx(i)
    enddo
    return
  end subroutine set__cSplineMatrix
  

end module CubSplineMod





  
  ! ! ====================================================== !
  ! ! === cubic spline 1D                                === !
  ! ! ====================================================== !
  ! subroutine cSpline1D( xD, yD, xI, yI, dydxI, coef, LD, LI )
  !   implicit none
  !   integer         , intent(in)    :: LD, LI
  !   double precision, intent(in)    :: xD(LD), yD(LD), xI(LI)
  !   double precision, intent(inout) :: yI(LI), dydxI(LI), coef(4,LD)
  !   integer                         :: i, j
  !   double precision                :: hx(LD), hxInv(LD), vx(LD), ux(LD)
  !   double precision                :: mat(LD,-1:+1)
  !   integer         , parameter     :: a_=1, b_=2, c_=3, d_=4
    
  !   ! ---------------------------------------------- !
  !   ! --  [ range ]                               -- !
  !   ! --       * xD :: x   :: LD   => (1,LD  )    -- !
  !   ! --       * hx :: f'  :: LD-1 => (1,LD-1)    -- !
  !   ! --       * vx :: f'' :: LD-2 => (2,LD-1)    -- !
  !   ! --       * ux :: f'' :: LD-2 => (2,LD-1)    -- !
  !   ! --                                          -- !
  !   ! --      x1  x2  x3       xN-2 xN-1 xN       -- !
  !   ! -- xD   *---*---*---*-....-*---*---*        -- !
  !   ! --        |   |   |          |   |          -- !
  !   ! -- hx     h1  h2  h3       hN-2 hN-1 hN     -- |
  !   ! --          *   *   * .... *   *   *        -- !
  !   ! -- ux       u2  u3  u3   uN-2 uN-1          -- |
  !   ! -- vx       v2  v3  v3   vN-2 vN-1          -- |
  !   ! --                                          -- !
  !   ! ---------------------------------------------- !

  !   ! ------------------------------------------------------ !
  !   ! --- [1] prepare :: h(x_j) : Delta ( x_j+1 - x_j )  --- !
  !   ! ------------------------------------------------------ !
  !   hx(:) = 0.d0
  !   do j=1, LD-1
  !      hx(j)    = xD(j+1) - xD(j)
  !      if ( hx(j) == 0.d0 ) then
  !         write(6,*) "[cSpline1D @ NewtonRaphson.f90] x-vector is not Monotonic... [ERROR] "
  !         do i=1, LD
  !            write(6,*) i, j, xD(i), yD(i)
  !         enddo
  !         stop
  !      else
  !         hxInv(j) = 1.d0 / hx(j)
  !      endif
  !   enddo
    
  !   ! ------------------------------------------------------ !
  !   ! --- [2] prepare :: v(j) : derivative = d2y/dx2     --- !
  !   ! ------------------------------------------------------ !
  !   vx(:) = 0.d0
  !   do j=2, LD-1
  !      vx(j) = 6.0d0*( ( yD(j+1)-yD(j  ) )*hxInv(j  ) - &
  !           &          ( yD(j  )-yD(j-1) )*hxInv(j-1)   )
  !   enddo
    
  !   ! ------------------------------------------------------ !
  !   ! --- [3] prepare :: u(j) : d2y / dx2 = f''          --- !
  !   ! ------------------------------------------------------ !
  !   !  -- [3-1] initialize                               --  !
  !   ux(:)  = 0.0d0
  !   !  -- [3-2] Boundary condition                       --  !
  !   !  --       ( natural spline )                       --  !
  !   ux( 1) = 0.0d0
  !   ux(LD) = 0.0d0

  !   ! ------------------------------------------------------ !
  !   ! --- [4] Solve equations for u(j)                   --- !
  !   ! ------------------------------------------------------ !
  !   !  -- [4-1] Preparation :: Make Matrix               --  !
  !   do i=2, LD-1
  !      mat(i,-1) =         hx(i-1)
  !      mat(i, 0) = 2.0d0*( hx(i-1) + hx(i) )
  !      mat(i,+1) =                   hx(i)
  !   enddo
  !   !  -- [4-2] Solve Equation                           --  !
  !   call solve__triDiagonalEqs( ux, vx, mat, LD )
  !   !  -- [4-3] Obtain Coefficients                      --  !
  !   call get__coefficients( ux, xD, yD, coef, LD )
    
  !   ! ------------------------------------------------------ !
  !   ! --- [5] Interpolation using coefficients           --- !
  !   ! ------------------------------------------------------ !
  !   do i=1, LI
  !      yI(i)      = interpolate__cSpline1D           ( xD, yD, xI(i), coef, LD )
  !      dydxI(i)   = interpolate__derivative_cSpline1D( xD,     xI(i), coef, LD )
  !   enddo
    
  !   return
  ! end subroutine cSpline1D


  ! ! ====================================================== !
  ! ! === cubic-spline 2D ( for grided 2D Data:fxy )     === !
  ! ! ====================================================== !
  ! subroutine cSpline2D( xD, yD, fxy, xk, yk, fk, dfdx, dfdy )
  !   implicit none
  !   double precision, intent(in)  :: xk, yk
  !   double precision, intent(in)  :: xD(:), yD(:), fxy(:,:)
  !   double precision, intent(out) :: fk, dfdx, dfdy
  !   integer                       :: i, j, Nx, Ny
  !   double precision              :: fnew, fnewprime
  !   double precision              :: xtoy, ytox
  !   double precision, allocatable :: xinterpd(:), yinterpd(:), xi_(:), yi_(:), coef(:,:)

  !   ! ------------------------------------------------------ !
  !   ! --- [1] Preparation                                --- !
  !   ! ------------------------------------------------------ !
  !   Nx = size( xD )
  !   Ny = size( yD )
  !   allocate( xI(1), fI_x(Ny), dfdxI(Ny), coef_x(4,Nx), source=0.d0 )
  !   allocate( yI(1), fI_y(Ny), dfdyI(Ny), coef_y(4,Nx), source=0.d0 )


  !   ! ------------------------------------------------------ !
  !   ! --- [2] calculate df/dy at (xk,yk)                 --- !
  !   ! ------------------------------------------------------ !
  !   xI(1) = xk
  !   yI(1) = yk
  !   do j=1, Ny
  !      call cSpline1D( xD, fxy(:,j), xI, fI_x(j), dfdxI(j), coef_x, Nx, 1 )
  !   enddo
  !   call cSpline1D( yD, fI_x, yI, fnew, dfdyI(j), coef_y, Ny, 1 )
  !   xtoy = fnew
  !   dfdy = fnewprime

  !   ! ! ------------------------------------------------------ !
  !   ! ! --- [3] calculate df/dx at  (xk,yk)                --- !
  !   ! ! ------------------------------------------------------ !
  !   ! do i=1, Nx
  !   !    call cSpline1D( yD, fxy(i,:), Ny, yk, fnew, fnewprime )
  !   !    yinterpd(i) = fnew
  !   ! enddo
  !   ! call cSpline1D( xD, yinterpd, Nx, xk, fnew, fnewprime )
  !   ! ytox = fnew
  !   ! dfdx = fnewprime

  !   ! ------------------------------------------------------ !
  !   ! --- [4] choose whether xtoy / ytox                 --- !
  !   ! ------------------------------------------------------ !
  !   fk   = xtoy
  !   ! fk   = ytox

  !   return
  ! end subroutine cSpline2D


  ! ! ====================================================== !
  ! ! === interpolate derivative using cSpline         === !
  ! ! ====================================================== !
  ! function return__deriv_cSpline( xD, xI, coef, LD )
  !   implicit none
  !   integer                      :: i, j, count
  !   double precision             :: return__deriv_cSpline
  !   double precision             :: yIprime, s_prime
  !   integer         , intent(in) :: LD
  !   double precision, intent(in) :: xD(LD), coef(4,LD)
  !   double precision, intent(in) :: xI
  !   integer         , parameter  :: a_=1, b_=2, c_=3, d_=4

  !   ! ------------------------------------------------------ !
  !   ! --- [1] interpolate using ax**3 + bx**2 + cx + d   --- !
  !   ! ------------------------------------------------------ !
  !   if ( xI .lt. xD(1) ) then
  !      ! yIprime =   3.d0*coef(a_,1)*( xD(1)-xD(1) )**2 &
  !      !      &    + 2.d0*coef(b_,1)*( xD(1)-xD(1) )    &
  !      !      &    +      coef(c_,1)
  !      yIprime = coef(c_,1)
       
  !   else if ( ( xI.ge.xD(1) ).and.( xI.lt.xD(LD) ) ) then
  !      do j=1, LD-1
  !         if ( xI.ge.xD(j).and.( xI.lt.xD(j+1)) ) then
  !            yIprime =   3.d0*coef(a_,j)*( xI-xD(j) )**2 &
  !                 &    + 2.d0*coef(b_,j)*( xI-xD(j) ) &
  !                 &    +      coef(c_,j)
  !            exit
  !         endif
  !      enddo
       
  !   else
  !      yIprime =   3.d0*coef(a_,LD-1)*( xD(LD)-xD(LD-1) )**2 &
  !           &    + 2.d0*coef(b_,LD-1)*( xD(LD)-xD(LD-1) )    &
  !           &    +      coef(c_,LD-1)
  !   endif

  !   ! ------------------------------------------------------ !
  !   ! --- [2] return value                               --- !
  !   ! ------------------------------------------------------ !
  !   return__deriv_cSpline = yIprime
    
  !   return
  ! end function return__deriv_cSpline


