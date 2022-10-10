module newtonRapMod
  use cubSplineMod
  implicit none
  logical, parameter :: verbose  = .false.
contains

  ! ====================================================== !
  ! === solve__linearNewtonRaphson1D                   === !
  ! ====================================================== !
  subroutine solve__linear_NewtonRaphson1D( fx, rix, flag, find, converge_factor_in )
    implicit none
    double precision, intent(in)    :: fx(:) ! function to be searched. fx(rix).
    double precision, intent(inout) :: rix   ! initial / return (real) index of extremum value.
    integer         , intent(inout) :: flag  ! flag to return success or fail
    character(3)    , intent(in)    :: find
    double precision, intent(in)    :: converge_factor_in
    double precision                :: converge_factor
    double precision, parameter     :: converge_factor_def = 1.e-4
    double precision, parameter     :: iterMax_factor      = 100.d0
    double precision, parameter     :: stepMax_factor      = 0.1d0
    double precision, parameter     :: smallstep_factor    = 0.05d0
    integer                         :: ix, Lx, iter, iterMax
    double precision                :: fk, dfdx, alpha, delta_x, stepMax
    double precision                :: rix_old, fk_old, dfMax, sign
    double precision                :: converge_x, converge_f, threshold_x, threshold_f

    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !
    if ( trim(find) == "min" ) then
       sign = - 1.d0
    else if ( trim(find) == "max" ) then
       sign = + 1.d0
    else
       write(6,*) "[solve__linear_NewtonRaphson1D] unknwon find keyword....[ERROR]"
       write(6,*) "[solve__linear_NewtonRaphson1D]   find = [ min,max ]....[ERROR]"
       stop
    end if
    if ( converge_factor_in == 0.d0 ) then
       converge_factor = converge_factor_def
    else
       converge_factor = converge_factor_in
    endif
    Lx          = size(    fx )
    dfMax       = abs( maxval( fx ) - minval( fx ) ) 
    stepMax     = dble(    Lx *  stepMax_factor )
    iterMax     = nint(    Lx *  iterMax_factor )
    threshold_x = dble(    Lx * converge_factor )
    threshold_f = dble( dfMax * converge_factor )
    flag        = 0
    ix          = ceiling( rix ) - 1
    if ( ( ix <= 0 ).or.( ix >= Lx+1 ) ) then
       write(6,*)
       write(6,*) "[solve__linear_NewtonRaphson1D]  rix ( 0 < rix < Lx ) is NOT given.... [CAUTION]"
       write(6,*) "[solve__linear_NewtonRaphson1D]  rix is set to 0.5*Lx                  [CAUTION]"
       write(6,*)
       rix = 0.1d0 * Lx
    endif
    rix_old     = 0.d0
    fk_old      = 0.d0
    
    ! ------------------------------------------------------ !
    ! --- [2] Main Loop                                  --- !
    ! ------------------------------------------------------ !
    do iter=1, itermax

       ! ------------------------------------------------------ !
       ! --- [2-1] calculate Derivative                     --- !
       ! ------------------------------------------------------ !
       ix     = ceiling( rix ) - 1
       alpha  = rix - dble( ix )
       fk     = alpha * fx(ix+1) + ( 1.d0 - alpha ) * fx(ix)
       dfdx   = fx(ix+1) - fx(ix)

       ! ------------------------------------------------------ !
       ! --- [2-2] calculate delta_x :: small move step     --- !
       ! ------------------------------------------------------ !
       if ( dfdx == 0.d0  ) then
          write(6,*)
          write(6,*) "[solve__linear_NewtonRaphson1D] dfdx == 0.d0  ==> converged at flat position.... [CAUTION]"
          write(6,*) "[solve__linear_NewtonRaphson1D]   failed to find extremum value....              [CAUTION]"
          write(6,*)
          flag = - 1
          exit
       else
          delta_x = sign * smallstep_factor * fk / dfdx
       endif
       do while ( abs(delta_x) >= stepMax )
          delta_x = smallstep_factor * delta_x
       enddo
       
       ! ------------------------------------------------------ !
       ! --- [2-3] step forward                             --- !
       ! ------------------------------------------------------ !
       rix    = rix + delta_x

       ! ------------------------------------------------------ !
       ! --- [2-4] check convergence                        --- !
       ! ------------------------------------------------------ !
       converge_x = sqrt( ( rix - rix_old )**2 )
       converge_f = sqrt( ( fk  -  fk_old )**2 )
       rix_old    = rix
       fk_old     = fk
       if ( iterMax >= 2 ) then
          if ( ( converge_x <= threshold_x ).and.( converge_f <= threshold_f ) ) then
             flag = 1
             exit
          endif
       endif
       if ( ( rix.le.0.d0 ).or.( rix.ge.dble(Lx) ) ) then
          flag = -1
          write(6,*) "[solve__linear_NewtonRaphson1D]  out of defined region of given function fx.... [CAUTION]"
          write(6,*) "[solve__linear_NewtonRaphson1D]    function fx seems to be monotonous.....      [CAUTION]"
          write(6,*) "[solve__linear_NewtonRaphson1D]    failed to find extremum value....            [CAUTION]"
          write(6,*) "[solve__linear_NewtonRaphson1D]    no solution was found                        [CAUTION]"
          exit
       end if
       
       ! ------------------------------------------------------ !
       ! --- [2-5] debug write out                          --- !
       ! ------------------------------------------------------ !
       if ( verbose ) then
          write(6,*) iter, rix, fk, delta_x, converge_x, converge_f
       endif
       
    enddo
    return
  end subroutine solve__linear_NewtonRaphson1D


  ! ====================================================== !
  ! === solve__cSpline_NewtonRaphson1D                 === !
  ! ====================================================== !
  subroutine solve__cSpline_NewtonRaphson1D( xD, fD, xk, fk, flag, find, converge_factor_in )
    implicit none
    double precision, intent(in)    :: xD(:), fD(:)        ! -- function to be searched
    double precision, intent(inout) :: xk, fk              ! -- initial/return (x,f) values
    integer         , intent(inout) :: flag                ! -- convergence state flag
    character(3)    , intent(in)    :: find                ! -- [min,max]
    double precision, intent(in)    :: converge_factor_in  ! -- convergence value factor (0:auto)
    double precision                :: converge_factor
    double precision, parameter     :: converge_factor_def = 1.e-4
    double precision, parameter     :: iterMax_factor      = 100.d0
    double precision, parameter     :: stepMax_factor      = 0.1d0
    double precision, parameter     :: smallstep_factor    = 0.05d0
    integer                         :: Lx, iter, iterMax
    double precision                :: dfdx, delta_x, stepMax
    double precision                :: xk_old, fk_old, dxMax, dfMax, sign, xI(1), fI(1), dfdxI(1)
    double precision                :: converge_x, converge_f, threshold_x, threshold_f


    ! ------------------------------------------------------ !
    ! --- [1] preparation                                --- !
    ! ------------------------------------------------------ !
    if ( trim(find) == "min" ) then
       sign = - 1.d0
    else if ( trim(find) == "max" ) then
       sign = + 1.d0
    else
       write(6,*) "[solve__cSpline_NewtonRaphson1D @newtonRapMod] unknwon find keyword....[ERROR]"
       write(6,*) "[solve__cSpline_NewtonRaphson1D @newtonRapMod]   find = [ min,max ]....[ERROR]"
       stop
    end if
    if ( converge_factor_in == 0.d0 ) then
       converge_factor = converge_factor_def
    else
       converge_factor = converge_factor_in
    endif
    Lx          = size( xD )
    dxMax       =  abs( xD(2) - xD(1) )
    dfMax       =  abs( maxval( fD ) - minval( fD ) ) 
    iterMax     = nint(    Lx *  iterMax_factor )
    stepMax     = dble( dxMax *  stepMax_factor )
    threshold_x = dble( dxMax * converge_factor )
    threshold_f = dble( dfMax * converge_factor )
    flag        = 0
    xk_old      = 0.d0
    fk_old      = 0.d0
    if ( size(fD).ne.Lx ) then
       write(6,*)
       write(6,*) "[solve__cSpline_NewtonRaphson1D] incompatible xD & fD....[ERROR]"
       write(6,*) "[solve__cSpline_NewtonRaphson1D]    STOP             ....[ERROR]"
       write(6,*)
       stop
    endif
    
    ! ------------------------------------------------------ !
    ! --- [2] Main Loop                                  --- !
    ! ------------------------------------------------------ !
    do iter=1, itermax

       ! ------------------------------------------------------ !
       ! --- [2-1] calculate Derivative                     --- !
       ! ------------------------------------------------------ !
       xI(1) = xk
       call interpolate__cubicSpline1D( xD, fD, xI, fI, dfdxI )
       fk    =    fI(1)
       dfdx  = dfdxI(1)
       
       ! ------------------------------------------------------ !
       ! --- [2-2] calculate delta_x :: small move step     --- !
       ! ------------------------------------------------------ !
       if ( dfdx == 0.d0  ) then
          write(6,*)
          write(6,*) "[solve__cSpline_NewtonRaphson1D] dfdx == 0.d0  ==> converged at flat position.... [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson1D]   failed to find extremum value....              [CAUTION]"
          write(6,*)
          flag = - 1
          exit
       else
          delta_x = sign * smallstep_factor * fk / dfdx
       endif
       do while ( abs(delta_x) >= stepMax )
          delta_x = smallstep_factor * delta_x
       enddo
       
       ! ------------------------------------------------------ !
       ! --- [2-3] step forward                             --- !
       ! ------------------------------------------------------ !
       xk    = xk + delta_x

       ! ------------------------------------------------------ !
       ! --- [2-4] check convergence                        --- !
       ! ------------------------------------------------------ !
       converge_x = sqrt( ( xk  - xk_old )**2 )
       converge_f = sqrt( ( fk  - fk_old )**2 )
       xk_old     = xk
       fk_old     = fk
       if ( iterMax >= 2 ) then
          if ( ( converge_x <= threshold_x ).and.( converge_f <= threshold_f ) ) then
             flag = 1
             exit
          endif
       endif
       if ( ( xk <= xD(1) ).or.( xk >= xD(Lx) ) ) then
          flag = -1
          write(6,*) "[solve__cSpline_NewtonRaphson1D]  out of defined region of given function fx.... [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson1D]    function fx seems to be monotonous.....      [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson1D]    failed to find extremum value....            [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson1D]    no solution was found                        [CAUTION]"
          exit
       end if
       
       ! ------------------------------------------------------ !
       ! --- [2-5] debug write out                          --- !
       ! ------------------------------------------------------ !
       if ( verbose ) then
          write(6,*) iter, xk, fk, delta_x, converge_x, converge_f
       endif
       
    enddo

    return
  end subroutine solve__cSpline_NewtonRaphson1D


  ! ====================================================== !
  ! === solve__cSpline_NewtonRaphson2D                 === !
  ! ====================================================== !
  subroutine solve__cSpline_NewtonRaphson2D( xD, yD, fD, xk, yk, fk, flag, find, converge_factor_in )
    implicit none
    double precision, intent(in)    :: xD(:), yD(:)        ! -- function to be searched
    double precision, intent(in)    :: fD(:,:)             ! -- function to be searched
    double precision, intent(inout) :: xk, yk, fk          ! -- initial/return (x,y,f) values
    integer         , intent(inout) :: flag                ! -- convergence state flag
    character(3)    , intent(in)    :: find                ! -- [min,max]
    double precision, intent(in)    :: converge_factor_in  ! -- convergence value factor (0:auto)
    double precision                :: converge_factor
    double precision, parameter     :: converge_factor_def = 1.e-4
    double precision, parameter     :: iterMax_factor      = 100.d0
    double precision, parameter     :: stepMax_factor      = 0.1d0
    double precision, parameter     :: smallstep_factor    = 0.05d0
    integer                         :: i, j, Lx, Ly, iter, iterMax
    double precision                :: det, detInv
    double precision                :: delta_x, delta_y, stepMax
    double precision                :: xk_old, yk_old, fDxk_old, fDyk_old
    double precision                :: dxMax, dfMax, sign
    double precision                :: converge_xy , converge_fg
    double precision                :: threshold_xy, threshold_fg
    double precision, allocatable   :: fDx(:,:), fDy(:,:)
    double precision                :: fDxk, fDyk, fDxx, fDxy, fDyx, fDyy

    ! flag :: -1 => does not converged
    ! flag ::  0 => saddle point ( det < 0 )
    ! flag :: +1 => minimum      ( det > 0 & fxx > 0 )
    ! flag :: +2 => maximum      ( det > 0 & fxx < 0 )

    ! ------------------------------------------------------ !
    ! --- [1] Arguments                                  --- !
    ! ------------------------------------------------------ !
    if      ( trim(find) == "min" ) then
       sign = - 1.d0
    else if ( trim(find) == "max" ) then
       sign = + 1.d0
    else
       write(6,*) "[solve__cSpline_NewtonRaphson1D @newtonRapMod] unknwon find keyword....[ERROR]"
       write(6,*) "[solve__cSpline_NewtonRaphson1D @newtonRapMod]   find = [ min,max ]....[ERROR]"
       stop
    end if
    if ( converge_factor_in == 0.d0 ) then
       converge_factor = converge_factor_def
    else
       converge_factor = converge_factor_in
    endif
    Lx           = size( xD )
    Ly           = size( yD )
    dxMax        = sqrt( ( xD(2) - xD(1) )**2 + ( yD(2) - yD(1) )**2 )
    dfMax        = sqrt( ( maxval( fD ) - minval( fD ) )**2 )
    iterMax      = nint( sqrt( dble( Lx**2 + Ly**2 ) ) *  iterMax_factor )
    stepMax      = dxMax *  stepMax_factor
    threshold_xy = dxMax * converge_factor
    threshold_fg = dfMax * converge_factor
    flag         = -1
    xk_old       = 0.d0
    yk_old       = 0.d0
    fDxk_old     = 0.d0
    fDyk_old     = 0.d0

    ! ------------------------------------------------------ !
    ! --- [2] partial derivative                         --- !
    ! ------------------------------------------------------ !
    allocate( fDx(Lx,Ly), fDy(Lx,Ly) )
    do j=1, Ly
       do i=2, Lx-1
          fDx(i,j) = ( fD(i+1,j) - fD(i-1,j) ) / ( xD(i+1) - xD(i-1) )
       enddo
       fDx( 1,j) = ( fD( 2,j) - fD(   1,j) ) / ( xD( 2) - xD(   1) )
       fDx(Lx,j) = ( fD(Lx,j) - fD(Lx-1,j) ) / ( xD(Lx) - xD(Lx-1) )
    enddo
    do i=1, Lx
       do j=2, Ly-1
          fDy(i,j) = ( fD(i,j+1) - fD(i,j-1) ) / ( yD(j+1) - yD(j-1) )
       enddo
       fDy(i, 1) = ( fD(i, 2) - fD(i,   1) ) / ( yD( 2) - yD(   1) )
       fDy(i,Ly) = ( fD(i,Ly) - fD(i,Ly-1) ) / ( yD(Ly) - yD(Ly-1) )
    enddo

    ! ------------------------------------------------------ !
    ! --- [2] Main Loop                                  --- !
    ! ------------------------------------------------------ !
    do iter=1, iterMax

       ! ------------------------------------------------------ !
       ! --- [2-1] calculate Derivative                     --- !
       ! ------------------------------------------------------ !
       call interpolate__cubicSpline2D_point( xD, yD, fDx, xk, yk, fDxk, fDxx, fDxy )
       call interpolate__cubicSpline2D_point( xD, yD, fDy, xk, yk, fDyk, fDyx, fDyy )
       det = ( fDxx * fDyy ) - ( fDxy * fDyx )
       
       ! ------------------------------------------------------ !
       ! --- [2-2] calculate delta_x :: small move step     --- !
       ! ------------------------------------------------------ !
       if ( det == 0.d0  ) then
          write(6,*)
          write(6,*) "[solve__cSpline_NewtonRaphson1D] det == 0.d0  ==> Hessian is 0. Cannot Find. [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson1D]   failed to find extremum value....         [CAUTION]"
          write(6,*)
          flag = - 1
          exit
       else
          detInv  = 1.d0 / det
          delta_x = ( ( fDxy * fDyk ) - ( fDyy * fDxk ) ) * detInv
          delta_y = ( ( fDyx * fDxk ) - ( fDxx * fDyk ) ) * detInv
       endif
       delta_x = smallstep_factor * delta_x
       delta_y = smallstep_factor * delta_y
       do while ( abs(delta_x) >= stepMax )
          delta_x = smallstep_factor * delta_x
       enddo
       do while ( abs(delta_y) >= stepMax )
          delta_y = smallstep_factor * delta_y
       enddo
       
       ! ------------------------------------------------------ !
       ! --- [2-3] step forward                             --- !
       ! ------------------------------------------------------ !
       xk = xk + delta_x
       yk = yk + delta_y

       ! ------------------------------------------------------ !
       ! --- [2-4] check convergence                        --- !
       ! ------------------------------------------------------ !
       converge_xy = sqrt( (   xk  -   xk_old )**2 + (   yk  -   yk_old )**2 )
       converge_fg = sqrt( ( fDxk  - fDxk_old )**2 + ( fDyk  - fDyk_old )**2 )
       xk_old      = xk
       yk_old      = yk
       fDxk_old    = fDxk
       fDyk_old    = fDyk
       if ( iterMax >= 2 ) then
          if ( ( converge_xy <= threshold_xy ).and.( converge_fg <= threshold_fg ) ) then
             if ( det > 0.d0 ) then
                if      ( fDxx > 0.d0 ) then
                   flag = 1
                else if ( fDxx < 0.d0 ) then
                   flag = 2
                else
                   if ( fDyy > 0.d0 ) then
                      flag = 1
                   else if ( fDyy < 0.d0 ) then
                      flag = 2
                   else
                      write(6,*) "[solve__cSpline_NewtonRaphson2D] fDxx, fDyy == 0.d0  ==> Hessian is 0. Cannot Find. [ERROR]"
                      stop
                   endif
                endif
             else if ( det < 0.d0 ) then
                flag = 0
             else
                write(6,*) "[solve__cSpline_NewtonRaphson2D] det == 0.d0  ==> Hessian is 0. Cannot Find. [ERROR]"
             endif
             exit                
          endif
       endif
       if ( ( xk <= xD(1) ).or.( xk >= xD(Lx) ).or.( yk <= yD(1) ).or.( yk >= yD(Ly) ) ) then
          flag = -1
          write(6,*) "[solve__cSpline_NewtonRaphson2D] out of defined region of given function fx.... [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson2D]   function fx seems to be monotonous.....      [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson2D]   failed to find extremum value....            [CAUTION]"
          write(6,*) "[solve__cSpline_NewtonRaphson2D]   no solution was found                        [CAUTION]"
          exit
       end if
       
       ! ------------------------------------------------------ !
       ! --- [2-5] debug write out                          --- !
       ! ------------------------------------------------------ !
       if ( verbose ) then
          if ( mod(iter,20) == 1 ) then
             write(6,"(6x, a)") "# iter xk yk fDxk fDyk delta_x delta_y converge_xy converge_fg"
          endif
          write(6,"(i8,1x,8(e12.5,1x))") iter, xk, yk, fDxk, fDyk, delta_x, delta_y, &
               &                         converge_xy, converge_fg
       endif
       
    enddo

    return
  end subroutine solve__cSpline_NewtonRaphson2D

  
end module newtonRapMod
