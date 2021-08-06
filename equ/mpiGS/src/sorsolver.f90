module sorsolver

  use constants, only : Nr, Nz
  implicit none
  integer            :: flag
  double precision   :: psiold(Nz,Nr), dpq(Nr, Nz)
  double precision   :: c(Nr,Nz), d(Nr,Nz), e(Nr,Nz)
  double precision   :: p(Nr,Nz), q(Nr,Nz), s(Nr,Nz)
  double precision   :: res
  ! SOR parameters
  integer         , parameter :: poissonmax    = 100000
  double precision, parameter :: omega         = 1.8d0
  double precision, parameter :: threshold     = 1.0d-3
contains

  subroutine solv_poisson
    implicit none
    integer :: k
    ! --- [1] Initialize   --- !
    call initset
    call set_bdrc
    
    ! --- [2] Main Routine --- !
    k     = 1
    flag  = 0
    do while( ( k.lt.poissonmax ).and.( flag.eq.0 ) )
       !  -- [2-1] Calculate dpsi and step forward by SOR -- !
       call SOR
       !  -- [2-2] Boundary Condition set                 -- !
       call set_bdrc
       !  -- [2-3] Check Convergence                      -- !
       call checkconvg
       k = k + 1
    enddo
    write(6,'(" SOR, iteration :: ",i8,2x,"flag :: ",i1,2x,"residual ::",1x,e12.5)') k, flag, res
    
    return
  end subroutine solv_poisson

  
  subroutine initset

    use variables, only : r, z, dr, dz
    use variables, only : psi
    implicit none
    integer          :: i, j
    double precision :: rrplus, rrminus
    double precision :: drsqinv, dzsqinv

    !  --- [1]  Initialize  --- !
    drsqinv     = 1.d0 / dr**2
    dzsqinv     = 1.d0 / dz**2
    do j=1, Nr
       do i=1, Nz
          psi(i,j)    = 0.d0
          psiold(i,j) = 0.d0
       enddo
    enddo
    
    !  --- [2] calculate c, d, e, p, q, r ( see Lab. note ) --- !
    do j=2, Nr-1
       rrplus = r(j) / ( 0.5d0*( r(j+1) + r(j  ) ) )
       rrminus= r(j) / ( 0.5d0*( r(j  ) + r(j-1) ) )
       do i=2, Nz-1
          c(i,j)   =   drsqinv *   rrplus
          d(i,j)   = - drsqinv * ( rrplus + rrminus )
          e(i,j)   =   drsqinv *            rrminus
          p(i,j)   =          dzsqinv
          q(i,j)   = - 2.d0 * dzsqinv
          s(i,j)   =          dzsqinv
          dpq(i,j) = 1.d0 / ( d(i,j) + q(i,j) )
       enddo
    enddo
    
    return
  end subroutine initset

  
  subroutine SOR

    use constants, only : Nr, Nz
    use variables, only : psi, rhs
    implicit none
    integer            :: i, j
    double precision   :: psinew( 0:Nr,0:Nz ), dpsi( 0:Nr,0:Nz )
    double precision   :: seidel

    !  --- [1] Calculate next psi(k) from psi(k-1) --- !
    do j=2, Nr-1
       do i=2, Nz-1
          seidel   =   c(i,j)*psi(i,j+1) + e(i,j)*psi(i,j-1) &
               &     + p(i,j)*psi(i+1,j) + s(i,j)*psi(i-1,j)
          psi(i,j) = ( rhs(i,j) - seidel ) * dpq(i,j)
       enddo
    enddo
    
    !  --- [2] Calculate dpsi  &&  Step Forward by SOR --- !
    do j=1, Nr
       do i=1, Nz
          psi(i,j)  = omega * psi(i,j) + ( 1.d0 - omega ) * psiold(i,j)
       enddo
    enddo
    
    return
  end subroutine SOR

  
  subroutine checkconvg
    
    use variables, only : psi, rhs
    implicit none
    integer            :: i, j
    double precision   :: lhs

    !  --- [1] Calculate residual of the equations --- !
    res  = 0.d0
    do j=1, Nr
       do i=1, Nz
          lhs =    c(i,j)*psi(i,j+1) + d(i,j)*psi(i,j) + e(i,j)*psi(i,j-1) &
               & + p(i,j)*psi(i+1,j) + q(i,j)*psi(i,j) + s(i,j)*psi(i-1,j)
          
          res = res + ( lhs - rhs(i,j) )**2
       enddo
    enddo
    res    = sqrt( res / ( Nr*Nz ) )
    !  --- [2] Judge if convergence or NOT --- !
    if ( res .lt. threshold ) then
       flag = 1
    endif
    !  --- [3] Update psiOld --- !
    do j=1, Nr
       do i=1, Nz
          psiold(i,j) = psi(i,j)
       enddo
    enddo
    
    return
  end subroutine checkconvg

  
  subroutine set_bdrc

    use constants, only : Nr, Nz
    use variables, only : psi, BdrcLR, BdrcTB
    implicit none
    integer            :: i, j, m

    do j=1, Nr
       psi( 1,j) = BdrcLR(j,1)
       psi(Nz,j) = BdrcLR(j,2)
       ! psi( 1,j) = psi(2,j)
    enddo
    do i=1, Nz
       psi(i, 1) = BdrcTB(i,2)
       psi(i,Nr) = BdrcTB(i,1)
    enddo
    
    return
  end subroutine set_bdrc

end module sorsolver
