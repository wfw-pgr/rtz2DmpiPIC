module myBiCGSTAB

  implicit none
  integer         , parameter   :: Ncomm = 5
  double precision, parameter   :: err   = 1.d-8

  integer                       :: NxM, itermax
  integer         , allocatable :: pt(:,:), Nelement(:)
  double precision, allocatable :: A(:,:) , Minv(:)

  integer                       :: MakeMatrix = 0
  
contains

  ! source :: include boundary condition of x.
  ! --- if x=phi_b at (i=1) --> source(1,j) = phi_b


  

  subroutine BiCGSTABCyl2D( x, source, dz, dr, N, M, rmin )
    
    
    implicit none
    integer         , intent(in)    :: N, M
    double precision, intent(in)    :: dr, dz, rmin
    double precision, intent(in)    :: source(N,M)
    double precision, intent(inout) :: x(N,M)
    

    if ( MakeMatrix .eq. 0 ) then
       NxM        = N*M
       itermax    = nint( NxM * 1.2d0 )
       allocate( A(NxM,Ncomm), pt(NxM,Ncomm), Nelement(NxM), Minv(NxM) )
       call GSLaplacian( dz, dr, N, M, rmin )
       MakeMatrix = 1
    endif
    
    call bicgstab1d( source, x )

    
    return
  end subroutine BiCGSTABCyl2D



  
  
  
  subroutine GSLaplacian( dz, dr, N, M, rmin )
    
    
    implicit none
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: dz, dr, rmin
    
    !                                    ---   dirichlet or periodicZ   ---
    character(9)     :: zbc = 'Dirichlet'
    ! character(9)     :: zbc = 'symetricZ'
    integer          :: i, j, k, kR, kL
    double precision :: rcoef1, rcoef2, coef(5), drsqinv, dzsqinv, diaginv


    !  --- [1] preparation --- !
    
    ! write(6,*) ' [making matrix] A is ', NxM, ' x ', NxM, 'matrix. [making matrix] ' 
    
    !    --    coef    --      !
    drsqinv      = 1.d0 / dr**2
    dzsqinv      = 1.d0 / dz**2
    coef(2)      = dzsqinv
    coef(4)      = dzsqinv
    
    !    -- initialize --      !
    A(:,:)  = 0.d0
    pt(:,:) = 0
    k = 1
    
    !  --- [2] Make Matrix --- !

    !    -- j=1 boundary
    do i=1, N
       A(k,1)       = 1.d0
       pt(k,1)      = k
       A(k,2)       = 1.d0
       pt(k,2)      = k+N
       Nelement(k)  = 2
       Minv(k)      = 1.d0
       k            = k+1
    enddo

    !    -- j = 2~M-1 region
    do j=2, M-1

       rcoef1         = (  dble(j-1)*dr + rmin  ) / ( ( dble(j-1)+0.5d0 )*dr + rmin )
       rcoef2         = (  dble(j-1)*dr + rmin  ) / ( ( dble(j-1)-0.5d0 )*dr + rmin )

       coef(1)        =   drsqinv * rcoef2
       coef(3)        = - drsqinv * ( rcoef1 + rcoef2 ) - 2.d0 * dzsqinv
       coef(5)        =   drsqinv * rcoef1
       diaginv        = 1.d0 / coef(3)
       
       do i=1,N
          A(k,1)      = coef(1)
          A(k,2)      = coef(2)
          A(k,3)      = coef(3)
          A(k,4)      = coef(4)
          A(k,5)      = coef(5)
          pt(k,1)     = k-N
          pt(k,2)     = k-1
          pt(k,3)     = k
          pt(k,4)     = k+1
          pt(k,5)     = k+N
          Nelement(k) = 5
          Minv(k)     = diaginv

          k = k+1          
       enddo

    enddo

    !    -- j = N+1  boundary
    do i=1, N
       A(k,1)      = 1.d0
       pt(k,1)     = k
       A(k,2)      = 1.d0
       pt(k,2)     = k-N
       Nelement(k) = 2
       Minv(k)     = 1.d0
       k           = k+1
    enddo


    !  ---  [2] z-boundary condition  --- !

    if ( zbc .eq. 'Dirichlet' ) then
       ! write(6,*) '  --- Boundary Condition is Dirichlet condition ---  '
       kL = 1
       kR = N
       do j=1, M
          
          ! --  left edge -- !
          A(kL,1)      = 1.d0
          pt(kL,1)     = kL
          A(kL,2)      = 1.d0
          pt(kL,2)     = kL+1
          Nelement(kL) = 2
          Minv(kL)     = 1.d0 
          kL           = kL + N
          
          ! -- right edge -- !
          A(kR,1)      = 1.d0
          pt(kR,1)     = kR
          A(kR,2)      = 1.d0
          pt(kR,2)     = kR+1
          Nelement(kR) = 2
          Minv(kR)     = 1.d0
          kR           = kR + N

       enddo

    else if ( zbc .eq. 'periodicZ' ) then
       ! write(6,*) '  --- Boundary Condition is Periodic condition ---  '
       kL = 1
       kR = N
       do j=1, M
          
          ! --  left edge -- !
          A(kL,1)      = + 1.d0
          A(kL,2)      = - 1.d0
          pt(kL,1)     = kL
          pt(kL,2)     = kL + (N-1)
          Nelement(kL) = 2
          kL           = kL + N
          
          ! -- right edge -- !
          A(kR,1)      = + 1.d0
          A(kR,2)      = - 1.d0
          pt(kR,1)     = kR
          pt(kR,2)     = kR - (N-1)
          Nelement(kR) = 2
          kR           = kR + N

       enddo

     
    else if ( zbc .eq. 'symetricZ' ) then
       ! write(6,*) '  --- Boundary Condition is Dirichlet condition ---  '
       kL = 1
       kR = N
       do j=1, M
          
          ! --  left edge -- !
          A(kL,1)      = +1.d0
          pt(kL,1)     = kL
          A(kL,2)      = -1.d0
          pt(kL,2)     = kL+1
          Nelement(kL) = 2
          kL           = kL + N
          
          ! -- right edge -- !
          A(kR,1)      = 1.d0
          pt(kR,1)     = kR
          Nelement(kR) = 1
          kR           = kR + N

       enddo

       Nelement(        1) = 1
       Nelement(N*(M-1)+1) = 1
       
    endif

    if ( NxM .ne.(k-1) ) write(6,*) ' [WARNING] N and k is incompatable [WARNING] '
        
    return
  end subroutine GSLaplacian
  




  
  function MatrixMultiply( Amat, x, Apt, ANelement, N, Nitem )
    

    implicit none

    integer         , intent(in) :: N, Nitem
    integer         , intent(in) :: ANelement(N), Apt(N,Nitem)
    double precision, intent(in) :: Amat(N,Nitem), x(N)

    integer                      :: i, j 
    double precision             :: MatrixMultiply(N)

    
    MatrixMultiply(:) = 0.d0

    do i=1, N
       do j=1, ANelement(i)
          MatrixMultiply(i) = MatrixMultiply(i) + Amat(i,j) * x( Apt(i,j) )
       enddo
    enddo

    return
  end function MatrixMultiply



  function DiagonalMultiply( u, v, N )
    
    implicit none
    integer         , intent(in) :: N
    double precision, intent(in) :: u(N), v(N)
    double precision             :: DiagonalMultiply(N)
    integer                      :: i

    !$omp parallel default(none) shared(N,DiagonalMultiply,u,v)
    !$omp do private(i)
    do i=1, N
       DiagonalMultiply(i) = u(i) * v(i)
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end function DiagonalMultiply


  
  function myDotProduct( r1, r2 )

    implicit none
    integer                      :: i, N
    double precision, intent(in) :: r1(:), r2(:)
    double precision             :: myDotProduct

    
    N = ubound( r1, 1 )
    if ( N .ne. ubound(r2,1) ) stop ' [ERROR] vector size is different ( myDotProduct ) [ERROR] '

    myDotProduct = 0.d0
    !$omp parallel default(none) shared(N,r1,r2,myDotProduct)
    !$omp do private(i) reduction(+:myDotProduct)
    do i=1, N
       myDotProduct = myDotProduct + r1(i) * r2(i)
    enddo
    !$omp end do
    !$omp end parallel

    return
  end function myDotProduct


  

  function daxpy( acoef, xvec, yvec, N )

    implicit none

    integer         , intent(in)    :: N
    double precision, intent(in)    :: acoef
    double precision, intent(in)    :: xvec(N), yvec(N)
    integer                         :: i

    double precision                :: daxpy(N)

    
    !$omp parallel default(none) shared(acoef,xvec,yvec,daxpy,N)
    !$omp do private(i)
    do i=1, N
       daxpy(i) = acoef * xvec(i) + yvec(i)       
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end function daxpy


  

  subroutine bicgstab1d( b, x )

    implicit none
    double precision, intent(in)  :: b(NxM)
    double precision, intent(out) :: x(NxM)

    integer          :: iter, i
    double precision :: alpha, beta, c1, c2, c3, rr
    double precision :: r(NxM), r0(NxM), p(NxM), p0(NxM), y(NxM), z(NxM), e(NxM), v(NxM), u(NxM)

      
    ! 1st step
    x(:)  = 0.d0
    r(:)  = b(:) - MatrixMultiply( A, x, pt, Nelement, NxM, Ncomm )
    r0(:) = r(:)
    p(:)  = r(:)
    p0(:) = r0(:)
    c1    = myDotProduct( r0, r0 )
    

    ! iteration step k-th
    do iter=1, itermax

       z(:)  = DiagonalMultiply( Minv, p, NxM )
       y(:)  = MatrixMultiply( A, z, pt, Nelement, NxM, Ncomm )
       c2    = myDotProduct( r0, y )
       alpha = c1 / c2

       ! e(:)  = r(:) - alpha * y(:)
       e(:)  = daxpy( -alpha, y, r, NxM )
       u(:)  = DiagonalMultiply( Minv, e, NxM )
       v(:)  = MatrixMultiply( A, u, pt, Nelement, NxM, Ncomm )
       
       c3    = myDotProduct( e, v ) / myDotProduct( v, v )
       
       ! x(:)  = x(:) + alpha * p(:) + c3 * e(:)
       !$omp parallel default(none) shared(x,alpha,z,c3,u,NxM)
       !$omp do private(i)
       do i=1, NxM
          x(i) = x(i) + alpha*z(i) + c3*u(i)
       enddo
       !$omp end do
       !$omp end parallel
       
       r(:)  = daxpy( -c3, v, e, NxM )
       ! r(:)  = e(:) - c3 * v(:)
       rr    = myDotProduct( r, r )
       if ( rr .lt. err ) exit
       c1    = myDotProduct( r0, r )
       beta  = c1 / ( c2 * c3 )

       
       ! p(:)  = r(:) + beta * ( p(:) - c3 * y(:) )
       ! y(:)  = MatrixMultiply( A, p, pt, Nelement, NxM, Ncomm )
       !$omp parallel default(none) shared(p,r,beta,c3,y,NxM)
       !$omp do private(i)
       do i=1, NxM
          p(i) = r(i) + beta * ( p(i) - c3 * y(i) )
       enddo
       !$omp end do
       !$omp end parallel

    enddo
    
    write(6,*) ' Bi-CGSTAB ended... in ', iter, '/',itermax
    write(6,*) ' Error === ', rr, '/' ,err
    
    return
  end subroutine bicgstab1d
  

  

end module myBiCGSTAB

