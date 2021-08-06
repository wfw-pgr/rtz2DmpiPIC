module myBiCGSTAB
contains
  ! source :: include boundary condition of x.
  ! --- if x=phi_b at (i=1) --> source(1,j) = phi_b
  subroutine BiCGSTABCrtXZ( x, source, dx, dy, LI, LJ )
    implicit none
    integer         , parameter     :: Ncomm = 5
    double precision, parameter     :: err   = 1.d-8
    integer         , intent(in)    :: LI, LJ
    double precision, intent(in)    :: dx, dy
    double precision, intent(in)    :: source(LI*LJ)
    double precision, intent(inout) :: x(LI*LJ)
    integer                         :: N, itermax
    integer         , allocatable   :: pt(:,:), Nelement(:)
    double precision, allocatable   ::  A(:,:), Minv(:)
    ! -- [1] parameter -- !
    N       = LI*LJ
    itermax = nint( N * 1.2d0 )
    allocate( A(N,Ncomm), pt(N,Ncomm), Nelement(N), Minv(N) )
    ! -- [2] Define XZ Poisson Operator ( Make Matrix ) -- !
    call XZLaplacian( A, pt, Nelement, dx, dy, LI, LJ, N, Minv )
    ! -- [3] BiCGStab solver -- !
    call bicgstab1d( A, source, x, pt, Nelement, N, itermax, err, Minv )
    return
  end subroutine BiCGSTABCrtXZ

  ! --- Cylindrical BiCGSolver --- #
  subroutine BiCGSTABCyl2D( x, source, dr, dz, N, M, rmin )
    implicit none
    integer         , parameter     :: Ncomm = 5
    double precision, parameter     :: err   = 1.d-8
    integer         , intent(in)    :: N, M
    double precision, intent(in)    :: dr, dz, rmin
    double precision, intent(in)    :: source((N+1)*M)
    double precision, intent(inout) :: x((N+1)*M)
    integer         , allocatable   :: pt(:,:), Nelement(:)
    double precision, allocatable   ::  A(:,:), Minv(:)
    integer                         :: NxM, itermax
    
    ! -- [1] parameter -- !
    NxM     = (N+1)*M
    itermax = nint( NxM * 1.2d0 )
    allocate( A(NxM,Ncomm), pt(NxM,Ncomm), Nelement(NxM), Minv(NxM)  )
    ! -- [2] Define RZ poisson operater ( Make Matrix ) -- !
    call RZLaplacian( A, pt, Nelement, dr, dz, N, M, NxM, rmin, Minv )
    ! -- [3] call BiCGStab solver -- !
    call bicgstab1d( A, source, x, pt, Nelement, NxM, itermax, err, Minv )
    
    return
  end subroutine BiCGSTABCyl2D


  subroutine XZLaplacian( A, pt, Nelement, dx, dy, LI, LJ, N, Minv )
    implicit none
    integer         , parameter   :: Ncomm = 5
    integer         , intent(in)  :: LI, LJ, N
    double precision, intent(in)  :: dx, dy
    integer         , intent(out) :: pt(N,Ncomm), Nelement(N)
    double precision, intent(out) ::  A(N,Ncomm), Minv(N)
    integer                       :: i, j, k
    double precision              :: coef(5)
    
    write(6,*)
    write(6,*) ' [making matrix] Amat is ', n, ' x ', n, 'matrix. [making matrix] ' 
    write(6,*)
    ! --- [1] Preparation --- !
    !  -- [1-1] Coefficients -- !
    coef(1) = 1.d0 / dy**2
    coef(2) = 1.d0 / dx**2
    coef(3) = - 2.d0 * ( 1.d0 / dx**2 + 1.d0 / dy**2 )
    coef(4) = 1.d0 / dx**2
    coef(5) = 1.d0 / dy**2
    ! -- [1-2] Initialize CRS variables -- !
    A(:,:)  = 0.d0
    pt(:,:) = 0
    k       = 1
    Minv(:) = 1.d0
    
    ! --- [2] Matrix Definition --- !
    !  -- [2-1] j=1 Boundary -- !
    do i=1, LI
       A(k,1)      = 1.d0
       pt(k,1)     = k
       Nelement(k) = 1
       k = k+1
    enddo
    !  -- [2-2] Main Region (j=2~LJ-1) -- !
    do j=2, LJ-1
       do i=1, LI
          if      ( i.eq. 1 ) then
             ! - i= 1 ( Dirichlet B.C. ) - !
             A(k,1)      =  1.d0
             pt(k,1)     = k
             Nelement(k) = 1
          else if ( i.eq.LI ) then
             ! - i=LI ( Dirichlet B.C. ) - !
             A(k,1)      =  1.d0
             pt(k,1)     = k
             Nelement(k) = 1
          else
             ! - Laplacian - !
             A(k,1)      = coef(1)
             A(k,2)      = coef(2)
             A(k,3)      = coef(3)
             A(k,4)      = coef(4)
             A(k,5)      = coef(5)
             pt(k,1)     = k-LI
             pt(k,2)     = k-1
             pt(k,3)     = k
             pt(k,4)     = k+1
             pt(k,5)     = k+LI
             Nelement(k) = 5
          endif
          k = k+1
       enddo
    enddo
    !  -- [2-3] j=LJ Boundary -- !
    do i=1, LI
       A(k,1) = 1.d0
       pt(k,1) = k
       Nelement(k) = 1
       k = k+1
    enddo
    !  -- END -- !
    if ( N.ne.(k-1) ) write(6,*) ' [WARNING] N and k is incompatable [WARNING] '
    return
  end subroutine XZLaplacian
  

  subroutine RZLaplacian( A, pt, Nelement, dr, dz, N, M, NxM, rmin, Minv )
    implicit none
    integer         , parameter   :: Ncomm = 5
    integer         , intent(in)  :: N, M, NxM
    double precision, intent(in)  :: dr, dz, rmin
    integer         , intent(out) :: pt(NxM,Ncomm), Nelement(NxM)
    double precision, intent(out) ::  A(NxM,Ncomm), Minv(NxM)
    character(9)                  :: zbc = 'Dirichlet'
    integer                       :: i, j, k, kR, kL
    double precision              :: rcoef, coef(5), coef_base(5), rsqinv, diaginv

    !  --- [1] preparation --- !
    write(6,*) ' [making matrix] A is ', NxM, ' x ', NxM, 'matrix. [making matrix] ' 
    !    --    coef    --    !
    coef_base(1) = 1.d0 / dr**2
    coef(2)      = 1.d0 / dz**2
    coef_base(3) = - 2.d0 * ( 1.d0 / dz**2 + 1.d0 / dr**2 )
    coef(4)      = 1.d0 / dz**2
    coef_base(5) = 1.d0 / dr**2
    !    -- initialize --    !
    A(:,:)  = 0.d0
    pt(:,:) = 0
    k       = 1
    !  -----------------------  !
    !  --- [2] Make Matrix ---  !
    !  -----------------------  !
    !    -- j=1 boundary -- !
    if ( rmin .eq. 0.d0 ) then
       coef(3) = coef_base(3) - 1.d0 / ( + 0.5d0 *dr )**2
       coef(5) = 2.d0 * coef_base(5)
       do i=1, M
          A(k,1)       = coef(2)
          A(k,2)       = coef(3)
          A(k,3)       = coef(4)
          A(k,4)       = coef(5)
          pt(k,1)      = k-1
          pt(k,2)      = k
          pt(k,3)      = k+1
          pt(k,4)      = k+M
          Nelement(k)  = 4
          Minv(k)      = 1.d0 / coef(3)
          k = k+1
       enddo
    else if ( rmin.gt.0.d0 ) then
       do i=1, M
          A(k,1)       = 1.d0
          pt(k,1)      = k
          Nelement(k)  = 1
          Minv(k)      = 1.d0
          k            = k+1
       enddo
    else
       stop ' [ERROR] Select rmin > 0 :: illigal rmin value. [ERROR] '
    endif
    !    -- [2-2] Main Region ( j=2~N ) -- !
    do j=2, N
       ! -- Coef(1,3,5) Definition -- !
       if ( rmin .eq. 0.d0 ) then
          rcoef  = 1.d0 / ( 2.d0 * dr * (  ( dble(j)-0.5d0 )*dr ) )
          rsqinv = 1.d0 / ( ( dble(j)-0.5d0 )*dr )**2
       else
          rcoef  = 1.d0 / ( 2.d0 * dr * (  dble(j-1)*dr + rmin   ) )
          rsqinv = 1.d0 / ( dble(j-1)*dr + rmin )**2
       endif
       coef(1) = coef_base(1) - rcoef
       coef(3) = coef_base(3) - rsqinv
       coef(5) = coef_base(5) + rcoef
       diaginv = 1.d0 / coef(3)
       ! -- Make Matrix -- !
       do i=1,M
          A(k,1)      = coef(1)
          A(k,2)      = coef(2)
          A(k,3)      = coef(3)
          A(k,4)      = coef(4)
          A(k,5)      = coef(5)
          pt(k,1)     = k-M
          pt(k,2)     = k-1
          pt(k,3)     = k
          pt(k,4)     = k+1
          pt(k,5)     = k+M
          Nelement(k) = 5
          Minv(k)     = diaginv
          k = k+1          
       enddo
    enddo
    !    -- [2-3] j=N+1 Boundary -- !
    do i=1, M
       A(k,1)      = 1.d0
       pt(k,1)     = k
       Nelement(k) = 1
       Minv(k)     = 1.d0
       k           = k+1
    enddo

    !  ---  [3] z-boundary condition  --- !
    if ( zbc .eq. 'Dirichlet' ) then
       write(6,*) '  --- Boundary Condition is Dirichlet condition ---  '
       kL = 1
       kR = M
       do j=1, N+1
          ! --  Left Edge -- !
          A(kL,1)      = 1.d0
          pt(kL,1)     = kL
          Nelement(kL) = 1
          Minv(kL)     = 1.d0
          kL           = kL + M
          ! -- Right Edge -- !
          A(kR,1)      = 1.d0
          pt(kR,1)     = kR
          Nelement(kR) = 1
          Minv(kR)     = 1.d0
          kR           = kR + M
       enddo
    endif
    if ( NxM .ne.(k-1) ) write(6,*) ' [WARNING] N and k is incompatable [WARNING] '
    if ( rmin.eq.0.d0  ) write(6,*) '  --- RZ-Coordinate Laplacian with    r=0 axis --- '
    if ( rmin.gt.0.d0  ) write(6,*) '  --- RZ-Coordinate Laplacian without r=0 axis --- '
          
    return
  end subroutine RZLaplacian
  

  function MatrixMultiply( Amat, x, Apt, ANelement, N, Nitem )
    implicit none
    integer         , intent(in) :: N, Nitem
    integer         , intent(in) :: ANelement(N), Apt(N,Nitem)
    double precision, intent(in) :: Amat(N,Nitem), x(N)
    integer                      :: i, j 
    double precision             :: MatrixMultiply(N), sumi

    !$omp parallel default(none) shared(N,Amat,ANelement,Apt,x,MatrixMultiply)
    !$omp do private(i,j,sumi)
    do i=1, N
       sumi = 0.d0
       do j=1, ANelement(i)
           sumi = sumi + Amat(i,j) * x( Apt(i,j) )
       enddo
       MatrixMultiply(i) = sumi
    enddo
    !$omp end do
    !$omp end parallel
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
    double precision             :: myDotProduct
    double precision, intent(in) :: r1(:), r2(:)
    
    N = ubound( r1, 1 )
    ! if ( N .ne. ubound(r2,1) ) stop ' [ERROR] vector size is different ( myDotProduct ) [ERROR] '
    myDotProduct = 0.d0
    !$omp parallel default(none) &
    !$omp shared(N,r1,r2,myDotProduct) private(i)
    !$omp do reduction(+:myDotProduct)
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
  

  subroutine bicgstab1d( A, b, x, pt, Nelement, NxM, itermax, err, Minv )
    implicit none
    integer         , parameter   :: Ncomm = 5
    integer         , intent(in)  :: NxM, itermax
    integer         , intent(in)  :: Nelement(NxM), pt(NxM,Ncomm)
    double precision, intent(in)  :: A(NxM,Ncomm), b(NxM), err, Minv(NxM)
    double precision, intent(out) :: x(NxM)
    integer                       :: iter, i
    double precision              :: alpha, beta, c1, c2, c3, rr
    double precision              :: r(NxM), r0(NxM), p(NxM), p0(NxM), y(NxM), e(NxM), v(NxM), z(NxM), u(NxM)
    double precision              :: sec0, sec1, sec2, sec3, sec4, sec5
    double precision              :: sec6, sec7, sec8, sec9, sec10, sec11
    double precision              :: sec12, sec13, sec14, sec15
    double precision              :: time_substitute, time_matmul, time_dotproduct
    ! 1st step
    x(:)  = 0.d0
    r(:)  = b(:) - MatrixMultiply( A, x, pt, Nelement, NxM, Ncomm )
    r0(:) = r(:)
    p(:)  = r(:)
    p0(:) = r0(:)
    c1 = myDotProduct( r0, r0 )

    ! iteration step k-th
    do iter=1, itermax
       z(:)  = DiagonalMultiply( Minv, p, NxM )
       y(:)  = MatrixMultiply( A, z, pt, Nelement, NxM, Ncomm )       
       c2    = myDotProduct( r0, y )
       alpha = c1 / c2
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

       r     = daxpy( -c3, v, e, NxM )
       rr    = myDotProduct( r, r )
       if ( rr .lt. err ) exit
       c1    = myDotProduct( r0, r )
       beta  = c1 / ( c2 * c3 )

       ! p(:)  = r(:) + beta * ( p(:) - c3 * y(:) )
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
  
end module MyBiCGSTAB

