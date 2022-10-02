module PBiCGStabMod
  implicit none
  ! -- Constants        -- !
  integer         , parameter   :: Ncomm         = 5
  double precision, parameter   :: err           = 1.d-8
  integer         , parameter   :: NeibPEtot_Max = 2
  ! -- Global Variables -- !
  integer                       :: myRank, PEtot, CommLen, NeibPEtot
  integer                       :: NxM, NxMw, itermax, Nloc, Mloc
  integer                       :: NeibPE(NeibPEtot_Max)
  ! -- Allocatable Var. -- !
  integer         , allocatable ::       FromTo(:,:)
  integer         , allocatable ::           pt(:,:),     Nelement(:)
  double precision, allocatable ::            A(:,:),         Minv(:)
  integer         , allocatable :: export_Index(:)  , import_Index(:)
  integer         , allocatable :: export_Addrs(:)  , import_Addrs(:)
  integer         , allocatable :: request_Send(:)  , request_Recv(:)
  integer         , allocatable ::  status_Send(:,:),  status_Recv(:,:)
  double precision, allocatable ::         xloc(:,:),         sloc(:,:)

contains
  ! ---------------   USAGE   --------------------- !
  ! source :: include boundary condition of x.
  ! --- if x=phi_b at (i=1) --> source(1,j) = phi_b
  ! ----------------------------------------------- !
  
  subroutine BiCGSTABCyl2D_MPI( x, source, dz, dr, N, M, rmin )
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: N, M
    double precision, intent(in)    :: dr, dz, rmin
    double precision, intent(in)    :: source(N,M)
    double precision, intent(inout) :: x(N,M)
    integer                         :: ierr, iPE, Nitem_comm
    integer, save                   :: MakeMatrix = 0

    ! --- [1] Make Matrix if No Matrix --- !
    if ( MakeMatrix.eq.0 ) then
       ! -- [1-1] Define Variables for Comm. -- !
       call DefineCommunication( N, M )
       ! -- [1-2] Make Matrix -- !
       allocate( xloc(Nloc,Mloc), sloc(Nloc,Mloc) )
       allocate( A(NxM,Ncomm)   , pt(NxM,Ncomm),  &
            &    Nelement(NxM)  , Minv(NxM)       )
       call GSLaplacian( dz, dr, Nloc, Mloc, rmin )
       MakeMatrix = 1
    endif
    
    ! --- [2] call PBiCGSTAB_Engine --- !
    sloc(:,:) = source( FromTo(myRank,1):FromTo(myRank,2), 1:Mloc )
    call PBiCGSTAB_Engine_MPI( sloc, xloc )

    ! --- [3] BroadCast Solution    --- !
    do iPE=0, PEtot-1
       ! -- [3-1] Prepare iPE's Solution -- !
       if ( myRank.eq.iPE ) sloc = xloc
       Nitem_comm = ( FromTo(iPE,2)-FromTo(iPE,1)+1 ) * Mloc
       ! -- [3-2] BroadCast              -- !
       call MPI_Bcast( sloc(1,1), Nitem_comm    , MPI_DOUBLE_PRECISION, &
            &          iPE      , MPI_COMM_WORLD, ierr )
       ! -- [3-3] Save in x(From:To)     -- !
       x( FromTo(iPE,1):FromTo(iPE,2), 1:Mloc ) = sloc(:,:)
    enddo
    
    return
  end subroutine BiCGSTABCyl2D_MPI

  
  subroutine DefineCommunication( N, M )

    implicit none
    include 'mpif.h'
    integer, intent(in) :: N, M
    integer             :: j, jM, ierr, surplus, NlocW, iPE

    ! --- [1] Define Grid Number and Partition --- !
    !  -- [1-1]  PE infomation       --- !
    call MPI_COMM_RANK( MPI_COMM_WORLD, myRank, ierr )
    call MPI_COMM_SIZE( MPI_COMM_WORLD, PEtot , ierr )
    !  -- [1-2]  FromTo :: Partition --- !
    allocate( FromTo( 0:PEtot-1, 2 ) )
    Nloc    = N / PEtot
    surplus = N - Nloc * PEtot
    do iPE=0, PEtot-1
       FromTo(iPE  ,2) = Nloc
    enddo
    do iPE=1, surplus
       FromTo(iPE-1,2) = FromTo(iPE-1,2) + 1
    enddo
    FromTo(0,1) = 1
    do iPE=1, PEtot-1
       FromTo(iPE  ,1) = FromTo(iPE-1,2) + 1
       FromTo(iPE  ,2) = FromTo(iPE-1,2) + FromTo(iPE,2)
    enddo
    !  -- [1-3]  Grid Definition     --- !
    Nloc    = FromTo(myRank,2) - FromTo(myRank,1) + 1
    Mloc    = M
    NxM     = Nloc * Mloc
    NxMw    = Nloc * Mloc + 2*Mloc
    itermax = nint( NxM * 1.2d0 )
    CommLen = Mloc * NeibPEtot_Max

    ! --- [2] Allocate Variables for MPI Communication--- !
    allocate( export_Index(0:NeibPEtot_Max), import_Index(0:NeibPEtot_Max) )
    allocate( export_Addrs(CommLen),         import_Addrs(CommLen)         )
    allocate( request_Send(NeibPEtot_Max),   request_Recv(NeibPEtot_Max)   )
    allocate(  status_Send(MPI_STATUS_SIZE,NeibPEtot_Max), &
         &     status_Recv(MPI_STATUS_SIZE,NeibPEtot_Max)  )

    ! --- [3] Generate MPI Communication Table --- !
    !  -- [3-0] Exception for PEtot=1 -- !
    if ( PEtot .eq.1       ) then
       NeibPEtot = 0
       return
    endif
    !  -- [3-1] Define Index/Address for Comm. -- !
    export_Index = 0
    import_Index = 0
    if ( myRank.eq.0       ) then
       NeibPEtot           = 1
       NeibPE(1)           = 1
       export_Index(1)     = Mloc
       import_Index(1)     = Mloc
       do j=1, Mloc
          export_Addrs(j)  = Nloc* j
          import_Addrs(j)  = NxM + j
       enddo
    endif
    if ( myRank.eq.PEtot-1 ) then
       NeibPEtot           = 1
       NeibPE(1)           = PEtot-2
       export_Index(1)     = Mloc
       import_Index(1)     = Mloc
       do j=1, Mloc
          export_Addrs(j)  = Nloc*(j-1) + 1
          import_Addrs(j)  = NxM + j
       enddo
    endif
    if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
       NeibPEtot           = 2
       NeibPE(1)           = myRank-1
       NeibPE(2)           = myRank+1
       export_Index(1)     = Mloc
       export_Index(2)     = Mloc*2
       import_Index(1)     = Mloc
       import_Index(2)     = Mloc*2
       do j=1, Mloc
          jM               = Mloc+ j
          export_Addrs(j)  = Nloc*(j-1) + 1
          export_Addrs(jM) = Nloc*(j  )
          import_Addrs(j)  = NxM + j
          import_Addrs(jM) = NxM + jM
       enddo
    endif

    return
  end subroutine DefineCommunication

  
  subroutine GSLaplacian( dz, dr, N, M, rmin )
    
    implicit none
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: dz, dr, rmin
    integer                       :: i, j, k, kR, kL
    double precision              :: rcoef1, rcoef2, coef(5), drsqinv, dzsqinv, diaginv

    !  --- [1] preparation --- !
    !   -- coef       --   !
    drsqinv = 1.d0 / dr**2
    dzsqinv = 1.d0 / dz**2
    coef(2) = dzsqinv
    coef(4) = dzsqinv    
    !   -- Initialize --   !
    A(:,:)  = 0.d0
    pt(:,:) = 0
    k       = 1
    ! write(6,*) ' [making matrix] A is ', NxM, ' x ', NxM, 'matrix. [making matrix] '
    
    !  --- [2]  Make Matrix   --- !
    !   -- [2-1] j=1 Boundary --  !
    do i=1, N
       A(k,1)         = 1.d0
       pt(k,1)        = k
       Nelement(k)    = 1
       Minv(k)        = 1.d0
       k              = k+1
    enddo
    !   -- [2-2] j=1~M-1      --  !
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
          k           = k+1          
       enddo
    enddo
    !   -- [2-3] j=M Boundary --  !
    do i=1, N
       A(k,1)         = 1.d0
       pt(k,1)        = k
       Nelement(k)    = 1
       Minv(k)        = 1.d0
       k              = k+1
    enddo

    !  --- [3] z-boundary condition   --- !
    !   -- Left  B.C. from Neibour PE --  !
    if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
       kL = 1
       do j=1, M
          pt(kL,2)     = k
          kL           = kL + N
          k            = k  + 1
       enddo
       !   -- Right B.C. from Neibour PE --  !
       kR = N
       do j=1, M
          pt(kR,4)     = k
          kR           = kR + N
          k            = k  + 1
       enddo
    endif
    if ( myRank.eq.0 ) then
       !  -- Left  B.C. for Dirichlet --  !
       kL = 1
       do j=1, M
          A(kL,1)      = 1.d0
          pt(kL,1)     = kL
          Nelement(kL) = 1
          Minv(kL)     = 1.d0 
          kL           = kL + N
       enddo
       !  -- Right B.C. from Neibour PE --  !
       kR = N
       do j=1, M
          pt(kR,4)     = k
          kR           = kR + N
          k            = k  + 1
       enddo
    endif
    if ( myRank.eq.PEtot-1 ) then
       !  -- Left  B.C. from Neibour PE --  !
       kL = 1
       do j=1, M
          pt(kL,2)     = k
          kL           = kL + N
          k            = k  + 1
       enddo
       !  -- Right B.C. for Dirichlet --  !
       kR = N
       do j=1, M
          A(kR,1)      = 1.d0
          pt(kR,1)     = kR
          Nelement(kR) = 1
          Minv(kR)     = 1.d0 
          kR           = kR + N
       enddo
    endif
    ! if ( NxM.ne.(k-1) ) write(6,*) ' [WARNING] N and k is incompatable [WARNING] '
        
    return
  end subroutine GSLaplacian
  
  
  function MatrixMultiply( Amat, x, Apt, ANelement, N, M, Nitem )
    implicit none
    integer         , intent(in) :: N, M, Nitem
    integer         , intent(in) :: ANelement(N), Apt(N,Nitem)
    double precision, intent(in) :: Amat(N,Nitem), x(M)
    integer                      :: i, j 
    double precision             :: MatrixMultiply(N), sumi

    !$omp parallel default(none) &
    !$omp shared(N,Amat,ANelement,Apt,x,MatrixMultiply) private(i,j,sumi)
    !$omp do
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

    !$omp parallel default(none) &
    !$omp shared(N,DiagonalMultiply,u,v) private(i)
    !$omp do
    do i=1, N
       DiagonalMultiply(i) = u(i) * v(i)
    enddo
    !$omp end do
    !$omp end parallel
    return
  end function DiagonalMultiply

  
  function myDotProduct( r1, r2, N )
    implicit none
    include 'mpif.h'
    integer                      :: i, ierr
    integer         , intent(in) :: N
    double precision, intent(in) :: r1(N), r2(N)
    double precision             :: sum
    double precision             :: myDotProduct
    
    sum = 0.d0
    !$omp parallel default(none) &
    !$omp shared(N,r1,r2,sum) private(i) 
    !$omp do reduction(+:sum)
    do i=1, N
       sum = sum + r1(i)*r2(i)
    enddo
    !$omp end do
    !$omp end parallel
    call MPI_AllReduce( sum, myDotProduct, 1, MPI_DOUBLE_PRECISION, &
         &              MPI_SUM, MPI_COMM_WORLD, ierr )
    
    return
  end function myDotProduct


  function daxpy( acoef, xvec, yvec, N )
    implicit none
    integer         , intent(in)    :: N
    double precision, intent(in)    :: acoef
    double precision, intent(in)    :: xvec(N), yvec(N)
    integer                         :: i
    double precision                :: daxpy(N)
    
    !$omp parallel default(none) &
    !$omp shared(acoef,xvec,yvec,daxpy,N) private(i)
    !$omp do 
    do i=1, N
       daxpy(i) = acoef * xvec(i) + yvec(i)       
    enddo
    !$omp end do
    !$omp end parallel
    return
  end function daxpy


  subroutine BoundaryCommunicate( xvc )
    implicit none
    include 'mpif.h'
    double precision, intent(inout) :: xvc(NxMw)
    integer                         :: j, neib, ierr
    integer                         :: is, ir, len_s, len_r
    double precision                :: sendBuff(CommLen)
    double precision                :: recvBuff(CommLen)
    if ( neibPEtot.eq.0 ) return
    
    ! --- [1] Preparation --- !
    do neib=1, NeibPEtot
       do j=export_Index(neib-1)+1, export_Index(neib)
          sendBuff(j) = xvc(export_Addrs(j))
       enddo
    enddo
    
    ! --- [2] ISend / IRecv --- !
    !  -- [2-1] ISend       --  !
    do neib=1, NeibPEtot
       is    = export_Index(neib-1) + 1
       len_s = export_Index(neib  ) - export_Index(neib-1)
       call MPI_Isend( sendBuff(is)      , len_s, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Send(neib), ierr )
    enddo
    !  -- [2-2] IRecv       --  !
    do neib=1, NeibPEtot
       ir    = import_Index(neib-1) + 1
       len_r = import_Index(neib  ) - import_Index(neib-1)
       call MPI_Irecv( recvBuff(ir)      , len_r, MPI_DOUBLE_PRECISION &
            &        , NeibPE(neib)      , 0    , MPI_COMM_WORLD &
            &        , request_Recv(neib), ierr )
    enddo
    
    ! --- [3] WaitAll / Update --- !
    !  -- [3-1] WaitAll Receive -- !
    call MPI_WaitAll( NeibPEtot, request_Recv, status_Recv, ierr )
    !  -- [3-2] Update  xVector -- !
    do neib=1, NeibPEtot
       do j=import_Index(neib-1)+1, import_Index(neib)
          xvc(import_Addrs(j)) = recvBuff(j)
       enddo
    enddo
    !  -- [3-3] WaitAll Send    -- !
    call MPI_WaitAll( NeibPEtot, request_Send, status_Send, ierr )
    
    return
  end subroutine BoundaryCommunicate

  
  subroutine PBiCGSTAB_Engine_MPI( b, x )
    implicit none
    double precision, intent(in)  :: b(NxM)
    double precision, intent(out) :: x(NxM)
    integer                       :: iter, i
    double precision              :: alpha, beta, c1, c2, c3, rr
    double precision              :: w(NxMw), r(NxMw), p(NxMw), r0(NxMw), p0(NxMw)
    double precision              :: y(NxMw), z(NxMw), e(NxMw),  v(NxMw),  u(NxMw)
    ! --- Data Structure ------------------------------------------!
    ! b(NxM ) :: { |NxM| }        = { |Data| }                     !
    ! w(NxMw) :: { |NxM|, |2*M| } = { |Data|, |Commu. Region.| }   !
    ! -------------------------------------------------------------!
    ! --- [0] 1st step --- !
    w (:)     = 0.d0
    r (1:NxM) =  b(1:NxM)
    ! r (1:NxM) =  b(1:NxM) - MatrixMultiply( A, w, pt, Nelement, NxM, NxMw, Ncomm )
    r0(1:NxM) =  r(1:NxM)
    p (1:NxM) =  r(1:NxM)
    p0(1:NxM) = r0(1:NxM)
    c1        = myDotProduct( r0(1:NxM), r(1:NxM), NxM )
    
    ! --- [1] Iteration Step --- !
    do iter=1, itermax

       z(1:NxM)  = DiagonalMultiply( Minv(1:NxM), p(1:NxM), NxM )
       call BoundaryCommunicate( z )
       y(1:NxM)  = MatrixMultiply( A, z, pt, Nelement, NxM, NxMw, Ncomm )
       c2        = myDotProduct( r0(1:NxM), y(1:NxM), NxM )
       alpha     = c1 / c2

       e(1:NxM)  = daxpy( -alpha, y(1:NxM), r(1:NxM), NxM )
       u(1:NxM)  = DiagonalMultiply( Minv(1:NxM), e(1:NxM), NxM )
       call BoundaryCommunicate( u )
       v(1:NxM)  = MatrixMultiply( A, u, pt, Nelement, NxM, NxMw, Ncomm )
       
       c3        =   myDotProduct( e(1:NxM), v(1:NxM), NxM ) &
            &      / myDotProduct( v(1:NxM), v(1:NxM), NxM )
       
       !$omp parallel default(none) &
       !$omp shared(w,alpha,z,c3,u,NxM) private(i)
       !$omp do
       do i=1, NxM
          w(i)   = w(i) + alpha*z(i) + c3*u(i)
       enddo
       !$omp end do
       !$omp end parallel
       
       r (1:NxM) = daxpy( -c3, v(1:NxM), e(1:NxM), NxM )
       rr        = myDotProduct( r(1:NxM), r(1:NxM), NxM )
       if ( rr .lt. err ) exit
       c1        = myDotProduct( r0(1:NxM), r(1:NxM), NxM )
       beta      = c1 / ( c2 * c3 )
       
       !$omp parallel default(none) &
       !$omp shared(p,r,beta,c3,y,NxM) private(i)
       !$omp do
       do i=1, NxM
          p(i)   = r(i) + beta * ( p(i) - c3 * y(i) )
       enddo
       !$omp end do
       !$omp end parallel

    enddo

    x(1:NxM) = w(1:NxM)

    if ( myRank.eq.0 ) then
       write(6,'(4x,a)') '[ P-BiCGSTAB Solver ]'
       write(6,'(8x,2(a,i16  ))') 'Iteration  = ', iter, '      / ', itermax
       write(6,'(8x,2(a,e16.9))') 'Error      = ', rr,   '      / ', err
    endif
    return
  end subroutine PBiCGSTAB_Engine_MPI

  
end module PBiCGStabMod
