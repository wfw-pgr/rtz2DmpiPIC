module cylPoisMod
  use mPBiCGSTAB
contains

  ! =================================================================== !
  ! ===  BiCGSTABCyl2D_MPI  :: BiCGSTAB solver for Cylindrical 2D   === !
  ! =================================================================== !
  subroutine BiCGSTABCyl2D_MPI( x, source, dz, dr, N, M, rMin )
    implicit none
    include 'mpif.h'
    integer         , intent(in)    :: N, M
    double precision, intent(in)    :: dr, dz, rMin
    double precision, intent(in)    :: source(N,M)
    double precision, intent(inout) :: x(N,M)
    integer                         :: ierr, iPE, Nbuf, Mbuf, Nitem_buf
    logical, save                   :: Flag__MakeMatrix = .true.
    double precision, allocatable   :: xbuf(:,:)

    ! ------------------------------------- !
    ! --- [1] Make Matrix if No Matrix  --- !
    ! ------------------------------------- !
    if ( Flag__MakeMatrix ) then
       ! -- [1-1] variables for Comm.   --  !
       call DefineCommunication( N, M )
       ! -- [1-2] Make Matrix           --  !
       allocate( xloc(Nloc,Mloc), sloc(Nloc,Mloc) )
       allocate( A(NxM,Ncomm)   , pt(NxM,Ncomm),  &
            &    Nelement(NxM)  , Minv(NxM)       )
       call RZLaplacian( dz, dr, Nloc, Mloc, rMin )
       Flag__MakeMatrix = .false.
    endif
    
    ! ------------------------------------- !
    ! --- [2] call PBiCGSTAB_Engine     --- !
    ! ------------------------------------- !
    sloc(:,:) = source( FromTo(myRank,1):FromTo(myRank,2), 1:Mloc )
    call PBiCGSTAB_Engine_MPI( sloc, xloc )

    ! ------------------------------------- !
    ! --- [3] BroadCast Solution        --- !
    ! ------------------------------------- !
    do iPE=0, PEtot-1
       ! -- [3-1] Prepare iPE's Solution -- !
       Mbuf       = Mloc
       Nbuf       = FromTo(iPE,2) - FromTo(iPE,1) + 1
       Nitem_buf  = Nbuf * Mbuf
       allocate( xbuf(Nbuf,Mbuf) )
       if ( myRank.eq.iPE ) xbuf = xloc
       ! -- [3-2] BroadCast              -- !
       call MPI_Bcast( xbuf, Nitem_buf    , MPI_DOUBLE_PRECISION, &
            &          iPE  , MPI_COMM_WORLD, ierr )
       ! -- [3-3] Save in x(From:To)     -- !
       x( FromTo(iPE,1):FromTo(iPE,2), 1:Mloc ) = xbuf(1:Nbuf,1:Mbuf)
       deallocate( xbuf )
    enddo
    
    return
  end subroutine BiCGSTABCyl2D_MPI


  ! =================================================================== !
  ! ===  RZLaplacian  ::  Calculate Laplacian ( cylindrical )       === !
  ! =================================================================== !
  subroutine RZLaplacian( dz, dr, N, M, rMin )
    implicit none
    include 'mpif.h'
    integer         , intent(in) :: N, M
    double precision, intent(in) :: dz, dr, rMin
    integer                      :: i, j, k, kR, kL
    logical                      :: Flag_SizeERROR
    double precision             :: rcoef, coef(5), coef_base(5), rsqinv, diaginv

    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    !  -- [1-1] Display                 --  !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'       ) '[ RZLaplacian  @ myPBiCGSTAB ]'
       write(6,'(5x,2(a,i12))') '* A-Matrix :: ', NxM, ' x ', NxM
    endif
    !  -- [1-2] Coefficient             --  !
    coef_base(1) =   1.d0 / dr**2
    coef(2)      =   1.d0 / dz**2
    coef_base(3) = - 2.d0 * ( 1.d0 / dz**2 + 1.d0 / dr**2 )
    coef(4)      =   1.d0 / dz**2
    coef_base(5) =   1.d0 / dr**2
    !  -- [1-3] Initialization          --  !
    A(:,:)  = 0.d0
    pt(:,:) = 0
    k       = 1
    Minv(:) = 0.d0
    
    ! ------------------------------------- !
    ! --- [2] Make Matrix               --- !
    ! ------------------------------------- !
    !  -- [2-1] j=1 boundary            --  !
    if ( rMin.eq.0.d0 ) then
       if ( myRank.eq.0 ) write(6,'(5x,a)') '* Geometry     includes r=0 :: rMin = 0.d0'
       coef(3) = coef_base(3) - 1.d0 / ( + 0.5d0 *dr )**2
       coef(5) = 2.d0 * coef_base(5)
       do i=1, N
          A(k,1)      = 0.d0
          A(k,2)      = coef(2)
          A(k,3)      = coef(3)
          A(k,4)      = coef(4)
          A(k,5)      = coef(5)
          pt(k,1)     = k   ! - arbitral - !
          pt(k,2)     = k-1
          pt(k,3)     = k
          pt(k,4)     = k+1
          pt(k,5)     = k+N
          Nelement(k) = 5
          Minv(k)     = 1.d0 / coef(3)
          k = k+1
       enddo
    else if ( rMin.gt.0.d0 ) then
       if ( myRank.eq.0 ) write(6,'(5x,a)') '* Geometry NOT includes r=0 :: rMin = 0.d0'
       do i=1, N
          A(k,1)      = 1.d0
          pt(k,1)     = k
          Nelement(k) = 1
          Minv(k)     = 1.d0
          k           = k+1
       enddo
    else
       stop ' [ERROR] Select rMin > 0 :: illigal rMin value. [ERROR] '
    endif
    !  -- [2-2] Main Region ( j=2~N )   --  !
    do j=2, M-1
       ! -- Coef(1,3,5) Definition -- !
       if ( rMin.eq.0.d0 ) then
          rcoef  = 1.d0 / ( 2.d0 * dr * (  ( dble(j)-0.5d0 )*dr ) )
          rsqinv = 1.d0 / ( ( dble(j)-0.5d0 )*dr )**2
       else
          rcoef  = 1.d0 / ( 2.d0 * dr * (  dble(j-1)*dr + rMin   ) )
          rsqinv = 1.d0 / ( dble(j-1)*dr + rMin )**2
       endif
       coef(1) = coef_base(1) - rcoef
       coef(3) = coef_base(3) - rsqinv
       coef(5) = coef_base(5) + rcoef
       diaginv = 1.d0 / coef(3)
       ! -- Make Matrix -- !
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
    !  -- [2-3] j=M boundary            --  !
    do i=1, N
       A(k,1)      = 1.d0
       pt(k,1)     = k
       Nelement(k) = 1
       Minv(k)     = 1.d0
       k           = k+1
    enddo

    ! ------------------------------------- !
    ! --- [3] Boundary ( Parallel )     --- !
    ! ------------------------------------- !
    if ( PEtot.gt.1 ) then
       if ( ( myRank.gt.0 ).and.( myRank.lt.PEtot-1 ) ) then
          !  -- (Left) Neibour PE       --  !
          kL       = 1
          do j=1, M
             pt(kL,2)     = k
             kL           = kL + N
             k            = k  + 1
          enddo
          !  -- (Right) Neibour PE      --  !
          kR       = N
          do j=1, M
             pt(kR,4)     = k
             kR           = kR + N
             k            = k  + 1
          enddo
       endif
       if ( myRank.eq.0 ) then
          !  -- (Left) Dirichlet        --  !
          kL = 1
          do j=1, M
             A(kL,1)      = 1.d0
             pt(kL,1)     = kL
             Nelement(kL) = 1
             Minv(kL)     = 1.d0 
             kL           = kL + N
          enddo
          !  -- (Right) Neibour PE      --  !
          kR = N
          do j=1, M
             pt(kR,4)     = k
             kR           = kR + N
             k            = k  + 1
          enddo
       endif
       if ( myRank.eq.PEtot-1 ) then
          !  -- (Left) Neibour PE       --  !
          kL = 1
          do j=1, M
             pt(kL,2)     = k
             kL           = kL + N
             k            = k  + 1
          enddo
          !  -- (Right) Dirichlet       --  !
          kR = N
          do j=1, M
             A(kR,1)      = 1.d0
             pt(kR,1)     = kR
             Nelement(kR) = 1
             Minv(kR)     = 1.d0 
             kR           = kR + N
          enddo
       endif
    endif
    
    ! ------------------------------------- !
    ! --- [4] Boundary ( Single Node )  --- !
    ! ------------------------------------- !
    if ( PEtot.eq.1 ) then
       !  -- Left  B.C. ( Dirichlet )   --  !
       kL = 1
       do j=1, M
          A(kL,1)      = 1.d0
          pt(kL,1)     = kL
          Nelement(kL) = 1
          Minv(kL)     = 1.d0 
          kL           = kL + N
       enddo
       !  -- Right B.C. ( Dirichlet )   --  !
       kR = N
       do j=1, M
          A(kR,1)      = 1.d0
          pt(kR,1)     = kR
          Nelement(kR) = 1
          Minv(kR)     = 1.d0 
          kR           = kR + N
       enddo
    endif
    
    ! ------------------------------------- !
    ! --- [5] Size Check                --- !
    ! ------------------------------------- !
    Flag_SizeERROR = .false.
    if ( PEtot.gt.1 ) then
       if ( ( myRank.eq.0 ).or.( myRank.eq.PEtot-1 ) ) then
          if ( (k-1).ne.(NxM+  M) ) Flag_SizeERROR = .true.
       else
          if ( (k-1).ne.(NxM+2*M) ) Flag_SizeERROR = .true.
       endif
    endif
    if ( PEtot.eq.1 ) then
       if    ( (k-1).ne.(NxM    ) ) Flag_SizeERROR = .true.
    endif
    if ( Flag_SizeERROR ) then
       write(6,*) ' [WARNING] Size and k_count are Incompatable.'
       write(6,*) "      myRank / PEtot :: ", myRank, " / ", PEtot
       write(6,*) "              NxM+2M :: ", NxM+2*M
       write(6,*) "              NxM+ M :: ", NxM+  M
       write(6,*) "              NxM    :: ", NxM
       write(6,*) "                 k-1 :: ", k-1
       stop
    endif
    return
  end subroutine RZLaplacian


end module cylPoisMod
