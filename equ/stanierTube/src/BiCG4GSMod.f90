module BiCGS4GSMod
  use mPBiCGSTAB
contains

  subroutine GSLaplacian( dz, dr, N, M, rMin )
    implicit none
    integer         , intent(in)  :: N, M
    double precision, intent(in)  :: dz, dr, rMin
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
       rcoef1         = (  dble(j-1)*dr + rMin  ) / ( ( dble(j-1)+0.5d0 )*dr + rMin )
       rcoef2         = (  dble(j-1)*dr + rMin  ) / ( ( dble(j-1)-0.5d0 )*dr + rMin )
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

  
end module BiCGS4GSMod
