module bdrc
contains
  
  subroutine updateBC
    use constants, only : bctype, Bv0, myRank
    implicit none
    !  --- [1]  solve poisson equation to obtain U  ---  !
    call BdrFeedBack
    if ( myRank.eq.0 ) then
       write(6,'(4x,a)') '[ Update Boundary   ]'
       write(6,'(8x,a,7x,a,e12.5)') &
            & 'Boundary   = Conducting Wall', '-- Bv0    = ', Bv0
    endif
    if ( bctype .eq. 'pcw' ) return

    return
  end subroutine updateBC

  
  subroutine BdrFeedBack
    use constants    , only : Nr, Nz
    use constants    , only : Bsensor_r1, Bsensor_r2, Bv0, alpha_star
    use variables    , only : r, z, psi, BdrcLR, BdrcTB
    use NewtonRaphson, only : cSpline2D
    implicit none
    integer                :: i, j
    double precision       :: psi1, psi2, dpsidr, dpsidz
    double precision       :: S1, S2, dBv0
    
    !  --- [1] Measure and Feedback  ---  !
    S1   = 0.5d0 * Bsensor_r1**2
    S2   = 0.5d0 * Bsensor_r2**2
    call cSpline2D( z, r, psi, 0.d0, Bsensor_r1, psi1, dpsidr, dpsidz )
    call cSpline2D( z, r, psi, 0.d0, Bsensor_r2, psi2, dpsidr, dpsidz )
    dBv0 = - alpha_star * ( psi1 - psi2 ) / ( S1 - S2 ) 
    Bv0  =   Bv0 + dBv0

    !  --- [2] Update Boundary  --- !
    do j=1, Nr
       BdrcLR(j,1) = 0.5d0 * r( j)**2 * Bv0
       BdrcLR(j,2) = 0.5d0 * r( j)**2 * Bv0
    enddo
    do i=1, Nz
       BdrcTB(i,1) = 0.5d0 * r(Nr)**2 * Bv0
       BdrcTB(i,2) = 0.5d0 * r( 1)**2 * Bv0
    enddo

    return
  end subroutine BdrFeedBack

  
end module bdrc
