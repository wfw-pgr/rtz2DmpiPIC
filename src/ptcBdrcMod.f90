module ptcBdrcMod
  implicit none
  logical, parameter :: Flag__BoundaryOutCheck = .true.
contains

  subroutine particleBoundary
    use constants,  only : Boundary1__pt, Boundary2__pt
    use ptCheckMod, only : ptBoundaryOutCheck
    implicit none
    ! -------------------------- !
    ! --- [1] Boundary - 2   --- !
    ! -------------------------- !
    !  -- No particle Boundary is needed for cylindrical Coordinate -- !
    !  -- Reflection @r=R0 in "ptPushMove" in particle.f90          -- !

    ! -------------------------- !
    ! --- [2] Boundary - 1   --- !
    ! -------------------------- !
    if ( trim(Boundary1__pt).eq.'reflect'  ) call ptEdgeReflection
    if ( trim(Boundary1__pt).eq.'periodic' ) call ptEdgeperiodic

    ! -------------------------------- !
    ! --- [3] Check ( Boundary1 )  --- !
    ! -------------------------------- !
    if ( Flag__BoundaryOutCheck ) call ptBoundaryOutCheck

    return
  end subroutine particleBoundary


  subroutine ptEdgeReflection
    use constants , only       : ns, np, myRank, PEtot
    use constants , only       : zp_, vz_, zo_
    use variables , only       : pxv, x1Lengloc
    use pMPIComMod, only       : ParticleExchange
    implicit none
    integer                   :: m, k
    double precision          :: x1Refl1, x1Refl2
    
    ! --------------------------------------- !
    ! --- [1]  Reflection  ( Boundary1 )  --- !
    ! --------------------------------------- !
    ! -- [3-1] Reflection ( Left  )    --  !
    if ( myRank.eq.0       ) then
       x1Refl1 = 0.d0
       !$omp parallel default(none) &
       !$omp shared(np,pxv,x1Refl1) private(m,k)
       !$omp do
       do k=1, ns
          do m=1, np(k)
             if ( pxv(zp_,m,k).le.x1Refl1 ) then
                pxv(zp_,m,k) = 2.d0*x1Refl1 - pxv(zp_,m,k)
                pxv(vz_,m,k) =              - pxv(vz_,m,k)
             endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif
    ! -- [3-2] Reflection ( Right )    --  !
    if ( myRank.eq.PEtot-1 ) then 
       x1Refl2 = x1Lengloc
       !$omp parallel default(none) &
       !$omp shared(np,pxv,x1Refl2) private(m,k)
       !$omp do
       do k=1, ns
          do m=1, np(k)
             if ( pxv(zp_,m,k).gt.x1Refl2 ) then
                pxv(zp_,m,k) = 2.d0*x1Refl2 - pxv(zp_,m,k)
                pxv(vz_,m,k) =              - pxv(vz_,m,k)
             endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif
    ! -- [3-3] Particle Exchange -- !
    call ParticleExchange

    return
  end subroutine ptEdgeReflection


  subroutine ptEdgePeriodic
    use pMPIComMod, only : ParticleExchange
    implicit none
    
    call ParticleExchange
    return
  end subroutine ptEdgePeriodic

  
  subroutine ReflecInXY( psuedxy, pxvrt, R0, sgn, rmin, gamma, dr, dz )
    implicit none
    double precision                :: vnrm2, vcrsx, vdotx, dtprm, theta, phi, costh, sinth
    double precision                :: vtmp(2), hdelta(2)
    double precision, intent(in)    :: R0, sgn, rmin, gamma, dr, dz
    double precision, intent(inout) :: psuedxy(3), pxvrt(8)
    
    !  --- [1] Delta t ' ---  ! :: dt' == > Actually including Gamma :: == Gamma dt'
    psuedxy(1) = pxvrt(1) + rmin
    psuedxy(2) = 0.d0
    vnrm2      = pxvrt(3)**2 + pxvrt(4)**2
    vcrsx      = pxvrt(3)*psuedxy(2) - pxvrt(4)*psuedxy(1)
    vdotx      = pxvrt(3)*psuedxy(1) + pxvrt(4)*psuedxy(2)
    dtprm      = (  - vdotx + sgn * sqrt( R0**2 * vnrm2 - vcrsx**2 ) ) / vnrm2
    
    !  --- [2] P(x,y), theta, phi --- !
    !   -- [2-1] update P(X,Y) --  !
    psuedxy(1) = psuedxy(1) + pxvrt(3) * dtprm
    psuedxy(2) =            + pxvrt(4) * dtprm
    !   -- [2-2] hDelta --  !
    hdelta(1)  = sign( 0.5d0*dz, psuedxy(1) )
    hdelta(2)  = sign( 0.5d0*dr, psuedxy(2) )
    psuedxy(1) = psuedxy(1) + hdelta(1)
    psuedxy(2) = psuedxy(2) + hdelta(2)
    !   -- [2-3] Angle Theta   --  !
    theta      = atan( psuedxy(2) / psuedxy(1) )
    phi        = atan( pxvrt(3)   / pxvrt(4)   )
    costh      = cos(  2.d0 * ( theta + phi )  )
    sinth      = sin(  2.d0 * ( theta + phi )  )

    !  --- [3] Rotation of vx, vy --- !
    vtmp(1)    = pxvrt(3)
    vtmp(2)    = pxvrt(4)
    pxvrt(3)   = costh*vtmp(1) - sinth*vtmp(2)
    pxvrt(4)   = sinth*vtmp(1) + costh*vtmp(2)

    !  --- [4] step forward : dt-dtprm  --- !
    psuedxy(1) = psuedxy(1) + pxvrt(3)*( gamma - dtprm ) - hdelta(1)
    psuedxy(2) = psuedxy(2) + pxvrt(4)*( gamma - dtprm ) - hdelta(2)
    psuedxy(3) = sqrt( psuedxy(1)**2 + psuedxy(2)**2 )

    return
  end subroutine ReflecInXY
  
end module ptcBdrcMod
