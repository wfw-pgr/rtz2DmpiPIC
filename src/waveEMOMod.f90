module waveEMOMod
  implicit none
  double precision            :: B1wave, E1wave, v1wave
  double precision            :: lambda, k_wNum, w_Freq
  double precision            :: kx, kz, wtheta, w_mhdt
  double precision, parameter :: B0wave = 1.d0
  double precision, parameter :: twopi  = 8.d0 * atan( 1.d0 )
contains
  
  
  subroutine setupEMO( wvMode, wvAmpl, wvAngl )
    use constants   , only        : myRank
    use displayMod  , only        : displaySectionTitle
    implicit none
    integer         , intent(in) :: wvMode
    double precision, intent(in) :: wvAmpl, wvAngl
    
    if ( myRank.eq.0 ) call displaySectionTitle( 'setup EMO', '-', 4, 4, 'subsection' )
    call setupEMO__DefwvConst( wvMode, wvAmpl, wvAngl )
    call setupEMO__setupField
    call setupEMO__ptPosition
    call setupEMO__DefNspRepr
    call setupEMO__ptVelocity
    
    return
  end subroutine setupEMO
  
  
  subroutine setupEMO__DefwvConst( wvMode, wvAmpl, wvAngl )
    use constants   , only        : cv, wp, qm, dt
    use variables   , only        : x1Leng
    implicit none
    integer         , intent(in) :: wvMode
    double precision, intent(in) :: wvAmpl, wvAngl

    ! ------------------------------ !
    ! --- [1] Parameter Settings --- !
    ! ------------------------------ !
    wtheta =   wvAngl / 360.d0 * twopi
    lambda =   x1Leng / dble( wvMode )
    k_wNum =    twopi / lambda
    w_Freq =   sqrt( wp(1)**2 + wp(2)**2 + cv**2 * k_wNum**2 )
    B1wave =   B0wave * wvAmpl
    E1wave =   w_Freq / k_wNum * B1wave
    v1wave =   qm(1)  / k_wNum * B1wave
    kz     =   k_wNum * cos( wtheta )
    kx     =   k_wNum * sin( wtheta )
    w_mhdt =   w_Freq * ( - 0.5d0*dt )
    
    return
  end subroutine setupEMO__DefwvConst

  
  subroutine setupEMO__setupField
    use constants, only : LIs, LJs
    use constants, only : ex_, ey_, ez_, bx_, by_, bz_
    use variables, only : EBf, x1gf, x1gh
    implicit none
    integer            :: i, j

    ! -------------------------- !
    ! --- [1] setup B-Field  --- !
    ! -------------------------- !
    do j=1, LJs
       do i=1, LIs
          EBf(bx_,i,j) = B0wave
          EBf(by_,i,j) = B1wave * cos( kz*x1gf(i) - w_mhdt )
          EBf(bz_,i,j) = 0.d0
       enddo
    enddo
    ! -------------------------- !
    ! --- [2] setup E-Field  --- !
    ! -------------------------- !
    do j=1, LJs
       do i=1, LIs
          EBf(ex_,i,j) = E1wave * cos( kz*x1gh(i)        )
          EBf(ey_,i,j) = 0.d0
          EBf(ez_,i,j) = 0.d0
       enddo
    enddo

    return
  end subroutine setupEMO__setupField


  subroutine setupEMO__ptPosition
    use constants , only : LI, LJ, ns, np, npc, npt, nptMax, ppcInit, myRank, PEtot
    use constants , only : xp_, zp_, xo_, zo_, wp_, el_, io_
    use variables , only : pxv, x1Lengloc, x2Lengloc, x1Leng
    use mt19937Mod, only : multi_grnd
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: k, m, iPE, ith, nptSum, ierr
    integer             :: npPE(2,0:PEtot-1)
    double precision    :: x1Minh, x1Maxh, x2Minh, x2Maxh

    ! ------------------------------ !
    ! --- [1] #.of Particle      --- !
    ! ------------------------------ !
    nptSum = (LI-1)*(LJ-1)*ppcInit*ns
    npt    = nint( ( x1Lengloc / x1Leng )*dble( nptSum ) )
    npc(1) = npt / 2
    npc(2) = npt / 2
    npt    = npc(1) + npc(2)
    np(1)  = npc(1)
    np(2)  = npc(2)
    call MPI_AllGather( npc , 2       , MPI_INTEGER, &
         &              npPE, 2       , MPI_INTEGER, &
         &              MPI_COMM_WORLD, ierr     )

    ! ------------------------------ !
    ! --- [2] Particle Showering --- !
    ! ------------------------------ !
    !  -- [2-1] Range Settings    -- !
    x1Minh = 0.d0
    x1Maxh = x1Lengloc
    x2Minh = 0.d0
    x2Maxh = x2Lengloc
    !  -- [2-2] Position ( ele. ) -- !
    !$omp parallel default(none) &
    !$omp shared(np,pxv,x1Minh,x1Maxh,x2Minh,x2Maxh) private(m,ith)
    !$ ith = omp_get_thread_num() + 1
    !$omp do
    do m=1, np(1)
       pxv(xp_,m,el_) = ( x2Maxh - x2Minh ) * multi_grnd( ith ) + x2Minh
       pxv(zp_,m,el_) = ( x1Maxh - x1Minh ) * multi_grnd( ith ) + x1Minh
       pxv(wp_,m,el_) = 1.d0
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [2-3] Position ( ion  ) -- !
    !$omp parallel default(none) &
    !$omp shared(np,pxv) private(m)
    !$omp do
    do m=1, np(2)
       pxv(xp_,m,io_) = pxv(xp_,m,el_)
       pxv(zp_,m,io_) = pxv(zp_,m,el_)
       pxv(wp_,m,io_) = 1.d0
    enddo
    !$omp end do
    !$omp end parallel
    
    ! -------------------------------------- !
    ! --- [3] Copy Old Particle Position --- !
    ! -------------------------------------- !
    !$omp parallel default(none) &
    !$omp shared(pxv,np) private(m,k)
    !$omp do
    do k=1, ns
       do m=1, np(k)
          pxv(xo_,m,k) = pxv(xp_,m,k)
          pxv(zo_,m,k) = pxv(zp_,m,k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! -------------------------------------- !
    ! --- [4] Display                    --- !
    ! -------------------------------------- !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'(2x,a)') '[ setupGS_ptPosition @setupGSMod ]'
       write(6,*)
       write(6,'(2x,6(a10,1x))') 'iPE', 'npc(e)', 'npc(i)', 'npc(1+2)', 'nptMax', 'nptSum'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       do iPE=0, PEtot-1
          write(6,'(2x,6(i10,1x))') iPE, npPE(1,iPE), npPE(2,iPE), &
               &                    npPE(1,iPE)+npPE(2,iPE), nptMax, nptSum
       enddo
       write(6,*)
    endif

    return
  end subroutine setupEMO__ptPosition


  subroutine setupEMO__DefNspRepr
    use constants , only : LIs, LJs, myRank, ppcInit, PEtot
    use variables , only : RhoPerNsp, NspRepres, volRepres, perCell, rwVh
    use momentsMod, only : nCountUp_SP
    implicit none
    include 'mpif.h'
    integer             :: NpMid, NpInner, NpOuter
    double precision    :: rhoRepres, wVhRepres

    ! ---------------------------------------- !
    ! --- [2] Determine Np0 representative --- !
    ! ---------------------------------------- !
    !   -- [2-1] Define Rho/Nsp & Nsp/Rho   -- !
    wVhRepres = rwVh( LJs/2 )
    rhoRepres = 1.d0
    NspRepres = ppcInit
    NspRepres = NspRepres / rhoRepres
    RhoPerNsp = 1.d0      / NspRepres
    !   -- [2-2] Redefine W_P -- !
    call nCountUp_SP( perCell )
    NpMid     = perCell( (LIs-2)/2,(LJs-2)/2,1 )
    NpInner   = perCell( (LIs-2)/2,        2,1 )
    NpOuter   = perCell( (LIs-2)/2,(LJs-3)  ,1 )

    ! ---------------------------- !
    ! --- [3]  Display  Info.  --- !
    ! ---------------------------- !
    if ( myRank.eq.PEtot/2 ) then
       !  -- [3-1] Display Info. PE  -- !
       write(6,'(2x,a)'      ) ' [ setupEMO__DefNspRepr @waveEMOMod ]'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(5x,a)') '---  Determination of Np0   ---'
       write(6,'(8x,a,f12.5)') '  RhoPerNsp  :: Phys.Pt. Density / SuperPaticle :: ', RhoPerNsp
       write(6,'(8x,a,f12.5)') '  NspRepres  :: Representative value of Nsp     :: ', NspRepres
       write(6,'(8x,a,f12.5)') '  rhoRepres  :: Representative value of Rho     :: ', rhoRepres
       write(6,'(8x,a,f12.5)') '  volRepres  :: Representative value of Volume  :: ', wVhRepres
       write(6,*)
       write(6,'(5x,a)') '---  SuperParticle @ Center ---'
       write(6,'(5x,a,i12  )') '  Np(j=   2) :: # of superparticle at  j =    2 :: ', NpInner
       write(6,'(5x,a,i12  )') '  Np(j=LJ/2) :: # of superparticle at  j = LJ/2 :: ', NpMid  
       write(6,'(5x,a,i12  )') '  Np(j=LJ-1) :: # of superparticle at  j = LJ-1 :: ', NpOuter
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,*)
    endif
    
    return
  end subroutine setupEMO__DefNspRepr


  subroutine setupEMO__ptVelocity
    use constants , only : ns , np , vthe, vthi, myRank
    use constants , only : zp_, vx_, vy_ , vz_ , el_, io_
    use variables , only : x1g, pxv
    use randGenMod, only : multi_gdev
    implicit none
    integer            :: m, k, ith
    double precision   :: zpp, u1x, vth(ns)
    integer, parameter :: x1Frm_ = 2

    ! --- [1]  Wave + thermal Velocity --- !
    !  -- [1-1] electron ( wave )      --  !
    ith    = 1
    do m=1, np(el_)
       zpp            = x1g(myRank,x1Frm_) + pxv(zp_,m,el_)
       u1x            = v1wave * ( - sin( kz*zpp - w_mhdt ) )
       pxv(vx_,m,el_) = u1x + vthe*multi_gdev( ith )
       pxv(vy_,m,el_) =       vthe*multi_gdev( ith )
       pxv(vz_,m,el_) =       vthe*multi_gdev( ith )
    enddo
    !  -- [1-2] ion ( stationary )     --  !
    ith    = 1
    do m=1, np(io_)
       pxv(vx_,m,io_) =       vthi*multi_gdev( ith )
       pxv(vy_,m,io_) =       vthi*multi_gdev( ith )
       pxv(vz_,m,io_) =       vthi*multi_gdev( ith )
    enddo
    
    return
  end subroutine setupEMO__ptVelocity

  
end module waveEMOMod
