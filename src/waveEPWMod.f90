module waveEPWMod
  implicit none
  double precision            :: B1wave, E1wave, v1wave, x1wave, n1wave
  double precision            :: lambda, k_wNum, w_Freq
  double precision            :: kx, kz, wtheta, w_mhdt
  double precision, parameter :: B0wave = 1.d0
  double precision, parameter :: twopi  = 8.d0 * atan( 1.d0 )
  double precision, allocatable :: rho(:,:)
contains
  
  
  subroutine setupEPW( wvMode, wvAmpl, wvAngl )
    use constants   , only        : myRank
    use displayMod  , only        : displaySectionTitle
    implicit none
    integer         , intent(in) :: wvMode
    double precision, intent(in) :: wvAmpl, wvAngl
    
    if ( myRank.eq.0 ) call displaySectionTitle( 'setup EPW', '-', 4, 4, 'subsection' )
    call setupEPW__DefwvConst( wvMode, wvAmpl, wvAngl )
    call setupEPW__setupField
    call setupEPW__ptPosition
    ! call setupEPW__uniformDst
    call setupEPW__DefNspRepr
    ! call setupEPW__ptDeviation
    call setupEPW__ptVelocity
    
    return
  end subroutine setupEPW
  
  
  subroutine setupEPW__DefwvConst( wvMode, wvAmpl, wvAngl )
    use constants   , only        : cv, wp, qm, dt, vthcv, valfe
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
    w_Freq =   sqrt( wp(1)**2 + 3.d0 * vthcv**2 * k_wNum**2 )
    B1wave =   0.d0
    E1wave =   wvAmpl * B0wave
    v1wave =   qm(1)  / w_Freq    * E1wave
    x1wave =   qm(1)  / w_Freq**2 * E1wave
    n1wave =   k_wNum * vAlfe**2  * E1wave
    kz     =   k_wNum * cos( wtheta )
    kx     =   k_wNum * sin( wtheta )
    w_mhdt =   w_Freq * ( - 0.5d0*dt )
    
    return
  end subroutine setupEPW__DefwvConst

  
  subroutine setupEPW__setupField
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
          EBf(bx_,i,j) = 0.d0
          EBf(by_,i,j) = 0.d0
          EBf(bz_,i,j) = B0wave
       enddo
    enddo
    ! -------------------------- !
    ! --- [2] setup E-Field  --- !
    ! -------------------------- !
    do j=1, LJs
       do i=1, LIs
          EBf(ex_,i,j) = 0.d0
          EBf(ey_,i,j) = 0.d0
          EBf(ez_,i,j) = E1wave * cos( kz*x1gf(i) )
       enddo
    enddo
    ! -------------------------- !
    ! --- [3] setup n-Field  --- !
    ! -------------------------- !
    allocate( rho(LIs,LJs) )
    do j=1, LJs
       do i=1, LIs
          rho(i,j) = 1.d0 + n1wave * sin( kz*x1gh(i) )
       enddo
    enddo
    return
  end subroutine setupEPW__setupField

  
  subroutine setupEPW__ptPosition
    use constants , only : LI, LJ, LIs, LJs, myRank, PEtot
    use constants , only : ns, nptTot, npt, np, npc, nptMax, ppcInit
    use constants , only : drInv, dzInv, dr, dz
    use constants , only : xp_, zp_, vx_, vy_, vz_, xo_, zo_, wp_, el_, io_
    use variables , only : pxv, x1Lengloc, x2Lengloc, rwVh
    use mt19937Mod, only : multi_grnd
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: k, m, i, j, ip, jp, ith, iPE, ierr, nptSum
    integer             :: npPE(2,0:PEtot-1)
    double precision    :: prob, sum, sumI, sumG, pdfMax
    double precision    :: pdf(LIs,LJs), x1Maxh, x1Minh, x2Maxh, x2Minh
    integer             :: LIpt, LJpt
    double precision    :: dxpt, dzpt, Nzp, Nxp



    ! --------------------------------------------------- !
    ! --- [1] Probability & Total Number of Particles --- !
    ! --------------------------------------------------- !
    !  -- [1-1] Probability Density Function           -- !
    pdf(:,:) = 0.d0
    sum      = 0.d0
    !$omp parallel default(none) &
    !$omp shared(pdf,rho,sum,LIs,LJs) private(i,j)
    !$omp do reduction(+:sum)
    do j=2, LJs-1
       do i=2, LIs-1
          pdf(i,j) =          rho(i,j)
          sum      =    sum + pdf(i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    call MPI_ALLREDUCE( sum, sumG, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )

    !  -- [1-2] Number of Particles in Each Processor  -- !
    nptSum   = (LI-1)*(LJ-1)*ppcInit*ns
    npt      = nint( ( sum / sumG ) * dble(nptSum) )
    npc(1)   = npt / 2
    npc(2)   = npt / 2
    npt      = npc(1) + npc(2)
    np(1)    = npc(1)
    np(2)    = npc(2)
    call MPI_AllGather( npc , 2       , MPI_INTEGER, &
         &              npPE, 2       , MPI_INTEGER, &
         &              MPI_COMM_WORLD, ierr     )
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
    !  -- [1-3] Redefinition of PDF -- !
    sumI     =      1.d0 / sum
    pdf(:,:) =  pdf(:,:) * sumI
    pdfMax   = maxval( pdf(2:LIs-1,2:LJs-1) )
    
    ! ------------------------------ !
    ! --- [2] Particle Showering --- !
    ! ------------------------------ !
    !  -- [2-1] Range Settings    -- !
    x1Minh = 0.d0
    x1Maxh = x1Lengloc
    x2Minh = 0.d0
    x2Maxh = x2Lengloc
    ! ith = 1
    ! if ( myRank.eq.0       ) x1Minh =           + 0.5d0*dz
    ! if ( myRank.eq.PEtot-1 ) x1Maxh = x1Lengloc - 0.5d0*dz
    !  -- [2-2] Put Particles     -- !
    !$omp parallel default(none) &
    !$omp shared (pxv,np,drInv,x1Maxh,x1Minh,x2Maxh,x2Minh,dzInv,pdfMax,pdf,rwVh,LIs,LJs) &
    !$omp private(m,ip,jp,prob,ith)
    !$ ith = omp_get_thread_num() + 1
    !$omp do
    do m=1, np(1)
100    continue
       !  - (1) Monte-Calro Position - !
       pxv(xp_,m,el_) = ( x2Maxh - x2Minh ) * multi_grnd( ith ) + x2Minh
       pxv(zp_,m,el_) = ( x1Maxh - x1Minh ) * multi_grnd( ith ) + x1Minh
       jp             = ceiling( pxv(xp_,m,el_)*drInv ) - 1 + 2
       ip             = ceiling( pxv(zp_,m,el_)*dzInv ) - 1 + 2
       !  - (2) Judge Compatibility  - !
       prob           = multi_grnd( ith ) * pdfMax
       if ( prob.gt.pdf(ip,jp) ) goto 100
       !  - (3) Store Data           - !
       pxv(wp_,m,el_) = 1.d0
    enddo
    !$omp end do
    !$omp end parallel

    Nzp  = sqrt( dble(ppcInit) )
    Nxp  = sqrt( dble(ppcInit) )
    dzpt =  dz / dble(Nzp)
    dxpt =  dr / dble(Nxp)
    LJpt = ( LJs-2 ) * Nxp
    LIpt = ( LIs-2 ) * Nzp
    !$omp parallel default(none) &
    !$omp shared(pxv,dzpt,dxpt,LIpt,LJpt) private(i,j,m)
    !$omp do
    do j=1, LJpt
       do i=1, LIpt
          m            = (j-1)*( LIpt ) + i
          pxv(zp_,m,io_) = dble( i-1 ) * dzpt
          pxv(xp_,m,io_) = dble( j-1 ) * dxpt
          pxv(zo_,m,io_) = dble( i-1 ) * dzpt
          pxv(xo_,m,io_) = dble( j-1 ) * dxpt
          pxv(wp_,m,io_) = 1.d0
       enddo
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

    return
  end subroutine setupEPW__ptPosition



  subroutine setupEPW__DefNspRepr
    use constants , only : LIs, LJs, myRank, ppcInit, PEtot, rhoRepres_0
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
    rhoRepres = rhoRepres_0
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
       write(6,'(2x,a)'      ) ' [ setupEPW__DefNspRepr @waveEPWMod ]'
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
  end subroutine setupEPW__DefNspRepr


  subroutine setupEPW__ptVelocity
    use constants , only : ns , np , vthe, vthi, myRank
    use constants , only : zp_, vx_, vy_ , vz_ , el_, io_
    use variables , only : x1g, pxv
    use randGenMod, only : multi_gdev
    implicit none
    integer            :: m, k, ith
    double precision   :: zpp, u1z, vth(ns)
    integer, parameter :: x1Frm_ = 2

    ! --- [1]  Wave + thermal Velocity --- !
    !  -- [1-1] electron ( wave )      --  !
    ith      = 1
    vth(el_) = vthe * 1.d0
    do m=1, np(el_)
       zpp            = x1g(myRank,x1Frm_) + pxv(zp_,m,el_)
       u1z            = v1wave * ( - sin( kz*zpp - w_mhdt ) )
       pxv(vx_,m,el_) =       vth(el_)*multi_gdev( ith )
       pxv(vy_,m,el_) =       vth(el_)*multi_gdev( ith )
       pxv(vz_,m,el_) = u1z + vth(el_)*multi_gdev( ith )
    enddo
    !  -- [1-2] ion ( stationary )     --  !
    ith      = 1
    vth(io_) = vthi * 1.d0
    do m=1, np(io_)
       pxv(vx_,m,io_) =       vth(io_)*multi_gdev( ith )
       pxv(vy_,m,io_) =       vth(io_)*multi_gdev( ith )
       pxv(vz_,m,io_) =       vth(io_)*multi_gdev( ith )
    enddo
    
    return
  end subroutine setupEPW__ptVelocity


  subroutine setupEPW__ptDeviation
    use constants , only : np , myRank
    use constants , only : zp_, el_
    use variables , only : x1g, pxv
    implicit none
    integer            :: m
    double precision   :: zpp
    integer, parameter :: x1Frm_ = 2

    ! --- [1]  Particle Deviation      --- !
    !  -- [1-1] electron ( wave )      --  !
    do m=1, np(el_)
       zpp            = pxv(zp_,m,el_) + x1g(myRank,x1Frm_)
       pxv(zp_,m,el_) = pxv(zp_,m,el_) + x1wave*( - cos( kz*zpp ) )
    enddo

    return
  end subroutine setupEPW__ptDeviation


  subroutine setupEPW__uniformDst
    use constants , only : LIs, LJs, ns, np, npc, npt, nptMax, ppcInit, myRank, PEtot
    use constants , only : xp_, zp_, xo_, zo_, wp_, dr, dz
    use variables , only : pxv
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: i, j, k, m, iPE, ierr
    integer             :: LIpt, LJpt, nptSum
    double precision    :: dxpt, dzpt, Nzp, Nxp
    integer             :: npPE(2,0:PEtot-1)

    ! ------------------------------ !
    ! --- [1] #.of Particle      --- !
    ! ------------------------------ !
    nptSum = 0
    npt    = (LIs-2) * ( LJs-2 ) * ppcInit
    npc(1) = npt
    npc(2) = npt
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
    Nzp  = sqrt( dble(ppcInit) )
    Nxp  = sqrt( dble(ppcInit) )
    dzpt =  dz / dble(Nzp)
    dxpt =  dr / dble(Nxp)
    LJpt = ( LJs-2 ) * Nxp
    LIpt = ( LIs-2 ) * Nzp
    !$omp parallel default(none) &
    !$omp shared(pxv,dzpt,dxpt,LIpt,LJpt) private(i,j,m,k)
    !$omp do
    do k=1, ns
       do j=1, LJpt
          do i=1, LIpt
             m            = (j-1)*( LIpt ) + i
             pxv(zp_,m,k) = dble( i-1 ) * dzpt
             pxv(xp_,m,k) = dble( j-1 ) * dxpt
             pxv(zo_,m,k) = dble( i-1 ) * dzpt
             pxv(xo_,m,k) = dble( j-1 ) * dxpt
             pxv(wp_,m,k) = 1.d0
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! -------------------------------------- !
    ! --- [3] Display                    --- !
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
  end subroutine setupEPW__uniformDst

  
end module waveEPWMod


  ! subroutine setupEPW__ptPosition
  !   use constants , only : LI, LJ, ns, np, npc, npt, nptMax, ppcInit, myRank, PEtot
  !   use constants , only : xp_, zp_, xo_, zo_, wp_, el_, io_
  !   use variables , only : pxv, x1Lengloc, x2Lengloc, x1Leng
  !   use mt19937Mod, only : multi_grnd
  !   !$ use omp_lib
  !   implicit none
  !   include 'mpif.h'
  !   integer             :: k, m, iPE, ith, nptSum, ierr
  !   integer             :: npPE(2,0:PEtot-1)
  !   double precision    :: x1Minh, x1Maxh, x2Minh, x2Maxh

  !   ! ------------------------------ !
  !   ! --- [1] #.of Particle      --- !
  !   ! ------------------------------ !
  !   nptSum = (LI-1)*(LJ-1)*ppcInit*ns
  !   npt    = nint( ( x1Lengloc / x1Leng )*dble( nptSum ) )
  !   npc(1) = npt / 2
  !   npc(2) = npt / 2
  !   npt    = npc(1) + npc(2)
  !   np(1)  = npc(1)
  !   np(2)  = npc(2)
  !   call MPI_AllGather( npc , 2       , MPI_INTEGER, &
  !        &              npPE, 2       , MPI_INTEGER, &
  !        &              MPI_COMM_WORLD, ierr     )

  !   ! ------------------------------ !
  !   ! --- [2] Particle Showering --- !
  !   ! ------------------------------ !
  !   !  -- [2-1] Range Settings    -- !
  !   x1Minh = 0.d0
  !   x1Maxh = x1Lengloc
  !   x2Minh = 0.d0
  !   x2Maxh = x2Lengloc
  !   !  -- [2-2] Position ( ele. ) -- !
  !   !$omp parallel default(none) &
  !   !$omp shared(np,pxv,x1Minh,x1Maxh,x2Minh,x2Maxh) private(m,ith)
  !   !$ ith = omp_get_thread_num() + 1
  !   !$omp do
  !   do m=1, np(1)
  !      pxv(xp_,m,el_) = ( x2Maxh - x2Minh ) * multi_grnd( ith ) + x2Minh
  !      pxv(zp_,m,el_) = ( x1Maxh - x1Minh ) * multi_grnd( ith ) + x1Minh
  !      pxv(wp_,m,el_) = 1.d0
  !   enddo
  !   !$omp end do
  !   !$omp end parallel
  !   !  -- [2-3] Position ( ion  ) -- !
  !   !$omp parallel default(none) &
  !   !$omp shared(np,pxv) private(m)
  !   !$omp do
  !   do m=1, np(2)
  !      pxv(xp_,m,io_) = pxv(xp_,m,el_)
  !      pxv(zp_,m,io_) = pxv(zp_,m,el_)
  !      pxv(wp_,m,io_) = 1.d0
  !   enddo
  !   !$omp end do
  !   !$omp end parallel
    
  !   ! -------------------------------------- !
  !   ! --- [3] Copy Old Particle Position --- !
  !   ! -------------------------------------- !
  !   !$omp parallel default(none) &
  !   !$omp shared(pxv,np) private(m,k)
  !   !$omp do
  !   do k=1, ns
  !      do m=1, np(k)
  !         pxv(xo_,m,k) = pxv(xp_,m,k)
  !         pxv(zo_,m,k) = pxv(zp_,m,k)
  !      enddo
  !   enddo
  !   !$omp end do
  !   !$omp end parallel

  !   ! -------------------------------------- !
  !   ! --- [4] Display                    --- !
  !   ! -------------------------------------- !
  !   if ( myRank.eq.0 ) then
  !      write(6,*)
  !      write(6,'(2x,a)') '[ setupGS_ptPosition @setupGSMod ]'
  !      write(6,*)
  !      write(6,'(2x,6(a10,1x))') 'iPE', 'npc(e)', 'npc(i)', 'npc(1+2)', 'nptMax', 'nptSum'
  !      write(6,'(2x,a)') '----------------------------------------------------------------------'
  !      do iPE=0, PEtot-1
  !         write(6,'(2x,6(i10,1x))') iPE, npPE(1,iPE), npPE(2,iPE), &
  !              &                    npPE(1,iPE)+npPE(2,iPE), nptMax, nptSum
  !      enddo
  !      write(6,*)
  !   endif

  !   return
  ! end subroutine setupEPW__ptPosition
