module setupGSMod
  implicit none
  integer         , parameter   :: uer_ = 1, uet_ = 2, uez_ = 3
  integer         , parameter   :: uir_ = 4, uit_ = 5, uiz_ = 6
  double precision, allocatable :: rho(:,:), psi(:,:)
  logical         , parameter   :: Flag__NullField       = .false.
  logical         , parameter   :: Flag__LightSpeedCheck = .false.
  logical         , parameter   :: Flag__NaNDetection    = .false.
  logical         , parameter   :: Flag__FlipBt          = .false.
contains
  
  subroutine setupGS
    use constants , only : LIs, LJs, myRank
    use displayMod, only : displaySectionTitle
    implicit none
    ! ----------------------- !
    ! --- [1]  Beginning  --- !
    ! ----------------------- !
    if ( myRank.eq.0 ) call displaySectionTitle( 'SETUP GS-FIELD', '-', 4, 4, 'subsection' )
    allocate( rho(LIs,LJs), psi(LIs,LJs) )
    ! ----------------------- !
    ! --- [2]  Routines   --- !
    ! ----------------------- !
    call setupGS_LoadFields
    ! call setupGS_NullField
    call setupGS_ptPosition
    call setupGS_DefNspRepr
    call setupGS_setJeiueiE
    call setupGS_ptVelocity
    call setupGS_BFadjustDT
    call setupGS_HalfDTStep
    ! ----------------------- !
    ! --- [3] PostProcess --- !
    ! ----------------------- !
    deallocate( rho, psi )
    if ( myRank.eq.0 ) call displaySectionTitle( ' END SETUP GS ', '-', 4, 4, 'subsection' )
    
    return
  end subroutine setupGS
 
  
  subroutine setupGS_LoadFields
    use constants , only : LIs, LJs, equDir, dr, cRank, myRank, PEtot
    use constants , only : br_, bt_, bz_, er_, et_, ez_, jr_, jt_, jz_
    use variables , only : rf , EBf, Jcr, EBo
    use pcwBdrcMod, only : conductingWall_divB
    use mxwBdrcMod, only : BoundaryB
    implicit none
    integer             :: i, j
    character(100)      :: BFileName, JFileName, FFileName
    double precision    :: EBl(3,LIs,LJs), Jcl(3,LIs,LJs), frp(3,LIs,LJs)
    integer, parameter  :: fl_ = 1, rh_ = 2,  pr_ = 3

    ! ------------------------------- !
    ! --- [1]   Load  GS-Fields   --- !
    ! ------------------------------- !
    BFileName = trim(equDir) // '/' // 'Bfd_' // cRank // '.bin'
    JFileName = trim(equDir) // '/' // 'Jcr_' // cRank // '.bin'
    FFileName = trim(equDir) // '/' // 'frp_' // cRank // '.bin'
    !   -- [1-1] B-Field           -- !
    open (50, file=trim(BFileName), status='old', form='unformatted', convert='LITTLE_ENDIAN' )
    read (50) EBl
    close(50)
    !   -- [1-2] Flux-Rho-Pressure -- !
    open (50, file=trim(FFileName), status='old', form='unformatted', convert='LITTLE_ENDIAN' )
    read (50) frp
    close(50)
    !   -- [1-3] Current           -- !
    open (50, file=trim(JFileName), status='old', form='unformatted', convert='LITTLE_ENDIAN' )
    read (50) Jcl
    close(50)
    !$omp parallel default(none) &
    !$omp shared(EBf,EBl,Jcr,Jcl,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBf(br_,i,j  ) = EBl(1,i,j)
          EBf(bt_,i,j  ) = EBl(2,i,j)
          EBf(bz_,i,j  ) = EBl(3,i,j)
          Jcr(jr_,i,j,0) = Jcl(1,i,j)
          Jcr(jt_,i,j,0) = Jcl(2,i,j)
          Jcr(jz_,i,j,0) = Jcl(3,i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)'            ) '[ setupGS_DefNspRepr @setupGSMod ]'
       write(6,'(2x,a10,a52,2x,a6)') "* SAVE :: ", trim(BFileName), '[ OK ]'
       write(6,'(2x,a10,a52,2x,a6)') "* SAVE :: ", trim(JFileName), '[ OK ]'
       write(6,'(2x,a10,a52,2x,a6)') "* SAVE :: ", trim(FFileName), '[ OK ]'
       write(6,*)
    endif
    
    ! ------------------------------- !
    ! --- [2] B.C. / Flux / Copy  --- !
    ! ------------------------------- !
    !  -- [2-1] Boundary Conditions / Flip Bt -- !
    !   -   B.C.  - !
    call conductingWall_divB( EBf )
    call BoundaryB
    !   - Flip Bt - !
    if ( Flag__FlipBt ) then
       write(6,'(2x,a)') '!!! -------------------------------------------------------------- !!!'
       write(6,'(2x,a)') ' [CAUTION] Toroidal Magnetic Field is flipped... => Case ? [CAUTION]'
       write(6,'(2x,a)') '!!! -------------------------------------------------------------- !!!'
       EBf(bt_,:,:)  = - EBf(bt_,:,:)
    endif
    !  -- [2-2] Magnetic Flux Function / Rho  -- !
    !   -   psi   - !
    psi(1:LIs,1:LJs) = 0.d0
    do j=2, LJs-1
       do i=1, LIs
          psi(i,j+1) = psi(i,j) - 0.5d0*( rf(j+1)*EBf(bz_,i,j+1) + rf(j)*EBf(bz_,i,j) ) * dr
       enddo
    enddo
    !   -   rho   - !
    rho(:,:) = frp(rh_,:,:)
    !  -- [2-3] E-Field, EBo   :: Copy Field  -- !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBo,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBf(er_,i,j) = 0.d0
          EBf(et_,i,j) = 0.d0
          EBf(ez_,i,j) = 0.d0
          EBo(er_,i,j) = EBf(er_,i,j)
          EBo(et_,i,j) = EBf(et_,i,j)
          EBo(ez_,i,j) = EBf(ez_,i,j)
          EBo(br_,i,j) = EBf(br_,i,j)
          EBo(bt_,i,j) = EBf(bt_,i,j)
          EBo(bz_,i,j) = EBf(bz_,i,j)          
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine setupGS_LoadFields


  subroutine setupGS_ptPosition
    use constants , only : LI, LJ, LIs, LJs, myRank, PEtot
    use constants , only : ns, nptTot, npt, np, npc, nptMax, ppcInit
    use constants , only : drInv, dzInv
    use constants , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_, el_, io_
    use variables , only : pxv, x1Lengloc, x2Lengloc, rwVh
    use mt19937Mod, only : multi_grnd
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: k, m, i, j, ip, jp, ith, iPE, ierr, nptSum
    integer             :: npPE(ns,0:PEtot-1)
    double precision    :: prob, sum, sumI, sumG, pdfMax
    double precision    :: pdf(LIs,LJs), x1Maxh, x1Minh, x2Maxh, x2Minh

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
          pdf(i,j) = 0.25d0 * ( rho(i,j  ) + rho(i+1,j  ) &
               &              + rho(i,j+1) + rho(i+1,j+1) )
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
       pxv(rp_,m,el_) = ( x2Maxh - x2Minh ) * multi_grnd( ith ) + x2Minh
       pxv(zp_,m,el_) = ( x1Maxh - x1Minh ) * multi_grnd( ith ) + x1Minh
       jp             = ceiling( pxv(rp_,m,el_)*drInv ) - 1 + 2
       ip             = ceiling( pxv(zp_,m,el_)*dzInv ) - 1 + 2
       !  - (2) Judge Compatibility  - !
       prob           = multi_grnd( ith ) * pdfMax
       if ( prob.gt.pdf(ip,jp) ) goto 100
       !  - (3) Store Data           - !
       pxv(wp_,m,el_) = rwVh(jp)
       pxv(rp_,m,io_) = pxv(rp_,m,el_)
       pxv(zp_,m,io_) = pxv(zp_,m,el_)
       pxv(wp_,m,io_) = pxv(wp_,m,el_)
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
          pxv(ro_,m,k) = pxv(rp_,m,k)
          pxv(zo_,m,k) = pxv(zp_,m,k)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine setupGS_ptPosition


  subroutine setupGS_DefNspRepr
    use constants, only : LI, LJ, LIs, LJs, nptTot, ppcMax, myRank, PEtot
    use variables, only : RhoPerNsp, NspRepres, rwVh, w_p, volRepres, perCell
    use momentsMod, only : nCountUp_CC, nCountUp_SP
    implicit none
    include 'mpif.h'
    integer            :: iMax, jMax, iPE, iPEMax, ierr
    integer            :: NpMid, NpInner, NpOuter
    double precision   :: NspMax, rhoMax, dnsMax, rwVMax, neMax, niMax
    double precision   :: max_sBuff(5), max_rBuff(5,0:PEtot-1)

    ! ---------------------------------------- !
    ! --- [1] Obtain Density Max in All PE --- !
    ! ---------------------------------------- !
    !  -- [1-1] Determine Reference     -- !
    NspRepres = 1.d0   ! ( Temporary )
    RhoPerNsp = 1.d0 / NspRepres
    call nCountUp_CC( w_p )
    !  -- [1-2] Standard Point          -- !
    !   - [CAUTION] jRepress = jMax ??  -  !
    iMax         = LIs / 2
    jMax         = LJs / 2
    rhoMax       =  rho(iMax,jMax  )
    dnsMax       =  w_p(iMax,jMax,1)
    rwVMax       = rwVh(jMax)
    !  -- [1-3] Max Info. Communicate   -- !
    max_sBuff(1) = dble(  iMax)
    max_sBuff(2) = dble(  jMax)
    max_sBuff(3) = dble(rhoMax)
    max_sBuff(4) = dble(dnsMax)
    max_sBuff(5) = dble(rwVMax)
    call MPI_ALLGATHER( max_sBuff,   5, MPI_DOUBLE_PRECISION, &
         &              max_rBuff,   5, MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, ierr )
    !  -- [1-4] Max among All Processor -- !
    iPEMax = 0
    rhoMax = max_rBuff(3,0)
    do iPE = 1, PEtot-1
       if ( max_rBuff(3,iPE).gt.rhoMax ) then
          rhoMax = max_rBuff(3,iPE)
          iPEMax = iPE
       endif
    enddo
    iMax   = nint( max_rBuff(1,iPEMax) )
    jMax   = nint( max_rBuff(2,iPEMax) )
    rhoMax = max_rBuff(3,iPEMax)
    dnsMax = max_rBuff(4,iPEMax)
    rwVMax = max_rBuff(5,iPEMax)

    ! ---------------------------------------- !
    ! --- [2] Determine Np0 representative --- !
    ! ---------------------------------------- !
    !   -- [2-1] Define Rho/Nsp & Nsp/Rho   -- !
    NspMax    = NspRepres * rwVMax * dnsMax
    NspRepres = NspMax / rhoMax
    RhoPerNsp = 1.d0   / NspRepres
    !   -- [2-2] Redefine W_P -- !
    call nCountUp_SP( perCell )
    NpMid     = perCell( (LIs-2)/2,(LJs-2)/2,1 )
    NpInner   = perCell( (LIs-2)/2,        2,1 )
    NpOuter   = perCell( (LIs-2)/2,(LJs-3)  ,1 )
    !   -- [2-3] Redefine W_P -- !
    call nCountUp_CC( w_p )
    neMax     = maxval( w_p(1:LIs,1:LJs,1) )
    niMax     = maxval( w_p(1:LIs,1:LJs,2) )
    call MPI_Bcast( neMax, 1, MPI_DOUBLE_PRECISION, iPEMax, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( niMax, 1, MPI_DOUBLE_PRECISION, iPEMax, MPI_COMM_WORLD, ierr )

    ! ---------------------------- !
    ! --- [3]  Display  Info.  --- !
    ! ---------------------------- !
    if ( myRank.eq.0 ) then
       !  -- [3-1] Display Info. PE  -- !
       write(6,'(2x,a)'      ) ' [ setupGS_DefNspRepr @setupGSMod ]'
       write(6,'(5x,a)'      ) '---  Density Max in Each PE ---'
       write(6,'(8x,a,i12)'  ) '  iPEMax  ==  ', iPEMax
       write(6,'(8x,a,f12.5)') '  rhoMax  ==  ', rhoMax
       write(6,*)
       write(6,'(2x,3(a8,1x),3(a12,1x))') 'iPE', 'iMax', 'jMax', 'rhoMax', 'dnsMax', 'rwVMax'
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       do iPE=0, PEtot-1
          write(6,'(2x,3(i8,1x),3(f12.5,1x))') iPE, nint( max_rBuff(1,iPE) ), nint( max_rBuff(2,iPE) ), &
               &                               max_rBuff(3,iPE), max_rBuff(4,iPE), max_rBuff(5,iPE)
       enddo
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,*)
       !  -- [3-2] Display Info.    --  !
       write(6,*)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(5x,a)') '---  Particle Distribution  ---'
       write(6,'(8x,a,i12  )') '  SuperParticle is Max at         :: iMax       :: ',   iMax
       write(6,'(8x,a,i12  )') '                                  :: jMax       :: ',   jMax
       write(6,'(8x,a,f12.5)') '                                  :: max(w*Np)  :: ', NspMax
       write(6,'(8x,a,f12.5)') '                   at the point   :: rhoMax     :: ', rhoMax
       write(6,*)
       write(6,'(5x,a)') '---  Determination of Np0   ---'
       write(6,'(8x,a,f12.5)') '  RhoPerNsp  :: Phys.Pt. Density / SuperPaticle :: ', RhoPerNsp
       write(6,'(8x,a,f12.5)') '  NspRepres  :: # of superparticle at n=1 in GS :: ', NspRepres
       write(6,'(8x,a,f12.5)') '  volRepres  :: Volume of superparticle at LJ/2 :: ', volRepres
       write(6,*)
       write(6,'(5x,a)') '---  Distributed Density    ---'
       write(6,'(8x,a,f12.5)') '  max(ne+ni)                                    :: ', neMax + niMax
       write(6,'(8x,a,f12.5)') '  max(ne   )                                    :: ', neMax
       write(6,'(8x,a,f12.5)') '  max(   ni)                                    :: ', niMax
       write(6,*)
    endif
    !   -- [3-3] Display Nsp example --  !
    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    if ( myRank.eq.PEtot/2 ) then
       write(6,*)
       write(6,'(5x,a)') '---  SuperParticle @ Center ---'
       write(6,'(5x,a,i12  )') '  Np(j=   2) :: # of superparticle at  j =    2 :: ', NpInner
       write(6,'(5x,a,i12  )') '  Np(j=LJ/2) :: # of superparticle at  j = LJ/2 :: ', NpMid  
       write(6,'(5x,a,i12  )') '  Np(j=LJ-1) :: # of superparticle at  j = LJ-1 :: ', NpOuter
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,*)
    endif
    
    return
  end subroutine setupGS_DefNspRepr


  subroutine setupGS_setJeiueiE
    use constants , only : LIs, LJs, ns, myRank, PEtot
    use constants , only : q, prtrb_coef, valfe, mr, InitJiJe
    use constants , only : br_, bt_, bz_, er_, et_, ez_, jr_, jt_, jz_, el_, io_
    use variables , only : Jcr, EBr, EBf, EBo, w_p
    use momentsMod, only : nCountUp_CC
    use mxwBdrcMod, only : BoundaryE
    implicit none
    include 'mpif.h'
    integer             :: i, j, k, ierr
    double precision    :: vprtrb, qInv(2)
    double precision    :: uetMin_s, uetMax_s, uitMin_s, uitMax_s
    double precision    :: uetMin  , uetMax  , uitMin  , uitMax
    
    !   --- [1] Electron / Ion Current  --- !
    !$omp parallel default(none) &
    !$omp shared(Jcr,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          ! - Electron - !
          Jcr(jr_,i,j,1) = ( 1.d0 - InitJiJe ) * Jcr(jr_,i,j,0)
          Jcr(jt_,i,j,1) = ( 1.d0 - InitJiJe ) * Jcr(jt_,i,j,0)
          Jcr(jz_,i,j,1) = ( 1.d0 - InitJiJe ) * Jcr(jz_,i,j,0)
          ! - Ion      - !
          Jcr(jr_,i,j,2) =          InitJiJe   * Jcr(jr_,i,j,0)
          Jcr(jt_,i,j,2) =          InitJiJe   * Jcr(jt_,i,j,0)
          Jcr(jz_,i,j,2) =          InitJiJe   * Jcr(jz_,i,j,0)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    !   --- [2] Flow Velocity  (u,v,w) ==> EBr --- !
    call nCountUp_CC( w_p )
    qInv(1) = 1.d0 / q(1)
    qInv(2) = 1.d0 / q(2)
    !$omp parallel default(none) &
    !$omp shared(q,Jcr,EBr,w_p,LIs,LJs,qInv) private(i,j,k)
    !$omp do
    do j=2, LJs
       do i=2, LIs
          EBr(uer_,i,j) = qInv(1) * Jcr(jr_,i,j,1) / ( ( w_p(i,j,1) + w_p(i-1,j  ,1) ) * 0.5d0 )
          EBr(uet_,i,j) = qInv(1) * Jcr(jt_,i,j,1) /     w_p(i,j,1)
          EBr(uez_,i,j) = qInv(1) * Jcr(jz_,i,j,1) / ( ( w_p(i,j,1) + w_p(i  ,j-1,1) ) * 0.5d0 )
          EBr(uir_,i,j) = qInv(2) * Jcr(jr_,i,j,2) / ( ( w_p(i,j,2) + w_p(i-1,j  ,2) ) * 0.5d0 )
          EBr(uit_,i,j) = qInv(2) * Jcr(jt_,i,j,2) /     w_p(i,j,2)
          EBr(uiz_,i,j) = qInv(2) * Jcr(jz_,i,j,2) / ( ( w_p(i,j,2) + w_p(i  ,j-1,2) ) * 0.5d0 )
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ++++++++++++++++++++++++++++++ !
    ! EBr = 0.d0
    ! ++++++++++++++++++++++++++++++ !
    
    !   --- [3] Add Purturbation in Flow  --- !
    vprtrb  = prtrb_coef * valfe / sqrt( mr )
    !$omp parallel default(none) &
    !$omp shared(EBr,psi,vprtrb,LIs,LJs) private(i,j,k)
    !$omp do
    do j=2, LJs
       do i=2, LIs
          if ( psi(i,j).le.0.d0 ) then
             if ( i.le.LIs/2 ) then
                EBr(uez_,i,j) = EBr(uez_,i,j) + vprtrb
                EBr(uiz_,i,j) = EBr(uiz_,i,j) + vprtrb
             else
                EBr(uez_,i,j) = EBr(uez_,i,j) - vprtrb
                EBr(uiz_,i,j) = EBr(uiz_,i,j) - vprtrb
             endif
          endif
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    EBr(:,    :,    1) = EBr(:,  :,  2)
    EBr(:,    :,LJs+1) = EBr(:,  :,LJs)
    EBr(:,    1,    :) = EBr(:,  2,  :)
    EBr(:,LIs+1,    :) = EBr(:,LIs,  :)
    uetMin_s           = minval( EBr(uet_,2:LIs-1,2:LJs-1) )
    uetMax_s           = maxval( EBr(uet_,2:LIs-1,2:LJs-1) )
    uitMin_s           = minval( EBr(uit_,2:LIs-1,2:LJs-1) )
    uitMax_s           = maxval( EBr(uit_,2:LIs-1,2:LJs-1) )
    call MPI_AllReduce( uetMin_s, uetMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( uetMax_s, uetMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( uitMin_s, uitMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
    call MPI_AllReduce( uitMax_s, uitMax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'( 2x,a)'      ) ' [ setupGS_setJeiueiE @setupGSMod ]'
       write(6,*)
       write(6,'(5x,a)'       ) '---  Perturbation Velocity  ---'
       write(6,'(8x,a,f12.5)') '  vprtrb                                        :: ', vprtrb
       write(6,*)
       write(6,'(5x,a)'       ) '---  Toroidal Drift Speed   ---'
       write(6,'(8x,2(a,1x,f12.6,3x))') 'vt(elc) :: min ==', uetMin, 'max ==', uetMax
       write(6,'(8x,2(a,1x,f12.6,3x))') 'vt(ion) :: min ==', uitMin, 'max ==', uitMax
       write(6,*)
    endif
    
    !   --- [4] E + ue x B = 0 ( For due/dt = 0 ) --- !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBr,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBf(er_,i,j) = - ( EBr(uet_,i,j)*EBf(bz_,i,j) - EBr(uez_,i,j)*EBf(bt_,i,j) )
          EBf(et_,i,j) = - ( EBr(uez_,i,j)*EBf(br_,i,j) - EBr(uer_,i,j)*EBf(bz_,i,j) )
          EBf(ez_,i,j) = - ( EBr(uer_,i,j)*EBf(bt_,i,j) - EBr(uet_,i,j)*EBf(br_,i,j) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    call BoundaryE
    
    !$omp parallel default(none) &
    !$omp shared(EBf,EBo,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBo(er_,i,j) = EBf(er_,i,j)
          EBo(et_,i,j) = EBf(et_,i,j)
          EBo(ez_,i,j) = EBf(ez_,i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine setupGS_setJeiueiE


  subroutine setupGS_ptVelocity
    use constants,     only : ns, np, npt, LIs, LJs
    use constants,     only : vthe, vthi, dz, dr, dzInv, drInv
    use constants,     only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables,     only : pxv, EBr
    use shapeFnMod,     only : shapef_r, shapef_z
    use randGenMod,    only : multi_gdev, gaussdev
    !$ use omp_lib
    implicit none
    integer                :: i, j, k, m, ith
    integer                :: iph, ipf, jph, jpf, iph_, ipf_, jph_, jpf_
    double precision       :: vth(ns), vdrf(3), sfrf(-2:2), sfzf(-2:2), sfrh(-2:2), sfzh(-2:2)

    !   --- [1]  Drift + Thermal + Bulk  ---  !
    vth(1)  = vthe
    vth(2)  = vthi
    do k=1, ns
       !$omp parallel default(none) &
       !$omp shared (k,np,pxv,dz,dzInv,dr,drInv,EBr,psi,vth,LIs,LJs ) &
       !$omp private(m,iph,jph,ipf,jpf,iph_,jph_,ipf_,jpf_,sfrf,sfzf,sfrh,sfzh,vdrf,ith)
       !$ ith = omp_get_thread_num() + 1
       !$omp do
       do m = 1, np(k)
          
          iph      = ceiling( pxv(zp_,m,k)*dzInv ) - 1
          jph      = ceiling( pxv(rp_,m,k)*drInv ) - 1
          ipf      =    nint( pxv(zp_,m,k)*dzInv )
          jpf      =    nint( pxv(rp_,m,k)*drInv )
          sfrf     = shapef_r( jpf, pxv(rp_,m,k)         , drInv )
          sfzf     = shapef_z( ipf, pxv(zp_,m,k)         , dzInv )
          sfrh     = shapef_r( jph, pxv(rp_,m,k)-0.5d0*dr, drInv )
          sfzh     = shapef_z( iph, pxv(zp_,m,k)-0.5d0*dz, dzInv )          
          iph      = iph + 2
          jph      = jph + 2
          ipf      = ipf + 2
          jpf      = jpf + 2
          if ( ( sfzf(-1).lt.0.d0 ).or.( sfzf(0).lt.0.d0 ).or.( sfzf(+1).lt.0.d0 ) ) then
             write(6,*) sfzf(:)
             stop
          endif
          if ( ( sfrf(-1).lt.0.d0 ).or.( sfrf(0).lt.0.d0 ).or.( sfrf(+1).lt.0.d0 ) ) then
             write(6,*) sfrf(:)
             stop
          endif
          if ( ( sfzh(-1).lt.0.d0 ).or.( sfzh(0).lt.0.d0 ).or.( sfzh(+1).lt.0.d0 ) ) then
             write(6,*) sfzh(:)
             stop
          endif
          if ( ( sfrh(-1).lt.0.d0 ).or.( sfrh(0).lt.0.d0 ).or.( sfrh(+1).lt.0.d0 ) ) then
             write(6,*) sfrh(:)
             stop
          endif
          vdrf(:) = 0.d0
          
          do j=-1,1
             do i=-1,1
                ipf_   = max( min( ipf+i, LIs ), 2 )
                jpf_   = max( min( jpf+j, LJs ), 2 )
                iph_   = max( min( iph+i, LIs ), 2 )
                jph_   = max( min( jph+j, LJs ), 2 )
                vdrf(1) = vdrf(1) + sfzh(i)*sfrf(j)*EBr(uer_+(k-1)*3,iph_,jpf_)
                vdrf(2) = vdrf(2) + sfzh(i)*sfrh(j)*EBr(uet_+(k-1)*3,iph_,jph_)
                vdrf(3) = vdrf(3) + sfzf(i)*sfrh(j)*EBr(uez_+(k-1)*3,ipf_,jph_)
             enddo
          enddo
          pxv(vr_,m,k) = vdrf(1) + vth(k)*multi_gdev( ith )
          pxv(vt_,m,k) = vdrf(2) + vth(k)*multi_gdev( ith )
          pxv(vz_,m,k) = vdrf(3) + vth(k)*multi_gdev( ith )
          ! -- Not Needed Gamma Correction -- !
          ! gamma    = 1.d0 / sqrt( 1.d0 - ( pxv(vr_,m,k)**2 + pxv(vt_,m,k)**2 + pxv(vz_,m,k)**2 ) )
          ! pxv(vr_,m) = gamma * pxv(vr_,m)
          ! pxv(vt_,m) = gamma * pxv(vt_,m)
          ! pxv(vz_,m) = gamma * pxv(vz_,m)
          ! -- -- !
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    !   --- [2]  Light Speed Check  --- !
    if ( flag__LightSpeedCheck ) then
       do k=1, ns
          do m=1, np(k)
             if ( ( pxv(vr_,m,k)**2 + pxv(vt_,m,k)**2 + pxv(vz_,m,k)**2 ).gt.1.d0 ) then
                write(6,*) ' [ERROR] Initial Velocity exceeds the light speed [ERROR] '
                write(6,'(a,1x,2(i8,1x))'   ) ' (ipf,jpf,iph,jph) = ', ipf, jpf, iph, jph
                write(6,'(a,1x,3(f10.5,1x))') ' vdrf    = ', vdrf
                write(6,'(a,1x,5(f10.5,1x))') ' sfzf    = ', sfzf
                write(6,'(a,1x,5(f10.5,1x))') ' sfrf    = ', sfrf
                write(6,'(a,1x,5(f10.5,1x))') ' sfzh    = ', sfzh
                write(6,'(a,1x,5(f10.5,1x))') ' sfrh    = ', sfrh
                write(6,'(a,1x,8(f10.5,1x))') ' r, z    = ', pxv(rp_:2,m,k)
                write(6,'(a,1x,8(f10.5,1x))') ' vrtz    = ', pxv(vr_:5,m,k)
                write(6,'(a,1x,8(f10.5,1x))') ' ro,zo,w = ', pxv(ro_:8,m,k)
                stop
             endif
          enddo
       enddo
    endif
    !   --- [3]  NaN detection :: ( x.eq.x ) => False for Nan  --- !
    !     - Nan Detection ( x.eq.x ) = False for Nan - !
    if ( flag__NaNDetection ) then
       do m=1, npt
          if ( ( pxv(rp_,m,k).ne.pxv(rp_,m,k) ).or.( pxv(zp_,m,k).ne.pxv(zp_,m,k) ) ) &
               write(6,*) ' [NaN] position', pxv(rp_:8,m,k)
          if ( ( pxv(vr_,m,k).ne.pxv(vr_,m,k) ).or.( pxv(vt_,m,k).ne.pxv(vt_,m,k) ).or.&
               & ( pxv(vz_,m,k).ne.pxv(vz_,m,k) ) ) &
               write(6,*) ' [NaN] velocity', pxv(rp_:8,m,k)
       enddo
    endif
    
    return
  end subroutine setupGS_ptVelocity


  subroutine setupGS_HalfDTStep
    use constants , only : dt, dtdr, drdt, dtdz, dzdt
    use ptAccelMod, only : position
    use maxwellMod, only : Efield
    use currentMod, only : DensityDecomposition
    use ptcBdrcMod, only : particleBoundary
    use qchargeMod, only : charge, Ecorrect
    use debugerMod, only : ampereconf
    implicit none
    
    ! --- [1] Half dt           --- !
    dt          = 0.5d0 * dt
    dtdr        = 0.5d0 * dtdr
    dtdz        = 0.5d0 * dtdz
    drdt        = 2.0d0 * drdt
    dzdt        = 2.0d0 * dzdt
    ! --- [2] PUSH Particles    --- !
    call position
    ! --- [3] Leap-Frog E-Field --- !
    call DensityDecomposition
    call particleBoundary
    ! call Efield
    ! call ampereconf
    ! --- [4] divE correction   --- !
    ! call charge
    ! call Ecorrect
    ! call charge
    ! --- [5] Restore Step Size --- !
    dt          = 2.0d0 * dt
    dtdr        = 2.0d0 * dtdr
    dtdz        = 2.0d0 * dtdz
    dzdt        = 0.5d0 * dzdt
    drdt        = 0.5d0 * drdt
    return
  end subroutine setupGS_HalfDTStep


  subroutine setupGS_BFadjustDT
    use constants, only : LIs, LJs, cv, wpewce, drInv, dzInv, myRank
    use constants, only : alpha_wce, alpha_wpe, alpha_CFL
    use constants, only : dt, dr, dz, dtdr, drdt, dtdz, dzdt
    use constants, only : br_, bt_, bz_
    use variables, only : EBf
    implicit none
    include 'mpif.h'
    integer            :: i, j, ierr
    double precision   :: dt_CFL, dt_wpe, dt_wce
    double precision   :: absB, absB2, DTcMax

    ! ------------------------------------ !
    ! --- [1] Obtain  Max( wce )       --- !
    ! ------------------------------------ !
    absB2 = - 1.d3
    do j=2, LJs-1
       do i=2, LIs-1
          absB2 = max( absB2, ( EBf(br_,i,j)**2 + EBf(bt_,i,j)**2 + EBf(bz_,i,j)**2 ) )
       enddo
    enddo
    call MPI_AllReduce( absB2, absB, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    absB    =  sqrt( absB )
    DTcMax  = 1.d0 / absB
    
    ! ------------------------------------ !
    ! --- [2] Redefine Time Step :: dt --- !
    ! ------------------------------------ !
    dt_wce  = alpha_wce          * DTcMax
    dt_wpe  = alpha_wpe / wpewce * DTcMax
    dt_CFL  = alpha_CFL / ( cv * sqrt( drInv**2 + dzInv**2 ) )
    dt      = min( dt_CFL, dt_wpe, dt_wce )
    dtdr    =   dt / dr
    drdt    =   dr / dt
    dtdz    =   dt / dz
    dzdt    =   dz / dt

    ! ------------------------------------ !
    ! --- [3] Display            :: dt --- !
    ! ------------------------------------ !
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'( 2x,a)') '----------------------------------------------------------------------'
       write(6,'(10x,a)'            ) '---    Time  Step  ( Re-Definition )   ---'
       write(6,'(10x,a,1(f12.6,1x))') '      B (   Max   )        :: ', absB
       write(6,'(10x,a,1(f12.6,1x))') '     dt (   CFL   )        :: ', dt_CFL
       write(6,'(10x,a,1(f12.6,1x))') '     dt (  1/wpe  )        :: ', dt_wpe
       write(6,'(10x,a,1(f12.6,1x))') '     dt (  1/wce  )        :: ', dt_wce
       write(6,'(10x,a)'            ) '---------------------------------------------'
       write(6,'(10x,a,1(f12.6,1x))') '       dt:PIC              :: ', dt
       write(6,*)
       write(6,'( 2x,a)') '----------------------------------------------------------------------'
    endif
    return
  end subroutine setupGS_BFadjustDT
  
  subroutine setupGS_NullField
    use variables, only : EBf, Jcr, EBo
    implicit none
    ! -- Set Null Field -- !
    EBf = 0.d0
    Jcr = 0.d0
    EBo = 0.d0
    rho = 1.d0
    EBf(6,:,:) = 1.d0
    ! -- End -- !
    return
  end subroutine setupGS_NullField
  
end module setupGSMod




!  -- [1-2] dnsMax / rhoMax / Index -- !
! iMax      = LIs / 2
! jMax      = LJs / 2
! dnsMax    = w_p(iMax,jMax,1)
! do j=1, LJs
!    do i=1, LIs
!       if ( w_p(i,j,1).gt.dnsMax ) then
!          dnsMax = w_p(i,j,1)
!          iMax   = i
!          jMax   = j
!       endif
!    enddo
! enddo
