module sFluxIpMod
  implicit none
  integer                       :: LIm, LJm, LIm_
  integer                       :: LatestFluxUpdate = -1
  integer                       :: ptOXO(3,2)
  logical                       :: Flag__sOpt, Flag__FixedXpt
  logical                       :: InitOXO   = .false.
  double precision              :: FluxSign  =  - 1.d0
  double precision              :: psiNX0    =  - 1.d0
  double precision              :: psiO1, psiO2  , psiFX
  double precision              :: rOpt1, rOpt2  , rXpt
  double precision              :: zOpt1, zOpt2  , zXpt
  double precision              :: psiNX, Iplasma, psiLim, psiMin, psiAmp, mergR, RateR
  integer         , parameter   :: id_=1, jd_=2
  double precision, allocatable :: psiN(:,:), psiF(:,:), psiNg(:,:), psiFg(:,:), recvBuff(:,:,:)
contains

  subroutine FluxCurrentAnalysis
    use diagnosMod, only : AngularMomentum
    implicit none

    call polFluxMap
    call AngularMomentum
    call xoPointAnalysis
    call sumFluxCurrent
    
    return
  end subroutine FluxCurrentAnalysis

  
  subroutine polFluxMap
    use constants , only         : LI , LJ , LIs, LJs, dr, myRank, PEtot, jobDir
    use constants , only         : br_, bt_, bz_, inn_, iFr_, isl_, LIs_
    use variables , only         : EBf, rf , kstep, ijDomain
    use getXOptMod, only         : solveXOpt
    implicit none
    include 'mpif.h'
    integer                     :: i, j, ir, ip, Nitem, iPE, ierr
    integer                     :: isOpt=0, iFixed=0
    integer                     :: request_Send(1), status_Send(MPI_STATUS_SIZE,1)
    integer                     :: request_Recv(PEtot), status_Recv(MPI_STATUS_SIZE,PEtot)
    character(100)              :: FileName
    double precision            :: coef, psiMin_s
    
    ! ----------------------- !
    ! --- [1] Preparation --- !
    ! ----------------------- !
    if ( .not.( InitOXO ) ) then
       ! -- [1-1] Obtain [ max(LI), max(LJ) ] -- !
       LIm = 0
       do iPE=0, PEtot-1
          LIm = max( ijDomain(iPE, LIs_), LIm )
       enddo
       LJm  = LJs
       LIm_ = LIm
       allocate( psiF(LIm,LJm), psiN(LIm,LJm), psiNg(LI,LJ), recvBuff(LIm,LJm,PEtot) )
       ! -- [1-2] Initial Estimation of O/X   -- !
       ptOXO(1,id_)   = nint( LI*0.33d0 )
       ptOXO(2,id_)   = nint( LI*0.50d0 )
       ptOXO(3,id_)   = nint( LI*0.67d0 )
       ptOXO(1,jd_)   = nint( LJ*0.50d0 )
       ptOXO(2,jd_)   = nint( LJ*0.50d0 )
       ptOXO(3,jd_)   = nint( LJ*0.50d0 )
       ! -- [1-3] Preparation of variables   -- !
       psiF           = 0.d0
       psiN           = 0.d0
       psiNg          = 0.d0
       recvBuff       = 0.d0
       Flag__sOpt     = .false.
       Flag__FixedXpt = .false.
       InitOXO        = .true.
    endif
    ! -- [1-4] update size if needed -- !
    LIm = 0
    do iPE=0, PEtot-1
       LIm = max( ijDomain( iPE, LIs_ ), LIm )
    enddo
    if ( LIm.ne.LIm_ ) then
       LIm_ = LIm
       deallocate( psiN, psiF, recvBuff )
       allocate  ( psiN(LIm,LJm), psiF(LIm,LJm), recvBuff(LIm,LJm,PEtot) )
       psiF           = 0.d0
       psiN           = 0.d0
       recvBuff       = 0.d0
    endif
    
    ! -------------------------------------------- !
    ! --- [2] Calculate Magnetic Flux Function --- !
    ! -------------------------------------------- !
    coef = FluxSign*dr
    !  -- [2-1] Magnetic Flux Function -- !
    psiF(:,:) = 0.d0
    do i=1, LIs
       psiF(i,2)    = 0.5d0*coef*EBf(bz_,i,2)
    enddo
    do j=3, LJs
       do i=1, LIs
          psiF(i,j) = psiF(i,j-1) + coef*rf(j)*EBf(bz_,i,j)
       enddo
    enddo
    !  -- [2-2] Magnetic Axis and LCFS -- !
    if      ( FluxSign.gt.0.d0 ) then
       psiMin_s = - 1.d8
       !$omp parallel default(none) shared(psiMin_s,psiF,LIs,LJs) private(i,j)
       !$omp do reduction(max:psiMin_s)
       do j=2, LJs
          do i=1, LIs
             psiMin_s = max( psiF(i,j), psiMin_s )
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       call MPI_AllReduce( psiMin_s, psiMin, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
    else if ( FluxSign.lt.0.d0 ) then
       psiMin_s = 1.d8
       !$omp parallel default(none) shared(psiMin_s,psiF,LIs,LJs) private(i,j)
       !$omp do reduction(min:psiMin_s)
       do j=2, LJs
          do i=1, LIs
             psiMin_s = min( psiF(i,j), psiMin_s )
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       call MPI_AllReduce( psiMin_s, psiMin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr )
    endif
    !  -- [2-3] Magnetic Axis and LCFS -- !
    psiLim = 0.d0
    psiAmp = psiLim - psiMin
    coef   = 1.d0 / psiAmp
    !$omp parallel default(none) &
    !$omp shared(coef,psiLim,psiF,psiN,LIs,LJs) private(i,j)
    !$omp do
    do j=2, LJs
       do i=1, LIs
          psiN(i,j)   = coef * ( psiLim - psiF(i,j) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! ----------------------------------- !
    ! --- [3] Communication  ( psi )  --- !
    ! ----------------------------------- !
    !  -- [3-1] Gather Data           --  !
    Nitem  = LIm*LJm
    call MPI_Isend( psiN,    Nitem, MPI_DOUBLE_PRECISION,  0, 0, &
         &          MPI_COMM_WORLD, request_Send        ,  ierr   )
    if ( myRank.eq.0 ) then
       do iPE=1, PEtot
          call MPI_Irecv( recvBuff(:,:,iPE), Nitem, MPI_DOUBLE_PRECISION, iPE-1, 0, &
               &          MPI_COMM_WORLD   , request_Recv(iPE)          , ierr )
          call MPI_WaitAll( 1, request_Recv(iPE), status_Recv(:,iPE), ierr )
       enddo
    endif
    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    call MPI_WaitAll( 1, request_Send, status_Send, ierr )
    !  -- [3-2] Storing Data          --  !
    psiNg(:,:) = 0.d0
    if ( myRank.eq.0 ) then
       !$omp parallel default(none) &
       !$omp shared(ijDomain,psiNg,recvBuff) private(iPE,j,i,ip,ir)
       !$omp do
       do iPE=0, PEtot-1
          do j=1, LJ
             do i=1, ijDomain(iPE,inn_)
                ip          = ijDomain(iPE,iFr_) + ( i-1 )
                ir          = ijDomain(iPE,isl_) + ( i-1 )
                psiNg(ip,j) = recvBuff(ir,j,iPE+1)
             enddo
          enddo
       enddo
       !$omp enddo
       !$omp end parallel
    endif
    
    ! ----------------------------------- !
    ! --- [4] Obtain (i,j) of Xpt,Opt --- !
    ! ----------------------------------- !
    if ( myRank.eq.0 ) call solveXOpt( psiNg, ptOXO, LI, LJ, Flag__sOpt, Flag__FixedXpt )
    
    ! ------------------------------ !
    ! --- [5]   Write in File    --- !
    ! ------------------------------ !
    !  -- X- and O-point Info.   --  !
    FileName  = trim(jobDir) // 'dat/' // 'XOpoint.dat'
    if ( Flag__sOpt     ) isOpt  = 1
    if ( Flag__FixedXpt ) iFixed = 1
    if ( myRank.eq.0 ) then
       open(50,file=trim(FileName),form='formatted',access='append')
       write(50,'(9(i10,1x))') kstep       , isOpt       , iFixed      , &
            &                  ptOXO(1,id_), ptOXO(2,id_), ptOXO(3,id_), &
            &                  ptOXO(1,jd_), ptOXO(2,jd_), ptOXO(3,jd_)
       close(50)
    endif
    
    return
  end subroutine polFluxMap

  
  subroutine xoPointAnalysis
    use constants , only         : dt, myRank, jobDir
    use constants , only         : jt_, bt_
    use variables , only         : kstep, rh, zh
    use recRateMod, only         : ReconnectionRate
    implicit none
    include 'mpif.h'
    integer                     :: ptOXO_(3,2), ierr
    double precision            :: isOpt=0.d0, iFixed=0.d0
    double precision            :: time, psiNX_, mergR_
    character(100)              :: FileName

    ! -------------------------- !
    ! --- [1]  Preparation   --- !
    ! -------------------------- !
    time    = dble(kstep) * dt
    ptOXO_  = ptOXO
    if ( Flag__sOpt     ) isOpt  = 1.d0
    if ( Flag__FixedXpt ) iFixed = 1.d0

    ! -------------------------- !
    ! --- [2]  X,O-point     --- !
    ! -------------------------- !
    if ( psiNX0.lt.0.d0 ) then
       ! --  =(undefined)  -- !
       if ( myRank.eq.0 ) psiNX0 = psiNg( ptOXO(2,id_), ptOXO(2,jd_) )
       call MPI_Bcast( psiNX0, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
       mergR    = 0.d0
    endif
    if ( myRank.eq.0 ) then
       ! -- try psiNX      -- !
       psiNX_   = psiNg( ptOXO(2,id_), ptOXO(2,jd_) )
       mergR_   = ( psiNX_ - psiNX0 ) / ( 1.d0 - psiNX0 )
       ! -- No-Update case -- !
       if ( abs( mergR_-mergR ).gt.0.1d0 ) ptOXO = ptOXO_
       if (      mergR_        .lt.0.0d0 ) ptOXO = ptOXO_
    endif
    call MPI_Bcast( ptOXO         , 6, MPI_INTEGER, 0,  MPI_COMM_WORLD, ierr )
    call MPI_Bcast( Flag__sOpt    , 1, MPI_LOGICAL, 0,  MPI_COMM_WORLD, ierr )
    call MPI_Bcast( Flag__FixedXpt, 1, MPI_LOGICAL, 0,  MPI_COMM_WORLD, ierr )
    psiO1    = psiNg( ptOXO(1,id_), ptOXO(1,jd_) ) * psiAmp
    psiFX    = psiNg( ptOXO(2,id_), ptOXO(2,jd_) ) * psiAmp
    psiO2    = psiNg( ptOXO(3,id_), ptOXO(3,jd_) ) * psiAmp
    psiNX    = psiNg( ptOXO(2,id_), ptOXO(2,jd_) )
    rXpt     =    rh( ptOXO(2,jd_) )
    rOpt1    =    rh( ptOXO(1,jd_) )
    rOpt2    =    rh( ptOXO(3,jd_) )
    zXpt     =    zh( ptOXO(2,id_) )
    zOpt1    =    zh( ptOXO(1,id_) )
    zOpt2    =    zh( ptOXO(3,id_) )
    mergR    = ( psiNX - psiNX0 ) / ( 1.d0 - psiNX0 )
    call ReconnectionRate( ptOXO, RateR, mergR, Flag__sOpt )

    ! -------------------------- !
    ! --- [3] Write in File  --- !
    ! -------------------------- !
    ! -- phi, psi WriteOut -- !
    FileName = trim(jobDir) // 'dat/' // 'MagneticFlux.dat'
    if ( myRank.eq.0 ) then
       open(50,file=trim(FileName),form='formatted',access='append')
       write(50,'(16(e15.8,1x))') time , psiO1, psiO2 , psiFX , psiNX   , psiMin, &
            &                     rOpt1, zOpt1, rOpt2 , zOpt2 , rXpt    , zXpt  , &
            &                     mergR, RateR, isOpt , iFixed
       close(50)
    endif
    
    return
  end subroutine xoPointAnalysis

  
  subroutine sumFluxCurrent
    use constants, only          : LIs, LJs, dt, dr, dz, jobDir, myRank
    use constants, only          : jt_, bt_
    use variables, only          : Jcr, EBf, w_p, x1s, kstep
    implicit none
    include 'mpif.h'
    integer                     :: i, j, ierr
    double precision            :: dS, time
    double precision            :: Ie       , Ii       ,            phi       , aMe       , aMi
    double precision            :: Ie95     , Ii95     , Ip95     , phi95     , aMe95     , aMi95
    double precision            :: IeSim    , IiSim    , IpSim    , phiSim    , aMeSim    , aMiSim
    double precision            :: IePriv1  , IiPriv1  , IpPriv1  , phiPriv1  , aMePriv1  , aMiPriv1
    double precision            :: IePriv2  , IiPriv2  , IpPriv2  , phiPriv2  , aMePriv2  , aMiPriv2
    double precision            :: Ie_s     , Ii_s     ,            phi_s     , aMe_s     , aMi_s
    double precision            :: Ie95_s   , Ii95_s   ,            phi95_s   , aMe95_s   , aMi95_s
    double precision            :: IeSim_s  , IiSim_s  ,            phiSim_s  , aMeSim_s  , aMiSim_s
    double precision            :: IePriv1_s, IiPriv1_s,            phiPriv1_s, aMePriv1_s, aMiPriv1_s
    double precision            :: IePriv2_s, IiPriv2_s,            phiPriv2_s, aMePriv2_s, aMiPriv2_s
    character(100)              :: FileName
    double precision, parameter :: LCFS95  = 0.05d0
    double precision, parameter :: LCFS    = 0.00d0
    
    ! ------------------------------------- !
    ! --- [1] Plasma Current ( < LCFS ) --- !
    ! ------------------------------------- !
    !  --  Ip, Ip95  --  !
    Ie95_s    = 0.d0
    Ii95_s    = 0.d0
    phi95_s   = 0.d0
    aMe95_s   = 0.d0
    aMi95_s   = 0.d0
    Ie_s      = 0.d0
    Ii_s      = 0.d0
    phi_s     = 0.d0
    aMe_s     = 0.d0
    aMi_s     = 0.d0
    dS        = dr*dz
    !$omp parallel default(none) &
    !$omp shared(LIs,LJs,psiN,Ie95_s,Ii95_s,phi95_s,aMe95_s,aMi95_s) &
    !$omp shared(Jcr,EBf,w_p,dS,Ie_s,Ii_s,phi_s,aMe_s,aMi_s) private(i,j)
    !$omp do reduction(+:Ie95_s,Ii95_s,phi95_s,aMe95_s,aMi95_s,Ie_s,Ii_s,phi_s,aMe_s,aMi_s)
    do j=2, LJs-1
       do i=2, LIs-1
          if ( psiN(i,j).ge.LCFS95 ) then
             Ie95_s  = Ie95_s  + Jcr(jt_,i,j,1) * dS
             Ii95_s  = Ii95_s  + Jcr(jt_,i,j,2) * dS
             phi95_s = phi95_s + EBf(bt_,i,j)   * dS
             aMe95_s = aMe95_s + w_p(i,j,1)
             aMi95_s = aMi95_s + w_p(i,j,2)             
          end if
          if ( psiN(i,j).ge.LCFS   ) then
             Ie_s    = Ie_s    + Jcr(jt_,i,j,1) * dS
             Ii_s    = Ii_s    + Jcr(jt_,i,j,2) * dS
             phi_s   = phi_s   + EBf(bt_,i,j)   * dS
             aMe_s   = aMe_s   + w_p(i,j,1)
             aMi_s   = aMi_s   + w_p(i,j,2)             
          end if
       enddo
    enddo
    !$omp enddo
    !$omp end parallel
    call MPI_Reduce(  Ie95_s,  Ie95, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce(  Ii95_s,  Ii95, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce( phi95_s, phi95, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce( aMe95_s, aMe95, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce( aMi95_s, aMi95, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce(    Ie_s,    Ie, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce(    Ii_s,    Ii, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce(   phi_s,   phi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce(   aMe_s,   aMe, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce(   aMi_s,   aMi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    Ip95    = Ie95  + Ii95
    Iplasma = Ie    + Ii

    ! ------------------------------------- !
    ! --- [2] All Current ( < Sim.Box ) --- !
    ! ------------------------------------- !
    !  --  Ipsim  --  !
    IeSim_s  = 0.d0
    IiSim_s  = 0.d0
    phiSim_s = 0.d0
    aMeSim_s = 0.d0
    aMiSim_s = 0.d0
    !$omp parallel default(none) &
    !$omp shared(LIs,LJs,Jcr,EBf,w_p,dS,IeSim_s,IiSim_s,phiSim_s,aMeSim_s,aMiSim_s)
    !$omp do reduction(+:IeSim_s,IiSim_s,phiSim_s,aMeSim_s,aMiSim_s)
    do j=2, LJs
       do i=2, LIs
          IeSim_s    = IeSim_s  + Jcr(jt_,i,j,1) * dS
          IiSim_s    = IiSim_s  + Jcr(jt_,i,j,2) * dS
          phiSim_s   = phiSim_s + EBf(bt_,i,j)   * dS
          aMeSim_s   = aMeSim_s + w_p(i,j,1)
          aMiSim_s   = aMiSim_s + w_p(i,j,2)             
       enddo
    enddo
    !$omp enddo
    !$omp end parallel
    call MPI_Reduce(  IeSim_s,  IeSim, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce(  IiSim_s,  IiSim, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce( phiSim_s, phiSim, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce( aMeSim_s, aMeSim, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_Reduce( aMiSim_s, aMiSim, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    IpSim = IeSim + IiSim
    
    ! ------------------------------------- !
    ! --- [3] Private Current ( < Xpt ) --- !
    ! ------------------------------------- !
    !  --  Ippriv  --  !
    if ( .not.(Flag__sOpt) ) then
       IePriv1_s  = 0.d0
       IiPriv1_s  = 0.d0
       phiPriv1_s = 0.d0
       aMePriv1_s = 0.d0
       aMiPriv1_s = 0.d0
       IePriv2_s  = 0.d0
       IiPriv2_s  = 0.d0
       phiPriv2_s = 0.d0
       aMePriv2_s = 0.d0
       aMiPriv2_s = 0.d0
       !$omp parallel default(none) &
       !$omp shared(LIs,LJs,Jcr,EBf,w_p,dS,IePriv1_s,IiPriv1_s,phiPriv1_s,aMePriv1_s,aMiPriv1_s) &
       !$omp shared(psiN,psiNX,x1s,zXpt) private(i,j)
       !$omp do reduction(+:IePriv1_s,IiPriv1_s,phiPriv1_s,aMePriv1_s,aMiPriv1_s)
       do j=2, LJs-1
          do i=2, LIs-1
             if ( ( psiN(i,j).ge.psiNX ).and.( x1s(i).lt.zXpt ) ) then
                IePriv1_s  = IePriv1_s  + Jcr(jt_,i,j,1) * dS
                IiPriv1_s  = IiPriv1_s  + Jcr(jt_,i,j,2) * dS
                phiPriv1_s = phiPriv1_s + EBf(bt_,i,j)   * dS
                aMePriv1_s = aMePriv1_s + w_p(i,j,1)
                aMiPriv1_s = aMiPriv1_s + w_p(i,j,2)             
             endif
          enddo
       enddo
       !$omp enddo
       !$omp end parallel
       !$omp parallel default(none) &
       !$omp shared(LIs,LJs,Jcr,EBf,w_p,dS,IePriv2_s,IiPriv2_s,phiPriv2_s,aMePriv2_s,aMiPriv2_s) &
       !$omp shared(psiN,psiNX,x1s,zXpt) private(i,j)
       !$omp do reduction(+:IePriv2_s,IiPriv2_s,phiPriv2_s,aMePriv2_s,aMiPriv2_s)
       do j=2, LJs-1
          do i=2, LIs-1
             if ( ( psiN(i,j).ge.psiNX ).and.( x1s(i).gt.zXpt ) ) then
                IePriv2_s  = IePriv2_s  + Jcr(jt_,i,j,1) * dS
                IiPriv2_s  = IiPriv2_s  + Jcr(jt_,i,j,2) * dS
                phiPriv2_s = phiPriv2_s + EBf(bt_,i,j)   * dS
                aMePriv2_s = aMePriv2_s + w_p(i,j,1)
                aMiPriv2_s = aMiPriv2_s + w_p(i,j,2)             
             endif
          enddo
       enddo
       !$omp enddo
       !$omp end parallel
       call MPI_Reduce(  IePriv1_s,  IePriv1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce(  IiPriv1_s,  IiPriv1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce( phiPriv1_s, phiPriv1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce( aMePriv1_s, aMePriv1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce( aMiPriv1_s, aMiPriv1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce(  IePriv2_s,  IePriv2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce(  IiPriv2_s,  IiPriv2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce( phiPriv2_s, phiPriv2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce( aMePriv2_s, aMePriv2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
       call MPI_Reduce( aMiPriv2_s, aMiPriv2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    else
       IePriv1  = 0.d0
       IiPriv1  = 0.d0
       IePriv2  = 0.d0
       IiPriv2  = 0.d0
       phiPriv1 = 0.d0
       phiPriv2 = 0.d0
       aMePriv1 = 0.d0
       aMePriv2 = 0.d0
       aMiPriv1 = 0.d0
       aMiPriv2 = 0.d0       
    endif
    !   - Ip -  !
    IpPriv1  = IePriv1 + IiPriv1
    IpPriv2  = IePriv2 + IiPriv2

    ! ------------------------------ !
    ! --- [4]   Write in File    --- !
    ! ------------------------------ !
    time    = dble(kstep) * dt
    if ( myRank.eq.0 ) then
       ! -- phi, psi WriteOut -- !
       FileName = trim(jobDir) // 'dat/' // 'Current.dat'
       open(50,file=trim(FileName),form='formatted',access='append')
       write(50,'(16(e15.8,1x))') time   , &
            &                     Ie     , Ie95, IePriv1, IePriv2, IeSim, &
            &                     Ii     , Ii95, IiPriv1, IiPriv2, IiSim, &
            &                     Iplasma, Ip95, IpPriv1, IpPriv2, IpSim
       close(50)
       ! -- phi, psi WriteOut -- !
       FileName = trim(jobDir) // 'dat/' // 'MagneticFlux.dat'
       open(50,file=trim(FileName),form='formatted',access='append')
       write(50,'(6(e15.8,1x))')  time , phi, phi95, phiPriv1, phiPriv2, phiSim
       close(50)
       ! -- Angular Momentum  -- !
       FileName = trim(jobDir) // 'dat/' // 'AngularMomentum.dat'
       open(50,file=trim(FileName),form='formatted',access='append')
       write(50,'(11(e15.8,1x))') time, &
            &                     aMe , aMe95, aMePriv1, aMePriv2, aMeSim, &
            &                     aMi , aMi95, aMiPriv1, aMiPriv2, aMiSim  
       close(50)

    endif
    
    return
  end subroutine sumFluxCurrent


end module sFluxIpMod


    ! call MPI_WaitAll( PEtot, request_Recv, status_Recv, ierr )
    ! write(6,*) myRank, 'here0'
    ! call MPI_Gather ( psiN(1,1)      , Nitem,  MPI_DOUBLE_PRECISION, &
    !      &            recvBuff(1,1,1), Nitem,  MPI_DOUBLE_PRECISION, &
    !      &            0              , MPI_COMM_WORLD, ierr )
    ! if ( myRank.eq.0 ) then
    !    write(6,*) 'here1'
    !    do iPE=0, PEtot-1
    !       write(rank_char,'(i2.2)') iPE
    !       open (50, file  ='psi'//rank_char//'.bin', form   ='unformatted',  &
    !            &    status='replace'               , convert='LITTLE_ENDIAN' )
    !       write(50) recvBuff(1:LIm,1:LJm,iPE+1)
    !       close(50)
    !    enddo
    ! endif
    ! if ( myRank.eq.0 ) then
    !    open (50, file  ='psi.bin', form   ='unformatted',  &
    !         &    status='replace', convert='LITTLE_ENDIAN'  )
    !    write(50) psiNg(1:LI,1:LJ)
    !    close(50)
    ! endif
    ! call MPI_Finalize( ierr )
    ! stop
    
    !     !  -- [3-1]  Preparation          --  !
    ! if      ( myRank.eq.0       ) then
    !    Nsend = ( ijDomain(myRank,4)-1 ) * ijDomain(myRank,5)
    ! else if ( myRank.eq.PEtot-1 ) then
    !    Nsend = ( ijDomain(myRank,4)-1 ) * ijDomain(myRank,5)
    ! else
    !    Nsend = ( ijDomain(myRank,4)-2 ) * ijDomain(myRank,5)
    ! endif
    ! !  -- [3-2]  NRecv & displs       --  !
    ! Nrecvs(0)       = ( ijDomain(      0,4) - 1 ) * ijDomain(0,5)
    ! displs(0)       = 1
    ! do iPE=1, PEtot-2
    !    Nrecvs(iPE)  = ( ijDomain(iPE,4)-2 )*( ijDomain(iPE,5) )
    !    displs(iPE)  = displs(iPE-1) + Nrecvs(iPE)
    ! enddo
    ! Nrecvs(PEtot-1) = ( ijDomain(PEtot-1,4) - 1 ) * ijDomain(PEtot-1,5)
    ! displs(PEtot-1) = displs(iPE-2) + Nrecvs(PEtot-1)
    ! !  -- [3-3] gatherV               --  !
    ! call MPI_GatherV( psiN    , Nsend ,         MPI_DOUBLE_PRECISION, &
    !      &            recvBuff, Nrecvs, displs, MPI_DOUBLE_PRECISION, &
    !      &            0       , MPI_COMM_WORLD, ierr )

    ! do j=1, LJs
    !    do i=1, LIs
    !       psiT(j,i) = psiN(i,j)
    !    enddo
    ! enddo
    ! Nsend     = ( ( ijDomain(myRank,3)-1 ) - (ijDomain(myRank,2)-1) + 1 ) * ijDomain(myRank,5)
    ! Nrecvs(0) = ( ( ijDomain(     0,3)   ) - (ijDomain(     0,2)-1) + 1 ) * ijDomain(myRank,5)

    ! do j=1, LJs
    !    do i=1, LIs-1
    !       psiNg(i,j) = psiN(i,j)
    !    enddo
    ! enddo
    ! do iPE=1, PEtot-1
    !    iFr   = ijDomain(iPE,2) + 1
    !    iTo   = ijDomain(iPE,3) - 1
    !    Nitem = (ijDomain(iPE,4)-1)*(ijDomain(iPE,5)-1)
    !    call MPI_iSend( psiN            , Nitem         , MPI_DOUBLE_PRECISION, &
    !         &                         0, MPI_COMM_WORLD, request_Send, ierr    )
    !    call MPI_iRecv( psiNg(iFr:iTo,:), Nitem         , MPI_DOUBLE_PRECISION, &
    !         &                         0, MPI_COMM_WORLD, request_Send, ierr    )
    !    call MPI_WaitAll( 1, request_Recv, status_Recv, ierr )
    !    call MPI_WaitAll( 1, request_Send, status_Send, ierr )
    ! enddo


  
  ! subroutine getXOpt( psi, LIp, LJp, iOXO, jOXO, Flag_OXO, Flag_sOpt )
  !   implicit none
  !   integer         , intent(in)    :: LIp, LJp
  !   double precision, intent(in)    :: psi(LIp,LJp)
  !   integer         , intent(inout) :: iOXO(3), jOXO(3)
  !   integer         , intent(inout) :: Flag_OXO(3), Flag_sOpt
  !   integer                         :: i, j, hLIp, iO1, jO1, iO2, jO2, iX, jX, check1, check2
  !   double precision                :: dpsidx1(LIp,LJp), dpsidx2(LIp,LJp)
  !   double precision                :: maxpsi, absBMin, absBxy
  !   logical                         :: Flag__NoXpt = .false.

  !   ! --------------- !
  !   ! --- O-point --- !
  !   ! --------------- !
  !   Flag_OXO(1:3) = 1
  !   if ( Flag_sOpt.eq.1 ) then
  !      ! -- Only 1 O-point exist -- !
  !      maxpsi = -1.d8
  !      iO1    = -1
  !      jO1    = -1
  !      do j=2, LJp-1
  !         do i=2, LIp-1
  !            if ( psi(i,j).gt.maxpsi ) then
  !               maxpsi = psi(i,j)
  !               iO1    = i
  !               jO1    = j
  !            endif
  !         enddo
  !      enddo
  !      iO2    = iO1
  !      jO2    = jO1
  !   endif
  !   if ( Flag_sOpt.eq.0 ) then
  !      ! -- 2 O-point exist -- !
  !      hLIp    = iOXO(2)
  !      !  -- search : Max( psi ) -- !
  !      !  --  (2~hLIp-1) -- !
  !      maxpsi = -1.d8
  !      iO1    = -1
  !      jO1    = -1
  !      do j=2, LJp-1
  !         do i=2, hLIp-1
  !            if ( psi(i,j).gt.maxpsi ) then
  !               maxpsi = psi(i,j)
  !               iO1    = i
  !               jO1    = j
  !            endif
  !         enddo
  !      enddo
  !      !  --  (hLIp~LIp-1) -- !
  !      maxpsi = -1.d8
  !      iO2    = -1
  !      jO2    = -1
  !      do j=2, LJp-1
  !         do i=hLIp,LIp-1
  !            if ( psi(i,j).gt.maxpsi ) then
  !               maxpsi = psi(i,j)
  !               iO2    = i
  !               jO2    = j
  !            endif
  !         enddo
  !      enddo
  !   endif
  !   !  -- check -- !
  !   if ( ( iO1.eq.-1 ).or.( jO1.eq.-1 ) ) then
  !      iO1 = nint( hLIp *0.5d0 )
  !      jO1 = nint(  LJp *0.5d0 )
  !      Flag_OXO(1) = 0
  !   endif
  !   if ( ( iO2.eq.-1 ).or.( jO2.eq.-1 ) ) then
  !      iO2 = nint( hLIp *0.5d0 ) + hLIp
  !      jO2 = nint(  LJp *0.5d0 )
  !      Flag_OXO(3) = 0
  !   endif
  !   !  -- whether Extremum Max. Value -- !
  !   !   - Left - !
  !   if ( Flag_OXO(1).eq.1 ) then
  !      check1 = 0
  !      if ( ( psi(iO1+1,jO1  ) - psi(iO1,jO1) )*( psi(iO1,jO1) - psi(iO1-1,jO1  ) ).gt.0.d0 ) check1=check1 + 1
  !      if ( ( psi(iO1  ,jO1+1) - psi(iO1,jO1) )*( psi(iO1,jO1) - psi(iO1  ,jO1-1) ).gt.0.d0 ) check1=check1 + 1
  !      if ( check1.ne.0 ) Flag_OXO(1) = 0
  !   endif
  !   !   - Right - !
  !   if ( Flag_OXO(3).eq.1 ) then
  !      check2 = 0
  !      if ( ( psi(iO2+1,jO2  ) - psi(iO2,jO2) )*( psi(iO2,jO2) - psi(iO2-1,jO2  ) ).gt.0.d0 ) check2=check2 + 1
  !      if ( ( psi(iO2  ,jO2+1) - psi(iO2,jO2) )*( psi(iO2,jO2) - psi(iO2  ,jO2-1) ).gt.0.d0 ) check2=check2 + 1
  !      if ( check2.ne.0 ) Flag_OXO(3) = 0
  !   endif
  !   !   - Too close O-point => Single O-point - !
  !   if ( abs( iO2-iO1 ).le.5 ) then
  !      Flag_OXO(1) = 1
  !      Flag_OXO(3) = 1
  !      Flag_sOpt   = 1
  !   endif
  !   !  -- O-point position -- !
  !   iOXO(1) = iO1
  !   jOXO(1) = jO1
  !   iOXO(3) = iO2
  !   jOXO(3) = jO2

  !   ! --------------- !
  !   ! --- X-point --- !
  !   ! --------------- !
  !   if ( Flag_sOpt.eq.1 ) Flag__NoXpt = .true.
  !   if ( Flag_sOpt.eq.0 ) then
  !      !  -- search : Min( |Bp| ) -- !
  !      do j=2, LJp-1
  !         do i=iO1, iO2
  !            dpsidx1(i,j) = psi(i+1,j  ) - psi(i-1,j  )
  !            dpsidx2(i,j) = psi(i  ,j+1) - psi(i  ,j-1)
  !         enddo
  !      enddo
  !      absBMin   = 1.d8
  !      iX        = -1
  !      jX        = -1
  !      do j=2, LJp-1
  !         do i=iO1+1, iO2-1
  !            absBxy    = sqrt( dpsidx1(i,j)**2 + dpsidx2(i,j)**2 )
  !            if ( absBxy.lt.absBMin ) then
  !               absBMin   = absBxy
  !               iX        = i
  !               jX        = j
  !            endif
  !         enddo
  !      enddo
  !      if ( ( iX.eq.-1 ).or.( jX.eq.-1 )                             ) Flag__NoXpt = .true.
  !      if ( ( abs( iX-iOXO(1) ).lt.2 ).or.( abs( iX-iOXO(3) ).lt.2 ) ) Flag__NoXpt = .true.

  !   endif
  !   if ( Flag__NoXpt ) then
  !      iX          = nint( (iOXO(1)+iOXO(3))*0.5d0 )
  !      jX          = nint( (jOXO(1)+jOXO(3))*0.5d0 )
  !      Flag_OXO(2) = 0
  !   endif
  !   iOXO(2) = iX
  !   jOXO(2) = jX

  !   return
  ! end subroutine getXOpt

  
  ! integer                       :: iOXO(3), jOXO(3), Flag_OXO(3)
  ! logical                       :: InitOXO   = .false.
