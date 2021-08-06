module ppcCtrlMod
contains

  subroutine ppcCtrl( action )
    use constants,          only : npt, LI, LJ, ns, np, jobdir, OMPNumThreads, myRank, PEtot
    use constants,          only : ppcMax, ppcMin, ppcMargin, drinv, dzinv
    use constants,          only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables,          only : pxv, goban, cellTop, perCell, goban, kstep, pCtrlDij
    use momentsMod,         only : nCountUp_SP
    use sortingMod,         only : sortptcl, ParallelSort
    use ptExchgMod,         only : prtclExchg
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    character(3),intent(in)     :: action ! 'APR' / 'COM'
    integer                     :: i, j, k, m, n2, ip, jp, iPE
    integer                     :: nwpt, count, ierr, nptLast, sumi
    integer                     :: ppcOpt
    integer                     :: nptAllMax, nptAllMin, nptAllAvg
    integer                     :: Npm(2,ns), perline(LI), NpmWhl(2,ns), npkWhl(ns)
    double precision            :: crit, nptRate(2)
    character(9),     parameter :: mode = 'ALLinCell'  ! ( 'RanSelect', 'ALLinCell', 'VDFSelect' )
    
    
    ! -- Diagnostic -- !
    ! if ( kstep.eq.1000 ) call VDFsample( 'out1.dat' )
    ! ---------------- !
    
    ! ===========     [1] Preparation    =========== !
    ppcOpt  = ppcMax  - ppcMargin
    if ( mode.eq.'VDFSelect' ) ppcOpt = 30
    ! --- (1-1) Categorize each particle into boxes :: GoBan --- !
    n2 = 0
    do k=1, ns
       !$omp parallel default(none) &
       !$omp shared(k,np,pxv,goban,dzinv,drinv) private(m,ip,jp)
       !$omp do
       do m=1, np(k)
          ip       = nint( pxv(zp_,m,k)*dzinv )
          jp       = nint( pxv(rp_,m,k)*drinv )
          goban(m) = LI*LJ*(k-1) + LI*jp + (ip-1) + 1
       enddo ! Simplicity -> jp=0 OK :: ip=0 xxx ( due to reflection )
       !$omp end do
       !$omp end parallel
    enddo
    ! --- (1-2) Load Balancing :: Determination of pCtrlDij :: --- !
    call nCountUp_SP( perCell )
    perline(:) = 0
    do i=1, LI
       sumi = 0
       !$omp parallel default(none) &
       !$omp shared(perCell,sumi,i) private(j,k)
       !$omp do reduction(+:sumi)
       do k=1, ns
          do j=1, LJ
             sumi = sumi + perCell(i,j,k)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       perline(i) = sumi
    enddo
    sumi  = sum( perline(:) )
    count = 0
    iPE   = 1
    pCtrlDij(1,1,1) = 1
    crit  = dble(sumi*iPE) / dble(PEtot)
    do i=1, LI
       count = count + perline(i)
       if ( count.ge.crit ) then
          pCtrlDij(2,1,iPE) = i
          iPE               = iPE + 1
          if ( iPE.eq.PEtot+1 ) exit
          pCtrlDij(1,1,iPE) = i+1
          crit              = dble(sumi*iPE) / dble(PEtot)
       endif
    enddo
    pCtrlDij(1,2,:) = 4
    pCtrlDij(2,2,:) = LJ-2
    
    ! ===========     [2] Particle Exchange    =========== !
    npm(:,:) = 0
    nptLast  = npt
    call prtclExchg(       1, np(1)      , nptLast, npm(:,1), 1 ) ! electron
    call prtclExchg( np(1)+1, np(1)+np(2), nptLast, npm(:,2), 2 ) ! ion
    np(1) = np(1) + ( npm(1,1) - npm(2,1) )
    np(2) = np(2) + ( npm(1,2) - npm(2,2) )
    npt   = np(1) + np(2)
    call MPI_REDUCE( npt, nptAllMin, 1, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_WORLD, ierr )
    call MPI_REDUCE( npt, nptAllMax, 1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierr )
    call MPI_REDUCE( npt, nptAllAvg, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    if ( myRank.eq.0 ) then
       nptAllAvg  = nptAllAvg / PEtot
       nptRate(1) = dble( nptAllAvg - nptAllMin ) / dble( nptAllAvg )
       nptRate(2) = dble( nptAllMax - nptAllAvg ) / dble( nptAllAvg )
       write(6,'(4x,a,i14,a,2(f6.3,a))') 'nptR:: Avg.', nptAllAvg, ' :: -', &
            & nptRate(1), '% -- +', nptRate(2), '%'
    endif
       

    ! ===========     [3] Sorting              =========== !
    ! ===== For MPI  ===== !
    if ( OMPNumThreads.gt.1 ) then
       call ParallelSort( nptLast, goban, 1, nptLast, pxv, OMPNumThreads )
    else
       call     sortptcl( nptLast, goban, 1, nptLast, pxv )
    endif
    if ( action.eq.'COM' ) return
    
    ! --- (1-1) Categorize each particle into boxes :: GoBan --- !
    n2 = 0
    do k=1, ns

       !$omp parallel default(none) &
       !$omp shared(k,np,pxv,goban,dzinv,drinv) private(m,ip,jp)
       !$omp do
       do m=1, np(k)
          ip       = nint( pxv(zp_,m,k)*dzinv )
          jp       = nint( pxv(rp_,m,k)*drinv )
          goban(m) = LI*LJ*(k-1) + LI*jp + (ip-1) + 1
       enddo ! Simplicity -> jp=0 OK :: ip=0 xxx ( due to reflection )
       !$omp end do
       !$omp end parallel
    enddo

    ! --- (1-2) Count up Number of SuperParticle in a cell. --- !
    call nCountUp_SP( perCell )
    n2 = 1
    do k=1, ns
       count = n2
       do j=1, LJ
          do i=1, LI
             cellTop(i,j,k) = count
             count          = count + perCell( i,j,k )
          enddo
       enddo
       n2 = n2 + np(k)
    enddo

    ! ======     [2] Particle Merge / Split     ====== !
    Npm(:,:) = 0
    call PrtclMerge( Npm(2,:), mode )
    call PrtclSplit( nwpt, Npm(1,:) )
    
    ! ======     [3] Post Process               ====== !
    ! --- (3-1) Update   npt / np   --- !
    npt = 0
    do k=1, ns
       np(k) = np(k) + Npm(1,k) - Npm(2,k)
       npt   = npt   + np(k)
    enddo
    call MPI_REDUCE( Npm  , NpmWhl, 4, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    call MPI_REDUCE( np(1), npkWhl, 2, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    if ( myRank.eq.0 ) then
       write(6,'(    12x,3(a14,4x))')           '  ne   '   , '  ni  '   , ' total '
       write(6,'(2x,a,2x,3(i14,4x))') ' #ptcl ', npkWhl(  1), npkWhl(  2), npkWhl(  1)+npkWhl(  2)
       write(6,'(2x,a,2x,3(i14,4x))') ' split ', NpmWhl(1,1), NpmWhl(1,2), NpmWhl(1,1)+NpmWhl(1,2)
       write(6,'(2x,a,2x,3(i14,4x))') ' merge ', NpmWhl(2,1), NpmWhl(2,2), NpmWhl(2,1)+NpmWhl(2,2)
       write(6,*)
       open( 50,file=trim(jobdir)//'dat/'//'Nparticles.dat',form='formatted',access='append')
       write(50,'(4(i14,1x),4(i8,1x))') kstep, npkWhl(  1) +npkWhl(  2), npkWhl(  1), npkWhl(  2) &
            &                                , NpmWhl(1,1), NpmWhl(1,2), NpmWhl(2,1), NpmWhl(2,2)
       close(50)
    endif
    
    ! --- (3-2)      Sorting        --- !
    if ( OMPNumThreads.gt.1 ) then
       call ParallelSort( nwpt-1, goban, 1, nwpt-1, pxv, OMPNumThreads )
    else
       call     sortptcl( nwpt-1, goban, 1, nwpt-1, pxv )
    endif
    
    ! ! ------- Diagnostic -------- !
    ! if ( kstep.eq.1000 ) call VDFsample( 'out2.dat' )

    return
  end subroutine ppcCtrl


  
  subroutine PrtclMerge( Npm, mode )
    
    use constants,  only          : LI, LJ, ns, ppcMax, ppcMargin, dr, dz, OMPNumThreads, myRank
    use constants,  only          : APRperCell, APRvdfResl
    use constants,  only          : drinv, dzinv
    use constants,  only          : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables,  only          : pxv, goban, cellTop, perCell, pCtrlDij
    use ptCoalsMod,  only          : VDFBinning, coalesptcl
    use randGenMod, only          : Knuth_Shuffle_OMP
    !$ use omp_lib
    implicit none
    integer,          parameter :: OutIdx       = LI*LJ*ns + 1
    integer,          parameter :: pBufSize     = 10*ppcMax     ! --     VDF Select -- !
    integer,          parameter :: InpThreshold = 40
    integer,          parameter :: ppcOptVDF    = 30
    integer,          parameter :: coalesBase   = 20           ! --  Random Select -- !
    double precision, parameter :: coalesBeta   = 0.6d0
    integer                     :: i, j, k, m, ivdf, ptidx, mythread, ierr
    integer                     :: pfirst, plast, coalesInp, coalesOut, ppcOpt
    integer                     :: ileft, irght, jbot, jtop
    integer                     :: Npmp(ns,OMPNumThreads), coalesInpList(APRperCell), idx(pBufSize,APRperCell)
    double precision            :: pxvr(8,pBufSize,APRperCell)
    integer ,     intent(inout) :: Npm(ns)
    character(9), intent(in)    :: mode
    integer :: i0, ip, j0, jp
    
    ! --- [1] Initialization --- !
    if ( mode.eq.'VDFSelect' ) ppcOpt = ppcOptVDF
    if ( mode.eq.'ALLinCell' ) ppcOpt = ppcMax  - ppcMargin
    if ( mode.eq.'RanSelect' ) ppcOpt = ppcMax  - ppcMargin
    Npmp(:,:) = 0
    ileft     = pCtrlDij( 1,1,myRank+1 )
    irght     = pCtrlDij( 2,1,myRank+1 )
    jbot      = pCtrlDij( 1,2,myRank+1 )
    jtop      = pCtrlDij( 2,2,myRank+1 )
    
    ! --- [2] Main Loop  --- !
    !!$omp parallel default(none) &
    !!$omp shared(pxv,goban,cellTop,perCell,dr,dz,Npmp,mode,ppcOpt,jbot,jtop,ileft,irght) &
    !!$omp private( mythread,ierr,pfirst,plast,coalesInp,coalesOut,coalesInpList,pxvr,idx,ptidx)
    !!$ mythread = omp_get_thread_num() + 1
    !!$omp do private(i,j,k,m,ivdf)
    mythread = 1
    do k=1, ns
       do j=jbot, jtop
          do i=ileft, irght
             ! -- Loop Begin -- !
             ierr   = 0
             if ( perCell(i,j,k).gt.ppcMax ) then
                ! ---   (1) Range for this cell   --- !
                pfirst = cellTop(i,j,k)
                plast  = cellTop(i,j,k) + perCell(i,j,k) - 1
                ! ---   (2) Mode Designation      --- !
                if      ( mode.eq.'VDFSelect' ) then  !   -- (1) VDF Select    -- !
                   call VDFBinning( pxv(1:8,pfirst:plast,k), perCell(i,j,k), APRperCell, APRvdfResl, &
                        &           coalesInpList(1:APRperCell), pxvr(1:8,1:perCell(i,j,k),1:APRperCell), &
                        &           idx(1:perCell(i,j,k),1:APRperCell) )
                   
                else if ( mode.eq.'RanSelect' ) then  !   -- (2) Random Select -- !
                   coalesInpList(1) = coalesBase + ( perCell(i,j,k) - ppcMax )
                   coalesOut        = ceiling( coalesBeta * coalesInpList(1) )
                   call Knuth_Shuffle_OMP( pfirst, plast )
                   plast            = cellTop(i,j,k) + coalesInpList(ivdf) - 1
                   do m=1, coalesInpList(1)
                      pxvr(1:8,m,1) = pxv(1:8,pfirst+m-1,k)
                      idx(m,1)      = m
                   enddo
                   
                else if ( mode.eq.'ALLinCell' ) then  !   -- (3) ALL Particle  -- !
                   coalesInpList(1) = perCell(i,j,k)
                   coalesOut        = ppcOpt
                   do m=1, perCell(i,j,k)
                      pxvr(1:8,m,1) = pxv(1:8,pfirst+m-1,k)
                      idx(m,1)      = m
                   enddo
                endif
                
                ! ---   (3) Coalessing            --- !
                do ivdf=1, APRperCell
                   if ( coalesInpList(ivdf).gt.InpThreshold ) then
                      coalesInp = coalesInpList(ivdf)
                      coalesOut = ppcOpt
                      call coalesptcl( pxvr(1:8,1:coalesInp,ivdf), dr       , dz  , &
                           &           coalesOut                 , coalesInp, ierr  )
                      if ( ierr.eq.0 ) then
                         Npmp(k,mythread) = Npmp(k,mythread) + ( coalesInp - coalesOut )
                         do m=1, coalesOut
                            ptidx          = pfirst + idx(m,ivdf) - 1
                            pxv(1:8,ptidx,k) = pxvr(1:8,m,ivdf)
                            goban(ptidx)   = LI*LJ*(k-1) + LI*j + (i-1) + 1
                         enddo
                         do m=coalesOut+1, coalesInp
                            ptidx          = pfirst + idx(m,ivdf) - 1
                            pxv(wp_,ptidx,k)   = 0.d0
                            goban(ptidx)   = OutIdx
                         enddo
                         i0       = nint( pxv(zp_,pfirst,k)*dzinv )
                         j0       = nint( pxv(rp_,pfirst,k)*drinv )
                         do m=1, coalesInp
                            ptidx          = pfirst + idx(m,ivdf) - 1
                            ip       = nint( pxv(zp_,ptidx,k)*dzinv )
                            jp       = nint( pxv(rp_,ptidx,k)*drinv )
                            if ( ( ip.ne.i0 ).or.( jp.ne.j0 ) ) then
                               write(6,'(a,4(i8,1x))') '[ERROR] ip, i0, jp, j0 ', ip, i0, jp, j0
                            endif
                         enddo
                      endif
                   endif
                enddo
                
             endif             
             ! -- Loop END  -- !
          enddo
       enddo
    enddo
    !!$omp end do
    !!$omp end parallel
    
    ! --- [3] Particle Income (+/-)  --- !
    Npm(1:ns)    = 0
    do m=1, OMPNumThreads
       do k=1, ns
          Npm(k) = Npm(k) + Npmp(k,m)
       enddo
    enddo

    return
  end subroutine PrtclMerge

  
  
  subroutine PrtclSplit( nwpt, Npm )

    use constants, only     : LI, LJ, ns, npt, nptMax, dr, dz, ppcMinThr, ppcMargin, myRank
    use constants, only     : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables, only     : pxv, goban, cellTop, perCell, pCtrlDij
    use ptSplitMod, only    : AssousSplit
    implicit none
    integer                :: i, j, k, pfirst, plast, ppcNewMax, ierr
    integer                :: ileft, irght, jbot, jtop
    integer, parameter     :: ppcSplMin = 30
    integer, intent(inout) :: nwpt, Npm(ns)

    Npm(:) = 0
    nwpt   = npt + 1
    ileft     = pCtrlDij( 1,1,myRank+1 )
    irght     = pCtrlDij( 2,1,myRank+1 )
    jbot      = pCtrlDij( 1,2,myRank+1 )
    jtop      = pCtrlDij( 2,2,myRank+1 )

    do k=1, ns
       do j=jbot, jtop
          do i=ileft, irght
             
             ! --- Main Loop Begin --- !
             ierr   = 0
             if ( ( perCell(i,j,k).lt.ppcMinThr(j) ).and.( perCell(i,j,k).gt.ppcSplMin ) ) then
                !  --   (1) Range for this cell   -- !
                pfirst    = cellTop(i,j,k)
                plast     = cellTop(i,j,k) + perCell(i,j,k) - 1
                ppcNewMax = ppcMinThr(j)   + ppcMargin
                !  --   (2) OverFlow Check        -- !
                if ( nwpt+(ppcNewMax-perCell(i,j,k)).ge.nptMax ) then
                   write(6,*) ' [ERROR] Split Memory is not sufficient nwpt [ERROR] '
                   write(6,*) '    ----- nwpt   == ', nwpt
                   write(6,*) '    ----- Nnew   == ', ppcNewMax
                   write(6,*) '    ----- nptMax == ', nptMax
                   exit
                endif
                !  --   (3) Assous's Splitting    -- !
                call AssousSplit( pxv(1:8,pfirst:plast,k), goban(nwpt:nptMax)  , dr  , dz    , &
                     &            ppcNewMax            , perCell(i,j,k)      , nwpt, nptMax, &
                     &            goban(pfirst)        , pxv(1:8,nwpt:nptMax,k), ierr          )
                if ( ierr.eq.0 ) then
                   Npm(k) = Npm(k) + ( ppcNewMax - perCell(i,j,k) )
                   ! else
                   ! write(6,*) ' [WARNING] No Split Assous in coalesptcl in (', i, ',', j, ') [WARNING] '
                   ! write(6,*) pfirst, plast, Npm(1,k)
                endif
             endif
             ! ---  Main Loop END  --- !
             
          enddo
       enddo
    enddo
    
    return
  end subroutine PrtclSplit
  
  
end module ppcCtrlMod








    !  -- [DEBUG] Leak Checker [DEBUG] -- !
    ! pLeak = 0
    ! do m=1, npt
    !    ip  = max( ceiling( pxv(zp_,m)*dzinv ), 1 )
    !    iPE = max( ceiling( ip * dble(PEtot) / dble(LI-1) ), 1 )
    !    if ( iPE.ne.myRank+1 ) then
    !       pLeak = pLeak + 1
    !    end if
    ! end do
    ! write(6,'(5(i8,1x))') myRank+1, pLeak, np(1), np(2), npt
    ! return
    !  -- [DEBUG] Leak Checker [DEBUG] -- !
    ! !  -- [DEBUG] npt  Checker [DEBUG] -- !
    ! count = 0
    ! do m=1, nptMax
    !    if ( pxv(wp_,m).gt.0.d0 ) then
    !       count = count + 1
    !    endif
    ! enddo
    ! write(6,*) npt, count
    ! write(6,*) 'between e-i'
    ! write(6,*) pxv(rp_:2,np(1)  ), goban(np(1)  )
    ! write(6,*) pxv(rp_:2,np(1)+1), goban(np(1)+1)
    ! write(6,*) 'between i-v'
    ! write(6,*) pxv(rp_:2,npt  ), goban(npt  ), pxv(wp_,npt  )
    ! write(6,*) pxv(rp_:2,npt+1), goban(npt+1), pxv(wp_,npt+1)
    ! write(6,*) 'outidx', LI*LJ*ns+1
    ! !  -- [DEBUG] npt  Checker [DEBUG] -- !
    

    ! if ( myRank.eq.0 ) then
    !    open(60,file='out.dat',form='formatted')
    !    do j=1, LJ
    !       write(60,*) ( perCell(i,j,1), i=1, LI )
    !    enddo
    !    close(60)
    !    stop
    ! endif

    ! if ( myRank.eq.1 ) then
    !    do iPE=1, PEtot
    !       write(6,*) iPE, pCtrlDij(1,1,iPE), pCtrlDij(2,1,iPE)
    !    enddo
    ! endif
    ! do iPE=1, PEtot
    !    pCtrlDij( 1,1,iPE ) = ceiling( dble(LI-1) / dble(PEtot) * (iPE-1) ) + 1
    !    pCtrlDij( 2,1,iPE ) = ceiling( dble(LI-1) / dble(PEtot) * (iPE  ) )
    ! enddo
    ! pCtrlDij( 1,1,1 ) = 1
    ! if ( myRank.eq.0 ) then
    !    do iPE=1, PEtot
    !       write(6,*) iPE, pCtrlDij(1,1,iPE), pCtrlDij(2,1,iPE)
    !    enddo
    ! endif
