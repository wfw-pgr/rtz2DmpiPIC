module diffDomMod
  use constants, only : LI, ns, PEtot
  implicit none
  integer              :: npLg(LI), npLk(LI,ns)
  integer, allocatable :: npLs(:,:)
  integer, parameter   :: neibPEtot_Max = 2
  integer, parameter   :: nCommLine_Max = 1
  integer, parameter   :: nBuffLine_Max = neibPEtot_Max * nCommLine_Max
  integer, parameter   :: ebCmp = 6
  integer, parameter   :: ptCmp = 8
  integer, parameter   :: rnk_=1, frm_=2, to_ =3, LIs_=4, LJs_=5
  integer, parameter   :: iFr_=6, iTo_=7, inn_=8, isl_=9, iel_=10
  integer              :: pBuffSize
  integer              :: nDst, nSrc, iFrontSpace, iRemnantFrm, niInc, niDec
  integer              :: export_Items(0:NeibPEtot_Max), export_Index(0:NeibPEtot_Max)
  integer              :: export_iFrom(  NeibPEtot_Max), export_iTo  (  NeibPEtot_Max)
  integer              :: export_Enemy(  NeibPEtot_Max)
  integer              :: import_Items(0:NeibPEtot_Max), import_Index(0:NeibPEtot_Max)
  integer              :: import_iFrom(  NeibPEtot_Max), import_iTo  (  NeibPEtot_Max)
  integer              :: import_Enemy(  NeibPEtot_Max), import_iSide(  NeibPEtot_Max)
  integer              :: import_gFrom(  NeibPEtot_Max)
  integer              :: nptPE(-1:PEtot)
  integer              :: export_ptItm(0:PEtot-1,ns), export_ptIdx(0:PEtot-1,ns)
  integer              :: import_ptItm(0:PEtot-1,ns), import_ptIdx(0:PEtot-1,ns)
  integer, allocatable :: request_Send(:), status_Send(:,:)
  integer, allocatable :: request_Recv(:), status_Recv(:,:)
  integer, allocatable :: ptrequest_Send(:), ptStatus_Send(:,:)
  integer, allocatable :: ptrequest_Recv(:), ptStatus_Recv(:,:)
  logical, parameter   :: Flag__Confirmation = .true.
  double precision, allocatable :: ptSendBuff(:,:,:), ptRecvBuff(:,:,:)
  double precision, parameter :: Crit = 0.1d0
  double precision     :: del1, del2
contains

  subroutine DiffusionDecomposition
    implicit none
    integer :: ierr
    call npcRedistribute
    ! call DefineExchangeTable
    ! call FieldExchange
    call ParticleExchange
    ! call RedefineFromTo

    call MPI_Finalize( ierr )
    stop
    return
  end subroutine DiffusionDecomposition


  subroutine npcRedistribute
    use constants , only : LI, LJ, LIs, ns, np, dzInv, OMPNumThreads, myRank, PEtot
    use constants , only : zp_, wp_
    use variables , only : pxv, FromTo
    use ptOrderMod, only : ptReOrdering
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: npLW(LIs+2,OMPNumThreads)
    integer             :: recvCnt(0:PEtot-1), displs(0:PEtot-1)
    integer             :: i, m, k, ip, ith, iPE, sendCnt, ierr
    integer             :: nComm, iDst, iSrc
    integer             :: isl, iel
    double precision    :: sum

    ! ----------------------------------------------- !
    ! --- [1] Sort Particle according to  i-Grid  --- !
    ! ----------------------------------------------- !
    call ptReOrdering( 'iGrid' )
    
    ! ----------------------------------------------- !
    ! --- [2] Count up #.of pt. on ith-Grid Line  --- !
    ! ----------------------------------------------- !
    !  -- [2-1] Count up k-th species :: npLs     --  !
    allocate( npLs(LIs,ns) )
    npLs = 0
    do k=1, ns
       if ( OMPNumThreads.eq.1 ) ith  = 1
       ! ith  = 1
       npLW = 0
       !$omp parallel default(none) &
       !$omp shared(k,np,pxv,dzInv,npLW) private(m,ip,ith)
       !$ ith = omp_get_thread_num() + 1
       !$omp do
       do m=1, np(k)
          ip           = ceiling( pxv(zp_,m,k)*dzInv ) + 1
          npLW(ip,ith) = npLW(   ip,ith) + 1
       enddo
       !$omp end do
       !$omp end parallel
       do ith=1, OMPNumThreads
          do i=1, LIs
             npLs(i,k) = npLs(i,k) + npLW(i,ith)
          enddo
       enddo
    enddo
    !  -- [2-3]   Communication  ( npL )          --  !
    isl     = FromTo(myRank,isl_)
    iel     = FromTo(myRank,iel_)
    sendCnt = FromTo(myRank,inn_)
    do iPE=0, PEtot-1
       recvCnt(iPE) = FromTo(iPE,inn_)
       displs (iPE) = FromTo(iPE,iFr_) - 1
    enddo
    npLk    = 0
    do k=1, ns
       call MPI_AllGatherV( npLs(isl:iel,k), sendCnt,         MPI_INTEGER, &
            &               npLk(      :,k), recvCnt, displs, MPI_INTEGER, &
            &               MPI_COMM_WORLD , ierr                          )
    enddo
    do i=1, LI
       npLg(i) = npLk(i,1) + npLk(i,2)
    enddo
    
    ! ----------------------------------------------- !
    ! --- [3] Judge whether exchange or NOT       --- !
    ! ----------------------------------------------- !
    !  -- [3-1] npt in all PE                     --  !
    do iPE=0, PEtot-1
       sum = 0
       do i = FromTo(iPE,iFr_), FromTo(iPE,iTo_)
          sum  = sum  + npLg(i)
       enddo
       nptPE(iPE) = sum
    enddo
    del1 = dble(         nptPE(myRank) - nptPE(myRank-1) ) &
         / dble( 0.5d0*( nptPE(myRank) + nptPE(myRank-1) ) )
    del2 = dble(         nptPE(myRank) - nptPE(myRank+1) ) &
         / dble( 0.5d0*( nptPE(myRank) + nptPE(myRank+1) ) )
    if ( myRank.eq.      0 ) del1 = 0.d0
    if ( myRank.eq.PEtot-1 ) del2 = 0.d0

    nDst         = 0
    nSrc         = 0
    iFrontSpace  = 0
    iRemnantFrm  = FromTo(myRank,isl_)
    export_Items = 0
    import_Items = 0
    if ( abs(del1).ge.Crit ) then
       if ( del1.gt.0.d0 ) then
          nDst               = nDst   + 1
          export_Enemy(nDst) = myRank - 1
          export_Items(nDst) = 1
          export_iFrom(nDst) = isl
          export_iTo  (nDst) = isl + export_Items(nDst) - 1
          iRemnantFrm        = export_iTo  (nDst)+1
       endif
       if ( del1.lt.0.d0 ) then
          nSrc               = nSrc   + 1
          import_Enemy(nSrc) = myRank - 1
          import_Items(nSrc) = 1
          import_iFrom(nSrc) = isl
          import_iTo  (nSrc) = isl + import_Items(nSrc) - 1
          import_gFrom(nSrc) = FromTo(myRank,iFr_) - 1 - import_Items(nSrc) + 1
          iFrontSpace        = import_Items(nSrc)
          import_iSide(nSrc) = 1
       endif
    endif
    if ( abs(del2).ge.Crit ) then
       if ( del2.gt.0.d0 ) then
          nDst               = nDst   + 1
          export_Enemy(nDst) = myRank + 1
          export_Items(nDst) = 1
          export_iFrom(nDst) = iel - export_Items(nDst) + 1
          export_iTo  (nDst) = iel
       endif
       if ( del2.lt.0.d0 ) then
          nSrc               = nSrc   + 1
          import_Enemy(nSrc) = myRank + 1
          import_Items(nSrc) = 1
          import_iFrom(nSrc) = iel - import_Items(nSrc) + 1
          import_gFrom(nSrc) = FromTo(myRank,iTo_) + 1
          import_iTo  (nSrc) = iel
          import_iSide(nSrc) = 2
       endif
    endif
    !  -- [3-2]  Export / Import --  !
    export_Index = 1
    do iDst=1, nDst
       export_Index(iDst) = export_Index(iDst-1) + export_Items(iDst-1)
    enddo
    import_Index = 1
    do iSrc=1, nSrc
       import_Index(iSrc) = import_Index(iSrc-1) + import_Items(iSrc-1)
    enddo
    nComm = nDst + nSrc
    niInc = 0
    niDec = 0
    do iDst=1, nDst
       niInc = niInc + import_Items(iDst)
    enddo
    do iSrc=1, nSrc
       niDec = niDec + export_Items(iSrc)
    enddo

    !  -- [2-3] Request / Status     --  !
    allocate( request_Send(nDst), status_Send(MPI_STATUS_SIZE,nDst) )
    allocate( request_Recv(nSrc), status_Recv(MPI_STATUS_SIZE,nSrc) )

    if ( myRank.eq.0 ) then
       open (40,file='send.dat',form='formatted',status='replace')
       write(40,'(6(a8,1x))') 'iPE', 'i', 'nDst', 'Items', 'Index', 'Enemy' 
       close(40)
    endif
    do iPE=0, PEtot-1
       if ( iPE.eq.myRank ) then
          open (40,file='send.dat',form='formatted',status='old',position='append')
          do i=1, nDst
             write(40,'(6(i8,1x))') iPE, i, nDst, export_Items(i), export_Index(i), export_Enemy(i)
          enddo
          close(40)
       endif
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
    enddo

    if ( myRank.eq.0 ) then
       open (40,file='recv.dat',form='formatted',status='replace')
       write(40,'(6(a8,1x))') 'iPE', 'i', 'nSrc', 'Items', 'Index', 'Enemy' 
       close(40)
    endif
    do iPE=0, PEtot-1
       if ( iPE.eq.myRank ) then
          open (40,file='recv.dat',form='formatted',status='old',position='append')
          do i=1, nSrc
             write(40,'(6(i8,1x))') iPE, i, nSrc, import_Items(i), import_Index(i), import_Enemy(i)
          enddo
          close(40)
       endif
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
    enddo

    
    return
  end subroutine npcRedistribute


  subroutine FieldExchange
    use constants, only : LIs, LJs, myRank, PEtot
    use variables, only : FromTo, EBf, EBr
    implicit none
    include 'mpif.h'
    integer :: ierr
    integer :: i, j, cmp, idx1, idx2, nRemnant
    integer :: iDst, iSrc, is, ir, len_s, len_r
    double precision, allocatable :: sendBuff(:,:,:), recvBuff(:,:,:)
    
    ! ------------------------- !
    ! --- [1]  Packing      --- !
    ! ------------------------- !
    !  -- [1-1] Allocation  --  !
    allocate( sendBuff(ebCmp,LJs,nBuffLine_Max) )
    allocate( recvBuff(ebCmp,LJs,nBuffLine_Max) )
    !  -- [1-2] sendBuff    --  !
    do iDst=1, nDst
       do i=1, export_Items(iDst)
          do j=1, LJs
             do cmp=1, ebCmp
                idx1                 = export_Index(iDst) + (i-1)
                idx2                 = export_iFrom(iDst) + (i-1)
                sendBuff(cmp,j,idx1) = EBf(cmp,idx2,j)
             enddo
          enddo
       enddo
    enddo
    
    ! ------------------------- !
    ! --- [2]  Exchange     --- !
    ! ------------------------- !
    !  -- [2-1] Send        --  !
    do iDst=1, nDst
       is    = export_Index(iDst)
       len_s = export_Items(iDst)*ebCmp*LJs
       call MPI_Isend( sendBuff(1,1,is)  , len_s, MPI_DOUBLE_PRECISION, &
            &          export_Enemy(iDst), 0    , MPI_COMM_WORLD      , &
            &          request_Send(iDst), ierr )
    enddo
    !  -- [2-2] Recv        --  !
    do iSrc=1, nSrc
       ir    = import_Index(iSrc)
       len_r = import_Items(iSrc)*ebCmp*LJs
       call MPI_Irecv( recvBuff(1,1,ir)  , len_r, MPI_DOUBLE_PRECISION, &
            &          import_Enemy(iSrc), 0    , MPI_COMM_WORLD      , &
            &          request_Recv(iSrc), ierr )
    enddo
    
    ! ------------------------- !
    ! --- [3]  Update       --- !
    ! ------------------------- !
    !  -- [3-1] WaitAll     --  !
    call MPI_WaitAll( nSrc, request_Recv, status_Recv, ierr )
    call MPI_WaitAll( nDst, request_Send, status_Send, ierr )
    !  -- [3-2] Re-Allocate --  !
    do j=1, LJs
       do i=1, LIs
          do cmp=1, ebCmp
             EBr(cmp,i,j) = EBf(cmp,i,j)
          enddo
       enddo
    enddo
    LIs = LIs + niInc - niDec
    deallocate( EBf )
    allocate( EBf(ebCmp,LIs,LJs) )
    do j=1, LJs
       do i=1, LIs
          do cmp=1, ebCmp
             EBf(cmp,i,j) = 0.d0
          enddo
       enddo
    enddo
    nRemnant = FromTo(myRank,inn_) - niDec
    do j=1, LJs
       do i=1, nRemnant
          do cmp=1, ebCmp
             idx1            = iFrontSpace + (i  )
             idx2            = iRemnantFrm + (i-1)
             EBf(cmp,idx1,j) = EBr(cmp,idx2,j)
          enddo
       enddo
    enddo
    !  -- [3-3] Update EBf  --  !
    do iSrc=1, nSrc
       do i=1, import_Items(iSrc)
          idx1            = import_iFrom(iSrc) + (i-1)
          idx2            = import_Index(iSrc) + (i-1)
          do j=1, LJs
             do cmp=1, ebCmp
                EBf(cmp,idx1,j) = recvBuff(cmp,j,idx2)
             enddo
          enddo
       enddo
    enddo
    return
  end subroutine FieldExchange


  subroutine ParticleExchange
    use constants, only : ns, np, myRank, PEtot
    use variables, only : pxv
    implicit none
    include 'mpif.h'
    integer            :: m, k, cmp, iFr, iTo, iDst, iSrc, idx1, idx2
    integer            :: is, ir, len_s, len_r, ierr, iPE
    integer            :: nptEx(ns), nptIm(ns)

    pBuffSize = nint( 1.0d0* max( np(1), np(2) ) )
    allocate( ptSendBuff(ptCmp,pBuffSize,ns) )
    allocate( ptRecvBuff(ptCmp,pBuffSize,ns) )
    !  -- [2-3] Request / Status     --  !
    allocate( ptRequest_Send(nDst), ptStatus_Send(MPI_STATUS_SIZE,nDst) )
    allocate( ptRequest_Recv(nSrc), ptStatus_Recv(MPI_STATUS_SIZE,nSrc) )
    ptSendBuff = 0
    ptRecvBuff = 0
    
    !  -- [2-1] Export Table         --  !
    export_ptItm = 0
    export_ptIdx = 1
    nptEx        = 0
    do k=1, ns
       do iDst=1, nDst
          iFr = export_iFrom(iDst)
          iTo = export_iFrom(iDst) + export_Items(iDst) - 1
          export_ptItm(iDst,k) = sum( npLs( iFr:iTo, k ) )
       enddo
    enddo
    do k=1, ns
       do iDst=1, nDst
          export_ptIdx(iDst,k) = export_ptIdx(iDst-1,k) + export_ptItm(iDst-1,k)
          nptEx(k)             = nptEx(k) + export_ptItm(iDst,k)
       enddo
    enddo
    !  -- [2-1] Import Table         --  !
    import_ptItm = 0
    import_ptIdx = 1
    nptIm        = 0
    do k=1, ns
       do iSrc=1, nSrc
          iFr = import_gFrom(iSrc)
          iTo = import_gFrom(iSrc) + import_Items(iSrc) - 1
          import_ptItm(iSrc,k) = sum( npLk( iFr:iTo, k ) )
       enddo
    enddo
    do k=1, ns
       do iSrc=1, nSrc
          import_ptIdx(iSrc,k) = import_ptIdx(iSrc-1,k) + import_ptItm(iSrc-1,k)
          nptIm(k)             = nptIm(k) + import_ptItm(iSrc,k)
       enddo
    enddo
    do k=1, ns
       ! --  Send  -- !
       do iDst=1, nDst
          is    = export_ptIdx(iDst,k)
          len_s = export_ptItm(iDst,k) * ptCmp
          ! is    = 1
          ! len_s = 1
          call MPI_iSend( ptSendBuff(1,is,k)  , len_s, MPI_DOUBLE_PRECISION, &
               &          export_enemy(iDst)  , 0    , MPI_COMM_WORLD      , &
               &          Request_Send(iDst), ierr )
       enddo
       ! --  Recv  -- !
       do iSrc=1, nSrc
          ir    = import_ptIdx(iSrc,k)
          len_r = import_ptItm(iSrc,k) * ptCmp
          ! ir    = 1
          ! len_r = 1
          call MPI_iRecv( ptRecvBuff(1,ir,k)  , len_r, MPI_DOUBLE_PRECISION, &
               &          import_enemy(iSrc)  , 0    , MPI_COMM_WORLD      , &
               &          Request_Recv(iSrc), ierr )
       enddo
       ! -- update -- !
       if ( nSrc.gt.0 ) call MPI_WaitAll( nSrc, Request_Recv, Status_Recv, ierr )
       if ( nDst.gt.0 ) call MPI_WaitAll( nDst, Request_Send, Status_Send, ierr )
       ! do iSrc=1, nSrc
       !    do m=1, import_ptItm(iSrc,k)
       !       idx1 = np(k)                + (m  )
       !       idx2 = import_ptIdx(iSrc,k) + (m-1)
       !       do cmp=1, ptCmp
       !          pxv(cmp,idx1,k) = ptRecvBuff(cmp,idx2,k)
       !       enddo
       !    enddo
       ! enddo
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
       ! -- update end -- !
    enddo
    
    ! if ( myRank.eq.0 ) then
    !    open (40,file='enemy.dat',form='formatted',status='replace')
    !    write(40,'(5(a8,1x))') 'iPE', 'iDst', 'nDst', 'Enemy' , 'inout'
    !    close(40)
    ! endif
    ! do iPE=0, PEtot-1
    !    if ( iPE.eq.myRank ) then
    !       open (40,file='enemy.dat',form='formatted',status='old',position='append')
    !       do iDst=1, nDst
    !          write(40,'(5(i8,1x))') iPE, iDst, nDst, export_Enemy(iDst), 1
    !       enddo
    !       do iSrc=1, nSrc
    !          write(40,'(5(i8,1x))') iPE, iSrc, nSrc, import_Enemy(iSrc), -1
    !       enddo
    !       close(40)
    !    endif
    !    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    ! enddo
    
    ! if ( myRank.eq.0 ) then
    !    open (40,file='send.dat',form='formatted',status='replace')
    !    write(40,'(7(a8,1x))') 'k', 'iPE', 'iDst', 'nDst', 'ptItm', 'ptIdx', 'Enemy' 
    !    close(40)
    ! endif
    ! do iPE=0, PEtot-1
    !    if ( iPE.eq.myRank ) then
    !       open (40,file='send.dat',form='formatted',status='old',position='append')
    !       do k=1, ns
    !          do iDst=1, nDst
    !             write(40,'(7(i8,1x))') k, iPE, iDst, nDst, export_ptItm(iDst,k), export_ptIdx(iDst,k), export_Enemy(iDst)
    !          enddo
    !       enddo
    !       close(40)
    !    endif
    !    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    ! enddo


    
    return
  end subroutine ParticleExchange

  

end module diffDomMod




    ! do k=1, 2
    !    nptEx(1) = 1
    !    nptEx(2) = 2
    !    nptIm    = 0
    !    len_s    = 2
    !    len_r    = 2
    !    do iDst=1, nDst
    !       call MPI_iSend( nptEx, len_s, MPI_INTEGER, &
    !            &          export_enemy(iDst), 0, MPI_COMM_WORLD, &
    !            &          Request_Send(iDst), ierr )
    !    enddo
    !    do iSrc=1, nSrc
    !       call MPI_iRecv( nptIm, len_r, MPI_INTEGER, &
    !            &          import_enemy(iSrc), 0, MPI_COMM_WORLD, &
    !            &          Request_Recv(iSrc), ierr )
    !    enddo
    !    if ( nSrc.gt.0 ) call MPI_WaitAll( nSrc, Request_Recv, Status_Recv, ierr )
    !    if ( nDst.gt.0 ) call MPI_WaitAll( nDst, Request_Send, Status_Send, ierr )
    !    ! call MPI_Barrier( MPI_COMM_WORLD, ierr )
    !    ! call MPI_Barrier( MPI_COMM_WORLD, ierr )
    ! enddo
    ! write(6,*) myRank, nptIm(1), nptIm(2)
    ! call MPI_Finalize( ierr )
    ! stop
