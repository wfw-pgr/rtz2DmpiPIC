module diffDomMod
  use constants, only             : LI, ns, PEtot
  implicit none
  integer                        :: npLg(LI), npLk(LI,ns)
  integer, allocatable           :: npLs(:,:)
  integer, parameter             :: neibPEtot_Max = 2
  integer, parameter             :: nCommLine_Max = 1
  integer, parameter             :: nBuffLine_Max = neibPEtot_Max * nCommLine_Max
  integer, parameter             :: ebCmp         = 6
  integer, parameter             :: ptCmp         = 8
  integer, parameter             :: rnk_=1, frm_=2, to_ =3, LIs_=4, LJs_=5
  integer, parameter             :: iFr_=6, iTo_=7, inn_=8, isl_=9, iel_=10
  integer                        :: pBuffSize
  integer                        :: nDst, nSrc, LIRemnant, NewLIs
  integer                        :: export_Items(0:NeibPEtot_Max), export_Index(0:NeibPEtot_Max)
  integer                        :: export_iFrom(  NeibPEtot_Max)
  integer                        :: export_Enemy(  NeibPEtot_Max)
  integer                        :: import_Items(0:NeibPEtot_Max), import_Index(0:NeibPEtot_Max)
  integer                        :: import_Enemy(  NeibPEtot_Max)
  integer                        :: import_gFrom(  NeibPEtot_Max)
  integer                        :: export_LItem, export_RItem, import_LItem, import_RItem
  integer                        :: nptPE(-1:PEtot)
  integer                        :: export_ptItm(0:PEtot-1,ns), export_ptIdx(0:PEtot-1,ns)
  integer                        :: import_ptItm(0:PEtot-1,ns), import_ptIdx(0:PEtot-1,ns)
  integer                        :: export_ptFrm(0:PEtot-1,ns)
  integer, allocatable           :: request_Send(:), status_Send(:,:)
  integer, allocatable           :: request_Recv(:), status_Recv(:,:)
  logical, parameter             :: Flag__Confirmation = .true.
  integer                        :: FromTo_(0:PEtot-1,10)
  double precision, allocatable  ::  ptSendBuff(:,:,:),  ptRecvBuff(:,:,:)
  double precision, allocatable  :: EBfSendBuff(:,:,:), EBfRecvBuff(:,:,:)
  double precision, allocatable  :: EBoSendBuff(:,:,:), EBoRecvBuff(:,:,:)
  double precision, allocatable  :: FieldStack(:,:,:,:)
  double precision, parameter    :: Crit = 0.1d0
  double precision               :: del1, del2
  double precision               :: zoffset(2)
contains

  
  subroutine DiffusionDecomposition
    ! use constants, only : LIs, LJs, cRank, myRank
    ! use variables, only : EBf, FromTo
    implicit none
    integer            :: ierr
    ! integer            :: i, j 
    ! character(100)     :: FileName
    
    call npcRedistribute
    call NewLIsSystemMaking
    ! call DefineExchangeTable

    ! do j=1, LJs
    !    do i=1, LIs
    !       EBf(5,i,j) = dble(i) + FromTo(myRank,frm_) - 1
    !    enddo
    ! enddo
    ! FileName = 'bin/out1_' // cRank // '.bin'
    ! open( 50, file=trim(FileName), form='unformatted' )
    ! write(50) EBf(5,:,:)
    ! close(50)
    
    call FieldExchange

    ! FileName = 'bin/out2_' // cRank // '.bin'
    ! open( 50, file=trim(FileName), form='unformatted' )
    ! write(50) EBf(5,:,:)
    ! close(50)

    call ParticleExchange
    call UpdateFromTo

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
    integer             :: npLW(LIs+2,OMPNumThreads), npLD(LIs,ns)
    integer             :: recvCnt(0:PEtot-1), displs(0:PEtot-1)
    integer             :: i, m, k, ip, ith, iPE, sendCnt, ierr
    integer             :: iDst, iSrc
    integer             :: isl, iel, iFr, iTo

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
       npLW = 0
       !$omp parallel default(none) &
       !$omp shared(k,np,pxv,dzInv,npLW) private(m,ip,ith)
       !$ ith = omp_get_thread_num() + 1
       !$omp do
       do m=1, np(k)
          ip           = ceiling( pxv(zp_,m,k)*dzInv ) + 1
          npLW(ip,ith) = npLW(ip,ith) + 1
       enddo
       !$omp end do
       !$omp end parallel
       do ith=1, OMPNumThreads
          do i=1, LIs
             npLs(i,k) = npLs(i,k) + npLW(i,ith)
          enddo
       enddo
    enddo
    !  -- [2-2]   Communication  ( npL )          --  !
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
    ! --- [3] Judge Exchange ( Import / Export )  --- !
    ! ----------------------------------------------- !
    !  -- [3-1] npt in all PE                     --  !
    do iPE=0, PEtot-1
       iFr        = FromTo(iPE,iFr_)
       iTo        = FromTo(iPE,iTo_)
       nptPE(iPE) = sum( npLg( iFr:iTo ) )
    enddo
    !  -- [3-2] del1 & del2 ( indicator of Exch.) --  !
    del1 = dble(         nptPE(myRank) - nptPE(myRank-1) ) &
         / dble( 0.5d0*( nptPE(myRank) + nptPE(myRank-1) ) )
    del2 = dble(         nptPE(myRank) - nptPE(myRank+1) ) &
         / dble( 0.5d0*( nptPE(myRank) + nptPE(myRank+1) ) )
    if ( myRank.eq.      0 ) del1 = 0.d0
    if ( myRank.eq.PEtot-1 ) del2 = 0.d0
    !  -- [3-3] Generalized Comm. Table           --  !
    nDst          = 0
    nSrc          = 0
    export_Items  = 0
    import_Items  = 0
    export_LItem  = 0
    export_RItem  = 0
    import_LItem  = 0
    import_RItem  = 0
    if ( abs(del1).ge.Crit ) then
       if ( del1.gt.0.d0 ) then
          nDst               = nDst   + 1
          export_Enemy(nDst) = myRank - 1
          export_LItem       = 1
          export_Items(nDst) = 1
          export_iFrom(nDst) = isl
       endif
       if ( del1.lt.0.d0 ) then
          nSrc               = nSrc   + 1
          import_Enemy(nSrc) = myRank - 1
          import_LItem       = 1
          import_Items(nSrc) = 1
          import_gFrom(nSrc) = FromTo(myRank,iFr_) - 1 - import_Items(nSrc) + 1
       endif
    endif
    if ( abs(del2).ge.Crit ) then
       if ( del2.gt.0.d0 ) then
          nDst               = nDst   + 1
          export_Enemy(nDst) = myRank + 1
          export_RItem       = 1
          export_Items(nDst) = 1
          export_iFrom(nDst) = iel - export_Items(nDst) + 1
       endif
       if ( del2.lt.0.d0 ) then
          nSrc               = nSrc   + 1
          import_Enemy(nSrc) = myRank + 1
          import_RItem       = 1
          import_Items(nSrc) = 1
          import_gFrom(nSrc) = FromTo(myRank,iTo_) + 1
       endif
    endif
    !  -- [3-4] Index & New LIs-Grid              --  !
    export_Index = 1
    do iDst=1, nDst
       export_Index(iDst) = export_Index(iDst-1) + export_Items(iDst-1)
    enddo
    import_Index = 1
    do iSrc=1, nSrc
       import_Index(iSrc) = import_Index(iSrc-1) + import_Items(iSrc-1)
    enddo
    LIRemnant = FromTo(myRank,inn_) - ( export_LItem + export_RItem )

    ! ----------------------------------------------- !
    ! --- [4] Generalized Comm.Table ( particle ) --- !
    ! ----------------------------------------------- !
    !  -- [4-1] Index for Particle  ( export )    --  !
    export_ptItm = 0
    export_ptIdx = 1
    do k=1, ns
       do iDst=1, nDst
          iFr = export_iFrom(iDst)
          iTo = export_iFrom(iDst) + export_Items(iDst) - 1
          export_ptItm(iDst,k) = sum( npLs( iFr:iTo, k ) )
       enddo
       do iDst=1, nDst
          export_ptIdx(iDst,k) = export_ptIdx(iDst-1,k) + export_ptItm(iDst-1,k)
       enddo
    enddo
    !  -- [4-2] Index for Particle  ( import )    --  !
    import_ptItm = 0
    import_ptIdx = 1
    do k=1, ns
       do iSrc=1, nSrc
          iFr = import_gFrom(iSrc)
          iTo = import_gFrom(iSrc) + import_Items(iSrc) - 1
          import_ptItm(iSrc,k) = sum( npLk( iFr:iTo, k ) )
       enddo
       do iSrc=1, nSrc
          import_ptIdx(iSrc,k) = import_ptIdx(iSrc-1,k) + import_ptItm(iSrc-1,k)
       enddo
    enddo
    !  -- [4-3] Export ptcl. Load From this point --  !
    npLD = 1
    do k=1, ns
       do i=2, LIs
          npLD(i,k) = npLD(i-1,k) + npLs(i-1,k)
       enddo
    enddo
    do k=1, ns
       do iDst=1, nDst
          export_ptFrm(iDst,k) = npLD( export_iFrom(iDst), k )
       enddo
    enddo

    ! ----------------------------------------------- !
    ! --- [5]    Request-Status Allocation        --- !
    ! ----------------------------------------------- !
    allocate( request_Send(nDst), status_Send(MPI_STATUS_SIZE,nDst) )
    allocate( request_Recv(nSrc), status_Recv(MPI_STATUS_SIZE,nSrc) )
    
    return
  end subroutine npcRedistribute


  subroutine FieldExchange
    use constants, only : LIs, LJs, myRank, PEtot
    use variables, only : FromTo, EBf, EBo
    implicit none
    include 'mpif.h'
    integer            :: i, j, cmp, idx1, idx2, isl
    integer            :: iDst, iSrc, is, ir, len_s, len_r, ierr
    
    ! ------------------------- !
    ! --- [1]  Packing      --- !
    ! ------------------------- !
    !  -- [1-1] Allocation  --  !
    allocate( EBfSendBuff(ebCmp,LJs,nBuffLine_Max) )
    allocate( EBfRecvBuff(ebCmp,LJs,nBuffLine_Max) )
    allocate( EBoSendBuff(ebCmp,LJs,nBuffLine_Max) )
    allocate( EBoRecvBuff(ebCmp,LJs,nBuffLine_Max) )
    allocate( FieldStack (ebCmp,LIs,LJs,2) )
    !  -- [1-2] ebSendBuff    --  !
    do iDst=1, nDst
       do i=1, export_Items(iDst)
          do j=1, LJs
             do cmp=1, ebCmp
                idx1                    = export_Index(iDst) + (i-1)
                idx2                    = export_iFrom(iDst) + (i-1)
                EBfSendBuff(cmp,j,idx1) = EBf(cmp,idx2,j)
                EBoSendBuff(cmp,j,idx1) = EBo(cmp,idx2,j)
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
       call MPI_Isend( EBfSendBuff(1,1,is), len_s, MPI_DOUBLE_PRECISION, &
            &          export_Enemy(iDst) , 0    , MPI_COMM_WORLD      , &
            &          request_Send(iDst) , ierr )
    enddo
    !  -- [2-2] Recv        --  !
    do iSrc=1, nSrc
       ir    = import_Index(iSrc)
       len_r = import_Items(iSrc)*ebCmp*LJs
       call MPI_Irecv( EBfRecvBuff(1,1,ir), len_r, MPI_DOUBLE_PRECISION, &
            &          import_Enemy(iSrc) , 0    , MPI_COMM_WORLD      , &
            &          request_Recv(iSrc) , ierr )
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
             FieldStack(cmp,i,j,1) = EBf(cmp,i,j)
          enddo
       enddo
    enddo
    LIs = NewLIs
    deallocate( EBf )
    allocate( EBf(ebCmp,LIs,LJs) )
    EBf = 0.d0
    !  -- [3-3] Restore EBf --  !
    isl = FromTo(myRank,isl_)
    do j=1, LJs
       do i=1, LIRemnant
          do cmp=1, ebCmp
             idx1            = isl + import_LItem + (i-1)
             idx2            = isl + export_LItem + (i-1)
             EBf(cmp,idx1,j) = FieldStack(cmp,idx2,j,1)
          enddo
       enddo
    enddo
    !  -- [3-3] Update EBf  --  !
    do iSrc=1, nSrc
       do i=1, import_Items(iSrc)
          idx1            =                isl + (i-1)
          idx2            = import_Index(iSrc) + (i-1)
          do j=1, LJs
             do cmp=1, ebCmp
                EBf(cmp,idx1,j) = EBfRecvBuff(cmp,j,idx2)
             enddo
          enddo
       enddo
    enddo
    
    return
  end subroutine FieldExchange


  subroutine ParticleExchange
    use constants , only : myRank, npc, PEtot
    use constants , only : ns, np, npt, wp_, zp_
    use variables , only : pxv
    use ptOrderMod, only : ptReOrdering
    implicit none
    include 'mpif.h'
    integer            :: m, k, cmp, iDst, iSrc, idx1, idx2
    integer            :: is, ir, len_s, len_r, ierr
    integer            :: npV(12), iPE

    ! ------------------------- !
    ! --- [1]  Packing      --- !
    ! ------------------------- !
    !  -- [1-1] Allocation  --  !
    pBuffSize = nint( 1.0d0* max( np(1), np(2) ) )
    allocate( ptSendBuff(ptCmp,pBuffSize,ns) )
    allocate( ptRecvBuff(ptCmp,pBuffSize,ns) )
    !  -- [1-2] Global zp   --  !
    do k=1, ns
       do m=1, np(k)
          pxv(zp_,m,k) = pxv(zp_,m,k) + zoffset(1)
       enddo
    enddo
       
    !  -- [1-3] ptSendBuff  --  !
    do k=1, ns
       do iDst=1, nDst
          do m=1, export_ptItm(iDst,k)
             idx1 = export_ptIdx(iDst,k) + (m-1)
             idx2 = export_ptFrm(iDst,k) + (m-1)
             do cmp=1, ptCmp
                ptSendBuff(cmp,idx1,k) = pxv(cmp,idx2,k)
                pxv(wp_,idx2,k)        = 0.d0
             enddo
          enddo
       enddo
    enddo
    
    ! ------------------------- !
    ! --- [2]  Exchange     --- !
    ! ------------------------- !
    do k=1, ns
       !  -- [2-1] Send     --  !
       do iDst=1, nDst
          is    = export_ptIdx(iDst,k)
          len_s = export_ptItm(iDst,k) * ptCmp
          call MPI_iSend( ptSendBuff(1,is,k)  , len_s, MPI_DOUBLE_PRECISION, &
               &          export_enemy(iDst)  , 0    , MPI_COMM_WORLD      , &
               &          Request_Send(iDst), ierr )
       enddo
       !  -- [2-2] Recv     --  !
       do iSrc=1, nSrc
          ir    = import_ptIdx(iSrc,k)
          len_r = import_ptItm(iSrc,k) * ptCmp
          call MPI_iRecv( ptRecvBuff(1,ir,k)  , len_r, MPI_DOUBLE_PRECISION, &
               &          import_enemy(iSrc)  , 0    , MPI_COMM_WORLD      , &
               &          Request_Recv(iSrc), ierr )
       enddo
       !  -- [2-3] Update   --  !
       call MPI_WaitAll( nSrc, Request_Recv, Status_Recv, ierr )
       call MPI_WaitAll( nDst, Request_Send, Status_Send, ierr )
       do iSrc=1, nSrc
          do m=1, import_ptItm(iSrc,k)
             idx1 = np(k)                + (m  )
             idx2 = import_ptIdx(iSrc,k) + (m-1)
             do cmp=1, ptCmp
                pxv(cmp,idx1,k) = ptRecvBuff(cmp,idx2,k)
             enddo
          enddo
       enddo
    enddo
    !  -- [1-2]  Local zp   --  !
    do k=1, ns
       do m=1, np(k)
          pxv(zp_,m,k) = pxv(zp_,m,k) - zoffset(2)
       enddo
    enddo

    !  -- [2-4] update np -- !
    npV(1)  = np (1)
    npV(4)  = np (2)
    npV(7)  = npc(1)
    npV(10) = npc(2)
    do k=1, ns
       do iSrc=1, nSrc
          np(k) = np(k) + import_ptItm(iSrc,k)
       enddo
    enddo
    npt = np(1) + np(2)
    npV(2)  = np (1)
    npV(5)  = np (2)
    npV(8)  = npc(1)
    npV(11) = npc(2)
    call ptReOrdering( 'iGrid' )
    npV(3)  = np (1)
    npV(6)  = np (2)
    npV(9)  = npc(1)
    npV(12) = npc(2)
    
    if ( myRank.eq.0 ) then
       open (40,file='npt.dat',form='formatted',status='replace')
       write(40,'(12(a8,1x))') 'np1(1)', 'np1(2)', 'npc1(1)', 'npc1(2)', &
            &                  'np2(1)', 'np2(2)', 'npc2(1)', 'npc2(2)', &
            &                  'np3(1)', 'np3(2)', 'npc3(1)', 'npc3(2)'
       close(40)
    endif
    do iPE=0, PEtot-1
       if ( iPE.eq.myRank ) then
          open (40,file='npt.dat',form='formatted',status='old',position='append')
          write(40,'(12(i8,1x))') npV(1), npV(2), npV(3), npV(4), &
               &                  npV(5), npV(6), npV(7), npV(8), &
               &                  npV(9), npV(10), npV(11), npV(12)
          close(40)
       endif
       call MPI_Barrier( MPI_COMM_WORLD, ierr )
    enddo
    
    return
  end subroutine ParticleExchange


  subroutine NewLIsSystemMaking
    use constants, only : LIs, LJs, dz, myRank, PEtot
    use variables, only : FromTo
    implicit none
    include 'mpif.h'
    integer            :: iPE, ierr
    integer            :: LIg(0:PEtot-1), LIi(0:PEtot-1)
    integer            :: iFrom(0:PEtot), iTo(0:PEtot)

    ! ------------------------- !
    ! --- [1] Update FromTo --- !
    ! ------------------------- !
    !  -- [1-1]  LIg & LIi  --  !
    NewLIs = LIs - ( export_LItem + export_RItem ) &
         &       + ( import_LItem + import_RItem )
    call MPI_AllGather( NewLIs, 1, MPI_INTEGER, &
         &              LIg   , 1, MPI_INTEGER, &
         &              MPI_COMM_WORLD , ierr )
    LIi(      0) = LIg(      0) - 1
    LIi(PEtot-1) = LIg(PEtot-1) - 1
    do iPE=1, PEtot-2
       LIi(iPE)  = LIg(iPE) - 2
    enddo
    !  -- [1-2] iFrom & iTo --  !
    iFrom(0) = 1
    do iPE=1, PEtot
       iFrom(iPE  ) = iFrom(iPE-1) + LIi(iPE-1)
    enddo
    do iPE=0, PEtot-1
       iTo  (iPE  ) = iFrom(iPE+1) - 1
    enddo
    !  -- [1-3] Update      --  !
    do iPE=0, PEtot-1
       FromTo_(iPE,rnk_)  = iPE
       FromTo_(iPE,frm_)  = iFrom(iPE) - 1
       FromTo_(iPE,to_)   = iTo  (iPE) + 1
       FromTo_(iPE,LIs_)  = LIg  (iPE)
       FromTo_(iPE,LJs_)  = LJs
       FromTo_(iPE,iFr_)  = iFrom(iPE)
       FromTo_(iPE,iTo_)  = iTo  (iPE)
       FromTo_(iPE,inn_)  = iTo  (iPE) - iFrom(iPE) + 1
       FromTo_(iPE,isl_)  = 2
       FromTo_(iPE,iel_)  = LIg(iPE) - 1
    enddo
    FromTo_(      0,isl_) = 1
    FromTo_(PEtot-1,iel_) = LIg(PEtot-1)
    FromTo_(      0,frm_) = iFrom(0)
    FromTo_(PEtot-1,to_ ) = iTo(PEtot-1)
    !  -- [3-8] zoffset -- !
    zoffset(1) = dble( FromTo (myRank,iFr_) - 2 )*dz
    zoffset(2) = dble( FromTo_(myRank,iFr_) - 2 )*dz
    
    return
  end subroutine NewLIsSystemMaking


  subroutine UpdateFromTo
    use constants, only : datDir, myRank, PEtot
    use variables, only : FromTo, kstep
    implicit none
    include 'mpif.h'
    integer            :: iPE
    character(  8)     :: cStep
    character(100)     :: FileName

    ! ------------------------- !
    ! --- [1] Update FromTo --- !
    ! ------------------------- !
    do iPE=0, PEtot-1
       FromTo(iPE,rnk_)  = FromTo_(iPE,rnk_)
       FromTo(iPE,frm_)  = FromTo_(iPE,rnk_)
       FromTo(iPE,to_)   = FromTo_(iPE,to_ )
       FromTo(iPE,LIs_)  = FromTo_(iPE,LIs_)
       FromTo(iPE,LJs_)  = FromTo_(iPE,LJs_)
       FromTo(iPE,iFr_)  = FromTo_(iPE,iFr_)
       FromTo(iPE,iTo_)  = FromTo_(iPE,iTo_)
       FromTo(iPE,inn_)  = FromTo_(iPE,inn_)
       FromTo(iPE,isl_)  = FromTo_(iPE,isl_)
       FromTo(iPE,iel_)  = FromTo_(iPE,iel_)
    enddo
    
    ! ------------------------- !
    ! --- [2]   Write out   --- !
    ! ------------------------- !
    if ( myRank.eq.0 ) then
       write(cStep,'(i8.8)') kstep
       FileName = trim(datDir) // 'PEinfo' // cStep //  '.dat'
       open( 50,file=trim(FileName),form='formatted',status='replace' )
       write(50,'(10(a9,1x))') '# PE', 'From', 'To', 'LIloc', 'LJloc', 'iFrom', 'iTo', 'inn', 'isl', 'iel'
       do iPE=0, PEtot-1
          write(50,'(10(i9,1x))') FromTo(iPE,rnk_), FromTo(iPE,frm_), FromTo(iPE,to_ ), &
               &                  FromTo(iPE,LIs_), FromTo(iPE,LJs_), FromTo(iPE,iFr_), FromTo(iPE,iTo_), &
               &                  FromTo(iPE,inn_), FromTo(iPE,isl_), FromTo(iPE,iel_)
       enddo
       close(50)
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       write(6,'(2x,a10,a45,2x,a6)') "* SAVE :: ", trim(FileName), '[ OK ]'
       write(6,*)
    endif

    return
  end subroutine UpdateFromTo



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


    ! if ( myRank.eq.0 ) then
    !    open (40,file='send.dat',form='formatted',status='replace')
    !    write(40,'(6(a8,1x))') 'iPE', 'i', 'nDst', 'Items', 'Index', 'Enemy' 
    !    close(40)
    ! endif
    ! do iPE=0, PEtot-1
    !    if ( iPE.eq.myRank ) then
    !       open (40,file='send.dat',form='formatted',status='old',position='append')
    !       do i=1, nDst
    !          write(40,'(6(i8,1x))') iPE, i, nDst, export_Items(i), export_Index(i), export_Enemy(i)
    !       enddo
    !       close(40)
    !    endif
    !    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    ! enddo

    ! if ( myRank.eq.0 ) then
    !    open (40,file='recv.dat',form='formatted',status='replace')
    !    write(40,'(6(a8,1x))') 'iPE', 'i', 'nSrc', 'Items', 'Index', 'Enemy' 
    !    close(40)
    ! endif
    ! do iPE=0, PEtot-1
    !    if ( iPE.eq.myRank ) then
    !       open (40,file='recv.dat',form='formatted',status='old',position='append')
    !       do i=1, nSrc
    !          write(40,'(6(i8,1x))') iPE, i, nSrc, import_Items(i), import_Index(i), import_Enemy(i)
    !       enddo
    !       close(40)
    !    endif
    !    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    ! enddo


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


    ! if ( myRank.eq.0 ) then
    !    open (40,file='ptLoad.dat',form='formatted',status='replace')
    !    write(40,'(8(a8,1x))') 'k', 'iPE', 'iDst', 'nDst', 'ptItm', 'ptIdx', 'Enemy' , 'ptFrm'
    !    close(40)
    ! endif
    ! do iPE=0, PEtot-1
    !    if ( iPE.eq.myRank ) then
    !       open (40,file='ptLoad.dat',form='formatted',status='old',position='append')
    !       do k=1, ns
    !          do iDst=1, nDst
    !             write(40,'(8(i8,1x))') k, iPE, iDst, nDst, &
    !                  &  export_ptItm(iDst,k), export_ptIdx(iDst,k), export_Enemy(iDst), export_ptFrm(iDst,k)
    !          enddo
    !       enddo
    !       close(40)
    !    endif
    !    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    ! enddo
    ! call MPI_Finalize( ierr )
    ! stop
