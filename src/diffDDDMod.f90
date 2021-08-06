module diffDDDMod
contains

  ! =================================================================== !
  ! ===  DiffuseDecomposition  :: Diffusion like DDD update         === !
  ! =================================================================== !
  subroutine DiffuseDecomposition
    use variables , only : EBf
    use fMPIComMod, only : FieldCommunicate
    use exchDDDMod, only : FieldExchange, ParticleExchange
    use com4DDDMod
    implicit none

    ! ------------------------------------- !
    ! --- [1] Compare npt between PE    --- !
    ! ------------------------------------- !
    call npcRedistribute

    ! ------------------------------------- !
    ! --- [2] Exchange Columns          --- !
    ! ------------------------------------- !
    if ( nColumnExport.gt.nColumnThreshold ) then
       call NewLIsSystemMaking
       call FieldExchange
       call ParticleExchange
       call UpdateijDomain
       call reAllocateVariables
       call reDefine_x1Geometry
       write(6,*) 'here 1'
       call reDefine_CommConfig
       write(6,*) 'here 2'
       call BcastDisplayNptInfo
       write(6,*) 'here 3'
       call deallocateVariables
       write(6,*) 'here 4'
       call FieldCommunicate( EBf, 'E' )
       write(6,*) 'here 5'
       call FieldCommunicate( EBf, 'B' )
       write(6,*) 'here 6'
    else
       deallocate( request_Send, status_Send )
       deallocate( request_Recv, status_Recv )
    endif

    return
  end subroutine DiffuseDecomposition


  ! =================================================================== !
  ! ===  npcRedistribute  ::  Judge Redistribution                  === !
  ! =================================================================== !
  subroutine npcRedistribute
    use constants , only : LI, LJ, LIs, ns, np, dzInv, OMPNumThreads, myRank, PEtot
    use constants , only : zp_, wp_
    use variables , only : pxv, ijDomain
    use ptOrderMod, only : ptReOrdering
    use com4DDDMod
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: i, m, k, ip, ith, iPE, sendCnt, ierr
    integer             :: exMaxLine, imMaxLine1, imMaxLine2 
    integer             :: iDst, iSrc, isl, iel, iFr, iTo, nCol1, nCol2
    integer             :: recvCnt(0:PEtot-1), displs(0:PEtot-1)
    integer             :: nptPE(-1:PEtot), npLW(LIs+2,OMPNumThreads)
    integer             :: npLD(LIs,ns), npLs(LIs,ns), npLg(LI), npLk(LI,ns)
    double precision    :: del1, del2

    ! ----------------------------------------------- !
    ! --- [1] Sort Particle according to  i-Grid  --- !
    ! ----------------------------------------------- !
    call ptReOrdering( 'iGrid' )

    ! ----------------------------------------------- !
    ! --- [2] Count up #.of pt. on ith-Grid Line  --- !
    ! ----------------------------------------------- !
    !  -- [2-1] Count up k-th species :: npLs     --  !
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
    isl     = ijDomain(myRank,isl_)
    iel     = ijDomain(myRank,iel_)
    sendCnt = ijDomain(myRank,inn_)
    do iPE=0, PEtot-1
       recvCnt(iPE) = ijDomain(iPE,inn_)
       displs (iPE) = ijDomain(iPE,iFr_) - 1
    enddo
    npLk    = 0
    do k=1, ns
       call MPI_AllGatherV( npLs(isl:iel,k), sendCnt,         MPI_INTEGER, &
            &               npLk(      :,k), recvCnt, displs, MPI_INTEGER, &
            &               MPI_COMM_WORLD , ierr                          )
    enddo
    !$omp parallel default(none) shared(npLg,npLk) private(i)
    !$omp do
    do i=1, LI
       npLg(i) = npLk(i,1) + npLk(i,2)
    enddo
    !$omp end do
    !$omp end parallel

    ! ----------------------------------------------- !
    ! --- [3] Judge Exchange ( Import / Export )  --- !
    ! ----------------------------------------------- !
    !  -- [3-1] npt in all PE                     --  !
    do iPE=0, PEtot-1
       iFr        = ijDomain(iPE,iFr_)
       iTo        = ijDomain(iPE,iTo_)
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
    nCol1         = int( del1 / DDDcritStep )
    nCol2         = int( del2 / DDDcritStep )
    exMaxLine     = min( DDDnLine_Max, ijDomain(myRank  ,inn_) / 4 )
    imMaxLine1    = min( DDDnLine_Max, ijDomain(myRank-1,inn_) / 4 )
    imMaxLine2    = min( DDDnLine_Max, ijDomain(myRank+1,inn_) / 4 )
    if ( nCol1.ge. 1 ) then
       nDst               = nDst   + 1
       export_Enemy(nDst) = myRank - 1
       export_LItem       = min( abs(nCol1), exMaxLine  )
       export_Items(nDst) = export_LItem
       export_iFrom(nDst) = isl
    endif
    if ( nCol1.le.-1 ) then
       nSrc               = nSrc   + 1
       import_Enemy(nSrc) = myRank - 1
       import_LItem       = min( abs(nCol1), imMaxLine1 )
       import_Items(nSrc) = import_LItem
       import_gFrom(nSrc) = ijDomain(myRank,iFr_) - 1 - import_Items(nSrc) + 1
    endif
    if ( nCol2.ge. 1 ) then
       nDst               = nDst   + 1
       export_Enemy(nDst) = myRank + 1
       export_RItem       = min( abs(nCol2), exMaxLine  )
       export_Items(nDst) = export_RItem
       export_iFrom(nDst) = iel - export_Items(nDst) + 1
    endif
    if ( nCol2.le.-1 ) then
       nSrc               = nSrc   + 1
       import_Enemy(nSrc) = myRank + 1
       import_RItem       = min( abs(nCol2), imMaxLine2 )
       import_Items(nSrc) = import_RItem
       import_gFrom(nSrc) = ijDomain(myRank,iTo_) + 1
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
    !  -- [3-5] #.of Column & Index               --  !
    LIRemnant       = ijDomain(myRank,inn_) - ( export_LItem + export_RItem )
    if      ( ( import_LItem.eq.0 ).and.( import_RItem.eq.0 ) ) then
       import_ebIdx(0) = ijDomain(myRank,isl_) + import_LItem
    else if ( ( import_LItem.eq.0 ).and.( import_RItem.gt.0 ) ) then
       import_ebIdx(0) = ijDomain(myRank,isl_) + import_LItem
       import_ebIdx(1) = ijDomain(myRank,isl_) + import_LItem + LIRemnant
    else if ( ( import_LItem.gt.0 ).and.( import_RItem.eq.0 ) ) then
       import_ebIdx(1) = ijDomain(myRank,isl_)
       import_ebIdx(0) = ijDomain(myRank,isl_) + import_LItem
    else if ( ( import_LItem.gt.0 ).and.( import_RItem.gt.0 ) ) then
       import_ebIdx(1) = ijDomain(myRank,isl_)
       import_ebIdx(0) = ijDomain(myRank,isl_) + import_LItem
       import_ebIdx(2) = ijDomain(myRank,isl_) + import_LItem + LIRemnant
    endif
    !  -- [3-6] #.of Column in exchange           --  !
    sendCnt = 0
    do iDst=1, nDst
       sendCnt = sendCnt + export_Items(iDst)
    enddo
    call MPI_AllReduce( sendCnt, nColumnExport , 1, MPI_INTEGER, &
         &              MPI_SUM, MPI_COMM_WORLD, ierr            )
    nColumnThreshold = PEtot * DDDcolThreshold
    
    ! ----------------------------------------------- !
    ! --- [4] Generalized Comm.Table ( particle ) --- !
    ! ----------------------------------------------- !
    !  -- [4-1] Index for Particle  ( export )    --  !
    export_ptItm = 0
    export_ptIdx = 1
    export_npkPE = 0
    do k=1, ns
       do iDst=1, nDst
          iFr = export_iFrom(iDst)
          iTo = export_iFrom(iDst) + export_Items(iDst) - 1
          export_ptItm(iDst,k) = sum( npLs( iFr:iTo, k ) )
       enddo
       do iDst=1, nDst
          export_ptIdx(iDst,k) = export_ptIdx(iDst-1,k) + export_ptItm(iDst-1,k)
          export_npkPE(k)      = export_npkPE(k)        + export_ptItm(iDst  ,k)
       enddo
    enddo
    !  -- [4-2] Index for Particle  ( import )    --  !
    import_ptItm = 0
    import_ptIdx = 1
    import_npkPE = 0
    do k=1, ns
       do iSrc=1, nSrc
          iFr = import_gFrom(iSrc)
          iTo = import_gFrom(iSrc) + import_Items(iSrc) - 1
          import_ptItm(iSrc,k) = sum( npLk( iFr:iTo, k ) )
       enddo
       do iSrc=1, nSrc
          import_ptIdx(iSrc,k) = import_ptIdx(iSrc-1,k) + import_ptItm(iSrc-1,k)
          import_npkPE(k)      = import_npkPE(k)        + import_ptItm(iSrc  ,k)
       enddo
    enddo
    !  -- [4-3] Pack  particle from this Index    --  !
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
    !  -- [4-4] Store particle from this Index    --  !
    do k=1, ns
       import_ptDst(0,k) = np(k)+1
       do iSrc=1, nSrc
          import_ptDst(iSrc,k) = import_ptDst(iSrc-1,k) + import_ptItm(iSrc-1,k)
       enddo
    enddo
    
    ! ----------------------------------------------- !
    ! --- [5]    Request-Status Allocation        --- !
    ! ----------------------------------------------- !
    allocate( request_Send(nDst), status_Send(MPI_STATUS_SIZE,nDst) )
    allocate( request_Recv(nSrc), status_Recv(MPI_STATUS_SIZE,nSrc) )

    return
  end subroutine npcRedistribute


  ! =================================================================== !
  ! ===  NewLIsSystemMaking  :: Rearrange #.of column in i-direc.   === !
  ! =================================================================== !
  subroutine NewLIsSystemMaking
    use constants, only : LIs, LJs, dz, myRank, PEtot
    use variables, only : ijDomain
    use com4DDDMod
    implicit none
    include 'mpif.h'
    integer            :: iPE, ierr
    integer            :: LIg(0:PEtot-1), LIi(0:PEtot-1)
    integer            :: iFrom(0:PEtot), iTo(0:PEtot)

    ! ------------------------------------- !
    ! --- [1] Update ijDomain           --- !
    ! ------------------------------------- !
    !  -- [1-1]  LIg & LIi              --  !
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
    !  -- [1-2] iFrom & iTo             --  !
    iFrom(0) = 1
    do iPE=1, PEtot
       iFrom(iPE  ) = iFrom(iPE-1) + LIi(iPE-1)
    enddo
    do iPE=0, PEtot-1
       iTo  (iPE  ) = iFrom(iPE+1) - 1
    enddo
    !  -- [1-3] Update                  --  !
    do iPE=0, PEtot-1
       ijDomain_(iPE,rnk_)  = iPE
       ijDomain_(iPE,frm_)  = iFrom(iPE) - 1
       ijDomain_(iPE,to_)   = iTo  (iPE) + 1
       ijDomain_(iPE,LIs_)  = LIg  (iPE)
       ijDomain_(iPE,LJs_)  = LJs
       ijDomain_(iPE,iFr_)  = iFrom(iPE)
       ijDomain_(iPE,iTo_)  = iTo  (iPE)
       ijDomain_(iPE,inn_)  = iTo  (iPE) - iFrom(iPE) + 1
       ijDomain_(iPE,isl_)  = 2
       ijDomain_(iPE,iel_)  = LIg(iPE) - 1
    enddo
    ijDomain_(      0,isl_) = 1
    ijDomain_(PEtot-1,iel_) = LIg(PEtot-1)
    ijDomain_(      0,frm_) = iFrom(0)
    ijDomain_(PEtot-1,to_ ) = iTo(PEtot-1)
    !  -- [3-8] zoffset -- !
    zoffset(1)  = dble( ijDomain (myRank,iFr_) - 2 )*dz
    zoffset(2)  = dble( ijDomain_(myRank,iFr_) - 2 )*dz
    if ( myRank.eq.0 ) zoffset(:) = 0.d0

    return
  end subroutine NewLIsSystemMaking

end module diffDDDMod

