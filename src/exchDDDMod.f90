module exchDDDMod
contains

  subroutine FieldExchange
    use constants , only : LIs, LJs, myRank
    use variables , only : EBf
    use com4DDDMod
    implicit none
    include 'mpif.h'
    integer             :: i, j, cmp, idx1, idx2, isl
    integer             :: iDst, iSrc, is, ir, len_s, len_r, ierr
    
    ! ------------------------- !
    ! --- [1]  Packing      --- !
    ! ------------------------- !
    !  -- [1-1] Allocation  --  !
    allocate( EBfSendBuff(ebCmp,LJs,nBuffLine_Max) )
    allocate( EBfRecvBuff(ebCmp,LJs,nBuffLine_Max) )
    allocate( FieldStack (ebCmp,LIs,LJs)           )
    !  -- [1-2] ebSendBuff  --  !
    !$omp parallel default(none) &
    !$omp shared(EBfSendBuff,nDst,export_Items,LJs,export_Index,export_iFrom,EBf) &
    !$omp private(iDst,i,j,cmp,idx1,idx2)
    !$omp do
    do iDst=1, nDst
       do i=1, export_Items(iDst)
          idx1                    = export_Index(iDst) + (i-1)
          idx2                    = export_iFrom(iDst) + (i-1)
          do j=1, LJs
             do cmp=1, ebCmp
                EBfSendBuff(cmp,j,idx1) = EBf(cmp,idx2,j)
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
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
    !$omp parallel default(none) &
    !$omp shared(LIs,LJs,FieldStack,EBf) private(i,j,cmp)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          do cmp=1, ebCmp
             FieldStack(cmp,i,j) = EBf(cmp,i,j)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  ------ LIs Update ------ !
    LIs = NewLIs
    deallocate( EBf )
    allocate( EBf(ebCmp,LIs,LJs) )
    EBf = 0.d0
    !  ------ LIs Update ------ !
    !  -- [3-3] Restore EBf --  !
    isl = ijDomain_(myRank,isl_)
    !$omp parallel default(none) &
    !$omp shared(LJs,LIRemnant,isl,import_LItem,export_LItem,EBf,FieldStack) &
    !$omp private(i,j,cmp,idx1,idx2)
    !$omp do
    do j=1, LJs
       do i=1, LIRemnant
          idx1            = isl + import_LItem + (i-1)
          idx2            = isl + export_LItem + (i-1)
          do cmp=1, ebCmp
             EBf(cmp,idx1,j) = FieldStack(cmp,idx2,j)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [3-3] Update EBf  --  !
    !$omp parallel default(none) &
    !$omp shared(nSrc,import_Items,LJs,import_ebIdx,import_Index,EBf,EBfRecvBuff) &
    !$omp private(iSrc,i,j,idx1,idx2,cmp)
    !$omp do
    do iSrc=1, nSrc
       do i=1, import_Items(iSrc)
          idx1            = import_ebIdx(iSrc) + (i-1)
          idx2            = import_Index(iSrc) + (i-1)
          do j=1, LJs
             do cmp=1, ebCmp
                EBf(cmp,idx1,j) = EBfRecvBuff(cmp,j,idx2)
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine FieldExchange


  subroutine ParticleExchange
    use constants , only : ns, np, npt, nptMax, ptExpandRatio, myRank
    use constants , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables , only : pxv
    use ptOrderMod, only : ptReOrdering
    use pMPIComMod, only : ResizeParticleMemory
    use com4DDDMod
    implicit none
    include 'mpif.h'
    integer             :: m, k, cmp, iDst, iSrc, idx1, idx2
    integer             :: is, ir, len_s, len_r, ierr
    integer             :: npkMax, npNew(ns), sBuffSize, rBuffSize

    ! ------------------------- !
    ! --- [1]  Packing      --- !
    ! ------------------------- !
    !  -- [1-1] Size check  --  !
    npNew(:) = np(:)
    do k=1, ns
       npNew(k) = npNew(k) + import_npkPE(k)
    enddo
    npkMax = max( npNew(1), npNew(2) )
    if ( nptMax.lt.npkMax ) then
       npkMax = nint( ptExpandRatio * dble(npkMax) )
       write(6,'(a,i6,2(i12,1x))') "[ResizeParticleMemory] ", myRank, nptMax, npkMax
       call ResizeParticleMemory( nptMax, npkMax )
    endif
    !  -- [1-2] Allocation  --  !
    sBuffSize = maxval( export_npkPE(:) )
    rBuffSize = maxval( import_npkPE(:) )
    allocate( ptSendBuff(ptCmp,sBuffSize,ns) )
    allocate( ptRecvBuff(ptCmp,rBuffSize,ns) )
    !  -- [1-3] Global zp   --  !
    !$omp parallel default(none) &
    !$omp shared(np,pxv,zoffset) private(k,m)
    !$omp do
    do k=1, ns
       do m=1, np(k)
          pxv(zp_,m,k) = pxv(zp_,m,k) + zoffset(1)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [1-4] ptSendBuff  --  !
    !$omp parallel default(none) &
    !$omp shared(nDst,export_ptItm,export_ptIdx,export_ptFrm,ptSendBuff,pxv) &
    !$omp private(k,iDst,m,idx1,idx2)
    !$omp do
    do k=1, ns
       do iDst=1, nDst
          do m=1, export_ptItm(iDst,k)
             idx1 = export_ptIdx(iDst,k) + (m-1)
             idx2 = export_ptFrm(iDst,k) + (m-1)
             ptSendBuff(rp_,idx1,k) = pxv(rp_,idx2,k)
             ptSendBuff(zp_,idx1,k) = pxv(zp_,idx2,k)
             ptSendBuff(vr_,idx1,k) = pxv(vr_,idx2,k)
             ptSendBuff(vt_,idx1,k) = pxv(vt_,idx2,k)
             ptSendBuff(vz_,idx1,k) = pxv(vz_,idx2,k)
             ptSendBuff(ro_,idx1,k) = pxv(ro_,idx2,k)
             ptSendBuff(zo_,idx1,k) = pxv(zo_,idx2,k)
             ptSendBuff(wp_,idx1,k) = pxv(wp_,idx2,k)
             pxv(wp_,idx2,k)        = 0.d0
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! ------------------------- !
    ! --- [2]  Exchange     --- !
    ! ------------------------- !
    do k=1, ns
       !  -- [2-1] Send     --  !
       do iDst=1, nDst
          is    = export_ptIdx(iDst,k)
          len_s = export_ptItm(iDst,k) * ptCmp
          call MPI_iSend( ptSendBuff(1,is,k), len_s, MPI_DOUBLE_PRECISION, &
               &          export_enemy(iDst), 0    , MPI_COMM_WORLD      , &
               &          Request_Send(iDst), ierr )
       enddo
       !  -- [2-2] Recv     --  !
       do iSrc=1, nSrc
          ir    = import_ptIdx(iSrc,k)
          len_r = import_ptItm(iSrc,k) * ptCmp
          call MPI_iRecv( ptRecvBuff(1,ir,k), len_r, MPI_DOUBLE_PRECISION, &
               &          import_enemy(iSrc), 0    , MPI_COMM_WORLD      , &
               &          Request_Recv(iSrc), ierr )
       enddo
       !  -- [2-3] Update   --  !
       call MPI_WaitAll( nSrc, Request_Recv, Status_Recv, ierr )
       call MPI_WaitAll( nDst, Request_Send, Status_Send, ierr )
       !$omp parallel default(none) &
       !$omp shared(k,nSrc,import_ptItm,import_ptDst,import_ptIdx,pxv,ptRecvBuff) &
       !$omp private(iSrc,m,idx1,idx2,cmp)
       !$omp do
       do iSrc=1, nSrc
          do m=1, import_ptItm(iSrc,k)
             idx1 = import_ptDst(iSrc,k) + (m-1)
             idx2 = import_ptIdx(iSrc,k) + (m-1)
             do cmp=1, ptCmp
                pxv(cmp,idx1,k) = ptRecvBuff(cmp,idx2,k)
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    
    ! ------------------------- !
    ! --- [3] Post Process  --- !
    ! ------------------------- !
    !  -- [3-1] update np   --  !
    do k=1, ns
       np(k) = npNew(k)
    enddo
    npt = np(1) + np(2)
    !  -- [3-2] Local zp    --  !
    !$omp parallel default(none) &
    !$omp shared(np,pxv,zoffset) private(k,m)
    !$omp do
    do k=1, ns
       do m=1, np(k)
          pxv(zp_,m,k) = pxv(zp_,m,k) - zoffset(2)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [3-3] Sorting     --  !
    call ptReOrdering( 'iGrid' )
    
    return
  end subroutine ParticleExchange


end module exchDDDMod
