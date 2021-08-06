module diffDDDMod
  use com4DDDMod
  integer :: npLk_(LI,ns)
contains


  subroutine DiffuseDecomposition
    use variables , only : EBf
    use fMPIComMod, only : FieldCommunicate
    use exchDDDMod, only : FieldExchange, ParticleExchange
    implicit none
    ! integer :: ierr
    
    call npcRedistribute
    call BcastDisplayNptInfo
    call NewLIsSystemMaking
    call FieldExchange
    call ParticleExchange
    call UpdateFromTo
    call reAllocateVariables
    call reDefine_x1Geometry
    call reDefine_CommConfig
    call BcastDisplayNptInfo
    call deallocateVariables
    call FieldCommunicate( EBf, 'E' )
    call FieldCommunicate( EBf, 'B' )

    return
  end subroutine DiffuseDecomposition


  subroutine npcRedistribute
    use constants , only : LI, LJ, LIs, ns, np, dzInv, OMPNumThreads, myRank, PEtot
    use constants , only : zp_, wp_
    use variables , only : pxv, FromTo
    use ptOrderMod, only : ptReOrdering
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer             :: i, m, k, ip, ith, iPE, sendCnt, ierr
    integer             :: iDst, iSrc, isl, iel, iFr, iTo
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
    !$omp parallel default(none) shared(npLg,npLk) private(i)
    !$omp do
    do i=1, LI
       npLg(i) = npLk(i,1) + npLk(i,2)
    enddo
    !$omp end do
    !$omp end parallel
    npLk_ = npLk

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
    !  -- [3-5] #.of Column & Index               --  !
    LIRemnant       = FromTo(myRank,inn_) - ( export_LItem + export_RItem )
    if      ( ( import_LItem.eq.0 ).and.( import_RItem.eq.0 ) ) then
       import_ebIdx(0) = FromTo(myRank,isl_) + import_LItem
    else if ( ( import_LItem.eq.0 ).and.( import_RItem.gt.0 ) ) then
       import_ebIdx(0) = FromTo(myRank,isl_) + import_LItem
       import_ebIdx(1) = FromTo(myRank,isl_) + import_LItem + LIRemnant
    else if ( ( import_LItem.gt.0 ).and.( import_RItem.eq.0 ) ) then
       import_ebIdx(1) = FromTo(myRank,isl_)
       import_ebIdx(0) = FromTo(myRank,isl_) + import_LItem
    else if ( ( import_LItem.gt.0 ).and.( import_RItem.gt.0 ) ) then
       import_ebIdx(1) = FromTo(myRank,isl_)
       import_ebIdx(0) = FromTo(myRank,isl_) + import_LItem
       import_ebIdx(2) = FromTo(myRank,isl_) + import_LItem + LIRemnant
    endif
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
    zoffset(1)  = dble( FromTo (myRank,iFr_) - 2 )*dz
    zoffset(2)  = dble( FromTo_(myRank,iFr_) - 2 )*dz
    if ( myRank.eq.0 ) zoffset(:) = 0.d0

    return
  end subroutine NewLIsSystemMaking

end module diffDDDMod


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


! if ( myRank.eq.3 ) then
!    do k=1, ns
!       do m=1, 10
!          write(6,*) myRank, k, m, ceiling( pxv(zp_,m,k) * dzInv ) - 1, pxv(zp_,m,k), pxv(wp_,m,k)
!       enddo
!       do m=np(k)-10, np(k)
!          write(6,*) myRank, k, m, ceiling( pxv(zp_,m,k) * dzInv ) - 1, pxv(zp_,m,k), pxv(wp_,m,k)
!       enddo
!       do m=1, np(k)
!          iPE = ceiling( pxv(zp_,m,k) * dzInv ) - 1
!          if ( iPE.gt.is ) is = iPE
!          if ( iPE.lt.ir ) ir = iPE
!       enddo
!       write(6,*) k, is, ir, LIs
!    enddo
! endif

! if ( myRank.eq.0 ) then
!    open (40,file='npt.dat',form='formatted',status='replace')
!    write(40,'(12(a8,1x))') 'np1(1)', 'np1(2)', 'npc1(1)', 'npc1(2)', &
!         &                  'np2(1)', 'np2(2)', 'npc2(1)', 'npc2(2)', &
!         &                  'np3(1)', 'np3(2)', 'npc3(1)', 'npc3(2)'
!    close(40)
! endif
! do iPE=0, PEtot-1
!    if ( iPE.eq.myRank ) then
!       open (40,file='npt.dat',form='formatted',status='old',position='append')
!       write(40,'(12(i8,1x))') npV(1), npV(2), npV(3), npV(4), &
!            &                  npV(5), npV(6), npV(7), npV(8), &
!            &                  npV(9), npV(10), npV(11), npV(12)
!       close(40)
!    endif
!    call MPI_Barrier( MPI_COMM_WORLD, ierr )
! enddo

! use constants, only : LIs, LJs, cRank, myRank
! use variables, only : EBf, FromTo
! integer            :: i, j 
! character(100)     :: FileName

! do j=1, LJs
!    do i=1, LIs
!       EBf(5,i,j) = dble(i) + FromTo(myRank,frm_) - 1
!    enddo
! enddo
! FileName = 'bin/out1_' // cRank // '.bin'
! open( 50, file=trim(FileName), form='unformatted' )
! write(50) EBf(5,:,:)
! close(50)

! FileName = 'bin/out2_' // cRank // '.bin'
! open( 50, file=trim(FileName), form='unformatted' )
! write(50) EBf(5,:,:)
! close(50)


! subroutine transferInvestigator
  !   use constants, only : myRank, PEtot, np, ns
  !   implicit none
  !   include 'mpif.h'
  !   integer            :: iPE, ierr, isrc, k, iFr, iTo, iDst
  !   integer            :: items(4), gitems(4,0:PEtot-1), npNew(2), npto, nptn, gnpto(0:PEtot-1), gnptn(0:PEtot-1)
  !   integer            :: sLnpt, sRnpt, rLnpt, rRnpt
  !   integer            :: gsLnpt(0:PEtot-1), gsRnpt(0:PEtot-1)
  !   integer            :: grLnpt(0:PEtot-1), grRnpt(0:PEtot-1)
  !   integer            :: import_ptItm_(0:PEtot-1,ns)   , import_ptIdx_(0:PEtot-1,ns)

  !   items(1) = export_LItem
  !   items(2) = export_RItem
  !   items(3) = import_LItem
  !   items(4) = import_RItem
  !   sLnpt    = 0
  !   sRnpt    = 0
  !   rLnpt    = 0
  !   rRnpt    = 0
  !   if      ( nDst.eq.1 ) then
  !      if      ( export_LItem.eq.1 ) then
  !         sLnpt = export_ptItm(1,1) + export_ptItm(1,2)
  !      else if ( export_RItem.eq.1 ) then
  !         sRnpt = export_ptItm(1,1) + export_ptItm(1,2)
  !      endif
  !   else if ( nDst.eq.2 ) then
  !      sLnpt = export_ptItm(1,1) + export_ptItm(1,2)
  !      sRnpt = export_ptItm(2,1) + export_ptItm(2,2)
  !   endif
  !   if      ( nSrc.eq.1 ) then
  !      if      ( import_LItem.eq.1 ) then
  !         rLnpt = import_ptItm(1,1) + import_ptItm(1,2)
  !      else if ( import_RItem.eq.1 ) then
  !         rRnpt = import_ptItm(1,1) + import_ptItm(1,2)
  !      endif
  !   else if ( nSrc.eq.2 ) then
  !      rLnpt = import_ptItm(1,1) + import_ptItm(1,2)
  !      rRnpt = import_ptItm(2,1) + import_ptItm(2,2)
  !   endif
  !   npNew(:) = np(:)
  !   do k=1, ns
  !      do iSrc=1, nSrc
  !         npNew(k) = npNew(k) + import_ptItm(iSrc,k)
  !      enddo
  !   enddo
  !   npto = np(1) + np(2)
  !   nptn = npNew(1) + npNew(2)

  !   import_ptItm_ = 0
  !   import_ptIdx_ = 1
  !   do k=1, ns
  !      do iSrc=1, nSrc
  !         iFr = import_gFrom(iSrc)
  !         iTo = import_gFrom(iSrc) + import_Items(iSrc) - 1
  !         import_ptItm_(iSrc,k) = sum( npLk_( iFr:iTo, k ) )
  !         ! write(6,'(6(i10,1x))') myRank, iFr, iTo, import_ptItm_(iSrc,k), npLK_(iFr,k), npLK_(iTo,k)
  !      enddo
  !      do iSrc=1, nSrc
  !         import_ptIdx_(iSrc,k) = import_ptIdx_(iSrc-1,k) + import_ptItm_(iSrc-1,k)
  !      enddo
  !   enddo

  !   ! if ( myRank.eq.0 ) then
  !   !    open (40,file='send.dat',form='formatted',status='replace')
  !   !    write(40,'(7(a8,1x))') 'k', 'iPE', 'iDst', 'nDst', 'ptItm', 'ptIdx', 'Enemy' 
  !   !    close(40)
  !   ! endif
  !   ! do iPE=0, PEtot-1
  !   !    if ( iPE.eq.myRank ) then
  !   !       open (40,file='send.dat',form='formatted',status='old',position='append')
  !   !       do k=1, ns
  !   !          do iDst=1, nDst
  !   !             write(40,'(7(i8,1x))') k, iPE, iDst, nDst, export_ptItm(iDst,k), &
  !   !                  &                 export_ptIdx(iDst,k), export_Enemy(iDst)
  !   !          enddo
  !   !       enddo
  !   !       close(40)
  !   !    endif
  !   !    call MPI_Barrier( MPI_COMM_WORLD, ierr )
  !   ! enddo
  !   ! if ( myRank.eq.0 ) then
  !   !    open (40,file='recv.dat',form='formatted',status='replace')
  !   !    write(40,'(8(a8,1x))') 'k', 'iPE', 'iDst', 'nDst', 'ptItm', 'ptIdx', 'Enemy' , 'ptFrm'
  !   !    close(40)
  !   ! endif
  !   ! do iPE=0, PEtot-1
  !   !    if ( iPE.eq.myRank ) then
  !   !       open (40,file='recv.dat',form='formatted',status='old',position='append')
  !   !       do k=1, ns
  !   !          do iSrc=1, nSrc
  !   !             write(40,'(8(i8,1x))') k, iPE, iSrc, nSrc, &
  !   !                  &  import_ptItm(iSrc,k), import_ptIdx(iSrc,k), import_Enemy(iSrc)
  !   !          enddo
  !   !       enddo
  !   !       close(40)
  !   !    endif
  !   !    call MPI_Barrier( MPI_COMM_WORLD, ierr )
  !   ! enddo
  !   ! call MPI_Finalize( ierr )
  !   ! stop

    

  !   call MPI_AllGather( items, 4, MPI_INTEGER, gitems, 4, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  !   call MPI_AllGather( sLnpt, 1, MPI_INTEGER, gsLnpt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  !   call MPI_AllGather( sRnpt, 1, MPI_INTEGER, gsRnpt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  !   call MPI_AllGather( rLnpt, 1, MPI_INTEGER, grLnpt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  !   call MPI_AllGather( rRnpt, 1, MPI_INTEGER, grRnpt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  !   call MPI_AllGather( npto , 1, MPI_INTEGER, gnpto , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  !   call MPI_AllGather( nptn , 1, MPI_INTEGER, gnptn , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
  !   if ( myRank.eq.0 ) then
  !      write(6,*) '----------------------------------------------------------------------'
  !      do iPE=0, PEtot-1
  !         write(6,'(5(i10,1x))') iPE, gitems(1,iPE), gitems(2,iPE), gitems(3,iPE), gitems(4,iPE)
  !      enddo
  !      write(6,*) '----------------------------------------------------------------------'
  !      do iPE=0, PEtot-1
  !         write(6,'(7(i10,1x))') iPE, gsLnpt(iPE), gsRnpt(iPE), grLnpt(iPE), grRnpt(iPE), gnpto(iPE), gnptn(iPE)
  !      enddo
  !      write(6,*) '----------------------------------------------------------------------'
  !   endif

  !   return
  ! end subroutine transferInvestigator
  
