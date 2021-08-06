module ptExchgMod
contains

  subroutine prtclExchg( n1, n2, n3, npm, ntype )
    use constants , only : LI, LJ, ns, nptMax
    use constants , only : myRank, PEtot, OMPNumThreads
    use constants , only : dzinv, drinv
    use constants , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables , only : pxv, goban, pCtrlDij
    use sortingMod, only : ParallelSort, sortptcl
    implicit none
    include 'mpif.h'
    integer, intent(in)           :: n1, n2, ntype
    integer, intent(inout)        :: n3, npm(2)
    integer                       :: k, m, iPE, ip, jp, myPE, ierr, idx, domainBdr(PEtot+1)
    integer                       :: nptAdd, nptSub, pBuffSize, sSize, rSize, nptRest, nptFill
    integer                       :: pCnt(PEtot), sCnt(PEtot), sdispls(PEtot), rCnt(PEtot), rdispls(PEtot)
    double precision              :: diPEblock
    integer, allocatable          :: ExchgIdx(:,:), ExchgIdx_(:)
    double precision, allocatable :: pxvS(:), pxvR(:)
    integer, parameter            :: OutIdx = LI*LJ*ns+1

    ! --- [1] Preparation --- !
    pBuffSize = nptMax / PEtot
    !$omp parallel default(none) &
    !$omp shared(pCnt,sCnt,rCnt) private(iPE)
    !$omp do
    do iPE=1, PEtot
       pCnt(iPE) = 0
       sCnt(iPE) = 0
       rCnt(iPE) = 0
    enddo
    !$omp end do
    !$omp end parallel
    diPEblock = dble( PEtot ) / dble( LI-1 )
    allocate( ExchgIdx(pBuffSize,PEtot) )
    do iPE=1, PEtot
       domainBdr(iPE)  = pCtrlDij(1,1,iPE)
    enddo
    domainBdr(PEtot+1) = pCtrlDij(2,1,PEtot)
    ! if ( myRank.eq.1 ) then
    !    do iPE=1, PEtot+1
    !       write(6,*) iPE, domainBdr(iPE)
    !    enddo
       
    !    do ip=1, LI-1
    !       iPE = binary_search( domainBdr, ip, 1, PEtot+1 )
    !       write(6,*) ip, iPE
    !    enddo
    ! endif
    ! stop !!!! --- [DEBUG] --- 
    
    ! --- [2] Packing of Particles -- !
    !  -- [2-1]   Classification   -- !
    myPE   = myRank + 1
    do m=n1, n2
       ip  = max( ceiling( pxv(zp_,m,ntype)*dzinv ), 1 )
       ! iPE = max( ceiling( ip * diPEblock ), 1 )
       iPE = binary_search( domainBdr, ip, 1, PEtot+1 )
       if ( iPE.ne.myPE ) then
          pCnt(iPE)               = pCnt(iPE) + 1
          ExchgIdx(pCnt(iPE),iPE) = m
       end if
    end do
    !  -- [2-2] Count & Displacement for SEND -- !
    nptSub = 0
    !$omp parallel default(none) &
    !$omp shared(sCnt,pCnt,nptSub) private(iPE)
    !$omp do reduction(+:nptSub)
    do iPE=1, PEtot
       sCnt(iPE) = pCnt(iPE)*8
       nptSub    = nptSub + pCnt(iPE)
    enddo
    !$omp end do
    !$omp end parallel
    sdispls(1) = 0
    do iPE=2, PEtot
       sdispls(iPE) = sdispls(iPE-1) + sCnt(iPE-1)
    enddo
    sSize = sdispls(PEtot) + sCnt(PEtot)
    allocate( pxvS( sSize ) )
    !  -- [2-3] Packing into pxvS Buffer -- !
    !$omp parallel default(none) &
    !$omp shared(pCnt,sdispls,pxv,pxvS,goban,ExchgIdx,ntype) private(m,iPE,idx)
    do iPE=1, PEtot
       !$omp do
       do m=1, pCnt(iPE)
          idx                       = sdispls(iPE)+8*(m-1)
          pxvS( 1+idx )             = pxv(rp_,ExchgIdx(m,iPE),ntype )
          pxvS( 2+idx )             = pxv(zp_,ExchgIdx(m,iPE),ntype )
          pxvS( 3+idx )             = pxv(vr_,ExchgIdx(m,iPE),ntype )
          pxvS( 4+idx )             = pxv(vt_,ExchgIdx(m,iPE),ntype )
          pxvS( 5+idx )             = pxv(vz_,ExchgIdx(m,iPE),ntype )
          pxvS( 6+idx )             = pxv(ro_,ExchgIdx(m,iPE),ntype )
          pxvS( 7+idx )             = pxv(zo_,ExchgIdx(m,iPE),ntype )
          pxvS( 8+idx )             = pxv(wp_,ExchgIdx(m,iPE),ntype )
          pxv ( 8,ExchgIdx(m,iPE),ntype ) = 0.d0
          goban(  ExchgIdx(m,iPE) ) = OutIdx
       enddo
       !$omp end do
    enddo
    !$omp end parallel
    
    ! --- [3] Communication --- !
    !  -- [3-1] Exchange the Size of Buffers   -- !
    call MPI_ALLtoALL( sCnt, 1, MPI_INTEGER, rCnt, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr )
    !  -- [3-2] Calculate RECIEVE Displacement -- !
    rdispls(1)   = 0
    do iPE=2, PEtot
       rdispls(iPE) = rdispls(iPE-1) + rCnt(iPE-1)
    enddo
    rSize  = rdispls(PEtot) + rCnt(PEtot)
    nptAdd = rSize   / 8
    allocate( pxvR( rSize ) )
    !  -- [3-3] Exchange Buffers -- !
    call MPI_ALLtoALLv( pxvS, sCnt, sdispls, MPI_DOUBLE_PRECISION, &
         &              pxvR, rCnt, rdispls, MPI_DOUBLE_PRECISION, &
         &              MPI_COMM_WORLD, ierr )
    
    ! --- [4] Storing & npt update --- !
    if ( nptAdd.ge.nptSub ) then
       ! -- OverWrite Region -- !
       !$omp parallel default(none) &
       !$omp shared(pCnt,sdispls,pxv,pxvR,goban,ExchgIdx,drinv,dzinv,ntype) &
       !$omp private(m,ip,jp,iPE,idx)
       do iPE=1, PEtot
          !$omp do
          do m=1, pCnt(iPE)
             idx = sdispls(iPE)+8*(m-1)
             pxv( 1,ExchgIdx(m,iPE),ntype ) = pxvR( 1+idx )
             pxv( 2,ExchgIdx(m,iPE),ntype ) = pxvR( 2+idx )
             pxv( 3,ExchgIdx(m,iPE),ntype ) = pxvR( 3+idx )
             pxv( 4,ExchgIdx(m,iPE),ntype ) = pxvR( 4+idx )
             pxv( 5,ExchgIdx(m,iPE),ntype ) = pxvR( 5+idx )
             pxv( 6,ExchgIdx(m,iPE),ntype ) = pxvR( 6+idx )
             pxv( 7,ExchgIdx(m,iPE),ntype ) = pxvR( 7+idx )
             pxv( 8,ExchgIdx(m,iPE),ntype ) = pxvR( 8+idx )
             ip                       = nint( pxvR( 2+idx )*dzinv )
             jp                       = nint( pxvR( 1+idx )*drinv )
             goban( ExchgIdx(m,iPE) ) = LI*LJ*(ntype-1) + LI*jp + (ip-1) + 1
          enddo
          !$omp end do
       enddo
       !$omp end parallel
       nptRest = nptAdd - nptSub
       nptFill = nptSub
       !  -- [4-1] Whether pxv can store whole new particles  -- !
       if ( nptRest.gt.nptMax-n3 ) then
          do m=1, n3
             if ( goban(m).eq.OutIdx ) then
                pxv(rp_,m,ntype) = pxvR( 1 + nptFill*8 )
                pxv(zp_,m,ntype) = pxvR( 2 + nptFill*8 )
                pxv(vr_,m,ntype) = pxvR( 3 + nptFill*8 )
                pxv(vt_,m,ntype) = pxvR( 4 + nptFill*8 )
                pxv(vz_,m,ntype) = pxvR( 5 + nptFill*8 )
                pxv(ro_,m,ntype) = pxvR( 6 + nptFill*8 )
                pxv(zo_,m,ntype) = pxvR( 7 + nptFill*8 )
                pxv(wp_,m,ntype) = pxvR( 8 + nptFill*8 )
                ip       = nint( pxv(zp_,m,ntype)*dzinv )
                jp       = nint( pxv(rp_,m,ntype)*drinv )
                goban(m) = LI*LJ*(ntype-1) + LI*jp + (ip-1) + 1
                nptFill  = nptFill + 1
                nptRest  = nptRest - 1
                if ( nptRest.eq.0 ) exit
             endif
          enddo
          if ( nptRest.gt.nptMax-n3 ) stop '[ERROR] Insufficient Memory [ERROR]'
       endif
       ! -- Vacant Region -- !
       !$omp parallel default(none) &
       !$omp shared(nptRest,nptFill,pxv,pxvR,goban,drinv,dzinv,n3,ntype) &
       !$omp private(idx,m,ip,jp)
       !$omp do 
       do m=1, nptRest
          idx = nptFill+m-1
          pxv(rp_,n3+m,ntype) = pxvR( 1 + ( idx )*8 )
          pxv(zp_,n3+m,ntype) = pxvR( 2 + ( idx )*8 )
          pxv(vr_,n3+m,ntype) = pxvR( 3 + ( idx )*8 )
          pxv(vt_,n3+m,ntype) = pxvR( 4 + ( idx )*8 )
          pxv(vz_,n3+m,ntype) = pxvR( 5 + ( idx )*8 )
          pxv(ro_,n3+m,ntype) = pxvR( 6 + ( idx )*8 )
          pxv(zo_,n3+m,ntype) = pxvR( 7 + ( idx )*8 )
          pxv(wp_,n3+m,ntype) = pxvR( 8 + ( idx )*8 )
          ip          = nint( pxv(zp_,n3+m,ntype)*dzinv )
          jp          = nint( pxv(rp_,n3+m,ntype)*drinv )
          goban(n3+m) = LI*LJ*(ntype-1) + LI*jp + (ip-1) + 1
       enddo
       !$omp end do
       !$omp end parallel
       n3 = n3 + ( nptAdd-nptSub )
    endif
    if ( nptAdd.lt.nptSub ) then
       k=0
       allocate( ExchgIdx_(nptSub) )
       do iPE=1, PEtot
          do m=1, pCnt(iPE)
             k = k+1
             ExchgIdx_(k) = ExchgIdx(m,iPE)
          enddo
       enddo
       !$omp parallel default(none) &
       !$omp shared(nptAdd,pxv,pxvR,goban,ExchgIdx_,drinv,dzinv,ntype) &
       !$omp private(idx,m,ip,jp) 
       !$omp do 
       do m=1, nptAdd
          idx                   = 8*(m-1)
          pxv( 1,ExchgIdx_(m),ntype ) = pxvR( 1+idx )
          pxv( 2,ExchgIdx_(m),ntype ) = pxvR( 2+idx )
          pxv( 3,ExchgIdx_(m),ntype ) = pxvR( 3+idx )
          pxv( 4,ExchgIdx_(m),ntype ) = pxvR( 4+idx )
          pxv( 5,ExchgIdx_(m),ntype ) = pxvR( 5+idx )
          pxv( 6,ExchgIdx_(m),ntype ) = pxvR( 6+idx )
          pxv( 7,ExchgIdx_(m),ntype ) = pxvR( 7+idx )
          pxv( 8,ExchgIdx_(m),ntype ) = pxvR( 8+idx )
          ip                    = nint( pxv(zp_,ExchgIdx_(m),ntype )*dzinv )
          jp                    = nint( pxv(rp_,ExchgIdx_(m),ntype )*drinv )
          goban( ExchgIdx_(m) ) = LI*LJ*(ntype-1) + LI*jp + (ip-1) + 1
       enddo
       !$omp end do
       !$omp end parallel
       deallocate( ExchgIdx_ )
    endif

    npm(1) = npm(1) + nptAdd
    npm(2) = npm(2) + nptSub

    return
  end subroutine PrtclExchg
  
  
  
  recursive function binary_search( array, key, imin, imax ) result( ret_binary_search )

    implicit none
    integer, intent(in) :: array(imax-imin+1)
    integer, intent(in) :: key, imin, imax
    integer             :: imid, ret_binary_search

    if ( imax.lt.imin ) then
       stop 'ERROR'
    else
       if ( imax-imin.le.1 ) then
          ret_binary_search = imin
          return
       else
          imid = imin + ceiling( ( imax - imin ) / 2.0 )
          if      ( key.ge.array(imid) ) then
             ret_binary_search = binary_search( array, key, imid, imax )
          else
             ret_binary_search = binary_search( array, key, imin, imid )
          end if
          return
       end if
    end if
    return
  end function binary_search

end module ptExchgMod



! !  -- [DEBUG] Leak Checker [DEBUG] -- !
! pLeak = 0
! do m=1, nptAdd
!    ip  = max( ceiling( pxvR( 2+8*(m-1) )*dzinv ), 1 )
!    iPE = max( ceiling( ip * diPEblock          ), 1 )
!    if ( iPE.ne.myPE ) then
!       pLeak = pLeak + 1
!    end if
! end do
! write(6,'(2(i8,1x))') myPE, pLeak
!  ---------------------------------- !
