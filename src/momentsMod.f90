module momentsMod
  implicit none
  logical         , parameter     :: FoldingBoundary = .true.
  integer         , parameter     :: nMoments        = 16
  integer         , parameter     :: nPrsQth         = 12
  integer         , parameter     :: nei_ = 1
  integer         , parameter     :: ur_  = 2,  ut_  = 3,  uz_  = 4
  integer         , parameter     :: pxx_ = 5,  pyy_ = 6,  pzz_ = 7
  integer         , parameter     :: pxy_ = 8,  pyz_ = 9,  pxz_ =10
  integer         , parameter     :: pt1_ =11,  pn1_ =12,  pn2_ =13
  integer         , parameter     :: qv1_ =14,  qv2_ =15,  qv3_ =16
  double precision, allocatable   :: fMoments(:,:,:,:)
contains
  
  subroutine MomentAnalysis
    use constants , only : LIs, LJs, ns, binDir, cRank, myRank
    use variables , only : kstep
    use p2ndMomMod, only : p2ndMoments
    implicit none
    include 'mpif.h'
    integer            :: ierr
    character(8)       :: cStep
    character(100)     :: FileName, binDirh
    ! --- [1] Measurements    --- !
    !  -- [1-1] Density       --  !
    call nCountUp_CC( fMoments(1:LIs,1:LJs,1:ns,nei_) )
    
    ! -- [1-2] Flow Velocity   -- !
    call u1stMoments
    ! -- [1-3] Pressure Tensor -- !
    call p2ndMoments( fMoments )
    
    ! --- [2] Write out -- !
    !  -- [2-1] FileName       -- !
    write(cStep,'(i8.8)') kstep
    binDirh  = trim(binDir)  // 'kstep'  // cStep // '/'
    FileName = trim(binDirh) // 'Moment' // cStep // '_' // cRank // '.bin'
    !  -- [2-2] Make Directory -- !
    if ( myRank.eq.0 ) call system( 'mkdir -p ' // trim(binDirh) )
    call MPI_Barrier( MPI_COMM_WORLD, ierr )
    !  -- [2-3] Write File     -- !
    open (50, file=trim(FileName),  form= 'unformatted',  &
         &    status= 'replace', convert= 'LITTLE_ENDIAN' )
    write(50) fMoments
    close(50)
    
    return
  end subroutine MomentAnalysis

  
  subroutine nCountUp_CC( sumup )
    use constants , only           : LIs, LJs, ns, np, OMPNumThreads
    use constants , only           : np, dr, dz, drInv, dzInv, rp_, zp_, wp_
    use variables , only           : pxv, RhoPerNsp, rwVhInv, sumupW
    use shapeFnMod, only           : shapef_z, shapef_r
    use depBdrcMod, only           : depFoldingBoundary
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer                       :: i, j, k, m, ip, jp, mythread
    double precision              :: sfrh(-2:2), sfzh(-2:2)
    double precision, intent(out) :: sumup  (LIs,LJs,ns)

    ! --- [1] Initialize  --- !
    !  -- [1-1] sumup     --  !
    !$omp parallel default(none) &
    !$omp shared(sumup,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             sumup(i,j,k) = 0.d0
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [1-2] sumupW    --  !
    !$omp parallel default(none) &
    !$omp shared(sumupW,LIs,LJs) private(i,j,k,m)
    !$omp do
    do m=1, OMPNumThreads
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                sumupW(i,j,k,m) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( OMPNumThreads.eq.1 ) mythread = 1
    
    ! --- [2] Count Up   --- !
    do k=1, ns
       !$omp parallel default(none) &
       !$omp shared (pxv,np,dr,dz,drInv,dzInv,k,sumupW) &
       !$omp private(mythread,m,ip,jp,i,j,sfrh,sfzh)
       !$ mythread = omp_get_thread_num() + 1
       !$omp do
       do m = 1, np(k)
          jp         =  ceiling( pxv(rp_,m,k)*drInv ) - 1
          ip         =  ceiling( pxv(zp_,m,k)*dzInv ) - 1
          sfrh(-2:2) = shapef_r( jp, pxv(rp_,m,k)-0.5d0*dr, drInv )
          sfzh(-2:2) = shapef_z( ip, pxv(zp_,m,k)-0.5d0*dz, dzInv )
          jp         = jp + 2
          ip         = ip + 2
          do j=-1, 1
             do i=-1, 1
                sumupW(ip+i,jp+j,k,mythread) = sumupW(ip+i,jp+j,k,mythread) + sfrh(j)*sfzh(i)*pxv(wp_,m,k)
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo

    ! --- [3] Reduction sumupW  => sumup --- !
    do m=2, OMPNumThreads
       !$omp parallel default(none) &
       !$omp shared(sumupW,m,LIs,LJs) private(i,j,k)
       !$omp do
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                sumupW(i,j,k,1) = sumupW(i,j,k,1) + sumupW(i,j,k,m)
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    
    ! --- [4] Boundary / Communication   --- !
    !  -- [ 2 -> 1 ] -- !
    call depFoldingBoundary( sumupW(:,:,:,1) )
    
    ! --- [5] Constants & Volume --- !
    !$omp parallel default(none) &
    !$omp shared(sumup,sumupW,RhoPerNsp,rwVhInv,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             sumup(i,j,k) = sumupW(i,j,k,1) * rwVhInv(j) * RhoPerNsp
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine nCountUp_CC

  
  subroutine u1stMoments
    use constants , only         : LIs, LJs, ns, q
    use constants , only         : jr_, jt_, jz_
    use variables , only         : Jcr
    use depBdrcMod, only         : FoldingCopy_CC
    implicit none
    integer                     :: i, j, k
    double precision            :: rhocInv
    double precision, parameter :: rhocMin = 1.d-5

    ! -- [1] u :: Flow Velocity ( u = J / qn ) -- !
    !$omp parallel default(none) &
    !$omp shared(fMoments,q,Jcr,LIs,LJs) private(i,j,k,rhocInv)
    !$omp do
    do k=1, ns
       do j=2, LJs-1
          do i=2, LIs-1
             rhocInv             =  1.d0 / ( q(k) * max( fMoments(i,j,k,nei_), rhocMin ) )
             fMoments(i,j,k,ur_) = 0.5d0*( Jcr(jr_,i,j,k) + Jcr(jr_,i  ,j+1,k) )*rhocInv
             fMoments(i,j,k,ut_) =       ( Jcr(jt_,i,j,k)                      )*rhocInv
             fMoments(i,j,k,uz_) = 0.5d0*( Jcr(jz_,i,j,k) + Jcr(jz_,i+1,j  ,k) )*rhocInv
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    do k=1, ns
       call FoldingCopy_CC( fMoments(1:LIs,1:LJs,k,ur_) )
       call FoldingCopy_CC( fMoments(1:LIs,1:LJs,k,ut_) )
       call FoldingCopy_CC( fMoments(1:LIs,1:LJs,k,uz_) )
    enddo

    return
  end subroutine u1stMoments


  subroutine nCountUp_SP( perCell )
    use constants , only         : LIs, LJs, ns, np, OMPNumThreads
    use constants , only         : drinv, dzinv
    use constants , only         : rp_, zp_, wp_
    use variables , only         : pxv, perCellW, sumupE
    use depBdrcMod, only         : depFoldingBoundary
    !$ use omp_lib
    implicit none
    integer                     :: i, j, k, m, ip, jp, mythread
    integer     , intent(inout) :: perCell(LIs,LJs,ns)
    
    ! --- [1] Initialize  --- !
    !$omp parallel default(none) &
    !$omp shared(perCell,perCellW,sumupE,LIs,LJs) &
    !$omp private(i,j,k,m)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             perCell(i,j,k) = 0
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do m=1, OMPNumThreads
       do k=1, ns
          do j=1, LJs
             do i=1, LIs
                perCellW(i,j,k,m) = 0
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp do
    do k=1, ns
       do j=0, LJs+2
          do i=0, LIs+2
             sumupE(i,j,k) = 0.d0
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( OMPNumThreads.eq.1 ) mythread = 1

    ! --- [2]   Count     --- !
    do k=1, ns
       !$omp parallel default(none) &
       !$omp shared(k,np,pxv,drinv,dzinv,perCellW) &
       !$omp private( mythread,m,ip,jp)
       !$ mythread = omp_get_thread_num() + 1
       !$omp do
       do m = 1, np(k)
          jp     = nint( pxv(rp_,m,k)*drinv ) + 2
          ip     = nint( pxv(zp_,m,k)*dzinv ) + 2
          perCellW(ip,jp,k,mythread) = perCellW(ip,jp,k,mythread) + min( 1, ceiling( pxv(wp_,m,k) ) )
       enddo
       !$omp end do
       !$omp end parallel
    enddo

    ! --- [3]  Reduction  --- !
    do m=1, OMPNumThreads
       !$omp parallel default(none) &
       !$omp shared(sumupE,perCellW,m,LIs,LJs) private(i,j,k)
       !$omp do
       do k=1, ns
          do j=1, LJs
             do i=1, LIs
                sumupE(i,j,k) = sumupE(i,j,k) + dble( perCellW(i,j,k,m) )
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    
    ! --- [4] Boundary Condition --- !
    call depFoldingBoundary( sumupE )
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             perCell(i,j,k) = nint( sumupE(i,j,k) )
          enddo
       enddo
    enddo

    return
  end subroutine nCountUp_SP

  
end module momentsMod
