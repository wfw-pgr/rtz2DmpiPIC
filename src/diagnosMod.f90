module diagnosMod
  implicit none
  integer :: LatestAMomUpdate = -1
contains

  subroutine AngularMomentum
    use constants , only : LIs, LJs, ns, np, rm, drInv, dzInv, OMPNumThreads
    use constants , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use constants , only : Boundary1__jc
    use variables , only : pxv, kstep, w_p, sumupW
    use depBdrcMod, only : depFoldingBoundary
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer            :: i, j, ip, jp, k, m, mythread
    double precision   :: sumup_s(LIs,LJs,ns)

    ! --- [1] Initialize  --- !
    !  -- [1-1] sumup_s     --  !
    !$omp parallel default(none) &
    !$omp shared(sumup_s,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             sumup_s(i,j,k) = 0.d0
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

    ! --- [2] Angular Momentum in Cell --- !
    do k=1, ns
       !$omp parallel default(none) &
       !$omp shared (pxv,np,drInv,dzInv,k,sumupW,LIs,LJs) &
       !$omp private(mythread,m,ip,jp)
       !$ mythread = omp_get_thread_num() + 1
       !$omp do
       do m = 1, np(k)
          jp = ceiling( pxv(rp_,m,k)*drInv ) + 1
          ip = ceiling( pxv(zp_,m,k)*dzInv ) + 1
          sumupW(ip,jp,k,mythread) = sumupW(ip,jp,k,mythread) + pxv(wp_,m,k)*pxv(rp_,m,k)*pxv(vt_,m,k)
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
    
    ! --- [4] MPI Communication --- !
    ! call MPI_ALLREDUCE( sumup_s, w_p, LIs*LJs*ns, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    if ( trim(Boundary1__jc).eq.'cWall'    ) call depFoldingBoundary( sumupW(0:LIs+2,0:LJs+2,1:ns,1) )
    
    ! --- [5] Copy sumupW --- !
    !$omp parallel default(none) &
    !$omp shared(LIs,LJs,sumupW,w_p,rm) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             w_p(i,j,k) = rm(k) * sumupW(i,j,k,1)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! --- [6] LatestUpdate  --- !
    LatestAMomUpdate = kstep
    
    return
  end subroutine AngularMomentum
  
end module diagnosMod
