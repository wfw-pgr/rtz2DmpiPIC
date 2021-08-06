module p2ndMomMod
  implicit none
  integer, parameter     :: ur_  = 2,  ut_  = 3,  uz_  = 4
  integer, parameter     :: pxx_ = 5,  pyy_ = 6,  pzz_ = 7
  integer, parameter     :: pxy_ = 8,  pyz_ = 9,  pxz_ =10
  integer, parameter     :: nPrsQth  = 12
  integer, parameter     :: nMoments = 16
contains

  ! =================================================================== !
  ! === [1] Pressure (2nd) Moments                                  === !
  ! =================================================================== !
  subroutine p2ndMoments( fMoments )
    use constants , only   : LIs, LJs, ns, np, OMPNumThreads
    use constants , only   : dr, dz, dzInv, drInv
    use constants , only   : rp_, zp_, vr_, vt_, vz_, wp_, br_, bt_, bz_
    use variables , only   : pxv, EBf, rwVhInv, RhoPerNsp, sumupP, sumupE
    use shapeFnMod, only   : shapef_z, shapef_r
    use depBdrcMod, only   : depFoldingBoundary
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer               :: i, j, k, m, cmp, pcmp, ip, jp, ipi, jpj, mythread
    double precision      :: sfrh(-2:2), sfzh(-2:2)
    double precision      :: up1, up2, up3, bf1, bf2, bf3, bpnrm, bfnrm
    double precision      :: plbf1, plbf2, plbf3, t1bf1, t1bf2, t1bf3
    double precision      :: n1bf1, n1bf2, n1bf3, n2bf1, n2bf2, n2bf3
    double precision      :: uf_t1, uf_n1, uf_n2, vp_t1, vp_n1, vp_n2
    double precision      :: delv1, delv2, delv3, delt1, deln1, deln2, delvS
    integer, parameter    :: d1_ =ur_, d2_ =ut_, d3_ =uz_
    integer, parameter    :: v1_ =vr_, v2_ =vt_, v3_ =vz_
    integer, parameter    :: p11_=1  , p22_=2  , p33_=3
    integer, parameter    :: p12_=4  , p23_=5  , p13_=6
    integer, parameter    :: pt1_=7  , pn1_=8  , pn2_=9
    integer, parameter    :: qv1_=10 , qv2_=11 , qv3_=12
    double precision, intent(inout) :: fMoments(LIs,LJs,ns,nMoments)

    ! ------------------------------------- !
    ! --- [1] Initialize                --- !
    ! ------------------------------------- !
    !$omp parallel default(none) &
    !$omp shared(sumupP,LIs,LJs) private(i,j,k,m,cmp)
    !$omp do
    do m=1, OMPNumThreads
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                do cmp=1, nPrsQth
                   sumupP(cmp,i,j,k,m) = 0.d0
                enddo
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( OMPNumThreads.eq.1 ) mythread = 1

    ! ------------------------------------- !
    ! --- [2] Count Up                  --- !
    ! ------------------------------------- !
    do k=1, ns
       !$omp parallel default(none) &
       !$omp shared (pxv,EBf,fMoments,np,dr,dz,drInv,dzInv,k,sumupP,LIs,LJs) &
       !$omp private(mythread,m,ip,jp,i,j,ipi,jpj,sfrh,sfzh) &
       !$omp private(up1,up2,up3,bf1,bf2,bf3,bpnrm,bfnrm,plbf1,plbf2,plbf3,t1bf1,t1bf2,t1bf3) &
       !$omp private(n1bf1,n1bf2,n1bf3,n2bf1,n2bf2,n2bf3,uf_t1,uf_n1,uf_n2,vp_t1,vp_n1,vp_n2) &
       !$omp private(delv1,delv2,delv3,delt1,deln1,deln2,delvS)
       !$ mythread = omp_get_thread_num() + 1
       !$omp do
       do m = 1, np(k)
          jp         = ceiling( pxv(rp_,m,k)*drInv ) - 1
          ip         = ceiling( pxv(zp_,m,k)*dzInv ) - 1
          sfrh(-2:2) = shapef_r( jp, pxv(rp_,m,k)-0.5d0*dr, drInv )
          sfzh(-2:2) = shapef_z( ip, pxv(zp_,m,k)-0.5d0*dz, dzInv )
          jp         = jp + 2
          ip         = ip + 2
          up1        = 0.d0
          up2        = 0.d0
          up3        = 0.d0
          bf1        = 0.d0
          bf2        = 0.d0
          bf3        = 0.d0
          do j=-1, 1
             do i=-1, 1
                ipi  = min( max( 2, ip+i ), LIs-1 )
                jpj  = min( max( 2, jp+j ), LJs-1 )
                up1  = up1 + fMoments(ipi,jpj,k,d1_)*sfrh(j)*sfzh(i)
                up2  = up2 + fMoments(ipi,jpj,k,d2_)*sfrh(j)*sfzh(i)
                up3  = up3 + fMoments(ipi,jpj,k,d3_)*sfrh(j)*sfzh(i)
                bf1  = bf1 +  EBf(br_,ipi,jpj)      *sfrh(j)*sfzh(i)
                bf2  = bf2 +  EBf(bt_,ipi,jpj)      *sfrh(j)*sfzh(i)
                bf3  = bf3 +  EBf(bz_,ipi,jpj)      *sfrh(j)*sfzh(i)
             enddo
          enddo
          bpnrm  = sqrt( bf1**2          + bf3**2 )
          bfnrm  = sqrt( bf1**2 + bf2**2 + bf3**2 )
          if ( bfnrm.eq.0.d0 ) then
             bfnrm = 0.d0
          else
             bfnrm = 1.d0 / bfnrm
          endif
          if ( bpnrm.eq.0.d0 ) then
             bpnrm = 0.d0
          else
             bpnrm = 1.d0 / bpnrm
          endif
          ! -- bp  -- !
          plbf1  = bf1 * bpnrm
          plbf2  = 0.d0
          plbf3  = bf3 * bpnrm
          ! -- t   -- !
          t1bf1  = bf1 * bfnrm
          t1bf2  = bf2 * bfnrm
          t1bf3  = bf3 * bfnrm
          ! -- n2  -- !
          n2bf1  = - plbf3
          n2bf2  =   0.d0
          n2bf3  =   plbf1
          ! -- n1  -- !
          n1bf1  = n2bf2*t1bf3 - n2bf3*t1bf2
          n1bf2  = n2bf3*t1bf1 - n2bf1*t1bf3
          n1bf3  = n2bf1*t1bf2 - n2bf2*t1bf1
          ! -- uth -- !
          uf_t1  =          up1*t1bf1 +          up2*t1bf2 +          up3*t1bf3
          uf_n1  =          up1*n1bf1 +          up2*n1bf2 +          up3*n1bf3
          uf_n2  =          up1*n2bf1 +          up2*n2bf2 +          up3*n2bf3
          vp_t1  = pxv(v1_,m,k)*t1bf1 + pxv(v2_,m,k)*t1bf2 + pxv(v3_,m,k)*t1bf3
          vp_n1  = pxv(v1_,m,k)*n1bf1 + pxv(v2_,m,k)*n1bf2 + pxv(v3_,m,k)*n1bf3
          vp_n2  = pxv(v1_,m,k)*n2bf1 + pxv(v2_,m,k)*n2bf2 + pxv(v3_,m,k)*n2bf3
          ! -- pressure -- !
          delv1  = pxv(v1_,m,k)-up1
          delv2  = pxv(v2_,m,k)-up2
          delv3  = pxv(v3_,m,k)-up3
          delt1  = vp_t1-uf_t1
          deln1  = vp_n1-uf_n1
          deln2  = vp_n2-uf_n2
          delvS  = delv1**2 + delv2**2 + delv3**2
          sumupP(p11_,ip,jp,k,mythread) = sumupP(p11_,ip,jp,k,mythread) + pxv(wp_,m,k)*delv1*delv1
          sumupP(p22_,ip,jp,k,mythread) = sumupP(p22_,ip,jp,k,mythread) + pxv(wp_,m,k)*delv2*delv2
          sumupP(p33_,ip,jp,k,mythread) = sumupP(p33_,ip,jp,k,mythread) + pxv(wp_,m,k)*delv3*delv3
          sumupP(p12_,ip,jp,k,mythread) = sumupP(p12_,ip,jp,k,mythread) + pxv(wp_,m,k)*delv1*delv2
          sumupP(p23_,ip,jp,k,mythread) = sumupP(p23_,ip,jp,k,mythread) + pxv(wp_,m,k)*delv2*delv3
          sumupP(p13_,ip,jp,k,mythread) = sumupP(p13_,ip,jp,k,mythread) + pxv(wp_,m,k)*delv3*delv1
          sumupP(pt1_,ip,jp,k,mythread) = sumupP(pt1_,ip,jp,k,mythread) + pxv(wp_,m,k)*delt1*delt1
          sumupP(pn1_,ip,jp,k,mythread) = sumupP(pn1_,ip,jp,k,mythread) + pxv(wp_,m,k)*deln1*deln1
          sumupP(pn2_,ip,jp,k,mythread) = sumupP(pn2_,ip,jp,k,mythread) + pxv(wp_,m,k)*deln2*deln2
          sumupP(qv1_,ip,jp,k,mythread) = sumupP(qv1_,ip,jp,k,mythread) + pxv(wp_,m,k)*delvS*delv1
          sumupP(qv2_,ip,jp,k,mythread) = sumupP(qv2_,ip,jp,k,mythread) + pxv(wp_,m,k)*delvS*delv2
          sumupP(qv3_,ip,jp,k,mythread) = sumupP(qv3_,ip,jp,k,mythread) + pxv(wp_,m,k)*delvS*delv3
       enddo
       !$omp end do
       !$omp end parallel
    enddo

    ! ------------------------------------------------ !
    ! --- [3] Reduction sumupP  => sumupP(:,:,:,1) --- !
    ! ------------------------------------------------ !
    do m=2, OMPNumThreads
       !$omp parallel default(none) &
       !$omp shared(sumupP,m,LIs,LJs) private(i,j,k)
       !$omp do
       do k=1, ns
          do j=1, LJs+1
             do i=1, LIs+1
                sumupP(p11_,i,j,k,1) = sumupP(p11_,i,j,k,1) + sumupP(p11_,i,j,k,m)
                sumupP(p22_,i,j,k,1) = sumupP(p22_,i,j,k,1) + sumupP(p22_,i,j,k,m)
                sumupP(p33_,i,j,k,1) = sumupP(p33_,i,j,k,1) + sumupP(p33_,i,j,k,m)
                sumupP(p12_,i,j,k,1) = sumupP(p12_,i,j,k,1) + sumupP(p12_,i,j,k,m)
                sumupP(p23_,i,j,k,1) = sumupP(p23_,i,j,k,1) + sumupP(p23_,i,j,k,m)
                sumupP(p13_,i,j,k,1) = sumupP(p13_,i,j,k,1) + sumupP(p13_,i,j,k,m)
                sumupP(pt1_,i,j,k,1) = sumupP(pt1_,i,j,k,1) + sumupP(pt1_,i,j,k,m)
                sumupP(pn1_,i,j,k,1) = sumupP(pn1_,i,j,k,1) + sumupP(pn1_,i,j,k,m)
                sumupP(pn2_,i,j,k,1) = sumupP(pn2_,i,j,k,1) + sumupP(pn2_,i,j,k,m)
                sumupP(qv1_,i,j,k,1) = sumupP(qv1_,i,j,k,1) + sumupP(qv1_,i,j,k,m)
                sumupP(qv2_,i,j,k,1) = sumupP(qv2_,i,j,k,1) + sumupP(qv2_,i,j,k,m)
                sumupP(qv3_,i,j,k,1) = sumupP(qv3_,i,j,k,1) + sumupP(qv3_,i,j,k,m)
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo

    ! ------------------------------------------------ !
    ! --- [4] x Const x Volume    &    Boundary    --- !
    ! ------------------------------------------------ !
    do cmp=1, nPrsQth
       pcmp = (cmp-1) + pxx_
       ! -- Constants & Volume -- !
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                sumupE(i,j,k) = sumupP(cmp,i,j,k,1) * rwVhInv(j)*RhoPerNsp
             enddo
          enddo
       enddo
       ! -- Boundary Condition -- !
       call depFoldingBoundary( sumupE )
       do k=1, ns
          do j=1, LJs
             do i=1, LIs
                fMoments(i,j,k,pcmp) = sumupE(i,j,k)
             enddo
          enddo
       enddo
    enddo
    return
  end subroutine p2ndMoments

  
end module p2ndMomMod


  ! subroutine p2ndMoments( fMoments )
  !   use constants , only   : LIs, LJs, ns, np, OMPNumThreads
  !   use constants , only   : dr, dz, dzInv, drInv
  !   use constants , only   : rp_, zp_, vr_, vt_, vz_, wp_
  !   use variables , only   : pxv, rwVhInv, RhoPerNsp, sumupE, sumupP
  !   use shapeFnMod, only   : shapef_z, shapef_r
  !   use depBdrcMod, only   : depFoldingBoundary
  !   !$ use omp_lib
  !   implicit none
  !   include 'mpif.h'
  !   integer               :: i, j, k, m, cmp, pcmp, ip, jp, ipi, jpj, mythread
  !   double precision      :: sfrh(-2:2), sfzh(-2:2)
  !   integer, parameter    :: d1_ =ur_, d2_ =ut_, d3_ =uz_
  !   integer, parameter    :: v1_ =vr_, v2_ =vt_, v3_ =vz_
  !   integer, parameter    :: p11_=1  , p22_=2  , p33_=3
  !   integer, parameter    :: p12_=4  , p23_=5  , p13_=6
  !   double precision, intent(inout) :: fMoments(LIs,LJs,ns,nMoments)

  !   ! --- [1] Initialize  --- !
  !   !$omp parallel default(none) &
  !   !$omp shared(sumupP,LIs,LJs) private(i,j,k,m,cmp)
  !   !$omp do
  !   do m=1, OMPNumThreads
  !      do k=1, ns
  !         do j=0, LJs+2
  !            do i=0, LIs+2
  !               do cmp=1, 6
  !                  sumupP(cmp,i,j,k,m) = 0.d0
  !               enddo
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   !$omp end do
  !   !$omp end parallel
  !   if ( OMPNumThreads.eq.1 ) mythread = 1
    
  !   ! --- [2] Count Up   --- !
  !   do k=1, ns
  !      !$omp parallel default(none) &
  !      !$omp shared (pxv,fMoments,np,dr,dz,drInv,dzInv,k,sumupP,LIs,LJs) &
  !      !$omp private(mythread,m,ip,jp,i,j,ipi,jpj,sfrh,sfzh)
  !      !$ mythread = omp_get_thread_num() + 1
  !      !$omp do
  !      do m = 1, np(k)
  !         jp         = ceiling( pxv(rp_,m,k)*drInv ) - 1
  !         ip         = ceiling( pxv(zp_,m,k)*dzInv ) - 1
  !         sfrh(-2:2) = shapef_r( jp, pxv(rp_,m,k)-0.5d0*dr, drInv )
  !         sfzh(-2:2) = shapef_z( ip, pxv(zp_,m,k)-0.5d0*dz, dzInv )
  !         jp         = jp + 2
  !         ip         = ip + 2
  !         do j=-1, 1
  !            do i=-1, 1
  !               ipi  = min( max( 2, ip+i ), LIs-1 )
  !               jpj  = min( max( 2, jp+j ), LJs-1 )
  !               sumupP(p11_,ipi,jpj,k,mythread) = + sumupP(p11_,ipi,jpj,k,mythread)            &
  !                    &                            + sfrh(j)*sfzh(i)*pxv(wp_,m,k)               &
  !                    &                              * ( pxv(v1_,m,k)-fMoments(ipi,jpj,k,d1_) ) &
  !                    &                              * ( pxv(v1_,m,k)-fMoments(ipi,jpj,k,d1_) )
  !               sumupP(p22_,ipi,jpj,k,mythread) = + sumupP(p22_,ipi,jpj,k,mythread)            &
  !                    &                            + sfrh(j)*sfzh(i)*pxv(wp_,m,k)               &
  !                    &                              * ( pxv(v2_,m,k)-fMoments(ipi,jpj,k,d2_) ) &
  !                    &                              * ( pxv(v2_,m,k)-fMoments(ipi,jpj,k,d2_) )
  !               sumupP(p33_,ipi,jpj,k,mythread) = + sumupP(p33_,ipi,jpj,k,mythread)            &
  !                    &                            + sfrh(j)*sfzh(i)*pxv(wp_,m,k)               &
  !                    &                              * ( pxv(v3_,m,k)-fMoments(ipi,jpj,k,d3_) ) &
  !                    &                              * ( pxv(v3_,m,k)-fMoments(ipi,jpj,k,d3_) )
  !               sumupP(p12_,ipi,jpj,k,mythread) = + sumupP(p12_,ipi,jpj,k,mythread)            &
  !                    &                            + sfrh(j)*sfzh(i)*pxv(wp_,m,k)               &
  !                    &                              * ( pxv(v1_,m,k)-fMoments(ipi,jpj,k,d1_) ) &
  !                    &                              * ( pxv(v2_,m,k)-fMoments(ipi,jpj,k,d2_) )
  !               sumupP(p13_,ipi,jpj,k,mythread) = + sumupP(p13_,ipi,jpj,k,mythread)            &
  !                    &                            + sfrh(j)*sfzh(i)*pxv(wp_,m,k)               &
  !                    &                              * ( pxv(v1_,m,k)-fMoments(ipi,jpj,k,d1_) ) &
  !                    &                              * ( pxv(v3_,m,k)-fMoments(ipi,jpj,k,d3_) )
  !               sumupP(p23_,ipi,jpj,k,mythread) = + sumupP(p23_,ipi,jpj,k,mythread)            &
  !                    &                            + sfrh(j)*sfzh(i)*pxv(wp_,m,k)               &
  !                    &                              * ( pxv(v2_,m,k)-fMoments(ipi,jpj,k,d2_) ) &
  !                    &                              * ( pxv(v3_,m,k)-fMoments(ipi,jpj,k,d3_) )
  !            enddo
  !         enddo
  !      enddo
  !      !$omp end do
  !      !$omp end parallel
  !   enddo
    
  !   ! --- [3] Reduction sumupP  => sumupP(:,:,:,1) --- !
  !   do m=2, OMPNumThreads
  !      !$omp parallel default(none) &
  !      !$omp shared(sumupP,m,LIs,LJs) private(i,j,k)
  !      !$omp do
  !      do k=1, ns
  !         do j=1, LJs+1
  !            do i=1, LIs+1
  !               sumupP(p11_,i,j,k,1) = sumupP(p11_,i,j,k,1) + sumupP(p11_,i,j,k,m)
  !               sumupP(p22_,i,j,k,1) = sumupP(p22_,i,j,k,1) + sumupP(p22_,i,j,k,m)
  !               sumupP(p33_,i,j,k,1) = sumupP(p33_,i,j,k,1) + sumupP(p33_,i,j,k,m)
  !               sumupP(p12_,i,j,k,1) = sumupP(p12_,i,j,k,1) + sumupP(p12_,i,j,k,m)
  !               sumupP(p23_,i,j,k,1) = sumupP(p23_,i,j,k,1) + sumupP(p23_,i,j,k,m)
  !               sumupP(p13_,i,j,k,1) = sumupP(p13_,i,j,k,1) + sumupP(p13_,i,j,k,m)
  !            enddo
  !         enddo
  !      enddo
  !      !$omp end do
  !      !$omp end parallel
  !   enddo

  !   do cmp=1, 6
  !      pcmp = ( cmp-1 ) + pxx_
  !      ! -- Constants & Volume -- !
  !      do k=1, ns
  !         do j=0, LJs+2
  !            do i=0, LIs+2
  !               sumupE(i,j,k) = sumupP(cmp,i,j,k,1) * rwVhInv(j)*RhoPerNsp
  !            enddo
  !         enddo
  !      enddo
  !      ! -- Boundary Condition -- !
  !      call depFoldingBoundary( sumupE )
  !      do k=1, ns
  !         do j=1, LJs
  !            do i=1, LIs
  !               fMoments(i,j,k,pcmp) = sumupE(i,j,k)
  !            enddo
  !         enddo
  !      enddo
  !   enddo
  !   return
  ! end subroutine p2ndMoments
  


          ! sumupP(p11_,ip,jp,k,mythread) = sumupP(p11_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( pxv(v1_,m,k)-up1  )*( pxv(v1_,m,k)-up1 )
          ! sumupP(p22_,ip,jp,k,mythread) = sumupP(p22_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( pxv(v2_,m,k)-up2  )*( pxv(v2_,m,k)-up2 )
          ! sumupP(p33_,ip,jp,k,mythread) = sumupP(p33_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( pxv(v3_,m,k)-up3  )*( pxv(v3_,m,k)-up3 )
          ! sumupP(p12_,ip,jp,k,mythread) = sumupP(p12_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( pxv(v1_,m,k)-up1  )*( pxv(v2_,m,k)-up2 )
          ! sumupP(p23_,ip,jp,k,mythread) = sumupP(p23_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( pxv(v2_,m,k)-up2  )*( pxv(v3_,m,k)-up3 )
          ! sumupP(p13_,ip,jp,k,mythread) = sumupP(p13_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( pxv(v3_,m,k)-up3  )*( pxv(v1_,m,k)-up1 )
          ! sumupP(pt1_,ip,jp,k,mythread) = sumupP(pt1_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( vp_t1-uf_t1       )*( vp_t1-uf_t1      )
          ! sumupP(pn1_,ip,jp,k,mythread) = sumupP(pn1_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( vp_n1-uf_n1       )*( vp_n1-uf_n1      )
          ! sumupP(pn2_,ip,jp,k,mythread) = sumupP(pn2_,ip,jp,k,mythread) &
          !      &                          + pxv(wp_,m,k)*( vp_n2-uf_n2       )*( vp_n2-uf_n2      )
