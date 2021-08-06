module currentMod
contains

  subroutine DensityDecomposition
    use constants , only         : LIs, LJs, ns, OMPNumThreads, Boundary1__jc
    use constants , only         : q, np, dr, dz, drinv, dzinv, drdt, dzdt
    use constants , only         : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_, jr_, jt_, jz_
    use variables , only         : Jcr, JcrW, pxv, RhoPerNsp, rwVfInv, rwVhInv
    use shapeFnMod, only         : shapef_r, shapef_z
    use jcrBdrcMod, only         : jcrFoldingBC, jcrPeriodicBC
    !$ use omp_lib
    implicit none
    include 'mpif.h'
    integer                     :: i, j, k, n, m, ip, jp
    integer                     :: i1, i2, j1, j2, iinc, jinc, mythread
    double precision            :: gamma, hdr, hdz, coefr(ns), coefz(ns), coeft(ns)
    double precision            :: s0(-2:2,2), s1(-2:2,2), ds(-2:2,2)
    double precision            :: jrp(-2:3,-2:3), jzp(-2:3,-2:3), jtp(-2:3,-2:3)
    double precision, parameter :: onethird = 1.d0 / 3.d0

    ! ----------------------------- !
    ! --- [1]  Initialization   --- !
    ! ----------------------------- !
    hdr = 0.5d0 * dr
    hdz = 0.5d0 * dz
    !$omp parallel default(none) &
    !$omp shared(Jcr,JcrW,LIs,LJs) private(i,j,k,n)
    !$omp do
    do k=0, ns
       do j=1, LJs
          do i=1, LIs
             Jcr (jr_,i,j,k) = 0.d0
             Jcr (jt_,i,j,k) = 0.d0
             Jcr (jz_,i,j,k) = 0.d0
          enddo
       enddo
    enddo
    !$omp enddo
    !$omp do
    do n=1, OMPNumThreads
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                JcrW(jr_,i,j,k,n) = 0.d0
                JcrW(jt_,i,j,k,n) = 0.d0
                JcrW(jz_,i,j,k,n) = 0.d0
             enddo
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( OMPNumThreads.eq.1 ) mythread = 1
    
    ! ----------------------------- !
    ! --- [2]  Sumup  Current   --- !
    ! ----------------------------- !
    do k=1, ns

       !$omp parallel default(none) &
       !$omp shared(k,np,pxv,dr,dz,drinv,dzinv,JcrW,hdr,hdz,LIs,LJs) &
       !$omp private(m,i,j,ip,jp,i1,j1,i2,j2,iinc,jinc,gamma,s0,s1,ds,jrp,jtp,jzp,mythread)
       !$ mythread = omp_get_thread_num() + 1
       !$omp do
       do m=1, np(k)
          ! -- (1) Initialize single particle current -- !
          do j=-2, 3
             do i=-2,3
                jrp(i,j) = 0.d0
                jtp(i,j) = 0.d0
                jzp(i,j) = 0.d0
             enddo
          enddo
          ! -- (2) particle Index ( 1:Old, 2:New )    -- !
          i1         = ceiling( pxv(zo_,m,k)*dzinv ) - 1
          j1         = ceiling( pxv(ro_,m,k)*drinv ) - 1
          i2         = ceiling( pxv(zp_,m,k)*dzinv ) - 1
          j2         = ceiling( pxv(rp_,m,k)*drinv ) - 1
          iinc       = i2-i1
          jinc       = j2-j1
          ! -- (3) Flux through cell-surface          -- !
          gamma      = 1.d0 / sqrt( 1.d0 + ( pxv(vr_,m,k)**2 + pxv(vt_,m,k)**2 + pxv(vz_,m,k)**2 ) )
          s0(-2:2,1) = shapef_z( i1, pxv(zo_,m,k)-hdz, dzinv )
          s0(-2:2,2) = shapef_r( j1, pxv(ro_,m,k)-hdr, drinv )
          s1(-2:2,1) = shapef_z( i2, pxv(zp_,m,k)-hdz, dzinv )
          s1(-2:2,2) = shapef_r( j2, pxv(rp_,m,k)-hdr, drinv )
          ds(-2:2,1) = cshift( s1(-2:2,1),-iinc,1 ) - s0(-2:2,1)
          ds(-2:2,2) = cshift( s1(-2:2,2),-jinc,1 ) - s0(-2:2,2)
          ! -- (4) Integration From Edge Region       -- !
          do j=-2,2
             do i=-2,2
                jzp(i+1,j) = jzp(i,j) - ds(i,1)*( s0(j,2) + 0.5d0 * ds(j,2) )*pxv(wp_,m,k)
                jrp(i,j+1) = jrp(i,j) - ds(j,2)*( s0(i,1) + 0.5d0 * ds(i,1) )*pxv(wp_,m,k)
                jtp(i,j  ) = gamma*pxv(vt_,m,k)*( s0(i,1)*s0(j,2) + 0.5d0*ds(i,1)*s0(j,2) &
                     &       + 0.5d0*s0(i,1)*ds(j,2) + onethird*ds(i,1)*ds(j,2) )*pxv(wp_,m,k)
             enddo
          enddo
          ! -- (5)  Sum up  -- !
          do j=-2, 2
             jp = j1 + 2 + j
             do i=-2, 2
                ip = i1 + 2 + i
                JcrW(jr_,ip,jp,k,mythread) = JcrW(jr_,ip,jp,k,mythread) + jrp(i,j)
                JcrW(jt_,ip,jp,k,mythread) = JcrW(jt_,ip,jp,k,mythread) + jtp(i,j)
                JcrW(jz_,ip,jp,k,mythread) = JcrW(jz_,ip,jp,k,mythread) + jzp(i,j)
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    
    ! --------------------------------- !
    ! --- [3]  Reduction in OpenMP  --- !
    ! --------------------------------- !
    do n=2, OMPNumThreads
       !$omp parallel default(none) &
       !$omp shared(JcrW,n,LIs,LJs) private(i,j,k)
       !$omp do
       do k=1, ns
          do j=0, LJs+2
             do i=0, LIs+2
                JcrW(jr_,i,j,k,1) = JcrW(jr_,i,j,k,1) + JcrW(jr_,i,j,k,n)
                JcrW(jt_,i,j,k,1) = JcrW(jt_,i,j,k,1) + JcrW(jt_,i,j,k,n)
                JcrW(jz_,i,j,k,1) = JcrW(jz_,i,j,k,1) + JcrW(jz_,i,j,k,n)
             enddo
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    enddo

    ! --------------------------------- !
    ! --- [4]  Overlay & B.C.       --- !
    ! --------------------------------- !
    if ( trim(Boundary1__jc).eq.'cWall'    ) call jcrFoldingBC ( 'Boundary_0', JcrW(1:3,0:LIs+2,0:LJs+2,1:ns,1) )
    if ( trim(Boundary1__jc).eq.'periodic' ) call jcrPeriodicBC( 'Boundary_0', JcrW(1:3,0:LIs+2,0:LJs+2,1:ns,1) )
    
    ! --------------------------------- !
    ! --- [5]  Ion/Electron Current --- !
    ! --------------------------------- !
    do k=1, ns
       coefr(k) = q(k) * drdt * RhoPerNsp
       coeft(k) = q(k)        * RhoPerNsp
       coefz(k) = q(k) * dzdt * RhoPerNsp
    enddo
    !$omp parallel default(none) &
    !$omp shared(Jcr,JcrW,rwVfInv,rwVhInv,coefr,coeft,coefz,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=1, LJs
          do i=1, LIs
             Jcr(jr_,i,j,k) = JcrW(jr_,i,j,k,1) * rwVfInv(j) * coefr(k)
             Jcr(jt_,i,j,k) = JcrW(jt_,i,j,k,1) * rwVhInv(j) * coeft(k)
             Jcr(jz_,i,j,k) = JcrW(jz_,i,j,k,1) * rwVhInv(j) * coefz(k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! --------------------------------- !
    ! --- [6] Net Current ( Je+Ji ) --- !
    ! --------------------------------- !
    !$omp parallel default(none) &
    !$omp shared(Jcr,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          Jcr(jr_,i,j,0) = Jcr(jr_,i,j,1) + Jcr(jr_,i,j,2) 
          Jcr(jt_,i,j,0) = Jcr(jt_,i,j,1) + Jcr(jt_,i,j,2)
          Jcr(jz_,i,j,0) = Jcr(jz_,i,j,1) + Jcr(jz_,i,j,2)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    return
  end subroutine DensityDecomposition

end module currentMod
