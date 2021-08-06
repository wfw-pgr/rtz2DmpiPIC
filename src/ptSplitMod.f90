module ptSplitMod
contains
  
  subroutine GaussEliminateOMP( A0, b, n, x, ierr )
    implicit none
    integer                         :: i, j, k, m
    integer         , intent(in)    :: n
    integer         , intent(inout) :: ierr
    double precision, intent(in)    :: A0(n,n), b(n)
    double precision, intent(out)   :: x(n)
    double precision                :: Ar, t, Am
    double precision                :: A(n,n), w(n)
    
    !  --- [0]  Preparation  ---  !
    A(:,:) = A0(:,:)
    x(:)   = b(:)
    
    !  --- [1]  Upper Triangular Matrization --- !
    do k=1,n
       
       ! -- [1-1] Select Pivot -- !
       m   = k
       Am  = abs( A(k,k) )
       do j=k+1, n
          if ( abs( A(k,j) ).gt.Am ) then
             Am = abs( A(k,j) )
             m  = j
          endif
       enddo
       if ( Am.eq.0.d0 ) then
          ierr = 1
          return
       endif
       if ( k.ne.m ) then
          w(k:n  ) = A(k:n,k)
          A(k:n,k) = A(k:n,m)
          A(k:n,m) = w(k:n  )
          t        = x(k)
          x(k)     = x(m)
          x(m)     = t
       endif
          
       ! -- [1-2] Diagonal Component -- !
       Ar      = 1.d0 / A(k,k)
       A(k,k)  = 1.d0

       ! -- [1-3] Non-Diagnal Component -- !
       if ( k.eq.n ) then  !     last column :: dividing only
          x(k) = Ar * x(k)
       else                ! not last column :: divide and subtracting under k
          A(k+1:n,k)    = Ar * A(k+1:n,k)
          x(k)          = Ar * x(k)
          do j=k+1,n
             A(k+1:n,j) = A(k+1:n,j) - A(k,j) * A(k+1:n,k)
             x(j)       = x(j)       - A(k,j) * x(k)
             A(k,j)     = 0.d0
          enddo
       endif

    enddo

    !  [2] Substitute From Back
    do k=n-1,1,-1
       do i=n,k+1,-1
          x(k) = x(k) - A(i,k)*x(i)
       enddo
    enddo

    return
  end subroutine GaussEliminateOMP

  
  
  subroutine getWeightOMP( lambda, rhs, weight, Nnode, Nnew, ierr )

    use constants, only              : OMPNumThreads
    !$ use omp_lib
    implicit none
    integer         , intent(in)    :: Nnode, Nnew
    integer         , intent(inout) :: ierr
    double precision, intent(in)    :: lambda(Nnode,Nnew), rhs(Nnode)
    double precision, intent(out)   :: weight(Nnew)
    integer                         :: i, j, m, mythread
    double precision                :: gj(Nnode)
    double precision                :: Lmat(Nnode,Nnode), LmatP(Nnode,Nnode,OMPNumThreads)
    
    ! --- [1] Calculate L-Matrix --- !    
    gj(:)        = 0.d0
    Lmat(:,:)    = 0.d0
    LmatP(:,:,:) = 0.d0
    !$omp parallel default(none) &
    !$omp shared(Nnew,Nnode,LmatP,lambda) private(m,i,j,mythread)
    !$ mythread = omp_get_thread_num() + 1
    !$omp do
    do m=1, Nnew
       do j=1, Nnode
          do i=1, Nnode
             LmatP(i,j,mythread) = LmatP(i,j,mythread) + lambda(i,m) * lambda(j,m)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    do mythread=1, OMPNumThreads
       do j=1, Nnode
          do i=1, Nnode
             Lmat(i,j) = Lmat(i,j) + LmatP(i,j,mythread)
          enddo
       enddo
    enddo
    
    !  --- [2] Get gj by Inverting L-matrix  ---  !
    call GaussEliminateOMP( Lmat, rhs, Nnode, gj, ierr )
    if ( ierr.eq.1 ) return

    !  --- [3] Calculation of each weight wm ---  !
    do m=1, Nnew
       weight(m) = 0.d0
       do j=1, Nnode
          weight(m) = weight(m) + gj(j) * lambda(j,m)
       enddo
    enddo

    return
  end subroutine getWeightOMP

  

  subroutine getLambdaOMP( lamb,pxvh,dr,dz,Nnew,i0,j0 )
    use constants, only : rp_, zp_, vr_, vt_, vz_, wp_
    use shapeFnMod, only : shapef_r, shapef_z
    implicit none
    integer                        :: i, j, ip, jp, m
    integer         , intent(in)   :: Nnew, i0, j0
    double precision, intent(in)   :: dr, dz, pxvh(8,Nnew)
    double precision, intent(out)  :: lamb(-2:2,-2:2,Nnew,4)
    double precision               :: drinv, dzinv
    double precision               :: sfrf(-2:2), sfrh(-2:2), sfzf(-2:2), sfzh(-2:2)
    
    drinv          = 1.d0 / dr
    dzinv          = 1.d0 / dz
    !$omp parallel default(none) &
    !$omp shared(Nnew,pxvh,i0,j0,dr,dz,drinv,dzinv,lamb) private(m,i,j,ip,jp,sfrf,sfrh,sfzf,sfzh)
    !$omp do
    do m=1, Nnew

       jp          = ceiling( pxvh(rp_,m)*drinv ) - 1
       ip          = ceiling( pxvh(zp_,m)*dzinv ) - 1
       sfrf(-2:2)  =         shapef_r( jp, pxvh(rp_,m)         , drinv )
       sfrh(-2:2)  = cshift( shapef_r( jp, pxvh(rp_,m)-0.5d0*dr, drinv ), j0-jp )
       sfzf(-2:2)  =         shapef_z( ip, pxvh(zp_,m)         , dzinv )
       sfzh(-2:2)  = cshift( shapef_z( ip, pxvh(zp_,m)-0.5d0*dz, dzinv ), i0-ip )
       do j=-2,2
          do i=-2, 2
             lamb( i,j,m,1 ) = sfrh(j) * sfzh(i)
             lamb( i,j,m,2 ) = sfrf(j) * sfzh(i)
             lamb( i,j,m,3 ) = sfrh(j) * sfzh(i)
             lamb( i,j,m,4 ) = sfrh(j) * sfzf(i)
          enddo
       enddo

    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine getLambdaOMP

  
  
  subroutine getMomentOMP( mmt, pxvh, dr, dz, Np, i0, j0, U0 )
    use constants, only : rp_, zp_, vr_, vt_, vz_, wp_
    use constants    , only : OMPNumThreads
    use shapeFnMod, only : shapef_r, shapef_z
    !$ use omp_lib
    implicit none
    integer                        :: i, j, ip, jp, m, mythread
    integer         , intent(in)   :: Np, i0, j0
    double precision, intent(in)   :: dr, dz, pxvh(8,Np)
    double precision, intent(out)  :: mmt(-2:2,-2:2,4)
    double precision, intent(out)  :: U0
    double precision               :: drinv, dzinv
    double precision               :: sfrf(-2:2), sfrh(-2:2), sfzf(-2:2), sfzh(-2:2)
    double precision               :: mmtP(-2:2,-2:2,4,OMPNumThreads)

    drinv          = 1.d0 / dr
    dzinv          = 1.d0 / dz
    mmt(:,:,:)     = 0.d0
    mmtP(:,:,:,:)  = 0.d0
    U0             = 0.d0
    !$omp parallel default(none) &
    !$omp shared(Np,mmtP,pxvh,dr,dz,drinv,dzinv,i0,j0,U0) &
    !$omp private(m,i,j,ip,jp,sfrf,sfrh,sfzf,sfzh,mythread)
    !$ mythread = omp_get_thread_num() + 1
    !$omp do reduction(+:U0)
    do m=1, Np

       jp          = ceiling( pxvh(rp_,m)*drinv ) - 1
       ip          = ceiling( pxvh(zp_,m)*dzinv ) - 1
       sfrf(-2:2)  =         shapef_r( jp, pxvh(rp_,m)         , drinv )
       sfrh(-2:2)  = cshift( shapef_r( jp, pxvh(rp_,m)-0.5d0*dr, drinv ), j0-jp )
       sfzf(-2:2)  =         shapef_z( ip, pxvh(zp_,m)         , dzinv )
       sfzh(-2:2)  = cshift( shapef_z( ip, pxvh(zp_,m)-0.5d0*dz, dzinv ), i0-ip )
       do j=-2, 2
          do i=-2, 2
             mmtP(i,j,1,mythread) = mmtP(i,j,1,mythread) + pxvh(wp_,m)          *sfrh(j)*sfzh(i)
             mmtP(i,j,2,mythread) = mmtP(i,j,2,mythread) + pxvh(wp_,m)*pxvh(vr_,m)*sfrf(j)*sfzh(i)
             mmtP(i,j,3,mythread) = mmtP(i,j,3,mythread) + pxvh(wp_,m)*pxvh(vt_,m)*sfrh(j)*sfzh(i)
             mmtP(i,j,4,mythread) = mmtP(i,j,4,mythread) + pxvh(wp_,m)*pxvh(vz_,m)*sfrh(j)*sfzf(i)
          enddo
       enddo
       U0          = U0 + pxvh(wp_,m) * sqrt( 1.d0 + pxvh(vr_,m)**2 + pxvh(vt_,m)**2 + pxvh(vz_,m)**2 )
       
    enddo
    !$omp end do
    !$omp end parallel
    do mythread=1, OMPNumThreads
       do j=-2,2
          do i=-2,2
             mmt(i,j,1) = mmt(i,j,1) + mmtP(i,j,1,mythread)
             mmt(i,j,2) = mmt(i,j,2) + mmtP(i,j,2,mythread)
             mmt(i,j,3) = mmt(i,j,3) + mmtP(i,j,3,mythread)
             mmt(i,j,4) = mmt(i,j,4) + mmtP(i,j,4,mythread)
          enddo
       enddo
    enddo
    
    return
  end subroutine getMomentOMP


  
  subroutine solveBetaOMP( HBx, U0, beta, Np, ierr )

    implicit none
    integer         , intent(in)    :: Np
    integer         , intent(inout) :: ierr
    double precision, intent(in)    :: HBx(Np,4,2)
    double precision, intent(inout) :: U0, beta
    integer                         :: m, k, try
    double precision                :: val, dif, res, winv
    double precision                :: coefb(3,Np), cbeta1, cbeta2, cbeta3
    integer         , parameter     :: itermax = 100
    double precision, parameter     :: cvg1    = 1.d-14
    double precision, parameter     :: cvg2    = 1.d-7

    !  ---  [1]  Calculation c1, c2, c3 of U( beta )  ---  !
    try = 1
    cbeta1 = 0.d0
    cbeta2 = 0.d0
    cbeta3 = 0.d0
    !$omp parallel default(none) &
    !$omp shared(Np,coefb,HBx,cbeta1,cbeta2,cbeta3) private(m,winv)
    !$omp do reduction(+:cbeta1,cbeta2,cbeta3)
    do m=1, Np
       coefb(1,m) = ( HBx(m,2,2)**2 + HBx(m,3,2)**2 + HBx(m,4,2)**2 )
       coefb(2,m) = 2.d0 * HBx(m,1,1) * ( HBx(m,2,1)*HBx(m,2,2) &
            &       + HBx(m,3,1)*HBx(m,3,2) + HBx(m,4,1)*HBx(m,4,2) )
       coefb(3,m) = ( HBx(m,2,1)**2 + HBx(m,3,1)**2 + HBx(m,4,1)**2 + 1.d0 ) * HBx(m,1,1)**2

       winv       = 1.d0   / HBx(m,1,1)
       cbeta1     = cbeta1 + coefb(1,m) * winv
       cbeta2     = cbeta2 + coefb(2,m) * winv
       cbeta3     = cbeta3 + coefb(3,m) * winv - winv
       
    enddo
    !$omp end do
    !$omp end parallel
    cbeta3 = cbeta3 - U0

    !   -- [1-2] Initial Guess of Beta  ---   !
    beta     = cbeta2**2 - 4.d0 * cbeta1 * cbeta3
       
    !  ---  [2] Newton Raphson Loop  ---  !
100 continue
    res = beta
    do k=1, itermax

       !  -- [2-1]  Caluculation U, dU/db ---  !
       val = 0.d0
       dif = 0.d0
       !$omp parallel default(none) &
       !$omp shared(Np,val,dif,coefb,beta) private(m)
       !$omp do reduction(+:val,dif)
       do m=1, Np
          val = val + sqrt( coefb(1,m)*beta**2 + coefb(2,m)*beta + coefb(3,m) )
          dif = dif + ( coefb(1,m) * beta + 0.5d0*coefb(2,m) ) &
               &        / sqrt( coefb(1,m)*beta**2 + coefb(2,m)*beta + coefb(3,m) )
       enddo
       !$omp end do
       !$omp end parallel
       val = val - U0
       
       !  -- [2-2]  Step   ---  !
       beta = beta - val / dif
       
       !  -- [2-3]  Beta Check   ---  !
       if( beta.ne.beta ) then
          ierr = 1
          return
       endif

       !  -- [2-4]  Judge  ---  !
       res  = abs( res - beta )
       if ( ( abs(beta).lt.cvg1 ).or.( res.lt.cvg2 ) ) return
       res  = beta
       
    enddo

    if ( try.eq.1 ) then
       try  = try + 1
       beta = ( - cbeta2 - sqrt( cbeta2**2 - 4.d0 * cbeta1 * cbeta3 ) ) / ( 2.d0 * cbeta1 )
       goto 100
    endif

    return
  end subroutine solveBetaOMP

  

  subroutine AssousSplit( pxvh, goban, dr, dz, Nnew, Nold, nwpt, nptMax, gobanRef, pxvOuter, ierr )
    use constants   , only           : rp_, zp_, vr_, vt_, vz_, wp_
    use mt19937Mod  , only           : grnd, multi_grnd
    !$ use omp_lib
    implicit none
    integer         , intent(in)    :: Nnew, Nold, gobanRef, nptMax
    integer         , intent(inout) :: nwpt, ierr, goban(nptMax-nwpt+1)
    double precision, intent(in)    :: dr, dz
    double precision, intent(inout) :: pxvh(8,Nold), pxvOuter(8,nptMax-nwpt+1)
    integer         , parameter     :: Nsub    = 6
    integer         , parameter     :: iterlim = 10
    double precision, parameter     :: vMax    = 0.5d0
    double precision, parameter     :: cLimit  = 0.99d0
    integer                         :: m, k, i, ip, jp, i0, j0, iter, idx, mythread
    double precision                :: drinv, dzinv, subrinv, subzinv, beta, coef
    double precision                :: mu, sd, gamma, tst, pdfmax, U0, U1, winv, wmin
    double precision                :: pxvhold(8,Nold), pxvhNew(8,Nnew)
    double precision                :: pdf(0:Nsub+1,0:Nsub+1), wcdf(0:Nold,2), psmp(8,Nnew)
    double precision                :: lamb(-2:2,-2:2,Nnew,4), HBx(Nnew,4,2), mmt(-2:2,-2:2,4), mmt2(-2:2,-2:2,4)
    
    ! --- [1] preparation  ---  !
    !  -- [1-1] parameters  --  !
    iter    = iterlim
    !$ mythread = omp_get_thread_num() + 1
    drinv   = 1.d0 / dr
    dzinv   = 1.d0 / dz
    subrinv = dble( Nsub ) * drinv
    subzinv = dble( Nsub ) * dzinv
    j0      = nint( pxvh(rp_,1)*drinv )
    i0      = nint( pxvh(zp_,1)*dzinv )
    !$omp parallel default(none) &
    !$omp shared(pxvhold,pxvh,Nold) private(k,m)
    !$omp do
    do m=1, Nold
       do k=1, 8
          pxvhold(k,m) = pxvh(k,m)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! -------------------------------------------- !
    ! --- [2] Begining of Assous-Welch Method  --- !
    ! -------------------------------------------- !
    !  -- [2-1] Sub-grid Particle Distribution  -- !
    pdf(:,:)  = 0.d0
    !$omp parallel default(none) &
    !$omp shared(pxvh,pdf,Nold,j0,i0,dr,dz,subrinv,subzinv) private(m,ip,jp)
    !$omp do
    do m=1, Nold
       jp         = ceiling( ( ( pxvh(rp_,m)+0.5d0*dr ) - dr*dble(j0) ) * subrinv )
       ip         = ceiling( ( ( pxvh(zp_,m)+0.5d0*dz ) - dz*dble(i0) ) * subzinv )
       pdf(ip,jp) = pdf(ip,jp) + 1.d0 / dble( Nold )
    enddo
    !$omp end do
    !$omp end parallel
    pdf(0,:)      = 0.d0
    pdf(:,0)      = 0.d0
    pdf(Nsub+1,:) = 0.d0
    pdf(:,Nsub+1) = 0.d0
    !  -- pdfMax --  !
    pdfmax = 0.d0
    do jp=1, Nsub
       do ip=1, Nsub
          pdfmax  = max( pdfmax, pdf(ip,jp) )
       enddo
    enddo

    !  -- [2-2] Get 0th, 1st, 2nd - Moment & Energy --  !
    call getMomentOMP( mmt, pxvh, dr, dz, Nold, i0, j0, U0 )

    ! ---------------------------- !
    ! === *** Return Point *** === !
    ! ---------------------------- !
300 continue    
    iter = iter - 1
    ierr = 0
    if ( iter .lt. 0 ) then
       !$omp parallel default(none) &
       !$omp shared(pxvh,pxvhold,Nold) private(m,k)
       !$omp do
       do m=1, Nold
          do k=1, 8
             pxvh(k,m) = pxvhold(k,m)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       ierr = 1
       return
    endif ! --- No Match --- !

    
    !  -- [2-3] Determine xm by sub-grid method  -- !
    !!!$omp parallel default(none) &
    !!!$omp shared(pxvhNew,dz,dr,i0,j0,pdfmax,pdf,Nnew,subrinv,subzinv) &
    !!!$omp private(mythread,jp,ip,coef,m)
    !!!$ mythread = omp_get_thread_num() + 1
    !!!$omp do
    do m=1, Nnew
310    continue
       pxvhNew(1,m) = dr * ( dble(j0) + multi_grnd( mythread ) - 0.5d0 )
       pxvhNew(2,m) = dz * ( dble(i0) + multi_grnd( mythread ) - 0.5d0 )
       jp           = ceiling( ( ( pxvhNew(1,m)+0.5d0*dr ) - dr*dble(j0) ) * subrinv )
       ip           = ceiling( ( ( pxvhNew(2,m)+0.5d0*dz ) - dz*dble(i0) ) * subzinv )       
       coef         = pdfmax * multi_grnd( mythread )
       if ( coef.ge.pdf(ip,jp) ) goto 310
    enddo
    !!!$omp end do
    !!!$omp end parallel

    !  -- [2-4] Get Lambda Matrix for m-particles Shape Functions -- !
    call getLambdaOMP( lamb, pxvhNew(:,1:Nnew), dr, dz, Nnew, i0, j0 )
    
    !  -- [2-5] Get Weight of charge / momentum conservation :: gj, hj --  !
    call getWeightOMP( lamb(-2:1,-2:1,:,1), mmt(-2:1,-2:1,1), HBx(:,1,1),16, Nnew, ierr )
    call getWeightOMP( lamb(-2:1,-1:1,:,2), mmt(-2:1,-1:1,2), HBx(:,2,1),12, Nnew, ierr )
    call getWeightOMP( lamb(-2:1,-2:1,:,3), mmt(-2:1,-2:1,3), HBx(:,3,1),16, Nnew, ierr )
    call getWeightOMP( lamb(-1:1,-2:1,:,4), mmt(-1:1,-2:1,4), HBx(:,4,1),12, Nnew, ierr )
    if ( ierr.eq.1 ) goto 300
    
    !  -- [2-6] Check Whether Moderate Weight --  !
    mu    = 0.d0
    sd    = 0.d0
    !$omp parallel default(none) &
    !$omp shared(HBx,Nnew,mu,sd) private(m)
    !$omp do reduction(+:mu,sd)
    do m=1, Nnew
       mu = mu + HBx(m,1,1)
       sd = sd + HBx(m,1,1)**2
    enddo
    !$omp end do
    !$omp end parallel
    coef  = 1.d0 / dble( Nnew )
    mu    =   mu * coef
    sd    = 1.d0 / sqrt( sd * coef - mu**2 )

    wmin  = 1.d6
    U1    = 0.d0
    tst   = 0.d0
    !$omp parallel default(none) &
    !$omp shared(Nnew,HBx,wmin,tst,U1,mu,sd) private(winv,m)
    !$omp do reduction(+:tst,U1) reduction(min:wmin)
    do m=1, Nnew
       winv       = 1.d0 / HBx(m,1,1)
       HBx(m,2,1) = HBx(m,2,1) * winv
       HBx(m,3,1) = HBx(m,3,1) * winv
       HBx(m,4,1) = HBx(m,4,1) * winv
       wmin       = min( HBx(m,1,1), wmin )
       tst        = tst + max( 0.d0, abs( ( HBx(m,1,1) - mu ) * sd ) - 3.d0 )
       U1         = U1  + HBx(m,1,1) * sqrt( 1.d0 + HBx(m,2,1)**2 + HBx(m,3,1)**2 + HBx(m,4,1)**2 )
    enddo
    !$omp end do
    !$omp end parallel
    if ( ( tst.gt.0.d0 ).or.( wmin.lt.0.d0 ).or.( U1.gt.U0 ) ) goto 300

    !  -- [2-7] Reconstruction of VDF  --  !
320 continue
    !    - [2-7-1] Probability of choose m-th particle as Refference -- !
    wcdf(:,:) = 0.d0
    do m=1, Nold
       wcdf(m,1) = pxvh(wp_,m)
       wcdf(m,2) = wcdf(m-1,2) + pxvh(wp_,m)
    enddo
    winv      = 1.d0 / wcdf(Nold,2)
    wcdf(:,2) = wcdf(:,2) * winv
    mythread  = 1
    !    - [2-7-2] Choose Reference particle randomly -- !
    do m=1, Nnew
330    continue
       coef        = multi_grnd( mythread )
       psmp(1:2,m) = pxvh(rp_:2,m)
       psmp(  7,m) = HBx(m,1,1)
       do k=1, Nold
          if ( coef.le.wcdf(k,2) ) then
             ! psmp(3,m) = pxvh(vr_,k) - HBx(m,2,1)
             ! psmp(4,m) = pxvh(vt_,k) - HBx(m,3,1)
             ! psmp(5,m) = pxvh(vz_,k) - HBx(m,4,1)
             psmp(3,m) = pxvh(vr_,k)
             psmp(4,m) = pxvh(vt_,k)
             psmp(5,m) = pxvh(vz_,k)
             psmp(8,m) = 1.d0
             wcdf(k,1) = 0.d0
             wcdf(k,2) = 0.d0
             exit
          endif
          if ( k.eq.Nold ) then
             coef      = 0.d0
             wcdf(0,2) = 0.d0
             do i=1, Nold
                wcdf(i,1) = pxvh(wp_,m)
                wcdf(i,2) = wcdf(i-1,2) + wcdf(i,1)
             enddo
             coef = 1.d0 / wcdf(Nold,2)
             do i=1, Nold
                wcdf(i,2) = wcdf(i,2) * coef
             enddo
             goto 330
          endif
       enddo
    enddo

    !  -- [2-8] Bloadening of Distribution :: Random Velocity --  !
    !   - Coefficients For Random Velocity --  !
    call getMomentOMP( mmt2, psmp, dr, dz, Nnew, i0, j0, U1 )
    call getWeightOMP( lamb(-2:1,-1:1,:,2), - mmt2(-2:1,-1:1,2), HBx(:,2,2),12, Nnew, ierr )
    call getWeightOMP( lamb(-2:1,-2:1,:,3), - mmt2(-2:1,-2:1,3), HBx(:,3,2),16, Nnew, ierr )
    call getWeightOMP( lamb(-1:1,-2:1,:,4), - mmt2(-1:1,-2:1,4), HBx(:,4,2),12, Nnew, ierr )
    if ( ierr.eq.1 ) goto 300
    !   - Random Velocity -  !
    !$omp parallel default(none) &
    !$omp shared(HBx,psmp,Nnew) private(m)
    !$omp do
    do m=1, Nnew
       HBx(m,2,2) = HBx(m,2,2) + psmp(3,m)
       HBx(m,3,2) = HBx(m,3,2) + psmp(4,m)
       HBx(m,4,2) = HBx(m,4,2) + psmp(5,m)
    enddo
    !$omp end do
    !$omp end parallel

    !  -- [2-9] Energy Adjustment --  !
    beta = 0.1d0
    call solveBetaOMP( HBx, U0, beta, Nnew, ierr )
    if ( ierr.eq.1 ) goto 300    
    !   - Update Particle Info. -  !
    !$omp parallel default(none) &
    !$omp shared(beta,HBx,Nnew,pxvhNew) private(m,coef)
    !$omp do
    do m=1, Nnew
       coef         = beta / HBx(m,1,1)
       pxvhNew(3,m) = HBx(m,2,1) + coef * HBx(m,2,2)
       pxvhNew(4,m) = HBx(m,3,1) + coef * HBx(m,3,2)
       pxvhNew(5,m) = HBx(m,4,1) + coef * HBx(m,4,2)
       pxvhNew(8,m) = HBx(m,1,1)
    enddo
    !$omp end do
    !$omp end parallel

    !  -- [2-10] c-Limit Check -- !
    tst = 0.d0
    !$omp parallel default(none) &
    !$omp shared(pxvhNew,Nnew,tst) private(gamma,m)
    !$omp do reduction(+:tst)
    do m=1, Nnew
       gamma  = pxvhNew(3,m)**2 + pxvhNew(4,m)**2 + pxvhNew(5,m)**2
       tst    = tst + max( 0.d0, gamma-vMax )
    enddo
    !$omp end do
    !$omp end parallel
    if ( tst.gt.1.d-8 ) then
       !$omp parallel default(none) &
       !$omp shared(pxvh,pxvhold,Nold) private(i,k)
       !$omp do
       do i=1, Nold
          do k=1, 8
             pxvh(k,i) = pxvhold(k,i)
          enddo
       enddo
       !$omp end do
       !$omp end parallel
       goto 300
    endif

    !  --- [3] Substitute back Particle Info.  ---  !
    !$omp parallel default(none) &
    !$omp shared(pxvh,pxvhNew,pxvOuter,goban,gobanRef,Nold,Nnew) private(k,m,idx)
    !$omp do
    do m=1, Nold
       do k=1, 8
          pxvh(k,m) = pxvhNew(k,m)
       enddo
    enddo
    !$omp end do
    !$omp barrier
    !$omp do
    do m=Nold+1, Nnew
       idx = m - ( Nold )
       do k=1, 8
          pxvOuter(k,idx) = pxvhNew(k,m)
       enddo
       goban(idx) = gobanRef
    enddo
    !$omp end do
    !$omp end parallel
    nwpt = nwpt + ( Nnew - Nold )
    
    return
  end subroutine AssousSplit

end module ptSplitMod
