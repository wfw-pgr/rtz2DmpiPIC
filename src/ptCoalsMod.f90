module ptCoalsMod
contains
  
  subroutine GaussEliminate( A0, b, n, x, ierr )

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
    do k   = 1, n
       
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
  end subroutine GaussEliminate

  
  
  subroutine getWeight( lambda, rhs, weight, Nnode, Nnew, ierr )

    implicit none
    integer         , intent(in)    :: Nnode, Nnew
    integer         , intent(inout) :: ierr
    double precision, intent(in)    :: lambda(Nnode,Nnew), rhs(Nnode)
    double precision, intent(out)   :: weight(Nnew)
    integer                         :: i, j, m
    double precision                :: gj(Nnode), Lmat(Nnode,Nnode)
    
    ! --- [1] Calculate L-Matrix --- !    
    gj(:)     = 0.d0
    Lmat(:,:) = 0.d0
    do m=1, Nnew
       do j=1, Nnode
          do i=1, Nnode
             Lmat(i,j) = Lmat(i,j) + lambda(i,m) * lambda(j,m)
          enddo
       enddo
    enddo
    
    !  --- [2] Get gj by Inverting L-matrix  ---  !
    call GaussEliminate( Lmat, rhs, Nnode, gj, ierr )
    if ( ierr.eq.1 ) return

    !  --- [3] Calculation of each weight wm ---  !
    do m=1, Nnew
       weight(m)    = 0.d0
       do j=1, Nnode
          weight(m) = weight(m) + gj(j) * lambda(j,m)
       enddo
    enddo

    return
  end subroutine getWeight

  

  subroutine getLambda( lamb,pxvh,dr,dz,Nnew,i0,j0 )
    use constants, only : rp_, zp_
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
    return
  end subroutine getLambda

  
  
  subroutine getMoment( mmt, pxvh, dr, dz, Np, i0, j0, U0 )
    use constants, only : rp_, zp_, vr_, vt_, vz_, wp_
    use shapeFnMod, only : shapef_r, shapef_z
    implicit none
    integer                        :: i, j, ip, jp, m
    integer         , intent(in)   :: Np, i0, j0
    double precision, intent(in)   :: dr, dz, pxvh(8,Np)
    double precision, intent(out)  :: mmt(-2:2,-2:2,4)
    double precision, intent(out)  :: U0
    double precision               :: drinv, dzinv
    double precision               :: sfrf(-2:2), sfrh(-2:2), sfzf(-2:2), sfzh(-2:2)

    drinv          = 1.d0 / dr
    dzinv          = 1.d0 / dz
    mmt(:,:,:)     = 0.d0
    U0             = 0.d0
    do m=1, Np

       jp          = ceiling( pxvh(rp_,m)*drinv ) - 1
       ip          = ceiling( pxvh(zp_,m)*dzinv ) - 1
       sfrf(-2:2)  =         shapef_r( jp, pxvh(rp_,m)         , drinv )
       sfrh(-2:2)  = cshift( shapef_r( jp, pxvh(rp_,m)-0.5d0*dr, drinv ), j0-jp )
       sfzf(-2:2)  =         shapef_z( ip, pxvh(zp_,m)         , dzinv )
       sfzh(-2:2)  = cshift( shapef_z( ip, pxvh(zp_,m)-0.5d0*dz, dzinv ), i0-ip )
       do j=-2, 2
          do i=-2, 2
             mmt(i,j,1) = mmt(i,j,1) + pxvh(wp_,m)          *sfrh(j)*sfzh(i)
             mmt(i,j,2) = mmt(i,j,2) + pxvh(wp_,m)*pxvh(vr_,m)*sfrf(j)*sfzh(i)
             mmt(i,j,3) = mmt(i,j,3) + pxvh(wp_,m)*pxvh(vt_,m)*sfrh(j)*sfzh(i)
             mmt(i,j,4) = mmt(i,j,4) + pxvh(wp_,m)*pxvh(vz_,m)*sfrh(j)*sfzf(i)
          enddo
       enddo
       U0          = U0 + pxvh(wp_,m) * sqrt( 1.d0 + pxvh(vr_,m)**2 + pxvh(vt_,m)**2 + pxvh(vz_,m)**2 )
       
    enddo

    return
  end subroutine getMoment


  
  subroutine solv_beta( HBx, U0, beta, Np, ierr )

    implicit none
    integer         , intent(in)    :: Np
    integer         , intent(inout) :: ierr
    double precision, intent(in)    :: HBx(Np,4,2)
    double precision, intent(inout) :: U0, beta
    integer                         :: m, k, try
    double precision                :: val, dif, res, winv
    double precision                :: coefb(3,Np), cbeta(3)
    integer         , parameter     :: itermax = 100
    double precision, parameter     :: cvg1    = 1.d-14
    double precision, parameter     :: cvg2    = 1.d-7

    try = 1
    !  ---  [1]  Calculation c1, c2, c3 of U( beta )  ---  !
    cbeta(:) = 0.d0
    do m=1, Np
       coefb(1,m) = ( HBx(m,2,2)**2 + HBx(m,3,2)**2 + HBx(m,4,2)**2 )
       coefb(2,m) = 2.d0 * HBx(m,1,1) * ( HBx(m,2,1)*HBx(m,2,2) &
            &       + HBx(m,3,1)*HBx(m,3,2) + HBx(m,4,1)*HBx(m,4,2) )
       coefb(3,m) = ( HBx(m,2,1)**2 + HBx(m,3,1)**2 + HBx(m,4,1)**2 + 1.d0 ) * HBx(m,1,1)**2

       winv       = 1.d0 / HBx(m,1,1)
       cbeta(1)   = cbeta(1) + coefb(1,m) * winv
       cbeta(2)   = cbeta(2) + coefb(2,m) * winv
       cbeta(3)   = cbeta(3) + coefb(3,m) * winv - winv
       
    enddo
    cbeta(3) = cbeta(3) - U0

    !   -- [1-2] Initial Guess of Beta  ---   !
    beta     = cbeta(2)**2 - 4.d0 * cbeta(1) * cbeta(3)
       
    !  ---  [2] Newton Raphson Loop  ---  !
100 continue
    res = beta
    do k=1, itermax

       !  -- [2-1]  Caluculation U, dU/db ---  !
       val = 0.d0
       dif = 0.d0
       do m=1, Np
          val = val + sqrt( coefb(1,m)*beta**2 + coefb(2,m)*beta + coefb(3,m) )
          dif = dif + ( coefb(1,m) * beta + 0.5d0*coefb(2,m) ) &
               &        / sqrt( coefb(1,m)*beta**2 + coefb(2,m)*beta + coefb(3,m) )
       enddo
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
       beta = ( - cbeta(2) - sqrt( cbeta(2)**2 - 4.d0 * cbeta(1) * cbeta(3) ) ) / ( 2.d0 * cbeta(1) )
       goto 100
    endif

    return
  end subroutine solv_beta

  
  
  subroutine coalesptcl( pxvh, dr, dz, Nnew, Nold, ierr )
    use constants , only : rp_, zp_, vr_, vt_, vz_, wp_
    use mt19937Mod, only : multi_grnd
    !$ use omp_lib
    implicit none
    integer         , intent(in)    :: Nnew, Nold
    double precision, intent(in)    :: dr, dz
    integer         , intent(inout) :: ierr
    double precision, intent(inout) :: pxvh(8,Nold)
    integer         , parameter     :: Nsub       = 8
    integer         , parameter     :: iterlim    = 10
    double precision, parameter     :: vMax       = 0.5d0
    character(7)    , parameter     :: prob       = 'uniform'
    integer                         :: m, k, i, ip, jp, i0, j0, iter, mythread
    double precision                :: drinv, dzinv, subrinv, subzinv, beta, coef
    double precision                :: mu(3,2), sd(3,2), gamma, tst, pdfmax, U0, U1, winv, wmin
    double precision                :: pdf(0:Nsub+1,0:Nsub+1), pxvhold(8,Nold), wcdf(0:Nold,2), psmp(8,Nnew)
    double precision                :: lamb(-2:2,-2:2,Nnew,4), HBx(Nnew,4,2)
    double precision                :: mmt(-2:2,-2:2,4), mmt2(-2:2,-2:2,4)
    
    ! --- [1] preparation  ---  !
    !  -- [1-1] parameters  --  !
    iter    = iterlim
    !$ mythread = omp_get_thread_num() + 1
    drinv   = 1.d0 / dr
    dzinv   = 1.d0 / dz
    subrinv = dble( Nsub    ) * drinv
    subzinv = dble( Nsub    ) * dzinv
    j0      = nint( pxvh(rp_,1) * drinv )
    i0      = nint( pxvh(zp_,1) * dzinv )
    do m=1, Nold
       do k=1, 8
          pxvhold(k,m) = pxvh(k,m)
       enddo
    enddo
    ! -- [1-2] Average, Deviation -- !
    do k=1, 3
       coef       = 0.d0
       mu(k,1)    = 0.d0
       sd(k,1)    = 0.d0
       do m=1, Nold
          coef    = coef    + pxvh(wp_,m)
          mu(k,1) = mu(k,1) + pxvh(wp_,m) * pxvh(k+2,m)
          sd(k,1) = sd(k,1) + pxvh(wp_,m) * pxvh(k+2,m)**2
       enddo
       coef       =    1.d0 / coef
       mu(k,1)    = mu(k,1) * coef
       sd(k,1)    = sqrt( sd(k,1)*coef - mu(k,1)**2 )
    enddo
    
    ! -------------------------------------------- !
    ! --- [2] Begining of Assous-Welch Method  --- !
    ! -------------------------------------------- !    
    !  -- [2-1] Sub-grid Particle Distribution  -- !
    pdf(:,:)      = 0.d0
    do m=1, Nold
       jp         = ceiling( ( ( pxvh(rp_,m)+0.5d0*dr ) - dr*dble(j0) ) * subrinv )
       ip         = ceiling( ( ( pxvh(zp_,m)+0.5d0*dz ) - dz*dble(i0) ) * subzinv )
       pdf(ip,jp) = pdf(ip,jp) + 1.d0 / dble( Nold )
    enddo
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
    call getMoment( mmt, pxvh, dr, dz, Nold, i0, j0, U0 )

    ! ---------------------------- !
    ! === *** Return Point *** === !
    ! ---------------------------- !
200 continue
    iter      = iter - 1
    ierr      = 0
    if ( iter .lt. 0 ) then
       do m=1, Nold
          do k=1, 8
             pxvh(k,m) = pxvhold(k,m)
          enddo
       enddo
       ierr = 1
       return
    endif ! --- No Match --- !
    
    !  -- [2-3] Determine xm by sub-grid method  -- !
    do m=1, Nnew
210    continue
       
       pxvh(rp_,m)  = dr * ( dble(j0) + multi_grnd(mythread) - 0.5d0 )
       pxvh(zp_,m)  = dz * ( dble(i0) + multi_grnd(mythread) - 0.5d0 )
       jp         = ceiling( ( ( pxvh(rp_,m)+0.5d0*dr ) - dr*dble(j0) ) * subrinv )
       ip         = ceiling( ( ( pxvh(zp_,m)+0.5d0*dz ) - dz*dble(i0) ) * subzinv )
       coef       = pdfmax * multi_grnd(mythread)
       if ( coef.ge.pdf(ip,jp)          ) goto 210
    enddo

    !  -- [2-4] Get Lambda Matrix for m-particles Shape Functions -- !
    call getLambda( lamb, pxvh(:,1:Nnew), dr, dz, Nnew, i0, j0 )
    
    !  -- [2-5] Get Weight of charge / momentum conservation :: gj, hj --  !
    call getWeight( lamb(-2:1,-2:1,:,1), mmt(-2:1,-2:1,1), HBx(:,1,1), 16, Nnew, ierr )
    call getWeight( lamb(-2:1,-1:1,:,2), mmt(-2:1,-1:1,2), HBx(:,2,1), 12, Nnew, ierr )
    call getWeight( lamb(-2:1,-2:1,:,3), mmt(-2:1,-2:1,3), HBx(:,3,1), 16, Nnew, ierr )
    call getWeight( lamb(-1:1,-2:1,:,4), mmt(-1:1,-2:1,4), HBx(:,4,1), 12, Nnew, ierr )
    ! write(6,*) ' Here 4 '

    if ( ierr.eq.1 ) goto 200
    ! write(6,*) ' Here 4-1 '
    
    !  -- [2-6] Check Whether Moderate Weight --  !
    mu(1,2)       = 0.d0
    sd(1,2)       = 0.d0
    do m=1, Nnew
       mu(1,2)    = mu(1,2) + HBx(m,1,1)
       sd(1,2)    = sd(1,2) + HBx(m,1,1)**2
    enddo
    coef          = 1.d0    / dble( Nnew )
    mu(1,2)       = mu(1,2) * coef
    sd(1,2)       = 1.d0    / sqrt( sd(1,2) * coef - mu(1,2)**2 )
    ! write(6,*) ' Here 5 '
    
    wmin          = 1.d6
    U1            = 0.d0
    tst           = 0.d0
    do m=1, Nnew
       winv       = 1.d0 / HBx(m,1,1)
       HBx(m,2,1) = HBx(m,2,1) * winv
       HBx(m,3,1) = HBx(m,3,1) * winv
       HBx(m,4,1) = HBx(m,4,1) * winv
       wmin       = min( HBx(m,1,1), wmin )
       tst        = tst + max( 0.d0, abs( ( HBx(m,1,1) - mu(1,2) ) * sd(1,2) ) - 3.d0 )
       U1         = U1  + HBx(m,1,1) * sqrt( 1.d0 + HBx(m,2,1)**2 + HBx(m,3,1)**2 + HBx(m,4,1)**2 )
    enddo
    ! write(6,*) ' Here 6 '

    if ( ( tst.gt.0.d0 ).or.( wmin.lt.0.d0 ).or.( U1.gt.U0 ) ) goto 200
    ! if ( ( tst.gt.0.d0 ).or.( wmin.lt.0.d0 ) ) goto 200

    !  -- [2-7] Reconstruction of VDF  --  !
    ! write(6,*) ' Here 7 '
220 continue
    
    ! -- Probability of choose m-th particle as Refference -- !
    wcdf(:,:)  = 0.d0
    if ( prob.eq.'uniform' ) then
       do m=1, Nold
          wcdf(m,1) = 1.d0
          wcdf(m,2) = wcdf(m-1,2) + wcdf(m,1)
       enddo
    endif
    if ( prob.eq.'weightd' ) then
       do m=1, Nold
          wcdf(m,1) = pxvh(wp_,m)
          wcdf(m,2) = wcdf(m-1,2) + wcdf(m,1)
       enddo
    endif
    winv      = 1.d0 / wcdf(Nold,2)
    wcdf(:,2) = wcdf(:,2) * winv

    ! -- Choose Reference particle randomly -- !
    do m=1, Nnew
230    continue
       coef = multi_grnd( mythread )
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
                wcdf(i,2) = wcdf(i-1,2) + wcdf(i,1)
             enddo
             coef = 1.d0 / wcdf(Nold,2)
             do i=1, Nold
                if ( wcdf(i,1).eq.0.d0 ) then
                   wcdf(i,2) = 0.d0
                else
                   wcdf(i,2) = wcdf(i,2) * coef
                endif
             enddo
             goto 230
          endif
       enddo
    enddo
    
    !  -- [2-8] Bloadening of Distribution :: Random Velocity --  !
    !   - Coefficients For Random Velocity --  !
    call getMoment( mmt2, psmp, dr, dz, Nnew, i0, j0, U1 )

    call getWeight( lamb(-2:1,-1:1,:,2), - mmt2(-2:1,-1:1,2), HBx(:,2,2),12, Nnew, ierr )
    call getWeight( lamb(-2:1,-2:1,:,3), - mmt2(-2:1,-2:1,3), HBx(:,3,2),16, Nnew, ierr )
    call getWeight( lamb(-1:1,-2:1,:,4), - mmt2(-1:1,-2:1,4), HBx(:,4,2),12, Nnew, ierr )
    if ( ierr.eq.1 ) goto 200
    !   - Random Velocity -  !
    do m=1, Nnew
       HBx(m,2,2) = HBx(m,2,2) + psmp(3,m)
       HBx(m,3,2) = HBx(m,3,2) + psmp(4,m)
       HBx(m,4,2) = HBx(m,4,2) + psmp(5,m)
    enddo ! px^ran(xm) = Sum( bjlj ) + px^sam @ Eq(1.13) Welch et.al.,

    !  -- [2-9]  Energy Adjustment --  !
    beta = 0.1d0
    call solv_beta( HBx, U0, beta, Nnew, ierr )
    if ( ierr.eq.1 ) goto 200
    ! write(6,*) beta

    !   - Update Particle Info. -  !
    do m=1, Nnew
       coef       = beta / HBx(m,1,1)
       pxvh(vr_,m)  = HBx(m,2,1) + coef * HBx(m,2,2)
       pxvh(vt_,m)  = HBx(m,3,1) + coef * HBx(m,3,2)
       pxvh(vz_,m)  = HBx(m,4,1) + coef * HBx(m,4,2)
       pxvh(wp_,m)  = HBx(m,1,1)
    enddo

    !  -- [2-10] c-Limit Check -- !
    do m=1, Nnew
       gamma = ( pxvh(vr_,m)**2 + pxvh(vt_,m)**2 + pxvh(vz_,m)**2 )
       if ( gamma .gt. vMax ) then
          do i=1, Nnew
             do k=1, 8
                pxvh(k,i) = pxvhold(k,i)
             enddo
          enddo
          goto 200
       endif
    enddo
    
    return
  end subroutine coalesptcl

  
  subroutine VDFBinning( pxvh, npth, APRperCell, APRvdfResl, Nold, pxvr, idx )
    use constants, only : rp_, zp_, vr_, vt_, vz_, wp_
    implicit none
    integer         , intent(in)  :: npth, APRperCell, APRvdfResl
    double precision, intent(in)  :: pxvh(8,npth)
    integer         , intent(out) :: idx(npth,APRperCell), Nold(APRperCell)
    double precision, intent(out) :: pxvr(8,npth,APRperCell)
    integer                       :: ip, jp, kp, k, m, count, minTopN, minIdx
    integer                       :: vdf(APRvdfResl,APRvdfResl,APRvdfResl), topN(2,APRperCell), goban(npth)
    double precision              :: dvrinv, dvtinv, dvzinv
    double precision              :: vMaxMin(3,2)

    ! --- [1] VDF Axis (vr,vt,vz) Making --- !
    vMaxMin(:,1)     =  1.d0
    vMaxMin(:,2)     = -1.d0
    do m=1, npth
       vMaxMin(1,1)  =  min( vMaxMin(1,1), pxvh(vr_,m) )
       vMaxMin(2,1)  =  min( vMaxMin(2,1), pxvh(vt_,m) )
       vMaxMin(3,1)  =  min( vMaxMin(3,1), pxvh(vz_,m) )
       vMaxMin(1,2)  =  max( vMaxMin(1,2), pxvh(vr_,m) )
       vMaxMin(2,2)  =  max( vMaxMin(2,2), pxvh(vt_,m) )
       vMaxMin(3,2)  =  max( vMaxMin(3,2), pxvh(vz_,m) )
    enddo ! :: vMaxMin = ( Min, Max )
    dvrinv           = dble( APRvdfResl ) / ( vMaxMin(1,2) - vMaxMin(1,1) )
    dvtinv           = dble( APRvdfResl ) / ( vMaxMin(2,2) - vMaxMin(2,1) )
    dvzinv           = dble( APRvdfResl ) / ( vMaxMin(3,2) - vMaxMin(3,1) )

    ! --- [2] Count VDF --- !
    vdf(:,:,:)       = 0
    do m=1, npth
       ip            = min( max( 1, ceiling( ( pxvh(vr_,m) - vMaxMin(1,1) ) * dvrinv ) ), APRvdfResl )
       jp            = min( max( 1, ceiling( ( pxvh(vt_,m) - vMaxMin(2,1) ) * dvtinv ) ), APRvdfResl )
       kp            = min( max( 1, ceiling( ( pxvh(vz_,m) - vMaxMin(3,1) ) * dvzinv ) ), APRvdfResl )
       vdf(ip,jp,kp) = vdf(ip,jp,kp) + 1
       goban(m)      = (ip-1) + (jp-1)*APRvdfResl + (kp-1)*APRvdfResl*APRvdfResl + 1
    enddo
    
    ! --- [3] Search for TopN VDF index --- !
    do k=1, APRperCell
       topN(1,k) = 0
       topN(2,k) = 0
    enddo
    minTopN = 0
    minIdx  = 1
    do kp=1, APRvdfResl
       do jp=1, APRvdfResl
          do ip=1, APRvdfResl
             if ( vdf(ip,jp,kp).gt.minTopN ) then
                ! -- Substitute into smallest position -- !
                topN(1,minIdx) = (ip-1) + (jp-1)*APRvdfResl + (kp-1)*APRvdfResl*APRvdfResl + 1
                topN(2,minIdx) = vdf(ip,jp,kp)
                ! -- Update minTopN & minIdx -- !
                minIdx         = 1
                minTopN        = topN(2,minIdx)
                do k=2, APRperCell
                   if ( topN(2,k).lt.minTopN ) then
                      minIdx   = k
                      minTopN  = topN(2,k)
                   endif
                enddo
             endif
          end do
       enddo
    enddo

    ! --- [4] Make New particle Set --- !
    do k=1, APRperCell
       count = 0
       do m=1, npth
          if ( goban(m).eq.topN(1,k) ) then
             count           = count + 1
             pxvr(:,count,k) = pxvh(:,m)
             idx(count,k)    = m
          endif
       enddo
       Nold(k) = count
    enddo
    
    return
  end subroutine VDFBinning

end module ptCoalsMod
