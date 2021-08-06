module collideMod
  implicit none
  integer              :: rtime
  integer              :: flag = 0
  double precision     :: ctime
  double precision     :: accTime(2,8) = 0.d0  
contains
  
  
  subroutine collision
    use constants , only : npt, ns, np, LI, LJ, OMPNumThreads
    use constants , only : q, rm, vthcv, dt, drinv, dzinv
    use constants , only : nuwce
    use constants , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables , only : pxv, goban, perCell, cellTop
    use variables , only : RhoPerNsp, w_p
    use sortingMod, only : ParallelSort
    use randGenMod, only : gaussdev, multi_gdev
    use mt19937Mod, only : grnd, multi_grnd
    use momentsMod, only : nCountUp_CC, nCountUp_SP
    !$ use omp_lib    
    implicit none
    integer             :: i, j, k, m, ip, jp, ic, jc, ii, n1, n2
    integer             :: maxcol, rstpt, idx, first, last, count, alloc_err
    integer             :: idxe, idxi, icat1, ecat1, icat2, ecat2
    integer             :: colend(3), triflag(3), rndindx(npt)
    integer             :: myTHnum
    integer, parameter  :: chunk = 7
    integer, parameter  :: Ncell = LI*LJ
    double precision    :: phi, delta, sinth, cosmod, sincos, coef
    double precision    :: uperp, uabsl, uxup, uyup, uzmod, uamod
    double precision    :: mu(3), delcoef(3), nlower(3), ulab(3), rmc(2,3), wt(2)
    double precision, parameter   :: twopi = 8.d0*atan(1.d0)
    integer         , allocatable :: collbuf(:,:)
    double precision, allocatable :: deltau(:,:)
    
    
    ! --- (0) preparation --- !
    ! call TATimeMeasure(0)
    
    delcoef(1) = 0.5d0 * nuwce * vthcv**3 * ( q(1)*q(1) )**2 * dt
    delcoef(2) = 0.5d0 * nuwce * vthcv**3 * ( q(2)*q(2) )**2 * dt
    delcoef(3) = 0.5d0 * nuwce * vthcv**3 * ( q(1)*q(2) )**2 * dt
    
    rmc(:,1)   = rm(1)
    rmc(:,2)   = rm(2)
    rmc(1,3)   = rm(1)
    rmc(2,3)   = rm(2)

    
    ! ===========     [1] Preparation    =========== !
    ! --- (1-1) Categorize each particle into boxes :: GoBan --- !
    
    n2 = 0
    do k=1, ns
       !$omp parallel default(none) shared(k,np,pxv,goban,dzinv,drinv) &
       !$omp private(m,ip,jp)
       !$omp do
       do m=1, np(k)
          ip       = nint( pxv(zp_,m,k)*dzinv )
          jp       = nint( pxv(rp_,m,k)*drinv )
          goban(m) = Ncell*(k-1) + LI*jp + (ip-1) + 1
       enddo ! Simplicity -> jp=0 OK :: ip=0 xxx ( due to reflection )
       !$omp end do
       !$omp end parallel

    enddo

    !  --- (1-2) Sort according to Cell Number ---  !
    ! call sortptcl( npt, goban, 1, npt, pxv )
    call ParallelSort( npt, goban, 1, npt, pxv, OMPNumThreads )


    !  --- (1-4) count up Number of SuperParticle in a cell.
    call nCountUp_CC( w_p )
    call nCountUp_SP( percell )
    n2 = 1
    do k=1, ns
       count = n2
       do jc=1, LJ
          do ic=1, LI
             cellTop(ic,jc,k) = count 
             count            = count + perCell( ic,jc,k )
          enddo
       enddo
       n2 = n2 + np(k)
    enddo

    ! pmm(:,:,:) = 0.d0
    ! Ek(:,:,:)  = 0.d0
    
    ! first = cellTop(LI/2,LJ/2,1)
    ! last  = cellTop(LI/2,LJ/2,1) + perCell(LI/2,LJ/2,1) - 1
    ! do m=first, last
    !    pmm(1,1,1) = pmm(1,1,1) + pxv(vr_,m)*rm(1)*pxv(wp_,m)
    !    pmm(2,1,1) = pmm(2,1,1) + pxv(vt_,m)*rm(1)*pxv(wp_,m)*pxv(rp_,m)
    !    pmm(3,1,1) = pmm(3,1,1) + pxv(vz_,m)*rm(1)*pxv(wp_,m)
    !    gamma= 1.d0 / ( 1.d0 + sqrt( 1.d0 + ( pxv(vr_,m)**2 + pxv(vt_,m)**2 + pxv(vz_,m)**2 ) ) )
    !    Ek(1,1,1) = Ek(1,1,1) + pxv(vr_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(2,1,1) = Ek(2,1,1) + pxv(vt_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(3,1,1) = Ek(3,1,1) + pxv(vz_,m)**2 * gamma * pxv(wp_,m)
    ! enddo
    ! first = cellTop(LI/2,LJ/2,2)
    ! last  = cellTop(LI/2,LJ/2,2) + perCell(LI/2,LJ/2,2) - 1
    ! do m=first, last
    !    pmm(1,2,1) = pmm(1,2,1) + pxv(vr_,m)*rm(2)*pxv(wp_,m)
    !    pmm(2,2,1) = pmm(2,2,1) + pxv(vt_,m)*rm(2)*pxv(wp_,m)*pxv(rp_,m)
    !    pmm(3,2,1) = pmm(3,2,1) + pxv(vz_,m)*rm(2)*pxv(wp_,m)
    !    gamma= 1.d0 / ( 1.d0 + sqrt( 1.d0 + ( pxv(vr_,m)**2 + pxv(vt_,m)**2 + pxv(vz_,m)**2 ) ) )
    !    Ek(1,2,1) = Ek(1,2,1) + pxv(vr_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(2,2,1) = Ek(2,2,1) + pxv(vt_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(3,2,1) = Ek(3,2,1) + pxv(vz_,m)**2 * gamma * pxv(wp_,m)
    ! enddo


    ! --- (3) shuffle --- !
    
    !$omp parallel shared(cellTop,perCell,rndindx)
    !$omp do private(i)
    do i=1, npt
       rndindx(i) = i
    enddo
    !$omp end do
    !$omp do private(k,ic,jc,first,last)
    do k=1, ns
       do jc=1, LJ
          do ic=1, LI
             first = cellTop(ic,jc,k)
             last  = cellTop(ic,jc,k) + perCell(ic,jc,k) - 1
             call KnuthShuffleIdx_OMP( rndindx(first:last), perCell(ic,jc,k) )
             ! call Knuth_Shuffle_OMP( first, last )
             ! call Knuth_Shuffle( first, last )
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    
    ! call TATimeMeasure(3)



    
    ! --- (4) make pairs --- !


    maxcol = 2 * maxval( perCell(:,:,:) ) + 3

    !$omp parallel default(none) &
    !$omp shared(perCell,cellTop,rndindx,RhoPerNsp,pxv,delcoef,maxcol,rmc,w_p) &
    !$omp private(count,triflag,k,rstpt,idx,collbuf,colend,nlower) &
    !$omp private(idxe,idxi,ii,ecat1,icat1,ecat2,icat2,i,j,myTHnum) &
    !$omp private(ulab,n1,n2,uperp,uabsl,phi,delta,sinth,cosmod) &
    !$omp private(sincos,uyup,uxup,uzmod,uamod,deltau,alloc_err,coef,mu,wt,ic,jc)
    allocate( collbuf( 2, maxcol ) , deltau( 3, maxcol ) )
    !$ myTHnum = omp_get_thread_num()
    !$omp do schedule(static,chunk)
    do jc=2, LJ-1
       do ic=2, LI-1

          ! electron-electron , ion-ion !
          count      = 1
          triflag(:) = 0
          do k=1, ns
             rstpt   = perCell(ic,jc,k)
             idx     = cellTop(ic,jc,k)

             if ( ( perCell(ic,jc,k).ge.3 ).and.( mod(perCell(ic,jc,k),2).eq.1 ) ) then
                collbuf( 1, count   ) = rndindx( idx     )
                collbuf( 2, count   ) = rndindx( idx + 1 )
                collbuf( 1, count+1 ) = rndindx( idx     )
                collbuf( 2, count+1 ) = rndindx( idx + 2 )
                collbuf( 1, count+2 ) = rndindx( idx + 1 )
                collbuf( 2, count+2 ) = rndindx( idx + 2 )

                count      = count + 3
                rstpt      = rstpt - 3
                idx        = idx   + 3
                triflag(k) = 3

             endif             ! -- 3-Body collision -- !


             do i=rstpt, 2, -2
                collbuf( 1, count ) = rndindx( idx     )
                collbuf( 2, count ) = rndindx( idx + 1 )
                count  = count + 1
                idx    = idx   + 2
             enddo             ! -- 2-Body collision -- !

             colend(k) = count - 1
             nlower(k) = w_p(ic,jc,k)

          enddo


          !      electron - ion       !
          if ( ( perCell(ic,jc,1).ne.0 ).and.( perCell(ic,jc,2).ne.0 ) ) then
             if ( perCell(ic,jc,1) .eq. perCell(ic,jc,2) ) then

                rstpt     = perCell(ic,jc,1)
                nlower(3) = w_p(ic,jc,1)

                idxe      = cellTop(ic,jc,1)
                idxi      = cellTop(ic,jc,2)
                do i=rstpt, 1, -1
                   collbuf( 1, count ) = rndindx( idxe )
                   collbuf( 2, count ) = rndindx( idxi )
                   idxi                = idxi  + 1
                   idxe                = idxe  + 1
                   count               = count + 1
                enddo

             else if ( perCell(ic,jc,1) .lt. perCell(ic,jc,2) ) then

                nlower(3) = w_p(ic,jc,1)

                ii        =   perCell(ic,jc,2) / perCell(ic,jc,1)
                ecat1     = ( perCell(ic,jc,2) - ii*perCell(ic,jc,1) )
                icat1     = ( ii+1 ) * ecat1
                ecat2     = perCell(ic,jc,1) - ecat1
                icat2     = perCell(ic,jc,2) - icat1

                idxi      = cellTop(ic,jc,2)
                do i=1, ii+1
                   idxe   = cellTop(ic,jc,1)
                   do j=1, ecat1
                      collbuf( 1, count ) = rndindx( idxe )
                      collbuf( 2, count ) = rndindx( idxi )
                      idxe                = idxe  + 1
                      idxi                = idxi  + 1
                      count               = count + 1                   
                   enddo
                enddo

                idx = idxe
                do i=1, ii
                   idxe = idx
                   do j=1, ecat2

                      collbuf( 1, count ) = rndindx( idxe )
                      collbuf( 2, count ) = rndindx( idxi )
                      idxe                = idxe  + 1
                      idxi                = idxi  + 1
                      count               = count + 1

                   enddo
                enddo


             else if ( perCell(ic,jc,1) .gt. perCell(ic,jc,2) )  then

                nlower(3) = w_p(ic,jc,2)

                ii        = perCell(ic,jc,1) / perCell(ic,jc,2)
                icat1     = ( perCell(ic,jc,1) - ii*perCell(ic,jc,2) )
                ecat1     = ( ii+1 ) * icat1
                ecat2     = perCell(ic,jc,1) - ecat1
                icat2     = perCell(ic,jc,2) - icat1

                idxe      = cellTop(ic,jc,1)
                do i=1, ii+1
                   idxi   = cellTop(ic,jc,2)
                   do j=1, icat1
                      collbuf( 1, count ) = rndindx( idxe )
                      collbuf( 2, count ) = rndindx( idxi )
                      idxe                = idxe  + 1
                      idxi                = idxi  + 1
                      count               = count + 1
                   enddo
                enddo

                idx = idxi
                do i=1, ii
                   idxi = idx
                   do j=1, icat2
                      collbuf( 1, count ) = rndindx( idxe )
                      collbuf( 2, count ) = rndindx( idxi )
                      idxe                = idxe  + 1
                      idxi                = idxi  + 1
                      count               = count + 1
                   enddo
                enddo

             endif
          end if
          
          colend(3) = count - 1


          ! --- (4) calculate delta u --- !

          n2 = 0
          do k=1, 3
             n1 = n2 + 1
             n2 = colend(k)
             do i=n1, n2

                ulab(1) = pxv(vr_, collbuf(1,i),k ) - pxv(vr_, collbuf(2,i),k )
                ulab(2) = pxv(vt_, collbuf(1,i),k ) - pxv(vt_, collbuf(2,i),k )
                ulab(3) = pxv(vz_, collbuf(1,i),k ) - pxv(vz_, collbuf(2,i),k )

                uperp   = sqrt( ulab(1)**2 + ulab(2)**2 )
                uabsl   = sqrt( ulab(1)**2 + ulab(2)**2 + ulab(3)**2 )



                ! -- reduced mass -- !
                coef    = 1.d0 / ( rmc(1,k)*pxv(wp_,collbuf(1,i),k) + rmc(2,k)*pxv(wp_,collbuf(2,i),k) )
                mu(1)   = coef * rmc(2,k) * pxv(wp_,collbuf(2,i),k) 
                mu(2)   = coef * rmc(1,k) * pxv(wp_,collbuf(1,i),k)
                mu(3)   = coef * rmc(1,k) * pxv(wp_,collbuf(1,i),k) * rmc(2,k)*pxv(wp_,collbuf(2,i),k)
                ! mu(3)   = rmc(1,k) * rmc(2,k) / ( rmc(1,k) + rmc(2,k) )
                coef    = 1.d0 / ( rmc(1,k)*pxv(wp_,collbuf(1,i),k)*pxv(rp_,collbuf(1,i),k) &
                     &           + rmc(2,k)*pxv(wp_,collbuf(2,i),k)*pxv(rp_,collbuf(2,i),k) )
                wt(1)   = coef   * rmc(2,k)*pxv(wp_,collbuf(2,i),k)*pxv(rp_,collbuf(2,i),k)
                wt(2)   = coef   * rmc(1,k)*pxv(wp_,collbuf(1,i),k)*pxv(rp_,collbuf(1,i),k)

                
                ! -- Scattering -- !
                
                ! phi   = twopi * grnd()
                ! delta = sqrt(delta) * gaussdev()
                phi     = twopi * multi_grnd( myTHnum )
                delta   = delcoef(k) * nlower(k) / ( mu(3)**2 * uabsl**3 )
                delta   = delta * ( pxv(wp_,collbuf(1,i),k) * pxv(wp_,collbuf(2,i),k) )**2
                if ( ( ( k.eq.1 ).or.( k.eq.2 ) ).and.( triflag(k).ge.1 ) ) then
                   delta      = 0.5d0 * delta
                   triflag(k) = triflag(k) - 1
                endif
                delta   = sqrt(delta) * multi_gdev( myTHnum )
                sinth   = 2.d0 * delta / ( 1.d0 + delta**2 )
                cosmod  = delta * sinth
                sincos  = sinth * cos( phi )

                uyup    =    1.d0 / uperp
                uxup    = ulab(1) * uyup
                uyup    = ulab(2) * uyup
                uzmod   = ulab(3) * sincos
                uamod   = uabsl   * sinth * sin(phi)

                if ( uperp.eq.0.d0 ) then
                   deltau(1,i) = uabsl * sincos 
                   deltau(2,i) = uamod 
                   deltau(3,i) = - uabsl*cosmod 
                   write(6,*) ' ########### uperp is zero ############ ', uperp
                   write(6,'(a,3(e10.3))') 'vx = ', pxv(vr_, collbuf(1,i),k ), pxv(vr_, collbuf(2,i),k ), ulab(1)
                   write(6,'(a,3(e10.3))') 'vy = ', pxv(vt_, collbuf(1,i),k ), pxv(vt_, collbuf(2,i),k ), ulab(2)
                   write(6,'(a,3(e10.3))') 'vz = ', pxv(vz_, collbuf(1,i),k ), pxv(vz_, collbuf(2,i),k ), ulab(3)
                else
                   deltau(1,i) = uxup*uzmod - uyup*uamod - ulab(1)*cosmod
                   deltau(2,i) = uyup*uzmod + uxup*uamod - ulab(2)*cosmod
                   deltau(3,i) = - uperp * sincos        - ulab(3)*cosmod
                end if

                ! if ( abs( pxv(rp_,collbuf(1,i)) - pxv(rp_,collbuf(2,i)) ).gt.dz ) then
                !    write(6,'(a,2(e10.3,1x),a,e10.3)') &
                !         & 'irregular grid x ::', pxv(rp_,collbuf(1,i)), pxv(rp_,collbuf(2,i)), 'dz=', dz
                !    write(6,*) ' i = ', i
                !    stop
                ! endif
                ! if ( abs( pxv(zp_,collbuf(1,i) ) - pxv(zp_, collbuf(2,i) ) ).gt.dy ) then                
                !    write(6,'(a,2(e10.3,1x),a,e10.3)') &
                !         & 'irregular grid y ::', pxv(zp_,collbuf(1,i)), pxv(zp_,collbuf(2,i)), 'dy=', dy
                !    write(6,*) ' i = ', i
                !    stop
                ! endif

                pxv(vr_, collbuf(1,i),k ) = pxv(vr_, collbuf(1,i),k ) + mu(1) * deltau(1,i)
                pxv(vr_, collbuf(2,i),k ) = pxv(vr_, collbuf(2,i),k ) - mu(2) * deltau(1,i)
                pxv(vt_, collbuf(1,i),k ) = pxv(vt_, collbuf(1,i),k ) + wt(1) * deltau(2,i)
                pxv(vt_, collbuf(2,i),k ) = pxv(vt_, collbuf(2,i),k ) - wt(2) * deltau(2,i)
                pxv(vz_, collbuf(1,i),k ) = pxv(vz_, collbuf(1,i),k ) + mu(1) * deltau(3,i)
                pxv(vz_, collbuf(2,i),k ) = pxv(vz_, collbuf(2,i),k ) - mu(2) * deltau(3,i)

             enddo
          enddo
       enddo
    enddo
    !$omp end do
    deallocate( collbuf, deltau, stat=alloc_err )
    if ( alloc_err .ne. 0 ) stop ' [ERROR] deallocate statement in grouping routine [ERROR] '
    !$omp end parallel
    
    
    ! call TATimeMeasure(4)
    ! call TestPtColl( cellTop, perCell )
    ! if ( mod( kstep, 20 ) .eq. 0 ) call TATimeMeasure(-1)

    ! --- momentum conservation --- !

    ! first = cellTop(LI/2,LJ/2,1)
    ! last  = cellTop(LI/2,LJ/2,1) + perCell(LI/2,LJ/2,1) - 1
    ! do m=first, last
    !    pmm(1,1,2) = pmm(1,1,2) + pxv(vr_,m)*rm(1)*pxv(wp_,m)
    !    pmm(2,1,2) = pmm(2,1,2) + pxv(vt_,m)*rm(1)*pxv(wp_,m)*pxv(rp_,m)
    !    pmm(3,1,2) = pmm(3,1,2) + pxv(vz_,m)*rm(1)*pxv(wp_,m)
    !    gamma= 1.d0 / ( 1.d0 + sqrt( 1.d0 + ( pxv(vr_,m)**2 + pxv(vt_,m)**2 + pxv(vz_,m)**2 ) ) )
    !    Ek(1,1,2) = Ek(1,1,2) + pxv(vr_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(2,1,2) = Ek(2,1,2) + pxv(vt_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(3,1,2) = Ek(3,1,2) + pxv(vz_,m)**2 * gamma * pxv(wp_,m)
    ! enddo
    ! first = cellTop(LI/2,LJ/2,2)
    ! last  = cellTop(LI/2,LJ/2,2) + perCell(LI/2,LJ/2,2) - 1
    ! do m=first, last
    !    pmm(1,2,2) = pmm(1,2,2) + pxv(vr_,m)*rm(2)*pxv(wp_,m)
    !    pmm(2,2,2) = pmm(2,2,2) + pxv(vt_,m)*rm(2)*pxv(wp_,m)*pxv(rp_,m)
    !    pmm(3,2,2) = pmm(3,2,2) + pxv(vz_,m)*rm(2)*pxv(wp_,m)
    !    gamma= 1.d0 / ( 1.d0 + sqrt( 1.d0 + ( pxv(vr_,m)**2 + pxv(vt_,m)**2 + pxv(vz_,m)**2 ) ) )
    !    Ek(1,2,2) = Ek(1,2,2) + pxv(vr_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(2,2,2) = Ek(2,2,2) + pxv(vt_,m)**2 * gamma * pxv(wp_,m)
    !    Ek(3,2,2) = Ek(3,2,2) + pxv(vz_,m)**2 * gamma * pxv(wp_,m)
    ! enddo

    ! write(6,*)
    ! do k=1, 3
    !    write(6,'(3(e12.5,1x))') pmm(k,1,1), pmm(k,2,1), pmm(k,1,1)+pmm(k,2,1)
    !    write(6,'(3(e12.5,1x))') pmm(k,1,2), pmm(k,2,2), pmm(k,1,2)+pmm(k,2,2)
    !    write(6,'(1(e12.5,1x))') ( pmm(k,1,1)+pmm(k,2,1) )-( pmm(k,1,2)+pmm(k,2,2) )
    !    write(6,*)
    !    write(6,'(3(e12.5,1x))') Ek(k,1,1), Ek(k,2,1), Ek(k,1,1)+Ek(k,2,1)
    !    write(6,'(3(e12.5,1x))') Ek(k,1,2), Ek(k,2,2), Ek(k,1,2)+Ek(k,2,2)
    !    write(6,'(1(e12.5,1x))') ( Ek(k,1,1)+Ek(k,2,1) )-( Ek(k,1,2)+Ek(k,2,2) )
    !    write(6,*)
    !    write(6,*)
    ! enddo
    
    return
  end subroutine collision
  
  
  
  
  subroutine Knuth_Shuffle_OMP( first, last )

    use variables , only : pxv
    use constants , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use mt19937Mod, only : grnd, multi_grnd
    !$ use omp_lib
    
    implicit none
    integer, intent(in) :: first, last
    integer             :: i, k, randpos, mythread
    double precision    :: r, tmp
    
    !$ mythread = omp_get_thread_num()

    do k=1, 2
       do i=last, first+1, -1

          r = multi_grnd( mythread )
          randpos = int( r*( i-first ) ) + first

          tmp = pxv(rp_,randpos,k); pxv(rp_,randpos,k) = pxv(rp_,i,k) ; pxv(rp_,i,k) = tmp
          tmp = pxv(zp_,randpos,k); pxv(zp_,randpos,k) = pxv(zp_,i,k) ; pxv(zp_,i,k) = tmp
          tmp = pxv(vr_,randpos,k); pxv(vr_,randpos,k) = pxv(vr_,i,k) ; pxv(vr_,i,k) = tmp
          tmp = pxv(vt_,randpos,k); pxv(vt_,randpos,k) = pxv(vt_,i,k) ; pxv(vt_,i,k) = tmp
          tmp = pxv(vz_,randpos,k); pxv(vz_,randpos,k) = pxv(vz_,i,k) ; pxv(vz_,i,k) = tmp

       enddo
    enddo
    return
  end subroutine Knuth_Shuffle_OMP



  
  subroutine KnuthShuffleIdx_OMP( idx, cellsize )

    !$ use omp_lib
    use mt19937Mod, only : multi_grnd
    
    implicit none
    integer                :: i, tmp, randpos, mythread
    double precision       :: r
    integer, intent(in)    :: cellsize
    integer, intent(inout) :: idx(:)
    
    
    !$ mythread = omp_get_thread_num()

    do i=cellsize, 2, -1
       ! -- Indexing -- !
       r            = multi_grnd( mythread )
       randpos      = int(r*i) + 1
       ! -- Exchange -- !
       tmp          = idx(randpos)
       idx(randpos) = idx(i)
       idx(i)       = tmp       
    enddo
    

    return
  end subroutine KnuthShuffleIdx_OMP
  


  
  subroutine Knuth_Shuffle( first, last )

    use variables, only : pxv
    use constants, only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use mt19937Mod, only : grnd

    implicit none
    integer, intent(in) :: first, last
    integer             :: i, k, randpos
    double precision    :: r, tmp

    do k=1, 2
       do i=last, first+1, -1

          r = grnd()
          randpos = int( r*( i-first ) ) + first

          tmp = pxv(rp_,randpos,k); pxv(rp_,randpos,k) = pxv(rp_,i,k) ; pxv(rp_,i,k) = tmp
          tmp = pxv(zp_,randpos,k); pxv(zp_,randpos,k) = pxv(zp_,i,k) ; pxv(zp_,i,k) = tmp
          tmp = pxv(vr_,randpos,k); pxv(vr_,randpos,k) = pxv(vr_,i,k) ; pxv(vr_,i,k) = tmp
          tmp = pxv(vt_,randpos,k); pxv(vt_,randpos,k) = pxv(vt_,i,k) ; pxv(vt_,i,k) = tmp
          tmp = pxv(vz_,randpos,k); pxv(vz_,randpos,k) = pxv(vz_,i,k) ; pxv(vz_,i,k) = tmp

       enddo
    enddo
    
    return
  end subroutine Knuth_Shuffle


  
  

  

  subroutine TATimeMeasure( subject )
    
    implicit none

    integer :: i, k, rtime_, trate, t_max
    double precision :: ctime_, hms(3)
    character(41) :: fmt

    integer, intent(in) :: subject
    
    if ( subject .eq. 0 ) then ! initialize
       
       call system_clock( rtime )
       call cpu_time( ctime )
       if ( flag .eq. 0 ) then
          accTime(:,:) = 0.d0
          flag = 1
       endif
          
    else if ( subject .eq. -1 ) then

       
       write(6,*)
       write(6,*) ' ---------------   Computation Time   --------------- '
       write(6,*)


       accTime(:,8) = 0.d0
       do k=1, 2
          do i=1, 4
             accTime(k,8) = accTime(k,8) + accTime(k,i)
          enddo

          if ( k.eq.1 ) write(6,*) ' ----       CPU    Time  ---- '
          if ( k.eq.2 ) write(6,*) ' ----       Real   Time  ---- '
          hms = sec2hms( accTime(k,1) )
          fmt = '(a,2x,e12.5,3x,2(i3,a),f7.3,a,3x,f7.3,a)'
          write(6,fmt) ' time pt(1)   :: ', accTime(k,1) &
               & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,1)/accTime(k,8)*100.0, ' % '
          hms = sec2hms( accTime(k,2) )
          write(6,fmt) ' time pt(2)   :: ', accTime(k,2) &
               & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,2)/accTime(k,8)*100.0, ' % '
          hms = sec2hms( accTime(k,3) )
          write(6,fmt) ' time pt(3)   :: ', accTime(k,3) &
               & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,3)/accTime(k,8)*100.0, ' % '
          hms = sec2hms( accTime(k,4) )
          write(6,fmt) ' time pt(4)   :: ', accTime(k,4) &
               & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,4)/accTime(k,8)*100.0, ' % '
         hms = sec2hms( accTime(k,8) )  
          write(6,fmt) ' ---- total ---- ', accTime(k,8) &
               & , int(hms(1)), ' h ', int(hms(2)), ' m ', hms(3), ' s ', accTime(k,8)/accTime(k,8)*100.0, ' % '
          write(6,*) ' ----------------------------- '
          write(6,*)
          
       enddo

    else if ( ( subject.gt.0 ).and.( subject.le.8 ) ) then
       
       ! CPU time 
       ctime_ = ctime
       call cpu_time( ctime )
       accTime(1,subject) = accTime(1,subject) + ( ctime - ctime_ )
       
       ! Real time
       rtime_ = rtime
       call system_clock( rtime, trate, t_max )
       if ( rtime .lt. rtime_ ) then
          accTime(2,subject) = accTime(2,subject) + ( ( t_max - rtime_ ) + rtime + 1 ) / dble( trate )
       else
          accTime(2,subject) = accTime(2,subject) + ( rtime - rtime_ ) / dble( trate )
       endif
    endif
       
    return
  end subroutine TATimeMeasure
  
  





  function sec2hms( sec )

    implicit none

    double precision, intent(in) :: sec
    real :: sec2hms(3)

    sec2hms(1) = real( int( sec ) / 3600 )
    sec2hms(2) = real( int( sec - 3600.0*sec2hms(1) ) / 60 )
    sec2hms(3) = real( sec - 3600.0*sec2hms(1) - 60.0*sec2hms(2) )

    return
  end function sec2hms

  
  

end module collideMod
