module randGenMod
  use constants, only : OMPNumThreads
  implicit none
  integer            :: flipflop
  double precision   :: devStock
  integer            :: multi_flipflop(OMPNumThreads)
  double precision   :: multi_devStock(OMPNumThreads)
contains

  
  subroutine InitRandomSeed
    use constants, only : OMPNumThreads, myRank, PEtot
    use mt19937Mod
    implicit none
    include 'mpif.h'
    integer            :: i, seeds(OMPNumThreads)
    integer, parameter :: InitSeed1 = 3141921
    integer, parameter :: InitSeed2 = 655351
    integer, parameter :: IntMax    = 2147483647

    ! ------------------------------------------------- !
    ! --- [1] Initialize Internal Module Variables. --- !
    ! ------------------------------------------------- !
    flipflop          = 0
    devStock          = 0.d0
    multi_flipflop(:) = 0
    multi_devStock(:) = 0.d0

    ! ---------------------------------------------------------------------------------- !
    ! --- [2] Call sgrnd for MT19937 Routines. :: arbitrary determined value... seed --- !
    ! ---------------------------------------------------------------------------------- !
    !  -- [2-1] For Every PE -- !
    call sgrnd( ( IntMax - InitSeed1 ) / PEtot * myRank + InitSeed1 )
    !  -- [2-2] For Every Thread, in Each PE      -- !
    do i=1, OMPNumThreads
       seeds(i)       = int( grnd() * IntMax )
    enddo
    call multi_sgrnd( seeds )
    !  -- [2-3] Initialize PE's Single Seed Again -- !
    call sgrnd( InitSeed2 / PEtot * ( myRank+1 ) )

    ! ------------------- !
    ! --- [3] Message --- !
    ! ------------------- !
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '[InitRandomSeed               ] RandGenMod           :: ( Initialized )'
       write(6,'(2x,a)') '[InitRandomSeed               ] mt19937Mod           :: ( Initialized )'
    endif
    return
  end subroutine InitRandomSeed

  
  function gaussdev()
    use mt19937Mod    
    implicit none
    double precision :: gaussdev
    double precision :: v1, v2, r, fac

    if ( flipflop.eq.0 ) then
10     continue
       v1 = 2.d0 * grnd() - 1.d0
       v2 = 2.d0 * grnd() - 1.d0
       r  = v1**2 + v2**2
       if ( ( r.ge.1.d0 ).or.( r.eq.0.d0 ) ) goto 10
       fac      = sqrt( -2.d0 * log(r) / r )
       gaussdev = v1 * fac
       devStock = v2 * fac
       flipflop = 1       
    else
       gaussdev = devStock
       flipflop = 0
    endif
    
    return
  end function gaussdev

  
  function multi_gdev( myTHnum )
    use mt19937Mod
    implicit none
    integer, intent(in) :: myTHnum
    double precision    :: multi_gdev
    double precision    :: v1, v2, r, fac
    
    if ( multi_flipflop(myTHnum).eq.0 ) then
20      continue
       v1 = 2.d0 * multi_grnd( myTHnum ) - 1.d0
       v2 = 2.d0 * multi_grnd( myTHnum ) - 1.d0
       r  = v1**2 + v2**2
       if ( ( r.ge.1.d0 ).or.( r.eq.0.d0 ) ) goto 20
       fac = sqrt( -2.d0 * log(r) / r )
       multi_gdev              = v1 * fac
       multi_devStock(myTHnum) = v2 * fac
       multi_flipflop(myTHnum) = 1
    else
       multi_gdev              = multi_devStock(myTHnum)
       multi_flipflop(myTHnum) = 0       
    endif
    
    return
  end function multi_gdev

  
  subroutine Knuth_Shuffle_OMP( first, last )
    use constants , only  : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_, ns
    use variables , only  : pxv
    use mt19937Mod, only  : grnd, multi_grnd
    !$ use omp_lib
    implicit none
    integer, intent(in)  :: first, last
    integer              :: i, k, randpos, mythread
    double precision     :: r, tmp

    !$ mythread = omp_get_thread_num() + 1
    do k=1, ns
       do i=last, first+1, -1
          r       = multi_grnd( mythread )
          randpos = int( r*( i-first ) ) + first
          tmp     = pxv(rp_,randpos,k); pxv(rp_,randpos,k) = pxv(rp_,i,k) ; pxv(rp_,i,k) = tmp
          tmp     = pxv(zp_,randpos,k); pxv(zp_,randpos,k) = pxv(zp_,i,k) ; pxv(zp_,i,k) = tmp
          tmp     = pxv(vr_,randpos,k); pxv(vr_,randpos,k) = pxv(vr_,i,k) ; pxv(vr_,i,k) = tmp
          tmp     = pxv(vt_,randpos,k); pxv(vt_,randpos,k) = pxv(vt_,i,k) ; pxv(vt_,i,k) = tmp
          tmp     = pxv(vz_,randpos,k); pxv(vz_,randpos,k) = pxv(vz_,i,k) ; pxv(vz_,i,k) = tmp
          tmp     = pxv(ro_,randpos,k); pxv(ro_,randpos,k) = pxv(ro_,i,k) ; pxv(ro_,i,k) = tmp
          tmp     = pxv(zo_,randpos,k); pxv(zo_,randpos,k) = pxv(zo_,i,k) ; pxv(zo_,i,k) = tmp
          tmp     = pxv(wp_,randpos,k); pxv(wp_,randpos,k) = pxv(wp_,i,k) ; pxv(wp_,i,k) = tmp
       enddo
    enddo
    
    return
  end subroutine Knuth_Shuffle_OMP

end module randGenMod
