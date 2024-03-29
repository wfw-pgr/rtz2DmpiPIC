module mt19937Mod
  use constants, only : OMPNumThreads
  implicit none
  integer   , parameter :: N         = 624
  integer   , parameter :: M         = 397
  integer(8), parameter :: matA      = -1727483681_8
  integer(8), parameter :: upperMask = -2147483648_8
  integer(8), parameter :: lowerMask = +2147483647_8
  integer(8), parameter :: maskB     = -1658038656_8
  integer(8), parameter :: maskC     = -272236544_8
  integer(8), parameter :: maskI     = 1812433253_8
  integer   , parameter :: maskF     = -1
  integer   , parameter :: mask1     = +1
  integer   , parameter :: s_shift   =  7
  integer   , parameter :: t_shift   =  15
  integer   , parameter :: u_shift   = -11
  integer   , parameter :: f_shift   = -18
  integer               :: mti       = N+1
  integer               :: mt(0:N-1)
  integer(8)            :: mag01      (0:1) = (/0_8,matA/)
  integer               :: multi_mt (0:N-1,OMPNumThreads)
  integer               :: multi_mti(OMPNumThreads) = N+1
contains

  ! ====================================================== !
  ! === seed generation                                === !
  ! ====================================================== !
  subroutine sgrnd( seed )
    implicit none
    integer, intent(in) :: seed

    mt(0) = iand( seed, maskF )
    do mti=1, N-1
       mt(mti) = int(maskI) * ( ieor( mt(mti-1), ( ishft( mt(mti-1), -30) ) ) ) + mti
       mt(mti) = iand( mt(mti), maskF )
    enddo
    ! -----    old initialization routine    ----- !
    ! mt(0) = iand(seed, -1)
    ! do mti = 1, N - 1
    !   mt(mti) = iand(69069 * mt(mti - 1), -1)
    ! end do
    !
    ! .... both are available...
    ! -------------------------------------------- !
    return
  end subroutine sgrnd


  ! ====================================================== !
  ! === random number generation                       === !
  ! ====================================================== !
  function grnd()
    implicit none
    integer               :: y, kk
    integer, parameter    :: seed_default = 4357
    double precision      :: grnd

    if ( mti >= N ) then
       
       if ( mti == N+1 ) then
          call sgrnd( seed_default )
       endif
       
       do kk=0, N-M-1
          y      =  ior( iand( mt(kk)  , int(upperMask) ), iand( mt(kk+1), int(lowerMask) ) )
          mt(kk) = ieor( ieor( mt(kk+M), ishft(y,maskF) ),  int( mag01( iand( y,mask1 ) ) ) )
       enddo
       do kk=N-M, N-2
          y      =  ior( iand( mt(kk), int(upperMask)), iand( mt(kk+1), int(lowerMask) ) )
          mt(kk) = ieor( ieor( mt(kk+(M-N) ), ishft( y,maskF ) ), int( mag01( iand(y,mask1) ) ) )
       enddo

       y         =  ior( iand( mt(N-1),int(upperMask) ), iand( mt(0), int(lowerMask) ) )
       mt(N-1)   = ieor( ieor( mt(M-1), ishft(y,maskF) ), int( mag01( iand(y,mask1) ) ) )
       mti       = 0
    endif

    y   = mt(mti)
    mti = mti + 1
    y   = ieor( y,       ishft(y, u_shift) )
    y   = ieor( y, iand( ishft(y, s_shift), int(maskB) ) )
    y   = ieor( y, iand( ishft(y, t_shift), int(maskC) ) )
    y   = ieor( y,       ishft(y, f_shift) )

    if ( y < 0 ) then
       grnd = ( dble(y) + 2.0d0**32 ) / ( 2.0d0**32 )
    else
       grnd =   dble(y)               / ( 2.0d0**32 )
    endif

    return
  end function grnd


  !
  ! -- written by N.K.  -- !
  !
  ! ====================================================== !
  ! === multithreaded ver. of MT19937  ( setting )     === !
  ! ====================================================== !
  subroutine multi_sgrnd( seeds )
    !$ use omp_lib
    implicit none
    integer,intent(in) :: seeds(OMPNumThreads)
    integer            :: k, kmti
    
    ! ------------------------------------------------------ !
    ! --- [1] iterate grnd for each seed                 --- !
    ! ------------------------------------------------------ !
    do k=1, OMPNumThreads
       multi_mt(0,k) = iand( seeds(k), maskF )
       do kmti=1, N-1
          multi_mt(kmti,k) = int(maskI) * ( ieor( multi_mt(kmti-1,k), ( ishft( multi_mt(kmti-1,k), -30) ) ) ) + kmti
          multi_mt(kmti,k) = iand( multi_mt(kmti,k), maskF )
       enddo
       multi_mti(k) = kmti
    enddo

    return
  end subroutine multi_sgrnd

  
  !
  ! -- written by N.K.  -- !
  !  
  ! ====================================================== !
  ! === random number generation                       === !
  ! ====================================================== !
  function multi_grnd( myth )
    implicit none
    integer                :: kk
    integer, parameter     :: seed_default = 4357
    integer                :: seeds(OMPNumThreads)
    double precision       :: multi_grnd
    integer                :: y
    integer, intent(in)    :: myth

    ! ------------------------------------------------------ !
    ! --- [1] iterate grnd for every thread              --- !
    ! ------------------------------------------------------ !
    
    if ( multi_mti(myth) >= N ) then
       
       if ( multi_mti(myth) == N+1 ) then
          do kk=1, OMPNumThreads
             seeds(kk) = seed_default + kk
          enddo
          call multi_sgrnd( seeds )
       endif
       
       do kk=0, N-M-1
          y                 =  ior( iand( multi_mt(kk,myth)  , int(upperMask) ), iand( multi_mt(kk+1,myth), int(lowerMask) ) )
          multi_mt(kk,myth) = ieor( ieor( multi_mt(kk+M,myth), ishft(y,maskF) ), int( mag01( iand( y,mask1 ) ) ) )
       enddo
       do kk=N-M, N-2
          y                 =  ior( iand( multi_mt(kk,myth), int(upperMask)), iand( multi_mt(kk+1,myth), int(lowerMask) ) )
          multi_mt(kk,myth) = ieor( ieor( multi_mt(kk+(M-N),myth ), ishft( y,maskF ) ), int( mag01( iand(y,mask1) ) ) )
       enddo

       y                    =  ior( iand( multi_mt(N-1,myth),int(upperMask) ), iand( multi_mt(0,myth), int(lowerMask) ) )
       multi_mt(N-1,myth)   = ieor( ieor( multi_mt(M-1,myth), ishft(y,maskF) ), int( mag01( iand(y,mask1) ) ) )
       multi_mti(myth)      = 0
       
    endif

    y               = multi_mt( multi_mti(myth), myth )
    multi_mti(myth) = multi_mti(myth) + 1
    y               = ieor( y,       ishft(y, u_shift) )
    y               = ieor( y, iand( ishft(y, s_shift), int(maskB) ) )
    y               = ieor( y, iand( ishft(y, t_shift), int(maskC) ) )
    y               = ieor( y,       ishft(y, f_shift) )

    if ( y < 0 ) then
       multi_grnd = ( dble(y) + 2.0d0**32 ) / ( 2.0d0**32 )
    else
       multi_grnd =   dble(y)               / ( 2.0d0**32 )
    endif

    return
  end function multi_grnd

end module mt19937Mod
