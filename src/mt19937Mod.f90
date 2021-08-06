! MT19937, Mersenne Twister Random Number Generator ([0,1) Real Number)
! in Fortran90
!
! Usage: 
!   1) When you use mt19937, add the sentence "use mt19937" 
!      above the implicit sentence.
!   2) To set an initial seed, call sgrnd(seed). (The "seed" is an integer.)
!      If you do not call this, the seed is 4357.
!   3) Use the function grnd().
!      (Do not declare "real(8) :: grnd".)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! A C-program for MT19937: Real number version
!   genrand() generates one pseudorandom real number (double)
! which is uniformly distributed on [0,1]-interval, for each
! call. sgenrand(seed) set initial values to the working area
! of 624 words. Before genrand(), sgenrand(seed) must be
! called once. (seed is any 32-bit integer except for 0).
! Integer generator is obtained by modifying two lines.
!   Coded by Takuji Nishimura, considering the suggestions by
! Topher Cooper and Marc Rieffel in July-Aug. 1997.
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(seed)
!                     integer seed
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fortran90 translation by Yusuke Konishi. May. 13, 2013.
!
! The interval is changed from [0, 1] to [0, 1).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





module mt19937Mod

  use constants, only : OMPNumThreads
  
  implicit none
  integer, parameter         :: lk = selected_int_kind (16)
  integer (kind=lk), private :: N, N1, M, MATA, UMASK, LMASK, TMASKB, TMASKC
  parameter( &
             & N = 624_lk, &
             & N1 = 625_lk, &
             & M = 397_lk, &
             & MATA = -1727483681_lk, &
             & UMASK = -2147483648_lk, &
             & LMASK = 2147483647_lk, &
             & TMASKB = -1658038656_lk, &
             & TMASKC = -272236544_lk &
             & )
  integer,private :: mti = N1, mt(0:N-1), mag01(0:1) = (/0_lk, MATA/)

  ! integer,private, parameter :: Nthread = 12
  integer,private :: multi_mt(0:N-1,OMPNumThreads), multi_mti(OMPNumThreads)

  contains

  subroutine sgrnd(seed)
    integer,intent(in) :: seed
!
!      setting initial seeds to mt[N] using
!      the generator Line 25 of Table 1 in
!      [KNUTH 1981, The Art of Computer Programming
!         Vol. 2 (2nd Ed.), pp102]
!
    mt(0) = iand(seed, -1)
    do mti = 1, N - 1
      mt(mti) = iand(69069 * mt(mti - 1), -1)
    end do
  end subroutine sgrnd

  real(8) function grnd()
    integer :: y, kk

    if(mti >= N) then
!                   generate N words at one time
      if(mti == N + 1) then
!                        if sgrnd() has not been called,
        call sgrnd(4357)
!                          a default initial seed is used
      endif

      do kk = 0, N - M - 1
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + M), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      do kk = N - M, N - 2
        y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
        mt(kk) = ieor(ieor(mt(kk + (M - N)), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      y = ior(iand(mt(N - 1), UMASK), iand(mt(0), LMASK))
      mt(N - 1) = ieor(ieor(mt(M - 1), ishft(y, -1)), mag01(iand(y, 1)))
      mti = 0
    endif

    y = mt(mti)
    mti = mti + 1
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7), TMASKB))
    y = ieor(y, iand(ishft(y, 15), TMASKC))
    y = ieor(y, ishft(y, -18))

    if(y < 0) then
      grnd = (dble(y) + 2.0d0 ** 32) / (2.0d0 ** 32)
    else
      grnd = dble(y) / (2.0d0 ** 32)
    endif
  end function grnd



  ! written by kent
  subroutine multi_sgrnd(seeds)

    !$ use omp_lib
    integer,intent(in) :: seeds(OMPNumThreads)
    integer :: k
    !
    !      setting initial seeds to mt[N] using
    !      the generator Line 25 of Table 1 in
    !      [KNUTH 1981, The Art of Computer Programming
    !         Vol. 2 (2nd Ed.), pp102]
    !
    do k=1, OMPNumThreads
       multi_mt(0,k) = iand(seeds(k), -1)
       mti=multi_mti(k)
       do mti = 1, N - 1
          multi_mt( mti,k ) = iand( 69069 * multi_mt( mti-1, k ), -1)
       end do
       multi_mti(k) = mti
    end do
  end subroutine multi_sgrnd


  
  real(8) function multi_grnd( mythread )
    integer :: y, kk
    integer, intent(in) :: mythread
    
    if(multi_mti(mythread) >= N) then
!                   generate N words at one time
      if(multi_mti(mythread) == N + 1) then
!                        if sgrnd() has not been called,
         stop ' use multi_sgrnd function stop at mt '
!                          a default initial seed is used
      endif

      do kk = 0, N - M - 1
        y = ior(iand(multi_mt(kk,mythread), UMASK), iand(multi_mt(kk + 1,mythread), LMASK))
        multi_mt(kk,mythread) = ieor(ieor(multi_mt(kk + M,mythread), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      do kk = N - M, N - 2
        y = ior(iand(multi_mt(kk,mythread), UMASK), iand(multi_mt(kk + 1,mythread), LMASK))
        multi_mt(kk,mythread) = ieor(ieor(multi_mt(kk + (M - N),mythread), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      y = ior(iand(multi_mt(N - 1,mythread), UMASK), iand(multi_mt(0,mythread), LMASK))
      multi_mt(N - 1,mythread) = ieor(ieor(multi_mt(M - 1,mythread), ishft(y, -1)), mag01(iand(y, 1)))
      multi_mti(mythread) = 0
    endif

    y = multi_mt(multi_mti(mythread),mythread)
    multi_mti(mythread) = multi_mti(mythread) + 1
    y = ieor(y, ishft(y, -11))
    y = ieor(y, iand(ishft(y, 7), TMASKB))
    y = ieor(y, iand(ishft(y, 15), TMASKC))
    y = ieor(y, ishft(y, -18))

    if(y < 0) then
      multi_grnd = (dble(y) + 2.0d0 ** 32) / (2.0d0 ** 32)
    else
      multi_grnd = dble(y) / (2.0d0 ** 32)
    endif
  end function multi_grnd

end module mt19937Mod
