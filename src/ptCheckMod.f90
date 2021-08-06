module ptCheckMod
contains

  subroutine ptBoundaryOutCheck
    use constants , only       : ns, np, myRank, zp_
    use variables , only       : pxv, x1Lengloc
    implicit none
    integer                   :: m, k
    double precision          :: x1End1, x1End2

    x1End1 = 0.d0
    x1End2 = x1Lengloc
    !$omp parallel default(none) &
    !$omp shared(np,pxv,x1End1,x1End2,myRank) private(m,k)
    !$omp do
    do k=1, ns
       do m=1, np(k)
          if ( pxv(zp_,m,k).lt.x1End1 ) then
             write(6,'(a,i6,2(a4,e12.5))') ' [ERROR] zp < x1End1 @ rank = ', myRank, ' :: ', &
                  &                           pxv(zp_,m,k), ' > ', x1End1
             write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
          endif
          if ( pxv(zp_,m,k).gt.x1End2 ) then
             write(6,'(a,i6,2(a4,e12.5))') ' [ERROR] zp > x1End2 @ rank = ', myRank, ' :: ', &
                  &                           pxv(zp_,m,k), ' > ', x1End2
             write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
          endif
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    return
  end subroutine ptBoundaryOutCheck


  subroutine ptErrorCheck( chk, cSpeedLimit, EBr )
    use constants, only           : LI, LJ, np, ns, dzInv, drInv, rmin, myRank
    use constants, only           : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use constants, only           : er_, et_, ez_, br_, bt_, bz_ 
    use variables, only           : kstep, pxv, rf
    use shapeFnMod, only           : shapef_r, shapef_z
    implicit none
    include 'mpif.h'
    integer                      :: k, m, i, j, ip, jp
    double precision             :: rr, rInv, EBp(6), sfr(-2:2), sfz(-2:2), sfrz
    double precision, intent(in) :: chk, cSpeedLimit, EBr(6,LI+1,LJ+1)
    character(len=100)           :: filename     = 'Error.log'
    logical         , parameter  :: flag__stop   = .false.
    logical         , parameter  :: flag__remove = .false.
    logical         , parameter  :: flag__report = .false.
    logical         , parameter  :: flag__chkOut = .true.
    double precision, parameter  :: alpha        =  0.1d0

    if ( flag__report ) then
       open (50,file=trim(filename),form='formatted',access='append')
       write(50,'(a,i10,a,e12.5,a)') '[ERROR] Particle ERROR at kstep == ', kstep, ' check == ', chk, ' [ERROR]'
       do k=1, ns
          do m=1, np(k)
             if ( ( pxv(vr_,m,k)**2 + pxv(vt_,m,k)**2 + pxv(vz_,m,k)**2 ).gt.cSpeedLimit ) then
                ip     = nint( pxv(zp_,m,k)*dzInv )
                jp     = nint( pxv(rp_,m,k)*drInv )
                sfr    = shapef_r( jp, pxv(rp_,m,k), drInv )
                sfz    = shapef_z( ip, pxv(zp_,m,k), dzInv )
                jp     = jp + 2
                ip     = ip + 2
                rInv   = 1.d0 / ( rmin + pxv(rp_,m,k) )
                EBp(:) = 0.d0
                do j=-1, 1
                   rr = rInv*abs(rf(jp+j))
                   do i=-1, 1
                      sfrz     = sfz(i) * sfr(j)
                      EBp(er_) = EBp(er_) + sfrz*EBr(1,ip+i,jp+j)*rr
                      EBp(et_) = EBp(et_) + sfrz*EBr(2,ip+i,jp+j)
                      EBp(ez_) = EBp(ez_) + sfrz*EBr(3,ip+i,jp+j)
                      EBp(br_) = EBp(br_) + sfrz*EBr(4,ip+i,jp+j)*rr
                      EBp(bt_) = EBp(bt_) + sfrz*EBr(5,ip+i,jp+j)
                      EBp(bz_) = EBp(bz_) + sfrz*EBr(6,ip+i,jp+j)
                   enddo
                enddo
                write(50,'(3x,a,1x,2(i8   ,1x))') '(ip,jp)     = ', ip, jp
                write(50,'(3x,a,1x,2(e10.3,1x))') ' (r,z)      = ', pxv(1:2,m,k)
                write(50,'(3x,a,1x,3(e10.3,1x))') '(vr,vt,vz)  = ', pxv(3:5,m,k)
                write(50,'(3x,a,1x,3(e10.3,1x))') '(ro,zo,wt)  = ', pxv(6:8,m,k)
                write(50,'(3x,a,1x,3(e10.3,1x))') '(Er,Et,Ez)p = ', EBp(1:3)
                write(50,'(3x,a,1x,3(e10.3,1x))') '(Br,Bt,Bz)p = ', EBp(4:6)
                if ( flag__remove ) then
                   pxv(vr_,m,k) = pxv(vr_,m,k) * alpha
                   pxv(vt_,m,k) = pxv(vt_,m,k) * alpha
                   pxv(vz_,m,k) = pxv(vz_,m,k) * alpha
                endif
             endif
          enddo
       enddo
       close(50)
    endif
    if ( flag__remove ) then
       !$omp parallel default(none) &
       !$omp shared(pxv,np,cSpeedLimit) private(m,k)
       !$omp do
       do k=1, ns
          do m=1, np(k)
             if ( ( pxv(vr_,m,k)**2+pxv(vt_,m,k)**2+pxv(vz_,m,k)**2 ).gt.cSpeedLimit ) then
                pxv(vr_,m,k) = pxv(vr_,m,k) * alpha
                pxv(vt_,m,k) = pxv(vt_,m,k) * alpha
                pxv(vz_,m,k) = pxv(vz_,m,k) * alpha
             endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    endif
    if ( flag__chkOut ) then
       if ( myRank.eq.0 ) then
          open (50,file=trim(filename),form='formatted',access='append')
          write(50,*) kstep, chk
          close(50)
       endif
    endif
    if ( flag__stop   ) stop ' [ERROR] Velocity overcome light speed ( Big Gamma ) [ERROR] '
    return
  end subroutine ptErrorCheck


end module ptCheckMod
