module pBoundaryMod
contains
  
  subroutine ptBoundary( ptBdrType )
    use constants , only       : ns, np, myRank, PEtot
    use constants , only       : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables , only       : pxv, x1Lengloc, x1Leng
    use pMPIComMod, only       : particleBoundaryCommunication
    implicit none
    integer                   :: m, k
    double precision          :: zRefl1, zRefl2
    character(10), intent(in) :: ptBdrType

    ! --- [1] Periodic Boundary - 1 --- !
    if ( trim(ptBdrType).eq.'periodic_1' ) then
       write(6,*) '[ERROR] NO periodic Boundary Condition in R-direction [ERROR]'
    endif
    ! --- [2] Periodic Boundary - 2 --- !
    if ( trim(ptBdrType).eq.'periodic_2' ) then
       if ( myRank.eq.0       ) then
          zRefl1 = 0.d0
          !$omp parallel default(none) &
          !$omp shared(np,pxv,zRefl1,zRefl2,x1Leng) private(m,k)
          !$omp do
          do k=1, ns
             do m=1, np(k)
                if ( pxv(zp_,m,k).le.zRefl1 ) pxv(zp_,m,k) = pxv(zp_,m,k) + x1Leng
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       if ( myRank.eq.PEtot-1 ) then
          zRefl2 = x1Lengloc
          !$omp parallel default(none) &
          !$omp shared(np,pxv,zRefl1,zRefl2,x1Leng) private(m,k)
          !$omp do
          do k=1, ns
             do m=1, np(k)
                if ( pxv(zp_,m,k).gt.zRefl2 ) pxv(zp_,m,k) = pxv(zp_,m,k) - x1Leng
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
    endif
    ! --- [3] Reflection Boundary - 2 --- !
    if ( trim(ptBdrType).eq.'reflect_2' ) then
       ! -- [3-1] Particle Communication -- !
       call ParticleBoundaryCommunication
       if ( myRank.eq.0       ) then
          zRefl1 = 0.d0
          !$omp parallel default(none) shared(np,pxv,zRefl1,zRefl2,myRank) private(m,k)
          !$omp do
          do k=1, ns
             do m=1, np(k)
                if ( pxv(zp_,m,k).le.zRefl1 ) then
                   pxv(zp_,m,k) = 2.d0*zRefl1 - pxv(zp_,m,k)
                   pxv(vz_,m,k) =             - pxv(vz_,m,k)
                   if ( pxv(zp_,m,k).lt.zRefl1 ) then
                      write(6,*) ' [ERROR] z < x1End1 @ (myRank,m,k)', myRank, m, k
                      write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
                   endif
                endif
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       if ( myRank.eq.PEtot-1 ) then 
          zRefl2 = x1Lengloc
          !$omp parallel default(none) shared(np,pxv,zRefl1,zRefl2,myRank) private(m,k)
          !$omp do
          do k=1, ns
             do m=1, np(k)
                if ( pxv(zp_,m,k).gt.zRefl2 ) then
                   pxv(zp_,m,k) = 2.d0*zRefl2 - pxv(zp_,m,k)
                   pxv(vz_,m,k) =             - pxv(vz_,m,k)
                   if ( pxv(zp_,m,k).gt.zRefl2 ) then
                      write(6,*) ' [ERROR] z > x1End2 @ (myRank,m,k)', myRank, m, k
                      write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
                   endif
                endif
             enddo
          enddo
          !$omp end do
          !$omp end parallel
       endif
       do k=1, ns
          do m=1, np(k)
             if ( pxv(zp_,m,k).gt.x1Lengloc ) then
                write(6,*) ' [ERROR] zp > x1End2 ', pxv(zp_,m,k), '  >  ', x1Lengloc
                write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
             endif
             if ( pxv(zp_,m,k).lt.0.d0      ) then
                write(6,*) ' [ERROR] zp < x1End1 ', pxv(zp_,m,k), '  >  ', 0.0
                write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
             endif
          enddo
       enddo
    end if
    ! --- [3] Reflection Boundary - 2 --- !
    if ( trim(ptBdrType).eq.'reflect_s' ) then
       ! -- [3-1] Particle Communication -- !
       zRefl1 = 0.d0
       zRefl2 = x1Lengloc
       !$omp parallel default(none) shared(np,pxv,zRefl1,zRefl2,myRank) private(m,k)
       !$omp do
       do k=1, ns
          do m=1, np(k)
             if ( pxv(zp_,m,k).le.zRefl1 ) then
                pxv(zp_,m,k) = 2.d0*zRefl1 - pxv(zp_,m,k)
                pxv(vz_,m,k) =             - pxv(vz_,m,k)
                if ( pxv(zp_,m,k).lt.zRefl1 ) then
                   write(6,*) ' [ERROR] z < x1End1 @ (myRank,m,k)', myRank, m, k
                   write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
                endif
             endif
             if ( pxv(zp_,m,k).gt.zRefl2 ) then
                pxv(zp_,m,k) = 2.d0*zRefl2 - pxv(zp_,m,k)
                pxv(vz_,m,k) =             - pxv(vz_,m,k)
                if ( pxv(zp_,m,k).gt.zRefl2 ) then
                   write(6,*) ' [ERROR] z > x1End2 @ (myRank,m,k)', myRank, m, k
                   write(6,'(8(e12.5,1x))') pxv(1:8,m,k)
                endif
             endif
          enddo
       enddo
       !$omp end do
       !$omp end parallel
    end if
    return
  end subroutine ptBoundary


  subroutine ReflecInXY( psuedxy, pxvrt, R0, sgn, rmin, gamma, dr, dz )
    implicit none
    double precision                :: vnrm2, vcrsx, vdotx, dtprm, theta, phi, costh, sinth
    double precision                :: vtmp(2), hdelta(2)
    double precision, intent(in)    :: R0, sgn, rmin, gamma, dr, dz
    double precision, intent(inout) :: psuedxy(3), pxvrt(8)
    
    !  --- [1] Delta t ' ---  ! :: dt' == > Actually including Gamma :: == Gamma dt'
    psuedxy(1) = pxvrt(1) + rmin
    psuedxy(2) = 0.d0
    vnrm2      = pxvrt(3)**2 + pxvrt(4)**2
    vcrsx      = pxvrt(3)*psuedxy(2) - pxvrt(4)*psuedxy(1)
    vdotx      = pxvrt(3)*psuedxy(1) + pxvrt(4)*psuedxy(2)
    dtprm      = (  - vdotx + sgn * sqrt( R0**2 * vnrm2 - vcrsx**2 ) ) / vnrm2
    
    !  --- [2] P(x,y), theta, phi --- !
    !   -- [2-1] update P(X,Y) --  !
    psuedxy(1) = psuedxy(1) + pxvrt(3) * dtprm
    psuedxy(2) =            + pxvrt(4) * dtprm
    !   -- [2-2] hDelta --  !
    hdelta(1)  = sign( 0.5d0*dz, psuedxy(1) )
    hdelta(2)  = sign( 0.5d0*dr, psuedxy(2) )
    psuedxy(1) = psuedxy(1) + hdelta(1)
    psuedxy(2) = psuedxy(2) + hdelta(2)
    !   -- [2-3] Angle Theta   --  !
    theta      = atan( psuedxy(2) / psuedxy(1) )
    phi        = atan( pxvrt(3)   / pxvrt(4)   )
    costh      = cos(  2.d0 * ( theta + phi )  )
    sinth      = sin(  2.d0 * ( theta + phi )  )

    !  --- [3] Rotation of vx, vy --- !
    vtmp(1)    = pxvrt(3)
    vtmp(2)    = pxvrt(4)
    pxvrt(3)   = costh*vtmp(1) - sinth*vtmp(2)
    pxvrt(4)   = sinth*vtmp(1) + costh*vtmp(2)

    !  --- [4] step forward : dt-dtprm  --- !
    psuedxy(1) = psuedxy(1) + pxvrt(3)*( gamma - dtprm ) - hdelta(1)
    psuedxy(2) = psuedxy(2) + pxvrt(4)*( gamma - dtprm ) - hdelta(2)
    psuedxy(3) = sqrt( psuedxy(1)**2 + psuedxy(2)**2 )

    return
  end subroutine ReflecInXY


  subroutine ptErrorCheck( chk, cSpeedLimit, EBr )
    use constants, only           : LI, LJ, np, ns, dzInv, drInv, rmin, myRank
    use constants, only           : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use constants, only           : er_, et_, ez_, br_, bt_, bz_ 
    use variables, only           : kstep, pxv, rf
    use shapefMod, only           : shapef_r, shapef_z
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
    double precision, parameter  :: alpha        = 0.1d0

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

end module pBoundaryMod
