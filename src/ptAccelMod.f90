module ptAccelMod
contains
  
  subroutine ptPushMove
    use constants ,   only       : ns, rMin, rMax
    use constants ,   only       : qm , np , dt , dr , dz , drInv, dzInv
    use constants ,   only       : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use constants ,   only       : er_, et_, ez_, br_, bt_, bz_
    use variables ,   only       : pxv, EBr, rh , rf
    use shapeFnMod,   only       : shapef_r, shapef_z
    use ptcBdrcMod,   only       : ReflecInXY
    use ptCheckMod,   only       : ptErrorCheck
    implicit none
    include 'mpif.h'
    integer                     :: i, j, ip, jp, k, m, cmp, ierr
    double precision            :: fac1, fac2, costh, sinth, chk, chkSum, gamma, vnrm2
    double precision            :: sfr(-2:2), sfz(-2:2), vtmp(2), psuedxy(3), qmhdt(2)
    double precision            :: EBp(6), Bts(6), vct(3:5)
    double precision            :: rInv, rr, sfrz, rRefl1, rRefl2
    double precision, parameter :: cSpeedLimit = 1.00d0

    ! ------------------------- !
    ! --- [1]  Preparation  --- !
    ! ------------------------- !
    rRefl1 = rMin + 0.5d0*sqrt(2.d0)*dr
    rRefl2 = rMax
    EBr    = 0.d0
    call Relocation( EBr )
    
    ! ------------------------- !
    ! --- [2]  Main  Loop   --- !
    ! ------------------------- !
    chk      = 0.d0
    qmhdt(:) = qm(:) * 0.5d0 * dt
    do k=1, ns  !  ns :: # of species
       !$omp parallel default( none ) &
       !$omp  shared( k,np,pxv,dr,dz,drInv,dzInv,EBr,dt,qmhdt,chk,rMin,rMax,rh,rf,rRefl1,rRefl2 ) &
       !$omp private( m,i,j,ip,jp,cmp,vnrm2,gamma,fac1,fac2,costh,sinth,sfr,sfz,EBp,Bts,vct,psuedxy) &
       !$omp private( vtmp,rInv,rr,sfrz )
       !$omp do reduction(+:chk)
       do m=1, np(k)

          ! -- (1) Grid -- !
          jp     = nint( pxv(rp_,m,k)*drInv )
          ip     = nint( pxv(zp_,m,k)*dzInv )
          sfr    = shapef_r( jp, pxv(rp_,m,k), drInv )
          sfz    = shapef_z( ip, pxv(zp_,m,k), dzInv )
          jp     = jp + 2
          ip     = ip + 2
          
          ! -- (2) Interpolation -- !
          rInv   = 1.d0 / max( ( rMin + pxv(rp_,m,k) ), 0.001d0*dr )
          EBp(:) = 0.d0
          do j=-1, 1
             rr  = rInv*abs( rf(jp+j) + rMin )
             do i=-1, 1
                sfrz     =   sfz(i) * sfr(j)
                EBp(er_) = EBp(er_) + sfrz*EBr(er_,ip+i,jp+j)*rr
                EBp(et_) = EBp(et_) + sfrz*EBr(et_,ip+i,jp+j)*rr
                EBp(ez_) = EBp(ez_) + sfrz*EBr(ez_,ip+i,jp+j)   ! *rr
                EBp(br_) = EBp(br_) + sfrz*EBr(br_,ip+i,jp+j)*rr
                EBp(bt_) = EBp(bt_) + sfrz*EBr(bt_,ip+i,jp+j)*rr
                EBp(bz_) = EBp(bz_) + sfrz*EBr(bz_,ip+i,jp+j)   ! *rr
             enddo
          enddo
          
          ! -- (3) Acceleration by E ( 1st-step ) -- !
          do cmp=1, 3
             EBp(cmp)  = qmhdt(k) * EBp(cmp)
          enddo
          pxv(vr_,m,k) = pxv(vr_,m,k) + EBp(er_)
          pxv(vt_,m,k) = pxv(vt_,m,k) + EBp(et_)
          pxv(vz_,m,k) = pxv(vz_,m,k) + EBp(ez_)
          
          ! -- (4) Rotation by B     ( 2nd-step ) -- !
          vnrm2        = pxv(vr_,m,k)**2 + pxv(vt_,m,k)**2 + pxv(vz_,m,k)**2
          chk          = chk  +  max( 0.d0,  ( vnrm2 - cSpeedLimit ) )
          gamma        = 1.d0 / sqrt( 1.d0 +   vnrm2 )
          fac1         = qmhdt(k) * gamma
          Bts(br_)     =     fac1 * EBp(br_)
          Bts(bt_)     =     fac1 * EBp(bt_)
          Bts(bz_)     =     fac1 * EBp(bz_)
          vct(vr_)     = pxv(vr_,m,k) + ( pxv(vt_,m,k)*Bts(bz_) - pxv(vz_,m,k)*Bts(bt_) )
          vct(vt_)     = pxv(vt_,m,k) + ( pxv(vz_,m,k)*Bts(br_) - pxv(vr_,m,k)*Bts(bz_) )
          vct(vz_)     = pxv(vz_,m,k) + ( pxv(vr_,m,k)*Bts(bt_) - pxv(vt_,m,k)*Bts(br_) )
          
          fac2         = 2.d0 / ( 1.d0 + ( Bts(br_)**2 + Bts(bt_)**2 + Bts(bz_)**2 ) )
          Bts(br_)     = fac2 * Bts(br_)
          Bts(bt_)     = fac2 * Bts(bt_)
          Bts(bz_)     = fac2 * Bts(bz_)

          pxv(vr_,m,k) = pxv(vr_,m,k) + ( vct(vt_)*Bts(bz_) - vct(vz_)*Bts(bt_) )
          pxv(vt_,m,k) = pxv(vt_,m,k) + ( vct(vz_)*Bts(br_) - vct(vr_)*Bts(bz_) )
          pxv(vz_,m,k) = pxv(vz_,m,k) + ( vct(vr_)*Bts(bt_) - vct(vt_)*Bts(br_) )
          
          ! -- (5) Acceleration by E ( 3rd-step ) -- !
          pxv(vr_,m,k) = pxv(vr_,m,k) + EBp(er_)
          pxv(vt_,m,k) = pxv(vt_,m,k) + EBp(et_)
          pxv(vz_,m,k) = pxv(vz_,m,k) + EBp(ez_)

          ! -- (6) Get psued Coordinate (xyz)     -- !
          vnrm2        = pxv(vr_,m,k)**2 + pxv(vt_,m,k)**2 + pxv(vz_,m,k)**2
          chk          = chk  +  max( 0.d0,  ( vnrm2 - cSpeedLimit ) )
          gamma        = dt   / sqrt( 1.d0 +   vnrm2 )
          psuedxy(1)   = pxv(rp_,m,k) + gamma*pxv(vr_,m,k) + rMin
          psuedxy(2)   =              + gamma*pxv(vt_,m,k)
          psuedxy(3)   = sqrt( psuedxy(1)**2 + psuedxy(2)**2 )
          
          ! -- (7) Copy to Old Position           -- !
          pxv(ro_,m,k) = pxv(rp_,m,k)
          pxv(zo_,m,k) = pxv(zp_,m,k)

          ! -- (8) Reflection :: boundary condition for RMax, RMin
          if ( psuedxy(3).gt.rRefl2 ) call ReflecInXY( psuedxy, pxv(1:8,m,k), rRefl2, +1.d0, rMin, gamma, dr, dz )
          ! if ( psuedxy(3) .lt. rRefl1 ) call ReflecInXY( psuedxy, pxv(:,m,k), rRefl1, -1.d0, rMin, gamma, dr, dz )
          pxv(rp_,m,k) = psuedxy(3) - rMin
          pxv(zp_,m,k) = pxv(zp_,m,k) + gamma*pxv(vz_,m,k)
          
          ! -- (8) calculate cos(a) & sin(a)      -- !
          psuedxy(3)   = 1.d0  / psuedxy(3)
          if ( pxv(rp_,m,k).ne.0.d0 ) then
             costh     = psuedxy(1) * psuedxy(3)
             sinth     = psuedxy(2) * psuedxy(3)
          else
             costh     = 1.d0
             sinth     = 0.d0
          endif

          ! -- (9) rotate velocity                -- !
          vtmp(1)      = pxv(vr_,m,k)
          vtmp(2)      = pxv(vt_,m,k)
          pxv(vr_,m,k) =   costh * vtmp(1) + sinth * vtmp(2)
          pxv(vt_,m,k) = - sinth * vtmp(1) + costh * vtmp(2)
          ! pxv(vt_,m,k) = vtmp(2) * ( pxv(ro_,m,k)+rMin ) * psuedxy(3) ! angular-momentum conservation

       enddo
       !$omp end do
       !$omp end parallel
    enddo

    !  ---  [3] Error Check ( gamma )  ---  !
    call MPI_REDUCE( chk, chkSum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    if ( chkSum.gt.0.d0 ) call ptErrorCheck( chkSum, cSpeedLimit, EBr )
    
    return
  end subroutine ptPushMove


  subroutine Relocation( EBr )
    use constants , only           : LIs, LJs
    use constants , only           : br_, bt_, bz_, er_, et_, ez_
    use variables , only           : rh, rfinv, EBf, EBo
    use pcwBdrcMod, only           : conductingWall_divB
    use pcwBdrcMod, only           : conductingWall_Epara
    use pcwBdrcMod, only           : conductingWall_EBr
    implicit none
    integer                       :: i, j
    double precision, intent(out) :: EBr(6,LIs+1,LJs+1)

    ! --- [1] Boundary Condition Before Relocation --- !
    ! call conductingWall_Epara
    call conductingWall_divB( EBf )
    
    ! --- [2] Relocation --- !
    ! Cylindrical consideration induced Relocation ::
    !   Er :: rEr => relocate on Vp :: electrical potential
    !   Bt :: rBt => relocate on I  :: current
    !$omp parallel default( none ) &
    !$omp shared( rh,rfinv,EBf,EBr,EBo,LIs,LJs ) private(i,j)
    !$omp do
    do j=2, LJs
       do i=2, LIs
          EBr(er_,i,j) = 0.50d0 *           ( EBf(er_,i,j  ) + EBf(er_,i-1,j  ) )
          EBr(et_,i,j) = 0.25d0 * ( rh(j  )*( EBf(et_,i,j  ) + EBf(et_,i-1,j  ) ) &
               &                  + rh(j-1)*( EBf(et_,i,j-1) + EBf(et_,i-1,j-1) ) )*rfinv(j)
          EBr(ez_,i,j) = 0.50d0 *           ( EBf(ez_,i,j  ) + EBf(ez_,i  ,j-1) )
          EBr(br_,i,j) = 0.25d0 * ( rh(j  )*( EBf(br_,i,j  ) + EBo(br_,i  ,j  ) ) &
               &                  + rh(j-1)*( EBf(br_,i,j-1) + EBo(br_,i  ,j-1) ) )*rfinv(j)
          EBr(bt_,i,j) = 0.50d0 *           ( EBf(bt_,i,j  ) + EBo(bt_,i  ,j  ) )
          EBr(bz_,i,j) = 0.25d0 * (         ( EBf(bz_,i,j  ) + EBo(bz_,i  ,j  ) ) &
               &                  +         ( EBf(bz_,i-1,j) + EBo(bz_,i-1,j  ) ) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! --- [3] Boundary Condition For EBr --- !
    call conductingWall_EBr
    
    return
  end subroutine Relocation


  subroutine position
    use constants   , only : ns, np, dt, dr, dz, rMin, rMax, Boundary1__pt
    use constants   , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_
    use variables   , only : pxv
    use ptcBdrcMod  , only : particleBoundary, ReflecInXY
    implicit none
    integer               :: m, k
    double precision      :: gamma, costh, sinth, rInv, rRefl1, rRefl2
    double precision      :: psuedxy(3), vtmp(2)
    logical, parameter    :: flag__Inner = .false.

    ! -- [1] Step Forward -- !
    rRefl1 = rMin + 0.5d0*sqrt( 2.d0 )*dr
    rRefl2 = rMax
    ! rRefl2 = rMax - 0.5d0*sqrt( 2.d0 )*dr
    do k=1, ns
       !$omp parallel default(none) &
       !$omp shared(k,np,dt,dr,dz,pxv,rMin,rMax,rRefl1,rRefl2) &
       !$omp private(m,gamma,costh,sinth,psuedxy,vtmp,rInv)
       !$omp do
       do m=1, np(k)
          ! -- (1) Copy Position in Old Register -- !
          pxv(ro_,m,k) = pxv(rp_,m,k)
          pxv(zo_,m,k) = pxv(zp_,m,k)

          ! -- (2) Step Forward dt in XY plain   -- !
          gamma        = dt / sqrt( 1.d0 + ( pxv(vr_,m,k)**2 + pxv(vt_,m,k)**2 + pxv(vz_,m,k)**2 ) )
          psuedxy(1)   = pxv(rp_,m,k) + gamma*pxv(vr_,m,k) + rMin
          psuedxy(2)   =              + gamma*pxv(vt_,m,k)
          psuedxy(3)   = sqrt( psuedxy(1)**2 + psuedxy(2)**2 )

          ! -- (3) Boundary :: Reflect @ r=R0    -- !
          ! if ( psuedxy(3).lt.rRefl1 ) call ReflecInXY( psuedxy, pxv(1:8,m,k), rRefl1, -1.d0, rMin, gamma, dr, dz )
          if ( psuedxy(3).gt.rRefl2 ) call ReflecInXY( psuedxy, pxv(1:8,m,k), rRefl2, +1.d0, rMin, gamma, dr, dz )

          ! -- (4) Update Position               -- !
          pxv(rp_,m,k) = psuedxy(3) - rMin
          pxv(zp_,m,k) = pxv(zp_,m,k) + gamma*pxv(vz_,m,k)

          ! -- (5) Rotate velocity               -- !
          if ( psuedxy(3).ne.0.d0 ) then
             rInv      = 1.d0 / psuedxy(3)
             costh     = psuedxy(1) * rInv
             sinth     = psuedxy(2) * rInv
          else
             costh     = 1.d0
             sinth     = 0.d0
          endif
          vtmp(1)      = pxv(vr_,m,k)
          vtmp(2)      = pxv(vt_,m,k)
          pxv(vr_,m,k) =   costh*vtmp(1) + sinth*vtmp(2)
          pxv(vt_,m,k) = - sinth*vtmp(1) + costh*vtmp(2)
       enddo
       !$omp end do
       !$omp end parallel
    enddo
    ! ! -- [2] z-Boundary -- !
    ! call particleBoundary
    
    return
  end subroutine position


end module ptAccelMod
