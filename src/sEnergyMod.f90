module sEnergyMod
contains

  subroutine CalcEnergy
    use constants , only : LIs, LJs, npt, jobdir, myRank
    use constants , only : np, rm
    use constants , only : dt, vAlfe, pi
    use constants , only : er_, et_, ez_, br_, bt_, bz_
    use constants , only : rp_, zp_, vr_, vt_, vz_, ro_, zo_, wp_, el_, io_
    use variables , only : pxv, EBf, EBo, neMax, neMin
    use variables , only : RhoPerNsp, kstep, volRepres, gVh, gVf
    use sFluxIpMod, only : Iplasma, psiMin, mergR, RateR
    implicit none
    include 'mpif.h'
    integer             :: i, j, m, ierr
    double precision    :: gamma, factor, wcendt
    double precision    :: Utot, Umxw
    double precision    :: UPK, UMF, UEF
    double precision    :: UPK_e, UPK_i 
    double precision    :: Ukex, Ukey, Ukez, Ukix, Ukiy, Ukiz
    double precision    :: UEFx, UEFy, UEFz, UMFx, UMFy, UMFz
    double precision    :: ExMid, EyMid, EzMid
    double precision    :: sendUbuff(3,4), recvUbuff(3,4)
    character(100)      :: FileName
    integer, parameter  :: ke_ = 1, ki_ = 2, ef_ = 3, bf_ = 4
    
    !  --- [0]  preparation    ---  !
    wcendt      = dble( kstep ) * dt
    
    !  ---------------------------  !
    !  --- [1]  kinetic energy ---  !    
    !  ---------------------------  !
    !   -- [1-1] Electron      --   !
    Ukex      = 0.d0
    Ukey      = 0.d0
    Ukez      = 0.d0
    !$omp parallel default(none) &
    !$omp shared(np,pxv,Ukex,Ukey,Ukez) private(m,gamma) 
    !$omp do reduction(+:Ukex,Ukey,Ukez)
    do m=1, np(el_)
       gamma  = 1.d0 / ( 1.d0 + sqrt( 1.d0 + ( pxv(vr_,m,el_)**2 + pxv(vt_,m,el_)**2 + pxv(vz_,m,el_)**2 ) ) )
       Ukex   = Ukex + pxv(vr_,m,el_)**2 * gamma * pxv(wp_,m,el_)
       Ukey   = Ukey + pxv(vt_,m,el_)**2 * gamma * pxv(wp_,m,el_)
       Ukez   = Ukez + pxv(vz_,m,el_)**2 * gamma * pxv(wp_,m,el_)
    enddo
    !$omp end do
    !$omp end parallel

    !   -- [1-2] Ion           --   !
    Ukix        = 0.d0
    Ukiy        = 0.d0
    Ukiz        = 0.d0
    !$omp parallel default(none) &
    !$omp shared(np,npt,pxv,Ukix,Ukiy,Ukiz) private(m,gamma) 
    !$omp do reduction(+:Ukix,Ukiy,Ukiz)
    do m=1, np(io_)
       gamma    = 1.d0 / ( 1.d0 + sqrt( 1.d0 + ( pxv(vr_,m,io_)**2 + pxv(vt_,m,io_)**2 + pxv(vz_,m,io_)**2 ) ) )
       Ukix     = Ukix + pxv(vr_,m,io_)**2 * gamma * pxv(wp_,m,io_)
       Ukiy     = Ukiy + pxv(vt_,m,io_)**2 * gamma * pxv(wp_,m,io_)
       Ukiz     = Ukiz + pxv(vz_,m,io_)**2 * gamma * pxv(wp_,m,io_)
    enddo
    !$omp end do
    !$omp end parallel

    !  -----------------------------------  !
    !  --- [2]  Electric Field Energy  ---  !
    !  -----------------------------------  !
    UEFx        = 0.d0
    UEFy        = 0.d0
    UEFz        = 0.d0
    !$omp parallel default(none) &
    !$omp shared(EBf,EBo,gVh,gVf,UEFx,UEFy,UEFz,LIs,LJs) private(i,j,ExMid,EyMid,EzMid) 
    !$omp do reduction(+:UEFx,UEFy,UEFz)
    do j=2, LJs
       do i=2, LIs-1
          ExMid = 0.5d0 * ( EBf(er_,i,j) + EBo(er_,i,j) )
          EyMid = 0.5d0 * ( EBf(et_,i,j) + EBo(et_,i,j) )
          EzMid = 0.5d0 * ( EBf(ez_,i,j) + EBo(ez_,i,j) )
          UEFx  = UEFx  + 0.5d0 * gVf(j) * ExMid**2
          UEFy  = UEFy  + 0.5d0 * gVh(j) * EyMid**2
          UEFz  = UEFz  + 0.5d0 * gVh(j) * EzMid**2 
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !  -----------------------------------  !
    !  --- [3] Magnetic Field Energy   ---  !
    !  -----------------------------------  !
    UMFx        = 0.d0
    UMFy        = 0.d0
    UMFz        = 0.d0
    !$omp parallel default(none) &
    !$omp shared(EBf,gVh,gVf,UMFx,UMFy,UMFz,LIs,LJs) private(i,j)
    !$omp do reduction(+:UMFx,UMFy,UMFz)
    do j=2, LJs
       do i=2, LIs-1
          UMFx  = UMFx + 0.5d0 * gVh(j) * EBf(br_,i,j)**2
          UMFy  = UMFy + 0.5d0 * gVf(j) * EBf(bt_,i,j)**2
          UMFz  = UMFz + 0.5d0 * gVf(j) * EBf(bz_,i,j)**2
       enddo
    enddo
    !$omp end do
    !$omp end parallel    
    
    !  -------------------------------------------  !
    !  --- [4]  Communication / Total Energy   ---  !
    !  -------------------------------------------  !
    !   -- [4-1]    Packing into buffer         --  !
    !    - [4-1-1]  Kinetic Energy ( electron ) -   !
    factor           = rm(1) * RhoPerNsp / ( vAlfe**2 ) * volRepres
    sendUbuff(1,ke_) = Ukex * factor
    sendUbuff(2,ke_) = Ukey * factor
    sendUbuff(3,ke_) = Ukez * factor
    !    - [4-1-2]  Kinetic Energy (   ion    ) -   !
    factor           = rm(2) * RhoPerNsp / ( vAlfe**2 ) * volRepres
    sendUbuff(1,ki_) = Ukix * factor
    sendUbuff(2,ki_) = Ukiy * factor
    sendUbuff(3,ki_) = Ukiz * factor
    !    - [4-1-3]  Electric Field Energy       -   !
    sendUbuff(1,ef_) = UEFx
    sendUbuff(2,ef_) = UEFy
    sendUbuff(3,ef_) = UEFz
    !    - [4-1-4]  Magnetic Field Energy       -   !
    sendUbuff(1,bf_) = UMFx
    sendUbuff(2,bf_) = UMFy
    sendUbuff(3,bf_) = UMFz
    
    !   -- [4-2]  Communication [Reduction]     --  !
    call MPI_AllReduce( sendUbuff, recvUbuff,  12, MPI_DOUBLE_PRECISION, &
         &              MPI_SUM  , MPI_COMM_WORLD, ierr                  )    

    !   -- [4-3]  Total Energy                  --  !
    !    - [Kinetic] -    !
    UPK_e            = recvUbuff(1,ke_) + recvUbuff(2,ke_) + recvUbuff(3,ke_)
    UPK_i            = recvUbuff(1,ki_) + recvUbuff(2,ki_) + recvUbuff(3,ki_)
    UPK              = UPK_e + UPK_i
    !    - [ Field ] -    !
    UEF              = recvUbuff(1,ef_) + recvUbuff(2,ef_) + recvUbuff(3,ef_)
    UMF              = recvUbuff(1,bf_) + recvUbuff(2,bf_) + recvUbuff(3,bf_)
    !    - [ Total ] -    !
    Umxw             =        UEF  + UMF
    Utot             = UPK  + UEF  + UMF
    
    !  --------------------------  !
    !  --- [5] Display Energy ---  !
    !  --------------------------  !
    if ( myRank.eq.0 ) then
       !    -- [5-1] Standard Output -- !
       write(6,'(2x,a,i15,1x,a,e15.8,1x)') 'kstep :: ', kstep,  ' wcet:: ', wcendt
       write(6,'(2x,3(a,e15.8,1x))'      ) '  Ek  :: ', UPK  ,  ' Ee  :: ', UEF   , ' Eb  :: ', UMF
       write(6,'(2x,2(a,e15.8,1x))'      ) '  Emxw:: ', Umxw ,  ' Etot:: ', Utot
       write(6,'(2x,2(a,e15.8,1x))'      ) '  Eele:: ', UPK_e,  ' Eion:: ', UPK_i
       write(6,'(2x,2(a,e15.8,1x))'      ) '  nemx:: ', neMax,  ' nemn:: ', neMin
       write(6,'(2x,2(a,e15.8,1x))'      ) '  Ip  :: ', Iplasma,' psiM:: ', psiMin
       write(6,'(2x,2(a,e15.8,1x))'      ) '  mrgR:: ', mergR,  ' recR:: ', RateR
       write(6,*)
       !    -- [5-2] write into File  --  !
       !     - [5-2-1] Total   Energy - !
       FileName = trim(jobdir) // 'dat/' // 'energy.dat'
       open( 50, file=trim(FileName), form='formatted', position='append' )
       write(50,'((i8,1x),5(e15.7,1x))') kstep, wcendt, UPK, UEF, UMF, Utot
       close(50)
       !     - [5-2-2] Kinetic Energy - !
       FileName = trim(jobdir) // 'dat/' // 'energy_kinetic.dat'
       open( 50, file=trim(FileName), form='formatted', position='append' )
       write(50,'((i8,1x),10(e15.7,1x))') kstep, wcendt, UPK , &
            &                             recvUbuff(1,ke_), recvUbuff(2,ke_), recvUbuff(3,ke_), UPK_e, &
            &                             recvUbuff(1,ki_), recvUbuff(2,ki_), recvUbuff(3,ki_), UPK_i 
       close(50)
       !     - [5-2-3] Field's Energy - !
       FileName = trim(jobdir) // 'dat/' // 'energy_field.dat'
       open( 50, file=trim(FileName), form='formatted', position='append' )
       write(50,'((i8,1x),10(e15.7,1x))') kstep, wcendt, Umxw,      &
            &                             recvUbuff(1,ef_), recvUbuff(2,ef_), recvUbuff(3,ef_), UEF, &
            &                             recvUbuff(1,bf_), recvUbuff(2,bf_), recvUbuff(3,bf_), UMF 
       close(50)       
    endif

    return
  end subroutine CalcEnergy

  

  subroutine Flow_Thermal_Energy
    use constants, only : LIs, LJs, rm, jobDir, dt, valfe, myRank
    use variables, only : kstep, gVf, gVh
    use momentsMod, only : fMoments
    use momentsMod, only : nei_, ur_, ut_, uz_, pxx_, pyy_, pzz_
    implicit none
    include 'mpif.h'
    integer            :: i, j, ierr
    double precision   :: wcendt, Ukt, valfe2Inv
    double precision   :: Uuer, Uuet, Uuez, Uper, Upet, Upez, Uue, Upe, Uke
    double precision   :: Uuir, Uuit, Uuiz, Upir, Upit, Upiz, Uui, Upi, Uki
    double precision   :: sendUbuff(3,4), recvUbuff(3,4)
    character(100)     :: FileName
    integer, parameter :: ue_ = 1, ui_ = 2, pe_ = 3, pi_ = 4

    ! -------------------------- !
    ! --- [0] Preparation    --- !
    ! -------------------------- !
    wcendt    = dble( kstep ) * dt
    valfe2Inv = 1.d0 / ( valfe**2 )
    
    ! -------------------------- !
    ! --- [1] Flow Energy    --- !
    ! -------------------------- !
    Uuer   = 0.d0
    Uuet   = 0.d0
    Uuez   = 0.d0
    Uuir   = 0.d0
    Uuit   = 0.d0
    Uuiz   = 0.d0
    !$omp parallel default(none) &
    !$omp shared(LIs,LJs,fMoments,gVf,gVh) &
    !$omp shared(Uuer,Uuet,Uuez,Uuir,Uuit,Uuiz) &
    !$omp private(i,j)
    !$omp do reduction(+:Uuer,Uuet,Uuez,Uuir,Uuit,Uuiz)
    do j=2, LJs
       do i=2, LIs
          Uuer = Uuer + fMoments(i,j,1,nei_)*fMoments(i,j,1,ur_ )**2 *gVf(j)
          Uuet = Uuet + fMoments(i,j,1,nei_)*fMoments(i,j,1,ut_ )**2 *gVh(j)
          Uuez = Uuez + fMoments(i,j,1,nei_)*fMoments(i,j,1,uz_ )**2 *gVh(j)
          Uuir = Uuir + fMoments(i,j,2,nei_)*fMoments(i,j,2,ur_ )**2 *gVf(j)
          Uuit = Uuit + fMoments(i,j,2,nei_)*fMoments(i,j,2,ut_ )**2 *gVh(j)
          Uuiz = Uuiz + fMoments(i,j,2,nei_)*fMoments(i,j,2,uz_ )**2 *gVh(j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    ! -------------------------- !
    ! --- [2] Thermal Energy --- !
    ! -------------------------- !
    Uper   = 0.d0
    Upet   = 0.d0
    Upez   = 0.d0
    Upir   = 0.d0
    Upit   = 0.d0
    Upiz   = 0.d0
    !$omp parallel default(none) &
    !$omp shared(LIs,LJs,fMoments,gVf,gVh) &
    !$omp shared(Uper,Upet,Upez,Upir,Upit,Upiz) &
    !$omp private(i,j)
    !$omp do reduction(+:Uper,Upet,Upez,Upir,Upit,Upiz)
    do j=2, LJs
       do i=2, LIs
          Uper  = Uper + fMoments(i,j,1,pxx_)*gVf(j)
          Upet  = Upet + fMoments(i,j,1,pyy_)*gVh(j)
          Upez  = Upez + fMoments(i,j,1,pzz_)*gVh(j)
          Upir  = Upir + fMoments(i,j,2,pxx_)*gVf(j)
          Upit  = Upit + fMoments(i,j,2,pyy_)*gVh(j)
          Upiz  = Upiz + fMoments(i,j,2,pzz_)*gVh(j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    !  -------------------------------------------  !
    !  --- [3]  Communication / Total Energy   ---  !
    !  -------------------------------------------  !
    !   -- [3-1]    Packing into buffer         --  !
    !    - Flow    (electron) -    !
    sendUbuff(1,ue_) = Uuer * 0.5d0*rm(1)*valfe2Inv
    sendUbuff(2,ue_) = Uuet * 0.5d0*rm(1)*valfe2Inv
    sendUbuff(3,ue_) = Uuez * 0.5d0*rm(1)*valfe2Inv
    !    - Flow    (ion)      -    !
    sendUbuff(1,ui_) = Uuir * 0.5d0*rm(2)*valfe2Inv
    sendUbuff(2,ui_) = Uuit * 0.5d0*rm(2)*valfe2Inv
    sendUbuff(3,ui_) = Uuiz * 0.5d0*rm(2)*valfe2Inv
    !    - Thermal (electron) -    !
    sendUbuff(1,pe_) = Uper * 0.5d0*rm(1)*valfe2Inv
    sendUbuff(2,pe_) = Upet * 0.5d0*rm(1)*valfe2Inv
    sendUbuff(3,pe_) = Upez * 0.5d0*rm(1)*valfe2Inv
    !    - Thermal (ion)      -    !
    sendUbuff(1,pi_) = Upir * 0.5d0*rm(2)*valfe2Inv
    sendUbuff(2,pi_) = Upit * 0.5d0*rm(2)*valfe2Inv
    sendUbuff(3,pi_) = Upiz * 0.5d0*rm(2)*valfe2Inv

    !   -- [3-2]  Communication [Reduction]     --  !
    call MPI_AllReduce( sendUbuff, recvUbuff,  12, MPI_DOUBLE_PRECISION, &
         &              MPI_SUM  , MPI_COMM_WORLD, ierr                  )    

    !   -- [3-3]  Total Kinetic Energy          --  !
    Uue    = recvUbuff(1,ue_) + recvUbuff(2,ue_) + recvUbuff(3,ue_) 
    Uui    = recvUbuff(1,ui_) + recvUbuff(2,ui_) + recvUbuff(3,ui_) 
    Upe    = recvUbuff(1,pe_) + recvUbuff(2,pe_) + recvUbuff(3,pe_) 
    Upi    = recvUbuff(1,pi_) + recvUbuff(2,pi_) + recvUbuff(3,pi_) 
    Uke    = Uue  + Upe
    Uki    = Uui  + Upi
    Ukt    = Uke  + Uki

    ! ---------------------- !
    ! --- [4]  Display   --- !
    ! ---------------------- !
    if ( myRank.eq.0 ) then
       !  -- Display -- !
       write(6,'(2x,a)') "-------------------    Flow & Thermal Energy     ---------------------"
       write(6,'(2x,a,i8,1x,a,e12.5,1x)') 'kstep :: ', kstep, ' wce_t :: ', wcendt
       write(6,'(2x,1(a,e15.8,1x))'     ) '  Uk  :: ', Ukt  
       write(6,'(2x,3(a,e15.8,1x))'     ) '  Uke :: ', Uke  , ' Uue :: '  , Uue, ' Upe :: ', Upe
       write(6,'(2x,3(a,e15.8,1x))'     ) '  Uki :: ', Uki  , ' Uui :: '  , Uui, ' Upi :: ', Upi
       write(6,'(2x,a)') '----------------------------------------------------------------------'
       !  --  File   -- !
       FileName = trim(jobdir) // 'dat/' // 'energy_FlowThermal.dat'
       open( 50, file=trim(FileName), form='formatted', position='append' )
       write(50,'((i8,1x),20(e15.7,1x))') kstep, wcendt, Ukt, &
            &                            Uuer, Uuet, Uuez, Uue, Uper, Upet, Upez, Upe, Uke, &
            &                            Uuir, Uuit, Uuiz, Uui, Upir, Upit, Upiz, Upi, Uki
       close(50)
    endif

    return
  end subroutine Flow_Thermal_Energy

       
end module sEnergyMod
