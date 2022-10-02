module parameterMod
contains

  ! ====================================================== !
  ! === load parameters                                === !
  ! ====================================================== !
  subroutine load__parameters
    use variablesMod
    implicit none

    ! ------------------------------------------------------ !
    ! --- [1] set namelists                              --- !
    ! ------------------------------------------------------ !
    namelist /equSettings/ job, solverType, equiliType, parityType, caseIOType, &
         &                 normalizationType, normalizationInGS, coordinateType, boundaryType, &
         &                 gridType, normalizeField
    namelist /grid/        LI, LJ, Lmid, zMin_mhd, zMax_mhd, rMin_mhd, rMax_mhd
    namelist /pic/         vthcv, wpewce, TiTe, mr, dr_Debye, dz_Debye
    namelist /equ/         iterMax1, iterMax2, convergence1, convergence2, Picard_alphaB, &
         &                 psisw, pTop, beta_pol, beta_tor, aspectRatio, Bmax, rhoFloor, &
         &                 epsilon, alpha, eta, c1, c2, beta_sep, Djh, lambda, &
         &                 alpha_star, Bsensor_r1, Bsensor_r2, &
         &                 limiter_rPos, limiter_zPos

    ! ------------------------------------------------------ !
    ! --- [2] load from file                             --- !
    ! ------------------------------------------------------ !
    limiter_rPos(:) = -1.d0
    limiter_zPos(:) = -1.d0
    open (lun,file=trim(parameterFile),status="old",form="formatted")
    read (lun,nml=equSettings)
    read (lun,nml=grid)
    read (lun,nml=pic)
    read (lun,nml=equ)
    close(lun)

    return
  end subroutine load__parameters
  

  ! ====================================================== !
  ! ===  set parameters                                === !
  ! ====================================================== !
  subroutine set__parameters
    use variablesMod
    implicit none
    integer          :: ik
    double precision :: lDebye

    ! ------------------------------------------------------ !
    ! --- [1] make Grid parameters for PIC simulation    --- !
    ! ------------------------------------------------------ !
    if ( trim(normalizationType) == "PIC" ) then
       ! -- [1-1] Resolution of Grid  -- ! 
       lDebye     =    vthcv / wpewce
       dr         = dr_Debye * lDebye
       dz         = dz_Debye * lDebye
       drInv      =     1.d0 / dr
       dzInv      =     1.d0 / dz
       rLeng      = dble( LJ-1 ) * dr
       zLeng      = dble( LI-1 ) * dz
       ! -- [1-2] Grid Min-Max Region -- ! 
       zMin       = ( - 0.5d0 ) * zLeng
       zMax       = ( + 0.5d0 ) * zLeng
       if ( aspectRatio > 1.d0 ) then
          rMin    = rLeng * ( aspectRatio-1.d0 ) * ( 0.5d0 * ( Bsensor_r1 + Bsensor_r2 ) )
       else
          rMin    = 0.d0
       endif
       rMax       = rLeng + rMin
    endif

    ! ------------------------------------------------------ !
    ! --- [2] make Grid parameters for MHD simulation    --- !
    ! ------------------------------------------------------ !
    if ( trim(normalizationType) == "MHD" ) then
       ! -- Resolution of Grid  -- !
       zMin       = zMin_mhd
       zMax       = zMax_mhd
       rMin       = rMin_mhd
       rMax       = rMax_mhd
       rLeng      = rMax_mhd - rMin_mhd
       zLeng      = zMax_mhd - zMin_mhd
       dz         = zLeng / dble( LI-1 )
       dr         = rLeng / dble( LJ-1 )
       drInv      = 1.d0 / dr
       dzInv      = 1.d0 / dz
    endif

    ! ------------------------------------------------------ !
    ! --- [3] set other parameters                       --- !
    ! ------------------------------------------------------ !
    !  -- [3-1] Alfven velocity for electron             --  !
    valfe         = 1.0d0 / wpewce
    psiMin        = 0.d0
    psiLim        = 0.d0
    magAxis(z_)   = 0.5d0 * ( zMax + zMin )
    magAxis(r_)   = 0.5d0 * ( rMax + rMin )

    !  -- [3-2] B_EF settings                            --  !
    Bsensor_r1    = Bsensor_r1 * rLeng + rMin
    Bsensor_r2    = Bsensor_r2 * rLeng + rMin
    nLimiter      = -1
    do ik=1, nLimiter_max
       if ( ( limiter_rPos(ik) == -1.d0 ).and.( limiter_zPos(ik) == -1.d0 ) ) then
          nLimiter = ik-1
          exit
       endif
       limiter_rPos(ik) = limiter_rPos(ik) * rLeng + rMin
       limiter_zPos(ik) = limiter_zPos(ik) * zLeng + zMin
    enddo
    if ( nLimiter == -1 ) then
       if ( myRank.eq.0 ) then
          write(6,*) "[initiatorMod.f90] nLimiter is illegal          [ERROR]"
          write(6,*) "[initiatorMod.f90]    * check limiter_rPos/zPos value..."
       endif
    endif
    
    return
  end subroutine set__parameters

  
end module parameterMod
