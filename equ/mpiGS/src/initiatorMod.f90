module initiatorMod
contains
  
  ! ====================================================== !
  ! === initialize__mpiGS                              === !
  ! ====================================================== !
  subroutine initialize__mpiGS
    use variablesMod
    use parameterMod, only : load__parameters, set__parameters
    use allocatorMod, only : allocate__variables
    use PBiCGStabMod, only : BiCGSTABCyl2D_MPI
    use currentFnMod, only : set__initial_currentDistribution, update__currentDistribution
    use boundaryCMod, only : set__initial_boundaryCondition
    use solveGSEqMod, only : determine__psiPotential
    use saveFieldMod
    implicit none
    
    ! ------------------------------------------------------ !
    ! --- [1] load parameters & initialize variables     --- !
    ! ------------------------------------------------------ !
    call load__parameters
    call  set__parameters
    call allocate__variables
    call make__coordinate
    
    ! ------------------------------------------------------ !
    ! --- [2] try to find approx. GS-Solution            --- !
    ! ------------------------------------------------------ !
    call set__initial_currentDistribution( "arbitral" )
    call set__initial_boundaryCondition

    call BiCGSTABCyl2D_MPI( psi, rhs, dz, dr, LI, LJ, rMin )
    call determine__psiPotential
    call update__currentDistribution

    if ( myRank == 0 ) then
       call save__map2D_asPointFile( psi    , zAxis, rAxis, LI, LJ, 1, "dat/psi.dat"     )
       call save__map2D_asPointFile( current, zAxis, rAxis, LI, LJ, 1, "dat/current.dat" )
    endif

    ! ------------------------------------------------------ !
    ! --- [3] set initial current & B.C.                 --- !
    ! ------------------------------------------------------ !
    call set__initial_currentDistribution( "pressure" )
    call set__initial_boundaryCondition

    if ( myRank == 0 ) then
       call save__map2D_asPointFile( current, zAxis, rAxis, LI, LJ, 1, "dat/current2.dat" )
    endif

    return
  end subroutine initialize__mpiGS


  ! ====================================================== !
  ! === make rAxis / zAxis                             === !
  ! ====================================================== !
  subroutine make__coordinate
    use variablesMod
    implicit none
    integer :: i, j

    ! ------------------------------------------------------ !
    ! --- [1] set rAxis / zAxis                          --- !
    ! ------------------------------------------------------ !
    do i=1, LI
       zAxis(i) = dz * dble(i-1) + zMin
    enddo
    do j=1, LJ
       rAxis(j) = dr * dble(j-1) + rMin
    enddo
    ! ------------------------------------------------------ !
    ! --- [2] set rInv                                   --- !
    ! ------------------------------------------------------ !
    do j=1, LJ
       if ( rAxis(j) == 0.d0 ) then
          rInv(j) = 0.d0
       else
          rInv(j) = 1.d0 / rAxis(j)
       endif
    enddo
    return
  end subroutine make__coordinate

  
end module initiatorMod
