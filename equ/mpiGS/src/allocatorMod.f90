module allocatorMod
contains
  
  ! ====================================================== !
  ! ===  allocate variables                            === !
  ! ====================================================== !
  subroutine allocate__variables
    use variablesMod
    implicit none

    ! ------------------------------------------------------ !
    ! --- [1] allocation of variables                    --- !
    ! ------------------------------------------------------ !
    allocate(         zAxis(LI), source=0.d0 )
    allocate(         rAxis(LJ), source=0.d0 )
    allocate(          rInv(LJ), source=0.d0 )
    allocate(        rhs(LI,LJ), source=0.d0 )
    allocate(        psi(LI,LJ), source=0.d0 )
    allocate( psi_past(LI,LJ,2), source=0.d0 )
    allocate(       psis(LI,LJ), source=0.d0 )
    allocate(    current(LI,LJ), source=0.d0 )
    allocate(   pressure(LI,LJ), source=0.d0 )
    allocate(        rho(LI,LJ), source=0.d0 )
    allocate(      gfunc(LI,LJ), source=0.d0 )
    allocate(     bdrcLR(LI,LJ), source=0.d0 )
    allocate(     bdrcTB(LI,LJ), source=0.d0 )
    return
  end subroutine allocate__variables

end module allocatorMod
