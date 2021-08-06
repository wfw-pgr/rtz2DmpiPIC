module prdBdrcMod
contains
  
  subroutine periodic_EBr
    use constants , only : LIs, LJs
    use variables , only : EBr
    use rMPIComMod, only : RelocatedExchange
    implicit none
    integer            :: i, cmp
    ! ----------------------------------- !
    ! --- [1]    Boundary_2 ( E & B ) --- !
    ! ----------------------------------- !
    !$omp parallel default(none) &
    !$omp shared(EBr,LIs,LJs) private(i,cmp)
    !$omp do
    do i=1, LIs
       do cmp=1, 6
          EBr(cmp,i,    1) = EBr(cmp,i,LJs-1)
          EBr(cmp,i,LJs  ) = EBr(cmp,i,    2)
          EBr(cmp,i,LJs+1) = EBr(cmp,i,    3)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    ! ----------------------------------- !
    ! --- [1]    Boundary_2 ( E )     --- !
    ! ----------------------------------- !
    call RelocatedExchange( EBr )

    return
  end subroutine periodic_EBr

end module prdBdrcMod
