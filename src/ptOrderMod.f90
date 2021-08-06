module ptOrderMod
  use constants, only  : ns
  implicit none
  integer             :: exptc(ns)
contains


  subroutine InquireIndex( order )
    use constants , only : LIs, LJs, drInv, dzInv
    use constants , only : ns , np , OMPNumThreads
    use constants , only : rp_, zp_, wp_
    use variables , only : pxv, pCellIdx
    !$ use omp_lib
    implicit none
    integer                    :: k, m, ip, jp, ith, OutOfIndex
    character(5)  , intent(in) :: order

    if      ( order.eq.'iCell' ) then
       OutOfIndex = (LIs-2)*(LJs-2)+1
       do k=1, ns
          !$omp parallel default(none) &
          !$omp shared(k,np,pxv,pCellIdx,dzInv,drInv,LIs,LJs,OutOfIndex) private(m,ip,jp,ith)
          !$ ith = omp_get_thread_num() + 1
          !$omp do
          do m=1, np(k)
             if ( pxv(wp_,m,k).gt.0.d0 ) then
                ip            = ceiling( pxv(zp_,m,k)*dzInv ) 
                jp            = ceiling( pxv(rp_,m,k)*drInv )
                pCellIdx(m,k) = (LIs-2)*(jp-1) + ip
             else
                pCellIdx(m,k) = OutOfIndex
             endif
          enddo
          !$omp end do
          !$omp end parallel
       enddo
    else if ( order.eq.'jCell' ) then
       OutOfIndex = (LIs-2)*(LJs-2)+1
       do k=1, ns
          !$omp parallel default(none) &
          !$omp shared(k,np,pxv,pCellIdx,dzInv,drInv,LIs,LJs,OutOfIndex) private(m,ip,jp,ith)
          !$ ith = omp_get_thread_num() + 1
          !$omp do
          do m=1, np(k)
             if ( pxv(wp_,m,k).gt.0.d0 ) then
                ip            = ceiling( pxv(zp_,m,k)*dzInv ) 
                jp            = ceiling( pxv(rp_,m,k)*drInv )
                pCellIdx(m,k) = (LJs-2)*(ip-1) + jp
             else
                pCellIdx(m,k) = OutOfIndex
             endif
          enddo
          !$omp end do
          !$omp end parallel
       enddo
    else if ( order.eq.'iGrid' ) then
       OutOfIndex = (LIs-2)*(LJs-2)+1
       do k=1, ns
          !$omp parallel default(none) &
          !$omp shared(k,np,pxv,pCellIdx,dzInv,LIs,OutOfIndex) private(m,ith)
          !$ ith = omp_get_thread_num() + 1
          !$omp do
          do m=1, np(k)
             if ( pxv(wp_,m,k).gt.0.d0 ) then
                pCellIdx(m,k) = ceiling( pxv(zp_,m,k)*dzInv ) + 1
             else
                pCellIdx(m,k) = OutOfIndex
             endif
          enddo
          !$omp end do
          !$omp end parallel
       enddo
    endif
    
    return
  end subroutine InquireIndex
  

  subroutine ptReOrdering( order )
    use constants , only        : LIs, LJs, ns, np, npc, npt, nptMax, OMPNumThreads
    use variables , only        : pxv, pCellIdx
    use sortingMod, only        : ParallelSort
    implicit none
    integer                    :: k, m, OutIndex
    character(5)  , intent(in) :: order

    ! --- [1] Acquire Cell Index  ( Sorting Order )  --- !
    call InquireIndex( order )
    
    ! --- [2] Particle Sorting according to pCellIdx --- !
    do k=1, ns
       call ParallelSort( nptMax, pCellIdx(1:nptMax,k), 1, np(k),     &
            &                     pxv (1:8,1:nptMax,k), OMPNumThreads )
    enddo
    ! --- [3] np & npc Management                    --- !
    select case( order )
    case( 'iCell' )
       OutIndex = (LIs-2)*(LJs-2)+1
    case( 'jCell' )
       OutIndex = (LIs-2)*(LJs-2)+1
    case( 'iGrid' )
       OutIndex = (LIs-2)*(LJs-2)+1
    end select
    
    do k=1, ns
       do m=np(k), 1, -1
          if ( pCellIdx(m,k).lt.OutIndex ) then
             np (k) = m
             npc(k) = m
             exit
          endif
       enddo
    enddo
    npt = 0
    do k=1, ns
       npt    = npt + npc(k)
    enddo
    
    return
  end subroutine ptReOrdering

end module ptOrderMod



  ! subroutine iCellIndex
  !   use constants , only : LIs, LJs, drInv, dzInv
  !   use constants , only : ns , np , OMPNumThreads
  !   use constants , only : rp_, zp_, wp_
  !   use variables , only : pxv, pCellIdx
  !   !$ use omp_lib
  !   implicit none
  !   integer             :: k, m, ip, jp, ith
  !   integer             :: exCount(OMPNumThreads)
  !   do k=1, ns
  !      exCount = 0
  !      !$omp parallel default(none) &
  !      !$omp shared(k,np,pxv,pCellIdx,dzInv,drInv,LIs,LJs,exCount) &
  !      !$omp private(m,ip,jp,ith)
  !      !$ ith = omp_get_thread_num() + 1
  !      !$omp do
  !      do m=1, np(k)
  !         if ( pxv(wp_,m,k).gt.0.d0 ) then
  !            ip            = ceiling( pxv(zp_,m,k)*dzInv ) 
  !            jp            = ceiling( pxv(rp_,m,k)*drInv )
  !            pCellIdx(m,k) = (LIs-2)*(jp-1) + ip
  !         else
  !            pCellIdx(m,k) = (LIs-2)*(LJs-2)+1
  !            exCount(ith)  = exCount(ith) + 1
  !         endif
  !      enddo
  !      !$omp end do
  !      !$omp end parallel
  !      exptc(k) = sum( exCount )
  !   enddo
    
  !   return
  ! end subroutine iCellIndex


  ! subroutine iGridIndex
  !   use constants , only : LIs, dzInv
  !   use constants , only : ns , np , OMPNumThreads
  !   use constants , only : zp_, wp_
  !   use variables , only : pxv, pCellIdx
  !   !$ use omp_lib
  !   implicit none
  !   integer             :: k, m, ith
  !   integer             :: exCount(OMPNumThreads)
  !   do k=1, ns
  !      exCount = 0
  !      !$omp parallel default(none) &
  !      !$omp shared(k,np,pxv,pCellIdx,dzInv,LIs,exCount) &
  !      !$omp private(m,ith)
  !      !$ ith = omp_get_thread_num() + 1
  !      !$omp do
  !      do m=1, np(k)
  !         if ( pxv(wp_,m,k).gt.0.d0 ) then
  !            pCellIdx(m,k) = ceiling( pxv(zp_,m,k)*dzInv ) + 1
  !         else
  !            pCellIdx(m,k) = LIs+2
  !            exCount(ith)  = exCount(ith) + 1
  !         endif
  !      enddo
  !      !$omp end do
  !      !$omp end parallel
  !      exptc(k) = sum( exCount )
  !   enddo
    
  !   return
  ! end subroutine iGridIndex
