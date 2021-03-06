module jxBchckMod
contains

  ! =================================================================== !
  ! ===  jxB_Force  ::  calculate jxB Force on each Grid            === !
  ! =================================================================== !
  subroutine jxB_Force( onGrid )
    use constants , only      : coordinate
    use RgridRZMod, only      : LIr, LJr, BgR, JgR
    use BgridRZMod, only      : LIg, LJg, BgB, JgB, rf, rhinv
    implicit none
    integer                  :: i, j, LI, LJ, a1_, a2_, a3_
    double precision         :: a1Min, a1Max, a2Min, a2Max, aMin, aMax, aAbs
    double precision         :: Bfc(3,LIg,LJg), Jfc(3,LIg,LJg)
    double precision         :: cFr(LIg,LJg), cFt(LIg,LJg), acc(3,LIg,LJg)
    character(1), intent(in) :: onGrid
    integer, parameter       :: b1_=1, b2_=2, b3_=3, j1_=1, j2_=2, j3_=3

    ! ------------------------------------- !
    ! --- [1] which variable to use     --- !
    ! ------------------------------------- !
    !  -- [1-1] on Rgrid                --  !
    if ( onGrid.eq.'R' ) then
       LI  = LIr
       LJ  = LJr
       Bfc = BgR
       Jfc = JgR
    endif
    !  -- [1-2] on Bgrid                --  !
    if ( onGrid.eq.'B' ) then
       LI  = LIg
       LJ  = LJg
       if ( coordinate.eq.'XZ' ) then
          do j=2, LJ-1
             do i=2, LI-1
                Bfc(b1_,i,j) = 0.50d0 * ( BgB(b1_,i  ,j  ) + BgB(b1_,i+1,j  ) )
                Bfc(b2_,i,j) = 0.25d0 * ( BgB(b2_,i  ,j  ) + BgB(b2_,i+1,j  ) &
                     &                  + BgB(b2_,i  ,j+1) + BgB(b2_,i+1,j+1) )
                Bfc(b3_,i,j) = 0.50d0 * ( BgB(b3_,i  ,j  ) &
                     &                  + BgB(b3_,i  ,j+1) )
                Jfc(j1_,i,j) = 0.50d0 * ( JgB(j1_,i  ,j  ) &
                     &                  + JgB(j1_,i  ,j+1) )
                Jfc(j2_,i,j) =            JgB(j2_,i  ,j  )
                Jfc(j3_,i,j) = 0.50d0 * ( JgB(j3_,i  ,j  ) + JgB(j3_,i+1,j  ) )
             enddo
          enddo
       endif
       if ( coordinate.eq.'RZ' ) then
          do j=2, LJ-1
             do i=2, LI-1
                Bfc(b1_,i,j) = 0.50d0 * (           BgB(b1_,i  ,j  ) + BgB(b1_,i+1,j  ) )
                Bfc(b2_,i,j) = 0.25d0 * ( rf(j  )*( BgB(b2_,i  ,j  ) + BgB(b2_,i+1,j  ) ) &
                     &                  + rf(j+1)*( BgB(b2_,i  ,j+1) + BgB(b2_,i+1,j+1) ) ) * rhInv(j)
                Bfc(b3_,i,j) = 0.50d0 * ( rf(j  )*( BgB(b3_,i  ,j  )                  ) &
                     &                  + rf(j+1)*( BgB(b3_,i  ,j+1)                  ) ) * rhInv(j)
                Jfc(j1_,i,j) = 0.50d0 * ( rf(j  )*( JgB(j1_,i  ,j  )                  ) &
                     &                  + rf(j+1)*( JgB(j1_,i  ,j+1)                  ) ) * rhInv(j)
                Jfc(j2_,i,j) =                      JgB(j2_,i  ,j  )
                Jfc(j3_,i,j) = 0.50d0 * (           JgB(j3_,i  ,j  ) + JgB(j3_,i+1,j  ) )
             enddo
          enddo
       endif
    endif

    ! ------------------------------------- !
    ! --- [2] Force Calculation         --- !
    ! ------------------------------------- !
    !  -- [2-1] Centrifugal Force       --  !
    cFr(:,:)   = 0.d0
    cFt(:,:)   = 0.d0
    acc(:,:,:) = 0.d0
    if ( coordinate.eq.'RZ' ) then
       do j=1, LJ
          do i=1, LI
             cFr(i,j) = 0.d0
             cFt(i,j) = 0.d0
          enddo
       enddo
    endif
    !  -- [2-2] Acceleration            --  !
    a1_=1; a2_=2; a3_=3
    do j=2, LJ-1
       do i=2, LI-1
          acc(a1_,i,j) = ( Jfc(j2_,i,j)*Bfc(b3_,i,j) - Jfc(j3_,i,j)*Bfc(b2_,i,j) ) + cFr(i,j)
          acc(a2_,i,j) = ( Jfc(j3_,i,j)*Bfc(b1_,i,j) - Jfc(j1_,i,j)*Bfc(b3_,i,j) ) - cFt(i,j)
          acc(a3_,i,j) = ( Jfc(j1_,i,j)*Bfc(b2_,i,j) - Jfc(j2_,i,j)*Bfc(b1_,i,j) )
       enddo
    enddo

    ! ------------------------------------- !
    ! --- [3] Display Max Min of Acc.   --- !
    ! ------------------------------------- !
    !  -- [3-1] coordinate settings     --  !
    if ( coordinate.eq.'XZ' ) then
       write(6,'(2x,    a)') '[ jxB_Force     @jxBcheckMod ]'
       write(6,'(5x,    a)') ' -- Cartesian   Coordinate --'
       write(6,'(10x,   a)') 'a1 :: Fx'
       write(6,'(10x,   a)') 'a2 :: Fz'
       a1_ = 3 ; a2_ = 1
    endif
    if ( coordinate.eq.'RZ' ) then
       write(6,'(2x,    a)') '[ jxB_Force     @jxBcheckMod ]'
       write(6,'(5x,    a)') ' -- Cylindrical Coordinate --'
       write(6,'(10x,   a)') 'a1 :: Fz'
       write(6,'(10x,   a)') 'a2 :: Fr'
       a1_ = 3 ; a2_ = 1
    endif

    !  -- [3-2] Max & Min               --  !
    a1Min = + 1.d8; a2Min = + 1.d8; aMin  = + 1.d8
    a1Max = - 1.d8; a2Max = - 1.d8; aMax  = - 1.d8
    do j=5, LJ-5
       do i=5, LI-5
          aAbs  = sqrt( acc(1,i,j)**2 + acc(2,i,j)**2 + acc(3,i,j)**2 )
          a1Min =  min( a1Min, acc(a1_,i,j) )
          a1Max =  max( a1Max, acc(a1_,i,j) )
          a2Min =  min( a2Min, acc(a2_,i,j) )
          a2Max =  max( a2Max, acc(a2_,i,j) )
          aMin  =  min(  aMin, aAbs         )
          aMax  =  max(  aMax, aAbs         )
       enddo
    enddo

    ! ------------------------------------------------------ !
    ! --- [ERROR] somehow Seg.Fault occur with belows    --- !
    ! ------------------------------------------------------ !
    !  -- error depends on compiler...                   --  !
    !  --            * macOS(gfortran) :: x              --  !
    ! -----------------------------------------------------  !
    write(6,*)
    write(6,*) " NO Display.... ( because of the unreasonable Seg.Fault... see jxBcheckMod.f90 )"
    write(6,*)
    write(6,*)
    ! write(6,*)
    ! write(6,'(5x,20x,a)') '( Min ) | (Max)'
    ! write(6,*) a1Min, a1Max
    ! write(6,'(10x,a,f15.8,a,f15.8)') ' a1 : ', a1Min, ' | ', a1Max
    ! write(6,'(10x,a,e12.5,a,e12.5)') ' a2 : ', a2Min, ' | ', a2Max
    ! write(6,'(10x,a,e12.5,a,e12.5)') '|a| : ',  aMin, ' | ',  aMax
    ! write(6,*)
    return
  end subroutine jxB_Force


  ! =================================================================== !
  ! ===  checkNan  ::  check NaN of Avp & Bfd                       === !
  ! =================================================================== !
  subroutine checkNan
    use constants, only : N1, N2
    use variables, only : Avp, Bfd
    implicit none
    integer            :: i, j
    
    do j=1, N2
       do i=1, N1
          if ( Avp(  i,j).ne.Avp(  i,j) ) write(6,*) "Nan Detec. :: ( i, j, Avp ) ", i, j, Avp(  i,j)
          if ( Bfd(3,i,j).ne.Bfd(3,i,j) ) write(6,*) "Nan Detec. :: ( i, j, Bfd ) ", i, j, Bfd(3,i,j)
       enddo
    enddo
    return
  end subroutine checkNan
  
  
end module jxBchckMod
