module smoothfMod
contains
  
  subroutine smooth_E_Linear
    use constants , only : LIs, LJs, SmCoef, myRank, PEtot
    use constants , only : er_, et_, ez_
    use variables , only : EBf, rf , rh, rfInv, rhInv
    use mxwBdrcMod, only : BoundaryE
    implicit none
    integer             :: i, j
    double precision    :: denom
    double precision    :: EBs(3,LIs,LJs)

    ! -------------------------- !
    ! --- [1]  Preparation   --- !
    ! -------------------------- !
    !  -- [1-1] denominator  --  !
    denom = 1.d0 / ( 1.d0 + 2.d0 * SmCoef )**2
    !  -- [1-2] Copy to EBs  --  !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBs(er_,i,j) = EBf(er_,i,j)
          EBs(et_,i,j) = EBf(et_,i,j)
          EBs(ez_,i,j) = EBf(ez_,i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! -------------------------- !
    ! --- [2]   Smooth  Er   --- !
    ! -------------------------- !
    if ( myRank.eq.0       ) EBs(er_,  1,:) = + EBs(er_,    2,:)
    if ( myRank.eq.PEtot-1 ) EBs(er_,LIs,:) = + EBs(er_,LIs-1,:)
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,rf,rfInv,denom,LIs,LJs) private(i,j)
    !$omp do
    do j=3, LJs-1
       do i=2, LIs-1
          EBf(er_,i,j) = denom * rfInv(j) * ( rf(j  )*EBs(er_,i  ,j  ) &
               &               + SmCoef   * ( rf(j-1)*EBs(er_,i  ,j-1) + rf(j  )*EBs(er_,i-1,j  ) &
               &                            + rf(j  )*EBs(er_,i+1,j  ) + rf(j+1)*EBs(er_,i  ,j+1) ) &
               &               + SmCoef**2* ( rf(j-1)*EBs(er_,i-1,j-1) + rf(j-1)*EBs(er_,i+1,j-1) &
               &                            + rf(j+1)*EBs(er_,i-1,j+1) + rf(j+1)*EBs(er_,i+1,j+1) ) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! -------------------------- !
    ! --- [3]   Smooth  Et   --- !
    ! -------------------------- !
    if ( myRank.eq.0       ) EBs(et_,  1,:) = + EBs(et_,    2,:)
    if ( myRank.eq.PEtot-1 ) EBs(et_,LIs,:) = + EBs(et_,LIs-1,:)
    EBs(et_,:,  1) = - EBs(et_,:,    2)
    EBs(et_,:,LJs) = + EBs(et_,:,LJs-1)
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,rh,rhInv,denom,LIs,LJs) private(i,j)
    !$omp do
    do j=2, LJs-1
       do i=2, LIs-1
          EBf(et_,i,j) = denom * rhInv(j) * ( rh(j  )*EBs(et_,i,j) &
               &               + SmCoef   * ( rh(j-1)*EBs(et_,i  ,j-1) + rh(j  )*EBs(et_,i-1,j  ) &
               &                            + rh(j  )*EBs(et_,i+1,j  ) + rh(j+1)*EBs(et_,i  ,j+1) ) &
               &               + SmCoef**2* ( rh(j-1)*EBs(et_,i-1,j-1) + rh(j-1)*EBs(et_,i+1,j-1) &
               &                            + rh(j+1)*EBs(et_,i-1,j+1) + rh(j+1)*EBs(et_,i+1,j+1) ) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! -------------------------- !
    ! --- [4]   Smooth  Ez   --- !
    ! -------------------------- !
    ! EBs(ez_,    1, :) = + EBs(ez_,    3,   :)
    ! EBs(ez_,LIs+1, :) = + EBs(ez_,LIs-1,   :)
    EBs(ez_,   :,  1) = - EBs(ez_,   :,    2)
    EBs(ez_,   :,LJs) = + EBs(ez_,   :,LJs-1) * rh(LJs-1) * rhInv(LJs)
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,rh,rhInv,denom,LIs,LJs) private(i,j)
    !$omp do
    do j=2, LJs-1
       do i=3, LIs-1
          EBf(ez_,i,j) = denom * rhInv(j) * ( rh(j  )*EBs(ez_,i  ,j  ) &
               &               + SmCoef   * ( rh(j-1)*EBs(ez_,i  ,j-1) + rh(j  )*EBs(ez_,i-1,j  ) &
               &                            + rh(j  )*EBs(ez_,i+1,j  ) + rh(j+1)*EBs(ez_,i  ,j+1) ) &
               &               + SmCoef**2* ( rh(j-1)*EBs(ez_,i-1,j-1) + rh(j-1)*EBs(ez_,i+1,j-1) &
               &                            + rh(j+1)*EBs(ez_,i-1,j+1) + rh(j+1)*EBs(ez_,i+1,j+1) ) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel

    ! -------------------------- !
    ! --- [5]   Boundary E   --- !
    ! -------------------------- !
    call BoundaryE
    
    return
  end subroutine smooth_E_Linear


  subroutine smooth_B_Linear
    use constants   , only : LIs, LJs, SmCoef, myRank
    use constants   , only : br_, bt_, bz_
    use variables   , only : EBf, rf , rh, rfInv, rhInv
    use mxwBdrcMod  , only : BoundaryB
    use pcwBdrcMod  , only : conductingWall_divB
    implicit none
    integer               :: i, j
    double precision      :: denom
    double precision      :: EBs(6,LIs,LJs)
    
    ! -------------------------- !
    ! --- [1]  Preparation   --- !
    ! -------------------------- !
    !  -- [1-1] denominator  --  !
    denom = 1.d0 / ( 1.d0 + 2.d0 * SmCoef )**2
    !  -- [1-2] Copy to EBs  --  !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,LIs,LJs) private(i,j)
    !$omp do
    do j=1, LJs
       do i=1, LIs
          EBs(br_,i,j) = EBf(br_,i,j)
          EBs(bt_,i,j) = EBf(bt_,i,j)
          EBs(bz_,i,j) = EBf(bz_,i,j)
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    !  -- [1-3] B.C.          -- !
    call conductingWall_divB( EBs )
    EBs(br_,:,1) = - EBs(br_,:,2)
    
    ! -------------------------- !
    ! --- [2]   Smooth  Br   --- !
    ! -------------------------- !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,rh,rhInv,denom,LIs,LJs) private(i,j)
    !$omp do
    do j=2, LJs-1
       do i=2, LIs-1
          EBf(br_,i,j) = denom * rhInv(j) * ( rh(j  )*EBs(br_,i  ,j  ) &
               &               + SmCoef   * ( rh(j-1)*EBs(br_,i  ,j-1) + rh(j  )*EBs(br_,i-1,j  ) &
               &                            + rh(j  )*EBs(br_,i+1,j  ) + rh(j+1)*EBs(br_,i  ,j+1) ) &
               &               + SmCoef**2* ( rh(j-1)*EBs(br_,i-1,j-1) + rh(j-1)*EBs(br_,i+1,j-1) &
               &                            + rh(j+1)*EBs(br_,i-1,j+1) + rh(j+1)*EBs(br_,i+1,j+1) ) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( myRank.eq.0 ) EBf(br_,2,:) = EBs(br_,2,:)
    
    ! -------------------------- !
    ! --- [3]   Smooth  Bt   --- !
    ! -------------------------- !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,rf,rfInv,denom,LIs,LJs) private(i,j)
    !$omp do
    do j=3, LJs-1
       do i=2, LIs-1
          EBf(bt_,i,j) = denom * rfInv(j) * ( rf(j  )*EBs(bt_,i  ,j  ) &
               &               + SmCoef   * ( rf(j-1)*EBs(bt_,i  ,j-1) + rf(j  )*EBs(bt_,i-1,j  ) &
               &                            + rf(j  )*EBs(bt_,i+1,j  ) + rf(j+1)*EBs(bt_,i  ,j+1) ) &
               &               + SmCoef**2* ( rf(j-1)*EBs(bt_,i-1,j-1) + rf(j-1)*EBs(bt_,i+1,j-1) &
               &                            + rf(j+1)*EBs(bt_,i-1,j+1) + rf(j+1)*EBs(bt_,i+1,j+1) ) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    if ( myRank.eq.0 ) EBf(bt_,2,:) = EBs(bt_,2,:)
    
    ! -------------------------- !
    ! --- [4]   Smooth  Bz   --- !
    ! -------------------------- !
    !$omp parallel default(none) &
    !$omp shared(EBf,EBs,rf,rfInv,denom,LIs,LJs) private(i,j)
    !$omp do
    do j=3, LJs-1
       do i=2, LIs-1
          EBf(bz_,i,j) = denom * rfInv(j) * ( rf(j  )*EBs(bz_,i  ,j  ) &
               &               + SmCoef   * ( rf(j-1)*EBs(bz_,i  ,j-1) + rf(j  )*EBs(bz_,i-1,j  ) &
               &                            + rf(j  )*EBs(bz_,i+1,j  ) + rf(j+1)*EBs(bz_,i  ,j+1) ) &
               &               + SmCoef**2* ( rf(j-1)*EBs(bz_,i-1,j-1) + rf(j-1)*EBs(bz_,i+1,j-1) &
               &                            + rf(j+1)*EBs(bz_,i-1,j+1) + rf(j+1)*EBs(bz_,i+1,j+1) ) )
       enddo
    enddo
    !$omp end do
    !$omp end parallel
   
    ! -------------------------- !
    ! --- [5]  Boundary B    --- !
    ! -------------------------- !
    call BoundaryB

    return
  end subroutine smooth_B_Linear

  
  subroutine smooth_A_Linear( Ain )
    use constants, only : LIs, LJs, ns, SmCoef
    implicit none
    integer                         :: i, j, k
    double precision                :: denom
    double precision                :: smA(LIs,LJs,ns)
    double precision, intent(inout) :: Ain(LIs,LJs,ns)

    !  --- [1] smooth Ain ---  !
    denom = 1.d0 / ( 1.d0 + 2.d0*SmCoef )**2
    !$omp parallel default(none) &
    !$omp shared(Ain,smA,denom,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1, ns
       do j=2, LJs-1
          do i=2, LIs-1
             smA(i,j,k) =  denom     * ( Ain(i,j,k) &
                  &      + SmCoef    * ( Ain(i  ,j-1,k) + Ain(i-1,j  ,k) + Ain(i+1,j  ,k) + Ain(i  ,j+1,k) ) &
                  &      + SmCoef**2 * ( Ain(i-1,j-1,k) + Ain(i+1,j-1,k) + Ain(i-1,j+1,k) + Ain(i+1,j+1,k) ) )
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    !  --- [2] Substitute back --- !
    !$omp parallel default(none) &
    !$omp shared(Ain,smA,LIs,LJs) private(i,j,k)
    !$omp do
    do k=1,ns
       do j=2, LJs-1
          do i=2, LIs-1
             Ain(i,j,k) = smA(i,j,k)
          enddo
       enddo
    enddo
    !$omp end do
    !$omp end parallel
    
    return
  end subroutine smooth_A_Linear
  

  subroutine smooth_J_Linear
    use constants, only : LIs, LJs, SmCoef, SmCnt, ns
    use variables, only : Jcr, rf, rh, rfInv, rhInv
    implicit none
    integer            :: i, j, k
    double precision   :: denom
    double precision   :: smJcr(3,LIs,LJs,ns)
    integer, parameter :: nsm = 0  ! ns or 0
    
    ! --- [0] preparation ---  !
    ! --  boundary copy region --  !
    do k=0, nsm
       do i=1, LIs
          Jcr(3,i,1 ,k) = Jcr(3,i,2   ,k)
          Jcr(3,i,LJs,k) = Jcr(3,i,LJs-1,k)
          Jcr(2,i,1 ,k) = Jcr(2,i,2   ,k)
          Jcr(2,i,LJs,k) = Jcr(2,i,LJs-1,k)
       enddo
       do j=1, LJs
          Jcr(1,1 ,j,k) = Jcr(1,2   ,j,k)
          Jcr(1,LIs,j,k) = Jcr(1,LIs-1,j,k)
          Jcr(2,1 ,j,k) = Jcr(2,2   ,j,k)
          Jcr(2,LIs,j,k) = Jcr(2,LIs-1,j,k)
       enddo
    enddo
    
    !  --- [1] smooth Jr,Jt,Jz ---  !
    denom = 1.d0 / ( 1.d0 + 2.d0 * SmCoef )**2
    do k=0, nsm
       do j=3, LJs-1
          do i=3, LIs-1
             smJcr(1,i,j,k) = denom * rfInv(j) * ( rf(j)*Jcr(1,i,j,k) &
                  &         + SmCoef    * ( rf(j-1)*Jcr(1,i  ,j-1,k) + rf(j  )*Jcr(1,i-1,j  ,k) &
                  &                       + rf(j  )*Jcr(1,i+1,j  ,k) + rf(j+1)*Jcr(1,i  ,j+1,k) ) &
                  &         + SmCoef**2 * ( rf(j-1)*Jcr(1,i-1,j-1,k) + rf(j-1)*Jcr(1,i+1,j-1,k) &
                  &                       + rf(j+1)*Jcr(1,i-1,j+1,k) + rf(j+1)*Jcr(1,i+1,j+1,k) ) )
             smJcr(2,i,j,k) = denom * rhInv(j) * ( rh(j)*Jcr(2,i,j,k) &
                  &         + SmCoef    * ( rh(j-1)*Jcr(2,i  ,j-1,k) + rh(j  )*Jcr(2,i-1,j  ,k) &
                  &                       + rh(j  )*Jcr(2,i+1,j  ,k) + rh(j+1)*Jcr(2,i  ,j+1,k) ) &
                  &         + SmCoef**2 * ( rh(j-1)*Jcr(2,i-1,j-1,k) + rh(j-1)*Jcr(2,i+1,j-1,k) &
                  &                       + rh(j+1)*Jcr(2,i-1,j+1,k) + rh(j+1)*Jcr(2,i+1,j+1,k) ) )
             smJcr(3,i,j,k) = denom * rhInv(j) * ( rh(j)*Jcr(3,i,j,k) &
                  &         + SmCoef    * ( rh(j-1)*Jcr(3,i  ,j-1,k) + rh(j  )*Jcr(3,i-1,j  ,k) &
                  &                       + rh(j  )*Jcr(3,i+1,j  ,k) + rh(j+1)*Jcr(3,i  ,j+1,k) ) &
                  &         + SmCoef**2 * ( rh(j-1)*Jcr(3,i-1,j-1,k) + rh(j-1)*Jcr(3,i+1,j-1,k) &
                  &                       + rh(j+1)*Jcr(3,i-1,j+1,k) + rh(j+1)*Jcr(3,i+1,j+1,k) ) )
          enddo
       enddo
    enddo

    ! --- [2] Boundary Region --- !
    !  --  left  side   --  !
    i = 2
    do k=0, nsm
       do j=3, LJs-1
          smJcr(1,i,j,k) = denom * rfInv(j) * ( rf(j)*Jcr(1,i,j,k) &
               &         + SmCoef    * ( rf(j-1)*Jcr(1,i  ,j-1,k) + rf(j  )*Jcr(1,i-1,j  ,k) &
               &                       + rf(j  )*Jcr(1,i+1,j  ,k) + rf(j+1)*Jcr(1,i  ,j+1,k) ) &
               &         + SmCoef**2 * ( rf(j-1)*Jcr(1,i-1,j-1,k) + rf(j-1)*Jcr(1,i+1,j-1,k) &
               &                       + rf(j+1)*Jcr(1,i-1,j+1,k) + rf(j+1)*Jcr(1,i+1,j+1,k) ) )
          smJcr(2,i,j,k) = denom * rhInv(j) * ( rh(j)*Jcr(2,i,j,k) &
               &         + SmCoef    * ( rh(j-1)*Jcr(2,i  ,j-1,k) + rh(j  )*Jcr(2,i-1,j  ,k) &
               &                       + rh(j  )*Jcr(2,i+1,j  ,k) + rh(j+1)*Jcr(2,i  ,j+1,k) ) &
               &         + SmCoef**2 * ( rh(j-1)*Jcr(2,i-1,j-1,k) + rh(j-1)*Jcr(2,i+1,j-1,k) &
               &                       + rh(j+1)*Jcr(2,i-1,j+1,k) + rh(j+1)*Jcr(2,i+1,j+1,k) ) )
       enddo
    enddo
    !  --  bottom side  --  !
    j = 2
    do k=0, nsm
       do i=3, LIs-1
          smJcr(2,i,j,k) = denom * rhInv(j) * ( rh(j)*Jcr(2,i,j,k) &
               &         + SmCoef    * ( rh(j-1)*Jcr(2,i  ,j-1,k) + rh(j  )*Jcr(2,i-1,j  ,k) &
               &                       + rh(j  )*Jcr(2,i+1,j  ,k) + rh(j+1)*Jcr(2,i  ,j+1,k) ) &
               &         + SmCoef**2 * ( rh(j-1)*Jcr(2,i-1,j-1,k) + rh(j-1)*Jcr(2,i+1,j-1,k) &
               &                       + rh(j+1)*Jcr(2,i-1,j+1,k) + rh(j+1)*Jcr(2,i+1,j+1,k) ) )
          smJcr(3,i,j,k) = denom * rhInv(j) * ( rh(j)*Jcr(3,i,j,k) &
               &         + SmCoef    * ( rh(j-1)*Jcr(3,i  ,j-1,k) + rh(j  )*Jcr(3,i-1,j  ,k) &
               &                       + rh(j  )*Jcr(3,i+1,j  ,k) + rh(j+1)*Jcr(3,i  ,j+1,k) ) &
               &         + SmCoef**2 * ( rh(j-1)*Jcr(3,i-1,j-1,k) + rh(j-1)*Jcr(3,i+1,j-1,k) &
               &                       + rh(j+1)*Jcr(3,i-1,j+1,k) + rh(j+1)*Jcr(3,i+1,j+1,k) ) )
       enddo
    enddo
    !  -- corner region --  !
    i = 2
    j = 2
    do k=0, nsm
       smJcr(2,i,j,k) = denom * rhInv(j) * ( rh(j)*Jcr(2,i,j,k) &
            &         + SmCoef    * ( rh(j-1)*Jcr(2,i  ,j-1,k) + rh(j  )*Jcr(2,i-1,j  ,k) &
            &                       + rh(j  )*Jcr(2,i+1,j  ,k) + rh(j+1)*Jcr(2,i  ,j+1,k) ) &
            &         + SmCoef**2 * ( rh(j-1)*Jcr(2,i-1,j-1,k) + rh(j-1)*Jcr(2,i+1,j-1,k) &
            &                       + rh(j+1)*Jcr(2,i-1,j+1,k) + rh(j+1)*Jcr(2,i+1,j+1,k) ) )
    enddo
    
    !  --- [3] Substitute back --- !
    do k=0, nsm
       do j=3, LJs-1
          do i=3, LIs-1
             Jcr(1,i,j,k) = smJcr(1,i,j,k)
             Jcr(2,i,j,k) = smJcr(2,i,j,k)
             Jcr(3,i,j,k) = smJcr(3,i,j,k)
          enddo
       enddo
    enddo
    do k=0, nsm
       do j=3, LJs-1
          Jcr(1,2,j,k) = smJcr(1,2,j,k)
          Jcr(2,2,j,k) = smJcr(2,2,j,k)
       enddo
       do i=3, LIs-1
          Jcr(2,i,2,k) = smJcr(2,i,2,k)
          Jcr(3,i,2,k) = smJcr(3,i,2,k)
       enddo
       Jcr(2,2,2,k) = smJcr(2,2,2,k)
    enddo

    SmCnt = SmCnt + 1
    
    return
  end subroutine smooth_J_Linear

  
end module smoothfMod
