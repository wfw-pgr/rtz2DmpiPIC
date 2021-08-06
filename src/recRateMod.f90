module recRateMod
contains

  subroutine ReconnectionRate( ptOXO, RateR, mergR, Flag__sOpt )
    use constants,          only   : LIs, LJs, valfe, mr, dt, myRank, jobDir
    use constants,          only   : br_, bt_, bz_, er_, et_, iFr_, iTo_, isl_, iel_
    use variables,          only   : EBf, kstep, ijDomain
    implicit none
    include 'mpif.h'
    integer         , intent(in)  :: ptOXO(3,2)
    double precision, intent(in)  :: mergR
    double precision, intent(out) :: RateR
    logical         , intent(in)  :: Flag__sOpt
    integer         , parameter   :: ihwidth = 2, jhwidth = 5
    integer         , parameter   :: id_     = 1, jd_     = 2
    integer                       :: i, j, iL, iR, jOave, iOave, iRange(2), jRange(2), ierr
    double precision              :: almost, time, B2Max, absB2
    double precision              :: sumEx, sumEy, sumEr, sumE_s(2), sumE_r(2)
    double precision              :: RateX, RateY, vAin , Bxin , vABx, vABxInv
    character(100)                :: FileName
    logical                       :: Flag__SingleMode

    ! ------------------------------ !
    ! --- [1] Erec ( Averaged )  --- !
    ! ------------------------------ !
    almost    = 0.d0
    iRange(1) =      ptOXO(2,id_)-ihwidth
    iRange(2) =      ptOXO(2,id_)+ihwidth
    jRange(1) = max( ptOXO(2,jd_)-jhwidth,     3 )
    jRange(2) = min( ptOXO(2,jd_)+jhwidth, LJs-1 )
    sumE_s(:) = 0.d0
    sumE_r(:) = 0.d0
    do i=iRange(1), iRange(2)
       if ( ( i.ge.ijDomain(myRank,iFr_) ).and.( i.le.ijDomain(myRank,iTo_) ) ) then
          iL = i - ijDomain(myRank,iFr_) + ijDomain(myRank,isl_)
          do j=jRange(1), jRange(2)
             sumE_s(1) = sumE_s(1) + EBf(er_,iL,j)
             sumE_s(2) = sumE_s(2) + EBf(et_,iL,j)
          enddo
       endif
    enddo
    call MPI_AllReduce( sumE_s, sumE_r, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    sumEx  = sumE_r(1) / dble( ( 2*ihwidth+1 )*( 2*jhwidth+1 ) )
    sumEy  = sumE_r(2) / dble( ( 2*ihwidth+1 )*( 2*jhwidth+1 ) )
    
    ! ------------------------------ !
    ! --- [2] BxIn / vAin        --- !
    ! ------------------------------ !
    if ( ( Flag__sOpt ).or.( ptOXO(3,id_)-ptOXO(1,id_).lt.4 ).or.( mergR.gt.0.95 ) ) then
       Flag__SingleMode = .true.
    else
       Flag__SingleMode = .false.
    endif
    
    if ( .not.( Flag__SingleMode ) ) then
       jOave  = max( min( ( ptOXO(1,jd_)+ptOXO(3,jd_) )/2, LJs-2 ), 3 )
       if      ( ( ptOXO(1,id_).ge.ijDomain(myRank,iFr_) ).and.( ptOXO(1,id_).le.ijDomain(myRank,iTo_) ) ) then
          ! --     | iO1 | case -- !
          iL  = ptOXO(1,id_) - ijDomain(myRank,iFr_) + ijDomain(myRank,isl_)
          iR  = min( ijDomain(myRank,iel_), ptOXO(3,id_)-ijDomain(myRank,iFr_)+ijDomain(myRank,isl_) )
       else if ( ( ptOXO(3,id_).ge.ijDomain(myRank,iFr_) ).and.( ptOXO(3,id_).le.ijDomain(myRank,iTo_) ) ) then
          ! --     | iO3 | case -- !
          iL  = max( ijDomain(myRank,isl_), ptOXO(1,id_)-ijDomain(myRank,iFr_)+ijDomain(myRank,isl_) )
          iR  = ptOXO(3,id_) - ijDomain(myRank,iFr_) + ijDomain(myRank,isl_)
       else if ( ( ptOXO(1,id_).lt.ijDomain(myRank,iFr_) ).and.( ptOXO(3,id_).gt.ijDomain(myRank,iTo_) ) ) then
          ! -- iO1 |     | iO3 case -- !
          iL  = ijDomain(myRank,isl_)
          iR  = ijDomain(myRank,iel_)
       else
          ! -- [ iO1 iO3 |    |  ] or [ |    | iO1 iO3 ] case -- !
          iL  = 2
          iR  = 1
       endif
       B2Max  = 0.d0
       absB2  = 0.d0
       do i=iL, iR
          absB2 = EBf(br_,i,jOave)**2 + EBf(bt_,i,jOave)**2
          if ( absB2.gt.B2Max ) then
             B2Max = absB2
          endif
       enddo
       call MPI_AllReduce( B2Max, absB2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
       if ( absB2.eq.0.d0 ) Flag__SingleMode = .true.
    endif
    if ( Flag__SingleMode ) then
       ! --- Single O-point Case --- !
       ! -- use Common Poloidal Flux as vAin -- !
       iOave = max( min( ( ptOXO(1,id_)+ptOXO(3,id_) )/2, LIs-2 ), 3 )
       B2Max = 0.d0
       absB2 = 0.d0
       if ( ( iOave.ge.ijDomain(myRank,iFr_) ).and.( iOave.le.ijDomain(myRank,iTo_) ) ) then
          do j=3, LJs-2
             absB2 = EBf(br_,iOave,j)**2 + EBf(bt_,iOave,j)**2
             if ( absB2.gt.B2Max ) then
                B2Max = absB2
             endif
          enddo
       endif
       call MPI_AllReduce( B2Max, absB2, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr )
       almost= 1.d0
    endif
    !  - ( assumption :: ni ~ 1 ) - !
    Bxin  = sqrt( absB2 )
    vAin  = valfe / sqrt( mr ) * Bxin

    ! ------------------------------ !
    ! --- [3] Reconnection Rate  --- !
    ! ------------------------------ !    
    vABx   = vAin*Bxin
    if ( vABx.eq.0.d0 ) then
       vABxInv = 0.d0
    else
       vABxInv = 1.d0 / vABx
    endif
    sumEr  = sqrt( sumEx**2+sumEy**2 )
    RateX  = sumEx * vABxInv
    RateY  = sumEy * vABxInv
    RateR  = sumEr * vABxInv
    
    ! ------------------------------ !
    ! --- [4]   Write in File    --- !
    ! ------------------------------ !
    if ( myRank.eq.0 ) then
       ! --   Reconnection Rate   -- !
       time     = dble(kstep)*dt
       FileName = trim(jobdir) // 'dat/' // 'ReconnectionRate.dat'
       open(50,file=trim(FileName),form='formatted',access='append')
       write(50,'(10(e15.8,1x))') time , &
            &                     RateR, RateX , RateY, &
            &                     vAin , Bxin  , sumEx, sumEy, &
            &                     mergR, almost
       close(50)
    endif
    
    return
  end subroutine ReconnectionRate
  
end module recRateMod
