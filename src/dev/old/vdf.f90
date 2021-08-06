module VDFModule
contains
  
  subroutine VDFsample
    use constants, only : LI, LJ, ns, np, jobdir
    use variables, only : pxv, kstep, x1Leng, x2Leng
    implicit none
    include 'mpif.h'
    integer, parameter :: Ndiv = 101, Nset=11
    integer            :: i, k, iset, ierr
    double precision   :: lrbt(4,Nset), vMaxMin(2,ns,Nset)
    double precision   :: vdf(0:Ndiv+1,0:Ndiv+1,3,ns,Nset), vdf_s(0:Ndiv+1,0:Ndiv+1,3,ns)
    character(len=  8) :: cstep
    character(len=100) :: filename
    logical, save      :: flag__parameter = .true.

    ! --- [1] Preparation --- !
    !  -- [1-1] Set 1 -- !
    lrbt(1,1)      = x1Leng * 0.48d0
    lrbt(2,1)      = x1Leng * 0.52d0
    lrbt(3,1)      = x2Leng * 0.48d0
    lrbt(4,1)      = x2Leng * 0.52d0
    vMaxMin(1,1,1) = - 0.5d0
    vMaxMin(2,1,1) = + 0.5d0
    vMaxMin(1,2,1) = - 0.2d0
    vMaxMin(2,2,1) = + 0.2d0
    !  -- [1-2] Set 2 -- !
    lrbt(1,2)      = x1Leng * 0.48d0
    lrbt(2,2)      = x1Leng * 0.52d0
    lrbt(3,2)      = x2Leng * 0.62d0
    lrbt(4,2)      = x2Leng * 0.66d0
    vMaxMin(1,1,2) = - 0.5d0
    vMaxMin(2,1,2) = + 0.5d0
    vMaxMin(1,2,2) = - 0.2d0
    vMaxMin(2,2,2) = + 0.2d0
    !  -- [1-3] Set 3 -- !
    lrbt(1,3)      = x1Leng * 0.40d0
    lrbt(2,3)      = x1Leng * 0.44d0
    lrbt(3,3)      = x2Leng * 0.48d0
    lrbt(4,3)      = x2Leng * 0.52d0
    vMaxMin(1,1,3) = - 0.5d0
    vMaxMin(2,1,3) = + 0.5d0
    vMaxMin(1,2,3) = - 0.2d0
    vMaxMin(2,2,3) = + 0.2d0
    !  -- [1-4] Set 4 -- !
    lrbt(1,4)      = x1Leng * 0.56d0
    lrbt(2,4)      = x1Leng * 0.60d0
    lrbt(3,4)      = x2Leng * 0.48d0
    lrbt(4,4)      = x2Leng * 0.52d0
    vMaxMin(1,1,4) = - 0.5d0
    vMaxMin(2,1,4) = + 0.5d0
    vMaxMin(1,2,4) = - 0.2d0
    vMaxMin(2,2,4) = + 0.2d0
    !  -- [1-5] Set 5 -- !
    lrbt(1,5)      = x1Leng * 0.55d0
    lrbt(2,5)      = x1Leng * 0.60d0
    lrbt(3,5)      = x2Leng * 0.55d0
    lrbt(4,5)      = x2Leng * 0.60d0
    vMaxMin(1,1,5) = - 0.5d0
    vMaxMin(2,1,5) = + 0.5d0
    vMaxMin(1,2,5) = - 0.2d0
    vMaxMin(2,2,5) = + 0.2d0
    !  -- [1-6] Set 6 -- !
    lrbt(1,6)      = x1Leng * 0.55d0
    lrbt(2,6)      = x1Leng * 0.60d0
    lrbt(3,6)      = x2Leng * 0.40d0
    lrbt(4,6)      = x2Leng * 0.45d0
    vMaxMin(1,1,6) = - 0.5d0
    vMaxMin(2,1,6) = + 0.5d0
    vMaxMin(1,2,6) = - 0.2d0
    vMaxMin(2,2,6) = + 0.2d0
    !  -- [1-7] Set 7 -- !
    lrbt(1,7)      = x1Leng * 0.45d0
    lrbt(2,7)      = x1Leng * 0.55d0
    lrbt(3,7)      = x2Leng * 0.72d0
    lrbt(4,7)      = x2Leng * 0.78d0
    vMaxMin(1,1,7) = - 0.5d0
    vMaxMin(2,1,7) = + 0.5d0
    vMaxMin(1,2,7) = - 0.2d0
    vMaxMin(2,2,7) = + 0.2d0
    !  -- [1-8] Set 8 -- !
    lrbt(1,8)      = x1Leng * 0.45d0
    lrbt(2,8)      = x1Leng * 0.55d0
    lrbt(3,8)      = x2Leng * 0.22d0
    lrbt(4,8)      = x2Leng * 0.28d0
    vMaxMin(1,1,8) = - 0.5d0
    vMaxMin(2,1,8) = + 0.5d0
    vMaxMin(1,2,8) = - 0.2d0
    vMaxMin(2,2,8) = + 0.2d0
    !  -- [1-9] Set 9 -- !
    lrbt(1,9)      = x1Leng * 0.60d0
    lrbt(2,9)      = x1Leng * 0.65d0
    lrbt(3,9)      = x2Leng * 0.65d0
    lrbt(4,9)      = x2Leng * 0.70d0
    vMaxMin(1,1,9) = - 0.5d0
    vMaxMin(2,1,9) = + 0.5d0
    vMaxMin(1,2,9) = - 0.2d0
    vMaxMin(2,2,9) = + 0.2d0
    !  -- [1-10] Set 10 -- !
    lrbt(1,10)      = x1Leng * 0.80d0
    lrbt(2,10)      = x1Leng * 0.85d0
    lrbt(3,10)      = x2Leng * 0.75d0
    lrbt(4,10)      = x2Leng * 0.80d0
    vMaxMin(1,1,10) = - 0.5d0
    vMaxMin(2,1,10) = + 0.5d0
    vMaxMin(1,2,10) = - 0.2d0
    vMaxMin(2,2,10) = + 0.2d0
    !  -- [1-11] Set 11 -- !
    lrbt(1,11)      = x1Leng * 0.72d0
    lrbt(2,11)      = x1Leng * 0.80d0
    lrbt(3,11)      = x2Leng * 0.46d0
    lrbt(4,11)      = x2Leng * 0.54d0
    vMaxMin(1,1,11) = - 0.5d0
    vMaxMin(2,1,11) = + 0.5d0
    vMaxMin(1,2,11) = - 0.2d0
    vMaxMin(2,2,11) = + 0.2d0

    ! --- [1] Sampling  --- !
    do iset=1, Nset
       do k=1, ns
          call VDFfromPart( pxv(:,1:np(k),k), np(k), lrbt(1:4,iset), &
               &            Ndiv, vMaxMin(1:2,k,iset), vdf_s(0:Ndiv+1,0:Ndiv+1,1:3,k) )
       enddo
       call MPI_ALLREDUCE( vdf_s, vdf(0:Ndiv+1,0:Ndiv+1,1:3,1:ns,iset), (Ndiv+2)*(Ndiv+2)*3*ns, &
            &              MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr )
    enddo

    ! --- [2] Write Out --- !
    write(cstep,'(i8.8)') kstep
    filename       = trim(jobdir)//'bin/'//'VDFsample_'//cstep//'.bin'
    open (50,file=trim(filename),form='unformatted')
    write(50) vdf
    close(50)
    if ( flag__parameter ) then
       filename       = trim(jobdir)//'dat/'//'VDFconstant.dat'
       open (50,file=trim(filename),form='formatted')
       write(50,'(a,1x,i12)') "Ndiv", Ndiv
       write(50,'(a,1x,i12)') "Nset", Nset
       do i=1, Nset
          write(50,'(8(e12.5,1x))') lrbt(1:4,i), vMaxMin(1,1,i), vMaxMin(2,1,i), vMaxMin(1,2,i), vMaxMin(2,2,i)
       enddo
       close(50)
       flag__parameter = .false.
    endif
    
    return
  end subroutine VDFsample

  
  subroutine VDFfromPart( pxvh, Ninp, lrbt, Ndiv, vMaxMin, vdf )
    implicit none
    integer         , intent(in)  :: Ninp, Ndiv
    double precision, intent(in)  :: vMaxMin(2), lrbt(4), pxvh(8,Ninp)
    double precision, intent(out) :: vdf(0:Ndiv+1,0:Ndiv+1,3) ! -- Return == vdf( Ndiv+2, 3 ) -- !
    integer                       :: retSize
    double precision, allocatable :: pxvs(:,:)

    ! --- [1]  Preparation  (  Get Size / Extract pt.  ) --- !
    retSize = CountInBox( pxvh, Ninp, lrbt )
    allocate( pxvs(8,retSize) )
    call ExtractFromBox( pxvh, Ninp, lrbt, pxvs, retSize )    
    
    ! --- [2]  Make VDF in LRBT-Box --- !
    call VDFfromALL( pxvs, retSize, Ndiv, vMaxMin, vdf )
    
    return
  end subroutine VDFfromPart
  
  subroutine VDFfromALL( pxvh, Ninp, Ndiv, vMaxMin, vdf )
    implicit none
    integer         , intent(in)  :: Ninp, Ndiv
    double precision, intent(in)  :: vMaxMin(2), pxvh(8,Ninp)
    double precision, intent(out) :: vdf(0:Ndiv+1,0:Ndiv+1,3) ! -- Return == vdf( Ndiv+2, 3 ) -- !
    integer                       :: iv1, iv2, m, drc, v1_, v2_
    double precision              :: dv
    integer, parameter            :: EdgeInclude = 1

    ! --- [1] Initialization --- !
    dv       = ( vMaxMin(2) - vMaxMin(1) ) / dble( Ndiv )
    do drc=1, 3
       do iv2=0, Ndiv+1
          do iv1=0, Ndiv+1
             vdf(iv1,iv2,drc) = 0.d0
          enddo
       enddo
    enddo
    ! --- [2]  VDF Counting  --- !
    !  -- [2-1] v1 -- !
    drc = 1
    v1_ = 4
    v2_ = 5
    do m=1, Ninp
       iv1              = ceiling( ( pxvh(v1_,m) - vMaxMin(1) ) / dv )
       iv2              = ceiling( ( pxvh(v2_,m) - vMaxMin(1) ) / dv )
       iv1              = max( min( iv1, Ndiv+1 ), 0 )
       iv2              = max( min( iv2, Ndiv+1 ), 0 )
       vdf(iv1,iv2,drc) = vdf(iv1,iv2,drc) + pxvh(8,m)
    enddo
    !  -- [2-2] v2 -- !
    drc = 2
    v1_ = 3
    v2_ = 5
    do m=1, Ninp
       iv1              = ceiling( ( pxvh(v1_,m) - vMaxMin(1) ) / dv )
       iv2              = ceiling( ( pxvh(v2_,m) - vMaxMin(1) ) / dv )
       iv1              = max( min( iv1, Ndiv+1 ), 0 )
       iv2              = max( min( iv2, Ndiv+1 ), 0 )
       vdf(iv1,iv2,drc) = vdf(iv1,iv2,drc) + pxvh(8,m)
    enddo
    !  -- [2-3] v3 -- !
    drc = 3
    v1_ = 3
    v2_ = 4
    do m=1, Ninp
       iv1              = ceiling( ( pxvh(v1_,m) - vMaxMin(1) ) / dv )
       iv2              = ceiling( ( pxvh(v2_,m) - vMaxMin(1) ) / dv )
       iv1              = max( min( iv1, Ndiv+1 ), 0 )
       iv2              = max( min( iv2, Ndiv+1 ), 0 )
       vdf(iv1,iv2,drc) = vdf(iv1,iv2,drc) + pxvh(8,m)
    enddo
    
    return
  end subroutine VDFfromALL

  
  Function CountInBox( pxvh, Ninp, lrbt )
    implicit none
    integer         , intent(in)  :: Ninp
    double precision, intent(in)  :: pxvh(8,Ninp), lrbt(4)
    integer                       :: m, CountInBox

    CountInBox = 0
    do m=1, Ninp
       if (   ( pxvh(2,m).ge.lrbt(1) ).and.( pxvh(2,m).lt.lrbt(2) ).and. &
            & ( pxvh(1,m).ge.lrbt(3) ).and.( pxvh(1,m).lt.lrbt(4) ) ) then
          CountInBox = CountInBox + 1
       endif
    enddo
    
    return
  end Function CountInBox

  
  subroutine ExtractFromBox( pxvh, Ninp, lrbt, pxvr, retSize )
    implicit none
    integer         , intent(in)  :: Ninp, retSize
    double precision, intent(in)  :: pxvh(8,Ninp), lrbt(4)
    double precision, intent(out) :: pxvr(8,retSize)
    integer                       :: i, m, count

    count = 0
    do m=1, Ninp
       if (   ( pxvh(2,m).ge.lrbt(1) ).and.( pxvh(2,m).lt.lrbt(2) ).and. &
            & ( pxvh(1,m).ge.lrbt(3) ).and.( pxvh(1,m).lt.lrbt(4) ) ) then

          count = count + 1
          do i=1, 8
             pxvr(i,count) = pxvh(i,m)
          enddo
       endif
    enddo
    if ( count.ne.retSize ) write(6,*) ' [WARNING] retSize and Total # of pt. in Box is incompatible [WARNING] '    
    return
  end subroutine ExtractFromBox

  
end module VDFModule
