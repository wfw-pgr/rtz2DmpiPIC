! =================================================================== !
! === makeVDFMod :: sample VDF from sim. particles                === !
! =================================================================== !
module makeVDFMod
  implicit none
  integer                         :: Ndiv, Nset, commLen
  character(len=100), parameter   :: SetupFile = 'utl/config/vdfSetup.dat'
  double precision  , allocatable :: lrbt(:,:), vMinMax(:,:,:)
  double precision  , allocatable :: vdf(:,:,:,:,:), vdf_s(:,:,:,:)
  logical           , parameter   :: LoadConstFile     = .false.
  logical           , parameter   :: LoadSetupFile     = .true.
  logical                         :: Flag__initialized = .false.
contains


  ! =================================================================== !
  ! === setupVDF :: Initialization for vdf sampling                 === !
  ! =================================================================== !
  subroutine setupVDF
    use constants, only : ns, datDir, myRank
    use variables, only : x1Leng, x2Leng
    implicit none
    include 'mpif.h'
    integer            :: m, ierr
    character(len=100) :: buff, FileName, ConstFile

    ! --------------------------------- !
    ! --- [1]   Load VDF Settings   --- !
    ! --------------------------------- !
    !  -- [1-1] set FileNames       --  !
    ConstFile = trim(datDir) // 'vdfConst.dat'
    if ( LoadConstFile ) FileName = ConstFile
    if ( LoadSetupFile ) FileName = SetupFile
    !  -- [1-2] Load Ndiv & Nset    --  !
    if ( myRank.eq.0   ) then
       open (50,file=trim(FileName),form='formatted')
       read (50,*) buff, Ndiv
       read (50,*) buff, Nset
       read (50,*) buff
    endif
    call MPI_Bcast( Ndiv, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( Nset, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
    !  -- [1-2] Allocation          --  !
    allocate( lrbt (4,Nset), vMinMax(2,ns,Nset)  )
    allocate( vdf  (0:Ndiv+1,0:Ndiv+1,3,ns,Nset) )
    allocate( vdf_s(0:Ndiv+1,0:Ndiv+1,3,ns)      )
    commLen = (Ndiv+2)*(Ndiv+2)*3*ns
    !  -- [1-3] Load LRBT / vMinMax --  !
    if ( myRank.eq.0   ) then
       do m=1, Nset
          read(50,*) lrbt(1:4,m), vMinMax(1,1,m), vMinMax(2,1,m), vMinMax(1,2,m), vMinMax(2,2,m)
       enddo
       close(50)
    end if
    call MPI_Bcast( lrbt   , 4*   Nset, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    call MPI_Bcast( vMinMax, 2*ns*Nset, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
    !  -- [1-4] LRBT Convert        --  !
    if ( LoadSetupFile ) then
       do m=1, Nset
          lrbt(1,m) = x1Leng * lrbt(1,m)
          lrbt(2,m) = x1Leng * lrbt(2,m)
          lrbt(3,m) = x2Leng * lrbt(3,m)
          lrbt(4,m) = x2Leng * lrbt(4,m)
       enddo
    endif
    ! --- [2]   Save Const File     --- !
    if ( LoadSetupFile ) then
       if ( myRank.eq.0 ) then
          open (50,file=trim(ConstFile),form='formatted')
          write(50,*) 'Ndiv', Ndiv
          write(50,*) 'Nset', Nset
          write(50,'(8(a12,1x))') '#LeftEdge', 'RightEdge', 'BotEdge', 'TopEdge',&
               &                  'veMin'    , 'veMax'    , 'viMin'  , 'viMax'
          do m=1, Nset
             write(50,'(8(e12.5,1x))') lrbt(1:4,m), vMinMax(1,1,m), vMinMax(2,1,m), vMinMax(1,2,m), vMinMax(2,2,m)
          enddo
          close(50)
       endif
    endif

    ! --- [3] Flag Management --- !
    Flag__initialized = .true.
    
    return
  end subroutine setupVDF
  

  ! =================================================================== !
  ! === VDFsample :: vdf sampling routine (controler)               === !
  ! =================================================================== !
  subroutine VDFSample
    use constants, only : ns, np, binDir, myRank
    use variables, only : pxv, kstep, x1g
    implicit none
    include 'mpif.h'
    integer            :: k, iset, ierr
    double precision   :: vMinMaxh(2), lrbt_loc(4)
    character(len=  8) :: cStep
    character(len=100) :: FileName, binDirh
    integer, parameter :: x1Len_=1, x1Frm_=2, x1To_=3
    
    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    if ( .not.( Flag__initialized ) ) call setupVDF
    if ( myRank.eq.0 ) then
       write(cStep,'(i8.8)') kstep
       binDirh  = trim(binDir ) // 'kstep'      // cStep // '/'
       FileName = trim(binDirh) // 'VDFsample_' // cStep // '.bin'
       call system( 'mkdir -p ' // trim(binDirh) )
    endif
    ! ------------------------------------- !
    ! --- [2] Sampling                  --- !
    ! ------------------------------------- !
    do iset=1, Nset
    enddo
    do iset=1, Nset
       lrbt_loc(1) = lrbt(1,iset) - x1g(myRank,x1Frm_)
       lrbt_loc(2) = lrbt(2,iset) - x1g(myRank,x1Frm_)
       lrbt_loc(3) = lrbt(3,iset)
       lrbt_loc(4) = lrbt(4,iset)
       do k=1, ns
          vMinMaxh = vMinMax(1:2,k,iset)
          call VDFfromPart( pxv(:,:,k), np(k), lrbt_loc, &
               &            vMinMaxh  , vdf_s(0:Ndiv+1,0:Ndiv+1,1:3,k)   )
       enddo
       call MPI_ALLREDUCE( vdf_s, vdf(0:Ndiv+1,0:Ndiv+1,1:3,1:ns,iset), commLen, &
            &              MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr   )
    enddo
    ! ------------------------------------- !
    ! --- [3] save as binary File       --- !
    ! ------------------------------------- !
    if ( myRank.eq.0 ) then
       open (50,file  =trim(FileName), form   ='unformatted',  &
       &        status='replace'     , convert='LITTLE_ENDIAN' )
       write(50) vdf
       close(50)
    endif
    return
  end subroutine VDFSample


  ! =================================================================== !
  ! === VDFfromPart :: Extract particle & sumup in velocity-space   === !
  ! =================================================================== !
  subroutine VDFfromPart( pxvh, Ninp, lrbt, vMinMax, vdfh )
    implicit none
    integer         , intent(in)  :: Ninp
    double precision, intent(in)  :: vMinMax(2), lrbt(4), pxvh(8,Ninp)
    double precision, intent(out) :: vdfh(0:Ndiv+1,0:Ndiv+1,3) ! -- Return == vdfh( Ndiv+2, 3 ) -- !
    integer                       :: retSize
    double precision, allocatable :: pxvs(:,:)
    
    ! ------------------------------------- !
    ! --- [1] Get Size / Extract pt.    --- !
    ! ------------------------------------- !
    retSize = CountInBox( pxvh, Ninp, lrbt )
    allocate( pxvs(8,retSize) )
    call ExtractFromBox( pxvh, Ninp, lrbt, pxvs, retSize )    
    
    ! ------------------------------------- !
    ! --- [2] Make VDF from pt.         --- !
    ! ------------------------------------- !
    call VDFfromALL( pxvs, retSize, vMinMax, vdfh )
    
    return
  end subroutine VDFfromPart
  
  
  ! =================================================================== !
  ! === VDFfromALL :: Extract particle & sumup in velocity-space   === !
  ! =================================================================== !
  subroutine VDFfromALL( pxvh, Ninp, vMinMax, vdfh )
    implicit none
    integer         , intent(in)  :: Ninp
    double precision, intent(in)  :: vMinMax(2), pxvh(8,Ninp)
    double precision, intent(out) :: vdfh(0:Ndiv+1,0:Ndiv+1,3) ! -- Return == vdf( Ndiv+2, 3 ) -- !
    integer                       :: iv1, iv2, m, drc, v1_, v2_
    double precision              :: dv
    integer, parameter            :: EdgeInclude = 1

    ! ------------------------------------- !
    ! --- [1] Initialization            --- !
    ! ------------------------------------- !
    dv       = ( vMinMax(2) - vMinMax(1) ) / dble( Ndiv )
    do drc=1, 3
       do iv2=0, Ndiv+1
          do iv1=0, Ndiv+1
             vdfh(iv1,iv2,drc) = 0.d0
          enddo
       enddo
    enddo
    ! ------------------------------------- !
    ! --- [2] VDF Counting              --- !
    ! ------------------------------------- !
    !  -- [2-1] v1                      --  !
    drc = 1
    v1_ = 4
    v2_ = 5
    do m=1, Ninp
       iv1               = ceiling( ( pxvh(v1_,m) - vMinMax(1) ) / dv )
       iv2               = ceiling( ( pxvh(v2_,m) - vMinMax(1) ) / dv )
       iv1               = max( min( iv1, Ndiv+1 ), 0 )
       iv2               = max( min( iv2, Ndiv+1 ), 0 )
       vdfh(iv1,iv2,drc) = vdfh(iv1,iv2,drc) + pxvh(8,m)
    enddo
    !  -- [2-2] v2 -- !
    drc = 2
    v1_ = 3
    v2_ = 5
    do m=1, Ninp
       iv1               = ceiling( ( pxvh(v1_,m) - vMinMax(1) ) / dv )
       iv2               = ceiling( ( pxvh(v2_,m) - vMinMax(1) ) / dv )
       iv1               = max( min( iv1, Ndiv+1 ), 0 )
       iv2               = max( min( iv2, Ndiv+1 ), 0 )
       vdfh(iv1,iv2,drc) = vdfh(iv1,iv2,drc) + pxvh(8,m)
    enddo
    !  -- [2-3] v3 -- !
    drc = 3
    v1_ = 3
    v2_ = 4
    do m=1, Ninp
       iv1               = ceiling( ( pxvh(v1_,m) - vMinMax(1) ) / dv )
       iv2               = ceiling( ( pxvh(v2_,m) - vMinMax(1) ) / dv )
       iv1               = max( min( iv1, Ndiv+1 ), 0 )
       iv2               = max( min( iv2, Ndiv+1 ), 0 )
       vdfh(iv1,iv2,drc) = vdfh(iv1,iv2,drc) + pxvh(8,m)
    enddo
    
    return
  end subroutine VDFfromALL


  ! =================================================================== !
  ! === CountInBox :: Count the number of particles inside LRBT     === !
  ! =================================================================== !
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


  ! =================================================================== !
  ! === ExtractFromBox :: pxvh => pxvr if ( pxvh is inside lrbt )   === !
  ! =================================================================== !
  subroutine ExtractFromBox( pxvh, Ninp, lrbt, pxvr, retSize )
    implicit none
    integer         , intent(in)  :: Ninp, retSize
    double precision, intent(in)  :: pxvh(8,Ninp), lrbt(4)
    double precision, intent(out) :: pxvr(8,retSize)
    integer                       :: i, m, count

    count = 0
    do m=1, Ninp
       ! -- IN / OUT -- !
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

  
end module MakeVDFMod
