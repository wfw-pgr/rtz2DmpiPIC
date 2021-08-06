module com4DDDMod
  use constants, only             : LI, ns, PEtot
  use constants, only             : DDDnLine_Max, DDDcritStep, DDDcolThreshold
  use constants, only             : rnk_, frm_, to_, LIs_, LJs_, iFr_, iTo_, inn_, isl_, iel_
  implicit none
  integer         , parameter    :: ebCmp              = 6
  integer         , parameter    :: ptCmp              = 8
  integer         , parameter    :: neibPEtot_Max      = 2
  integer         , parameter    :: nBuffLine_Max = neibPEtot_Max * DDDnLine_Max
  logical         , parameter    :: Flag__Confirmation = .true.
  integer                        :: ijDomain_(0:PEtot-1,10)
  integer                        :: nDst, nSrc, LIRemnant, NewLIs, nColumnExport, nColumnThreshold
  integer                        :: export_LItem, export_RItem, import_LItem, import_RItem
  integer                        :: export_Items(0:NeibPEtot_Max)   , export_Index(0:NeibPEtot_Max)
  integer                        :: export_iFrom(  NeibPEtot_Max)   , export_Enemy(  NeibPEtot_Max)
  integer                        :: import_Items(0:NeibPEtot_Max)   , import_Index(0:NeibPEtot_Max)
  integer                        :: import_gFrom(  NeibPEtot_Max)   , import_Enemy(  NeibPEtot_Max)
  integer                        :: import_ebIdx(0:NeibPEtot_Max)
  integer                        :: export_ptItm(0:NeibPEtot_Max,ns), export_ptIdx(0:NeibPEtot_Max,ns)
  integer                        :: import_ptItm(0:NeibPEtot_Max,ns), import_ptIdx(0:NeibPEtot_Max,ns)
  integer                        :: export_ptFrm(0:NeibPEtot_Max,ns), import_ptDst(0:NeibPEtot_Max,ns)
  integer                        :: export_npkPE(ns)                , import_npkPE(ns)
  integer         , allocatable  :: request_Send(:), status_Send(:,:)
  integer         , allocatable  :: request_Recv(:), status_Recv(:,:)
  double precision, allocatable  ::  ptSendBuff(:,:,:),  ptRecvBuff(:,:,:)
  double precision, allocatable  :: EBfSendBuff(:,:,:), EBfRecvBuff(:,:,:)
  double precision, allocatable  ::  FieldStack(:,:,:)
  double precision               :: zoffset(2)
contains

  ! =================================================================== !
  ! ===  reAllocateVariables  :: Allocate & reInitialize variables  === !
  ! =================================================================== !
  subroutine reAllocateVariables
    use variables , only : EBr, EBo, Jcr, JcrW
    use variables , only : w_p, cellTop, perCell, sumupW, sumupP, sumupE, perCellW
    use momentsMod, only : fMoments
    use allocatMod, only : allocateEBrEBo, allocateCurrent
    use allocatMod, only : allocateOthers, allocateMoments
    implicit none

    ! --- [1] Re-Allocation      --- !
    !  -- [1-1]    Field         --  !
    deallocate( EBo, EBr             )
    call allocateEBrEBo
    !  -- [1-2]    Current       --  !
    deallocate( Jcr, JcrW            )
    call allocateCurrent
    !  -- [1-3]    Moment        --  !
    deallocate( fMoments, sumupW, sumupP, sumupE )
    call allocateMoments
    !  -- [1-4]    Others        --  !
    deallocate( cellTop, perCell, perCellW, w_p  )
    call allocateOthers
    
    return
  end subroutine reAllocateVariables


  ! =================================================================== !
  ! ===  deAllocateVariables  :: deAllocate DDD variables           === !
  ! =================================================================== !
  subroutine deallocateVariables
    implicit none
    deallocate(   ptSendBuff, ptRecvBuff  )
    deallocate( request_Send, status_Send )
    deallocate( request_Recv, status_Recv )
    deallocate(  EBfSendBuff, EBfRecvBuff )
    deallocate(  FieldStack               )
    return
  end subroutine deallocateVariables


  ! =================================================================== !
  ! ===  reDefine_x1Geometry  :: reDefine  x1s, x1g                 === !
  ! =================================================================== !
  subroutine reDefine_x1Geometry
    use constants, only : LIs, dz , myRank, PEtot
    use variables, only : x1g, x1s, ijDomain, x1Lengloc
    implicit none
    integer            :: i, iPE
    integer, parameter :: LIs_  =4, iFr_  =5
    integer, parameter :: x1Len_=1, x1Frm_=2, x1To_=3

    ! --- [1] x1Lengloc  --- !
    x1Lengloc = dble(LIs-2) * dz

    ! --- [2] x1g / x1s  --- !
    do iPE=0, PEtot-1
       x1g(iPE,x1Len_) = dble( ijDomain(iPE,LIs_) - 2 ) * dz
       x1g(iPE,x1Frm_) = dble( ijDomain(iPE,iFr_) - 2 ) * dz
       x1g(iPE,x1To_ ) = x1g(iPE,x1Frm_) + x1g(iPE,x1Len_)
    enddo
    x1g(0,x1Frm_) = 0.d0
    deallocate( x1s      )
    allocate  ( x1s(LIs) )
    do i=1, LIs
       x1s(i) = dz * dble(i-2) + x1g(myRank,x1Frm_) + 0.5d0*dz
    enddo
    
    return
  end subroutine reDefine_x1Geometry


  ! =================================================================== !
  ! ===  reDefine_CommConfig  :: reDefine communication config      === !
  ! =================================================================== !
  subroutine reDefine_CommConfig
    use constants , only : LIs, LJs
    use constants , only : Boundary1__EB, Boundary1__jc
    use pMPIComMod, only : setExchangeGeometry
    use fMPIComMod, only : setupExchangeTable_Field
    use jMPIComMod, only : setupExchangeTable_Current
    use rMPIComMod, only : setupExchangeTable_Relocated
    implicit none

    ! ---  Re-setup MPI config  --- !
    call setExchangeGeometry
    call setupExchangeTable_Field    ( LIs, LJs, Boundary1__EB )
    call setupExchangeTable_Current  ( LIs, LJs, Boundary1__jc )
    call setupExchangeTable_Relocated( LIs, LJs, Boundary1__EB )

    return
  end subroutine reDefine_CommConfig
  
  
  subroutine BcastDisplayNptInfo
    use constants, only : LIs, np, npc, myRank, PEtot
    use variables, only : npgPE
    implicit none
    include 'mpif.h'
    integer            :: j, iPE, ierr
    integer            :: sBuff(9), rBuff(9,0:PEtot-1), sumBuff(9)

    ! --- [1]  BroadCast np --- !
    sBuff(1) = export_LItem
    sBuff(2) = export_RItem
    sBuff(3) = import_LItem
    sBuff(4) = import_RItem
    sBuff(5) = np (1)
    sBuff(6) = np (2)
    sBuff(7) = npc(1)
    sBuff(8) = npc(2)
    sBuff(9) = LIs
    call MPI_AllGather( sBuff, 9, MPI_INTEGER, &
         &              rBuff, 9, MPI_INTEGER, &
         &              MPI_COMM_WORLD, ierr )
    sumBuff  = 0
    do iPE=0, PEtot-1
       do j=1, 4
          npgPE(iPE,j) = rBuff(j+4,iPE)
       enddo
       do j=1, 8
          sumBuff(j)   = sumBuff(j) + rBuff(j,iPE)
       enddo
       npgPE(iPE,5) = npgPE(iPE,3) + npgPE(iPE,4)
       sumBuff(9)   = sumBuff(9)   + npgPE(iPE,5)
    enddo
    ! --- [2]  Display      --- !
    if ( ( myRank.eq.0 ).and.( Flag__Confirmation ) ) then
       write(6,*)
       write(6,'(2x,a)'                 ) '[ BcastDisplayNptInfo @ diffDDDMod ] '
       write(6,'(2x,6(a4,1x),3(a12,1x))') 'PE#.', 'ex.L', 'ex.R', 'im.L', 'im.R', 'LIs', 'np(el_)', 'np(io_)', 'npt'
       write(6,'(2x,a)'                 ) '----------------------------------------------------------------------'
       do iPE=0, PEtot-1
          write(6,'(2x,6(i4,1x),3(i12,1x))') iPE, &
               &        rBuff(1,iPE), rBuff(2,iPE), rBuff(3,iPE), rBuff(4,iPE), rBuff(9,iPE), &
               &        npgPE(iPE,1), npgPE(iPE,2), npgPE(iPE,5)
       enddo
       write(6,'(2x,a)'               ) '----------------------------------------------------------------------'
       write(6,'(2x,a7,23x,3(i12,1x))') '[Total]', sumBuff(7), sumBuff(8), sumBuff(9)
       write(6,'(2x,a)'               ) '----------------------------------------------------------------------'
       write(6,*)
    endif
    
    return
  end subroutine BcastDisplayNptInfo


  ! =================================================================== !
  ! ===  UpdateijDomain  :: ijDomain_ => ijDomain                   === !
  ! =================================================================== !
  subroutine UpdateijDomain
    use constants, only : datDir, myRank, PEtot
    use variables, only : ijDomain, kstep
    implicit none
    include 'mpif.h'
    integer            :: iPE
    character(  8)     :: cStep
    character(100)     :: FileName

    ! ------------------------------------- !
    ! --- [1] Update ijDomain           --- !
    ! ------------------------------------- !
    do iPE=0, PEtot-1
       ijDomain(iPE,rnk_)  = ijDomain_(iPE,rnk_)
       ijDomain(iPE,frm_)  = ijDomain_(iPE,frm_)
       ijDomain(iPE,to_)   = ijDomain_(iPE,to_ )
       ijDomain(iPE,LIs_)  = ijDomain_(iPE,LIs_)
       ijDomain(iPE,LJs_)  = ijDomain_(iPE,LJs_)
       ijDomain(iPE,iFr_)  = ijDomain_(iPE,iFr_)
       ijDomain(iPE,iTo_)  = ijDomain_(iPE,iTo_)
       ijDomain(iPE,inn_)  = ijDomain_(iPE,inn_)
       ijDomain(iPE,isl_)  = ijDomain_(iPE,isl_)
       ijDomain(iPE,iel_)  = ijDomain_(iPE,iel_)
    enddo
    ! ------------------------------------- !
    ! --- [2] Write out                 --- !
    ! ------------------------------------- !
    if ( myRank.eq.0 ) then
       write(cStep,'(i8.8)') kstep
       FileName = trim(datDir) // 'ijDomain/' // 'ijDomain' // cStep //  '.dat'
       open( 50,file=trim(FileName),form='formatted',status='replace' )
       write(50,'(10(a12,1x))') '# PE', 'From', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
       do iPE=0, PEtot-1
          write(50,'(10(i12,1x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
               &                   ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
               &                   ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
       enddo
       close(50)
       if ( Flag__Confirmation ) then
          write(6,'(2x,a)') '----------------------------------------------------------------------'
          write(6,'(2x,a10,a45,2x,a6)') "* SAVE :: ", trim(FileName), '[ OK ]'
          write(6,*)
       endif
    endif
    return
  end subroutine UpdateijDomain
  

end module com4DDDMod
