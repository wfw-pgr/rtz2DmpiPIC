module writingMod
contains

  subroutine WriteField( onGrid )
    use RgridRZMod , only            : LIr, LJr, BgR, fgR, JgR, ugR
    use BgridRZMod , only            : LIg, LJg, BgB, fgB, JgB, ugB
    use constants, only            : job, Bv0
    implicit none
    character(1), intent(in)      :: onGrid
    character(100)                :: FileName, outDir
    double precision, allocatable :: Bwr(:,:,:), Jwr(:,:,:)
    double precision, allocatable :: Pwr(:,:,:), Uwr(:,:,:)

    !  --- [1] Prepare Grid Setting ( we have Only BgR, BgB option, now )   --- !
    if ( onGrid.eq.'R' ) then
       allocate( Bwr(3,LIr,LJr), Jwr(3,LIr,LJr), Pwr(LIr,LJr,3), Uwr(LIr,LJr,3) )
       Bwr(:,:,:) =   BgR(:,:,:)
       Jwr(:,:,:) =   JgR(:,:,:)
       Pwr(:,:,:) =   fgR(:,:,:)
       Uwr(:,:,:) =   ugR(:,:,:)
       outDir     = 'job/'//trim(job)//'/RgD/'
    endif
    if ( onGrid.eq.'B' ) then
       allocate( Bwr(3,LIg,LJg), Jwr(3,LIg,LJg), Pwr(3,LIg,LJg), Uwr(3,LIg,LJg) )
       Bwr(:,:,:) =   BgB(:,:,:)
       Jwr(:,:,:) =   JgB(:,:,:)
       Pwr(:,:,:) =   fgB(:,:,:)
       Uwr(:,:,:) =   ugB(:,:,:)
       outDir     = 'job/'//trim(job)//'/BgD/'
    endif
    write(6,*)
    write(6,'(2x,a)') "----------------------------------------------------------------------"
    write(6,'(2x,a)') '[ WriteField    @ writingMod ]'
    write(6,*)
    

    !  --- [2] Write Variables  --- !
    !   -- [2-1] Bfd.bin        --  !
    FileName = trim(outDir) // 'Bfd.bin'
    open( 50,file  =trim(FileName), form   ='unformatted' &
         &  ,status='replace'     , convert='LITTLE_ENDIAN' )
    write(50) Bwr
    close(50)
    write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
    !   -- [2-2] frp.bin    -- !
    FileName = trim(outDir) // 'frp.bin'
    open( 50,file  =trim(FileName), form   ='unformatted' &
         &  ,status='replace'     , convert='LITTLE_ENDIAN' )
    write(50) Pwr
    close(50)
    write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
    !   -- [2-3] Jcr.bin    -- !
    FileName = trim(outDir) // 'Jcr.bin'
    open( 50,file  =trim(FileName), form   ='unformatted' &
         &  ,status='replace'     , convert='LITTLE_ENDIAN' )
    write(50) Jwr
    close(50)
    write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
    !   -- [2-4] ugB.bin    -- !
    FileName = trim(outDir) // 'uvc.bin'
    open( 50,file  =trim(FileName), form   ='unformatted' &
         &  ,status='replace'     , convert='LITTLE_ENDIAN' )
    write(50) Uwr
    close(50)
    write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
    !   -- [2-5] Bv0.dat    -- !
    FileName = trim(outDir) // 'Bv0.dat'
    open( 50,file=trim(FileName),form='formatted',status='replace' )
    write(50,'(e15.8)') Bv0
    close(50)
    write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )

    deallocate( Bwr, Pwr, Jwr, Uwr )
    write(6,*)
    write(6,'(2x,a)') "----------------------------------------------------------------------"
    write(6,*)

    return
  end subroutine WriteField


  subroutine WriteParam( onGrid )
    use RgridRZMod, only      : LIr, LJr
    use BgridRZMod, only      : LIg, LJg
    use constants , only      : wpewce, vthcv, valfe, TiTe, rhofloor
    use constants , only      : job, N1, N2, PEpic
    use constants , only      : dx1_Debye, dx2_Debye, x1Min, x1Max, x2Min, x2Max
    use variables , only      : dx1, dx2
    implicit none
    integer                  :: LIh, LJh
    character(100)           :: FileName
    character(1), intent(in) :: onGrid

    if ( onGrid.eq.'R' ) then
       LIh   = LIr
       LJh   = LJr
       FileName = 'job/'//trim(job)//'/RgD/parameter.dat'
    endif
    if ( onGrid.eq.'B' ) then
       LIh   = LIg
       LJh   = LJg
       FileName = 'job/'//trim(job)//'/BgD/parameter.dat'
    endif
    
    open (30, file=trim(FileName), status='replace', form='formatted' )
    write(30,'(2(a12,2x),i15)'  ) 'LI'      ,'long'  , LIh
    write(30,'(2(a12,2x),i15)'  ) 'LJ'      ,'long'  , LJh
    write(30,'(2(a12,2x),e15.8)') 'dz'      ,'double', dx1
    write(30,'(2(a12,2x),e15.8)') 'dr'      ,'double', dx2
    write(30,'(2(a12,2x),e15.8)') 'dz_Debye','double', dx1_Debye
    write(30,'(2(a12,2x),e15.8)') 'dr_Debye','double', dx2_Debye
    write(30,'(2(a12,2x),e15.8)') 'wpewce'  ,'double', wpewce
    write(30,'(2(a12,2x),e15.8)') 'vthcv'   ,'double', vthcv
    write(30,'(2(a12,2x),e15.8)') 'valfe'   ,'double', valfe
    write(30,'(2(a12,2x),e15.8)') 'TiTe'    ,'double', TiTe
    write(30,'(2(a12,2x),e15.8)') 'rhofloor','double', rhofloor
    write(30,'(2(a12,2x),e15.8)') 'x1Min'   ,'double', x1Min
    write(30,'(2(a12,2x),e15.8)') 'x1Max'   ,'double', x1Max
    write(30,'(2(a12,2x),e15.8)') 'x2Min'   ,'double', x2Min
    write(30,'(2(a12,2x),e15.8)') 'x2Max'   ,'double', x2Max
    write(30,'(2(a12,2x),i15)'  ) 'PEtot'   ,'long'  , PEpic
    close(30)
    write(6,'(2x,a)') '[ WriteParam    @ writingMod ]'
    write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
    write(6,*)
    
    return
  end subroutine WriteParam
  

  ! =================================================================== !
  ! ===  writeF_MPI :: write Domain Decomposited Field ( for MPI )  === !
  ! =================================================================== !
  subroutine WriteF_MPI( onGrid )
    use RgridRZMod, only           : LIr, LJr, BgR, fgR, JgR, ugR
    use BgridRZMod, only           : LIg, LJg, BgB, fgB, JgB, ugB
    use constants , only           : job, Bv0, PEpic
    implicit none
    integer                       :: i, j, ip, jp, iPE, dir
    integer                       :: LIh, LJh, LIloc, LJloc, surplus
    integer         , parameter   :: rnk_=1, frm_ = 2, to_ = 3, LIs_ = 4, LJs_ = 5
    integer         , parameter   :: iFr_=6, iTo_ = 7, inn_= 8, isl_ = 9, iel_ = 10
    character(  1), intent(in)    :: onGrid
    character(100)                :: FileName, outDir
    character(  6)                :: cRank
    integer         , allocatable :: ijDomain(:,:)
    double precision, allocatable :: Bwr(:,:,:), Jwr(:,:,:)
    double precision, allocatable :: Pwr(:,:,:), Uwr(:,:,:)
    double precision, allocatable :: Bpt(:,:,:), Jpt(:,:,:)
    double precision, allocatable :: Ppt(:,:,:), Upt(:,:,:)
    

    ! -------------------------- !
    ! --- [1]  Grid Select   --- !
    ! -------------------------- !
    !  -- [1-1] Prepare Grid Setting ( we have Only BgR, BgB option, now ) -- !
    if ( onGrid.eq.'R' ) then
       allocate( Bwr(3,LIr,LJr), Jwr(3,LIr,LJr), Pwr(LIr,LJr,3), Uwr(LIr,LJr,3) )
       Bwr(:,:,:) =   BgR(:,:,:)
       Jwr(:,:,:) =   JgR(:,:,:)
       Pwr(:,:,:) =   fgR(:,:,:)
       Uwr(:,:,:) =   ugR(:,:,:)
       outDir     = 'job/'//trim(job)//'/RgD/'
       LIh        = LIr
       LJh        = LJr
    endif
    if ( onGrid.eq.'B' ) then
       allocate( Bwr(3,LIg,LJg), Jwr(3,LIg,LJg), Pwr(3,LIg,LJg), Uwr(3,LIg,LJg) )
       Bwr(:,:,:) =   BgB(:,:,:)
       Jwr(:,:,:) =   JgB(:,:,:)
       Pwr(:,:,:) =   fgB(:,:,:)
       Uwr(:,:,:) =   ugB(:,:,:)
       outDir     = 'job/'//trim(job)//'/BgD/'
       LIh        = LIg
       LJh        = LJg
    endif
    write(6,*)
    write(6,'(2x,a)') "----------------------------------------------------------------------"
    write(6,'(2x,a)') '[ WriteF_MPI    @ writingMod ]'
    write(6,*)

    ! -------------------------- !
    ! --- [2]   Partition    --- !
    ! -------------------------- !
    !  -- [2-1] Allocation                               -- !
    allocate( ijDomain( 0:PEpic-1, 10 ) )
    !  -- [2-2] LIs partition (Actual Data zone [2:LIh]) -- !
    LIloc   = LIh / PEpic
    surplus = LIh - LIloc * PEpic
    do iPE=0, PEpic-1
       ijDomain(iPE,rnk_)  = iPE
       ijDomain(iPE,inn_)  = LIloc
    enddo
    do iPE=0, surplus-1
       ijDomain(iPE,inn_)  = ijDomain(iPE,inn_) + 1
    enddo
    do iPE=0, PEpic-1
       ijDomain(iPE,LIs_)  = ijDomain(iPE,inn_) + 2
       ijDomain(iPE,LJs_)  = LJh
    enddo
    ijDomain(      0,LIs_) = ijDomain(      0,inn_) + 1
    ijDomain(PEpic-1,LIs_) = ijDomain(PEpic-1,inn_) + 1
    !  -- [2-3] End position :: From & To -- !
    ijDomain(0,iTo_) = ijDomain(0,inn_)
    do iPE=1, PEpic-1
       ijDomain(iPE,iTo_)  = ijDomain(iPE-1,iTo_) + ijDomain(iPE,inn_)
    enddo
    ijDomain(0,iFr_) = 1
    do iPE=1, PEpic-1
       ijDomain(iPE,iFr_)  = ijDomain(iPE-1,iTo_) + 1
    enddo
    do iPE=0, PEpic-1
       ijDomain(iPE,frm_)  = ijDomain(iPE,iFr_) - 1
       ijDomain(iPE,to_ )  = ijDomain(iPE,iTo_) + 1
    enddo
    ijDomain(      0,frm_) = 1
    ijDomain(PEpic-1,to_ ) = ijDomain(PEpic-1,iTo_)
    do iPE=0, PEpic-1
       ijDomain(iPE,isl_)  = 2
       ijDomain(iPE,iel_)  = ijDomain(iPE,LIs_) - 1
    enddo
    ijDomain(      0,isl_) = 1
    ijDomain(PEpic-1,iel_) = ijDomain(PEpic-1,LIs_)
    
    ! ---------------------------- !
    ! --- [3] Write Variables  --- !
    ! ---------------------------- !
    do iPE=0, PEpic-1
       !   -- [3-0] Preparation -- !
       write(cRank,'(i6.6)') iPE
       LIloc = ijDomain(iPE,LIs_)
       LJloc = ijDomain(iPE,LJs_)
       allocate( Bpt(3,LIloc,LJloc), Ppt(3,LIloc,LJloc), &
            &    Jpt(3,LIloc,LJloc), Upt(3,LIloc,LJloc)  )
       Bpt(:,:,:) = 0.d0
       Ppt(:,:,:) = 0.d0
       Jpt(:,:,:) = 0.d0
       Upt(:,:,:) = 0.d0
       do dir=1, 3
          do j=1, LJloc
             do i=1, LIloc
                ! - [ ip=From(iPE):@i=2, ip=To(iPE):@i=LIloc ] - !
                ip = ijDomain(iPE,frm_) + ( i-1 )
                jp = j
                Bpt(dir,i,j) = Bwr(dir,ip,jp)
                Jpt(dir,i,j) = Jwr(dir,ip,jp)
                Ppt(dir,i,j) = Pwr(dir,ip,jp)
                Upt(dir,i,j) = Uwr(dir,ip,jp)
             enddo
          enddo
       enddo
       !   -- [3-1] Bfd.bin     -- !
       FileName = trim(outDir) // "Bfd_" // cRank // ".bin"
       open( 50,file  =trim(FileName), form   ='unformatted' &
            &  ,status='replace'     , convert='LITTLE_ENDIAN' )
       write(50) Bpt
       close(50)
       write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
       !   -- [3-2] frp.bin     -- !
       FileName = trim(outDir) // "frp_" // cRank // ".bin"
       open( 50,file  =trim(FileName), form   ='unformatted' &
            &  ,status='replace'     , convert='LITTLE_ENDIAN' )
       write(50) Ppt
       close(50)
       write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
       !   -- [3-3] Jcr.bin     -- !
       FileName = trim(outDir) // "Jcr_" // cRank // ".bin"
       open( 50,file  =trim(FileName), form   ='unformatted' &
            &  ,status='replace'     , convert='LITTLE_ENDIAN' )
       write(50) Jpt
       close(50)
       write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
       !   -- [3-4] ugB.bin     -- !
       FileName = trim(outDir) // "uvc_" // cRank // ".bin"
       open( 50,file  =trim(FileName), form   ='unformatted' &
            &  ,status='replace'     , convert='LITTLE_ENDIAN' )
       write(50) Upt
       close(50)
       write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
       !   -- [3-5] PostProcess -- !
       deallocate( Bpt, Ppt, Jpt, Upt )
    enddo
    ! ---------------------------- !
    ! --- [4]  Write Bv0.dat   --- !
    ! ---------------------------- !
    FileName = trim(outDir) // 'Bv0.dat'
    open( 50,file=trim(FileName),form='formatted',status='replace' )
    write(50,'(e15.8)') Bv0
    close(50)
    write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
    
    ! ---------------------------- !
    ! --- [5] MPI config File  --- !
    ! ---------------------------- !
    FileName = trim(outDir) // 'ijDomain.dat'
    open( 50,file=trim(FileName),form='formatted',status='replace' )
    write(50,'(10(a12,1x))') '#.PE', 'From', 'To', 'LIs', 'LJs', 'iFr', 'iTo', 'inn', 'isl', 'iel'
    do iPE=0, PEpic-1
       write(50,'(10(i12,1x))') ijDomain(iPE,rnk_), ijDomain(iPE,frm_), ijDomain(iPE,to_ ), &
            &                   ijDomain(iPE,LIs_), ijDomain(iPE,LJs_), ijDomain(iPE,iFr_), ijDomain(iPE,iTo_), &
            &                   ijDomain(iPE,inn_), ijDomain(iPE,isl_), ijDomain(iPE,iel_)
    enddo
    close(50)
    write(6,'(a35,1x,a)') trim( FileName ), '   has been outputted...'
    write(6,*)
    write(6,*)

    deallocate( Bwr, Pwr, Jwr, Uwr )
    write(6,'(2x,a)') "----------------------------------------------------------------------"
    return
  end subroutine WriteF_MPI

  
end module writingMod


    ! allocate( FromTo( 0:PEpic-1, 4 ) )
    ! !  -- [2-2] LIs partition (Actual Data zone [2:LIh]) -- !
    ! LIloc   = (LIh-2) / PEpic
    ! surplus = (LIh-2) - LIloc * PEpic
    ! do iPE=0, PEpic-1
    !    FromTo(iPE,LIs_) = LIloc
    ! enddo
    ! do iPE=0, surplus-1
    !    FromTo(iPE,LIs_) = FromTo(iPE,LIs_) + 1
    ! enddo
    ! do iPE=0, PEpic-1
    !    FromTo(iPE,LIs_) = FromTo(iPE,LIs_) + 2
    !    FromTo(iPE,LJs_) = LJh
    ! enddo
    ! !  -- [2-3] End position :: From & To -- !
    ! FromTo(0,frm_) = 1
    ! FromTo(0,to_ ) = FromTo(0,frm_) + FromTo(0,LIs_) - 1
    ! do iPE=1, PEpic-1
    !    FromTo(iPE,frm_) = FromTo(iPE-1,to_ ) - 1
    !    FromTo(iPE,to_ ) = FromTo(iPE  ,frm_) + FromTo(iPE,LIs_) - 1
    ! enddo
    ! !   -- [5] MPI config File -- !
    ! FileName = trim(outDir) // 'PEinfo.dat'
    ! open( 50,file=trim(FileName),form='formatted',status='replace' )
    ! write(50,'(5(a12,2x))') '# PE', 'From', 'To', 'LIloc', 'LJloc'
    ! do iPE=0, PEpic-1
    !    write(50,'(5(i12,2x))') iPE, FromTo(iPE,frm_), FromTo(iPE,to_), FromTo(iPE,LIs_), FromTo(iPE,LJs_)
    ! enddo
    ! close(50)
    ! write(6,'(5x,a,2x,a)') '* SAVE', trim( FileName )
    ! write(6,*)
    ! write(6,*)
