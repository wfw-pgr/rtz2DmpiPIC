module saveFieldMod
contains
  
  ! ====================================================== !
  ! === save__map2D_asPointFile                        === !
  ! ====================================================== !
  subroutine save__map2D_asPointFile( Data, axis1, axis2, LI, LJ, nComponents, outFile )
    implicit none
    integer         , parameter  :: cLen = 300
    integer         , parameter  :: lun  = 50
    integer         , intent(in) :: LI, LJ, nComponents
    double precision, intent(in) :: Data(nComponents,LI,LJ), axis1(LI), axis2(LJ)
    character(*)    , intent(in) :: outFile
    integer                      :: i, j, k
    character(cLen)              :: num, names, shape1, shape2, fmt

    ! ------------------------------------------------------ !
    ! --- [1] prepare header                             --- !
    ! ------------------------------------------------------ !
    names  = "#"
    do k=1, nComponents+2
       write(num,"(i10)") k
       names = trim(names) // " x" // adjustL(num)
    enddo
    write(shape1,"((i10) (i10))"      ) LI*LJ,  nComponents+2
    write(shape2,"((i10) (i10) (i10))") LJ, LI, nComponents+2
    shape1 = "# " // trim(shape1)
    shape2 = "# " // trim(shape2)
    ! ------------------------------------------------------ !
    ! --- [2] prepare format                             --- !
    ! ------------------------------------------------------ !
    write(num,"(i10)") nComponents+2
    fmt    = "(" // trim(num) // "(e15.8,1x))"
    ! ------------------------------------------------------ !
    ! --- [2] save point Data file                       --- !
    ! ------------------------------------------------------ !
    open( lun,file=trim(outFile),form='formatted', status="replace" )
    write(lun,"(a)") trim(names)
    write(lun,"(a)") trim(shape1)
    write(lun,"(a)") trim(shape2)
    do j=1, LJ
       do i=1, LI
          write(lun,trim(fmt)) axis1(i), axis2(j), ( Data(k,i,j), k=1, nComponents )
       enddo
    enddo
    close(lun)
    write(6,'(a20,a25,1x,a)') "[save__pointFile] ", trim( outFile ), '  has been outputted...'

    return
  end subroutine save__map2D_asPointFile


end module saveFieldMod
