module utilityMod
contains

  ! =================================================================== !
  ! ===  makeJobFile  ::  arange job and its sub Directories        === !
  ! =================================================================== !
  subroutine makeJobFile
    use constants, only : job, myRank
    implicit none
    character(200)     :: command, jobdir

    ! ------------------------------------- !
    ! --- [1] Make job Directory        --- !
    ! ------------------------------------- !
    jobdir  = 'job/' // trim(job) // '/'
    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,'(2x,a)'     ) "----------------------------------------------------------------------"
       write(6,'(2x,a)'     ) '[ makeJobFile    @initialMod ]'
       write(6,'(5x,a,1x,a)') '* SAVE Directory :: ', trim(jobDir)
       write(6,*)
       ! --- [2] make job dir    --- !
       command = 'if [ ! -d ' // trim(jobdir) // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)//  ' ); fi'
       call system( trim(command) )
       !  -- [2-1] Bgrid dir     --  !
       command = 'if [ ! -d '//trim(jobdir)// 'BgD/' // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// 'BgD/' // ' ); fi'
       call system( trim(command) )
       !  -- [2-2] Egrid dir     --  !
       command = 'if [ ! -d '//trim(jobdir)// 'Egs/' // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// 'Egs/' // ' ); fi'
       call system( trim(command) )
       !  -- [2-3] Rgrid dir     --  !
       command = 'if [ ! -d '//trim(jobdir)// 'RgD/' // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// 'RgD/' // ' ); fi'
       call system( trim(command) )
       !  -- [2-3] PDF   dir     --  !
       command = 'if [ ! -d '//trim(jobdir)// 'pdf/' // ' ]; then'
       command = trim(command) // ' ( mkdir ' //trim(jobdir)// 'pdf/' // ' ); fi'
       call system( trim(command) )
       !  -- [2-4] PDF   dir     --  !
       command = 'cp -r '//'src '//trim(jobdir)
       call system( trim(command) )
       !  -- [2-5] src Copy      --  !
       call system( "cp -r src " // trim(jobdir) // "BgD/" )
    endif
    return
  end subroutine makeJobFile

end module utilityMod

