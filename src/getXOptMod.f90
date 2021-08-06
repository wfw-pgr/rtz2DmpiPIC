module getXOptMod
contains
  
  subroutine solveXOpt( psi, ptOXO, LIp, LJp, Flag__sOpt, Flag__FixedXpt )
    implicit none
    integer         , intent(in)    :: LIp, LJp
    double precision, intent(in)    :: psi(LIp,LJp)
    integer         , intent(inout) :: ptOXO(3,2)
    logical         , intent(inout) :: Flag__FixedXpt, Flag__sOpt
    integer                         :: i, j, iX, jX, iO1, jO1, iO2, jO2
    double precision                :: psiX, psiO1, psiO2
    integer         , parameter     :: id_=1, jd_=2

    ! ---------------------- !
    ! --- Search O-point --- !
    ! ---------------------- !
    iX    = ptOXO(2,id_)
    psiX  = psi( ptOXO(2,id_), ptOXO(2,jd_) )
    psiO1 = -1.d8 ; psiO2 = -1.d8
    iO1   = -1    ; jO1   = -1
    iO2   = -1    ; jO2   = -1
    do j=2, LJp
       do i=2, iX-1
          if ( psi(i,j).gt.psiO1 ) then
             psiO1 = psi(i,j)
             iO1   = i
             jO1   = j
          endif
       enddo
    enddo
    do j=2, LJp
       do i=iX+1, LIp
          if ( psi(i,j).gt.psiO2 ) then
             psiO2 = psi(i,j)
             iO2   = i
             jO2   = j
          endif
       enddo
    enddo
    ! ------------------------ !
    ! --- Xpt exist or Not --- !
    ! ------------------------ !
    if ( ( psiO1.gt.psiX ).and.( psiO2.gt.psiX ) ) then
       ! -- Xpt Exist case -- !
       Flag__FixedXpt  = .false.
       Flag__sOpt      = .false.
       ptOXO(1,id_)    = iO1
       ptOXO(1,jd_)    = jO1
       ptOXO(3,id_)    = iO2
       ptOXO(3,jd_)    = jO2
       call solveXpt_ij( psi, ptOXO, LIp, LJp, Flag__FixedXpt )
    else
       ! --  No Xpt case  -- !
       Flag__FixedXpt  = .false.
       Flag__sOpt      = .true.
       if ( ( psiO1.gt.psiX ).and.( psiO2.le.psiX ) ) then
          ptOXO(:,id_) = iO1
          ptOXO(:,jd_) = jO1
       endif
       if ( ( psiO1.le.psiX ).and.( psiO2.gt.psiX ) ) then
          ptOXO(:,id_) = iO2
          ptOXO(:,jd_) = jO2
       endif
       if ( ( psiO1.le.psiX ).and.( psiO2.le.psiX ) ) then
          call solveXpt_j( psi, iX, jX, jO1, jO2, LIp, LJp )
          ptOXO(:,id_) = iX
          ptOXO(:,jd_) = jX
       endif
    endif
    ! ------------------------ !
    ! ---   Error Check    --- !
    ! ------------------------ !
    !  -- Error in i-Range --  !
    if ( ( ptOXO(1,id_).lt.3 ).or.( ptOXO(1,id_).gt.LIp-2 ) ) then
       write(6,*) "[solveXOpt] ptOXO Error :: iO1 == ", ptOXO(1,id_)
       ptOXO(1,id_) = LIp / 2
    endif
    if ( ( ptOXO(2,id_).lt.3 ).or.( ptOXO(2,id_).gt.LIp-2 ) ) then
       write(6,*) "[solveXOpt] ptOXO Error :: iX  == ", ptOXO(2,id_)
       ptOXO(2,id_) = LIp / 2
    endif
    if ( ( ptOXO(3,id_).lt.3 ).or.( ptOXO(3,id_).gt.LIp-2 ) ) then
       write(6,*) "[solveXOpt] ptOXO Error :: iO2 == ", ptOXO(3,id_)
       ptOXO(3,id_) = LIp / 2
    endif
    !  -- Error in j-Range --  !
    if ( ( ptOXO(1,jd_).lt.3 ).or.( ptOXO(1,jd_).gt.LJp-2 ) ) then
       write(6,*) "[solveXOpt] ptOXO Error :: jO1 == ", ptOXO(1,jd_)
       ptOXO(1,jd_) = LJp / 2
    endif
    if ( ( ptOXO(2,jd_).lt.3 ).or.( ptOXO(2,jd_).gt.LJp-2 ) ) then
       write(6,*) "[solveXOpt] ptOXO Error :: jX  == ", ptOXO(2,jd_)
       ptOXO(2,jd_) = LJp / 2
    endif
    if ( ( ptOXO(3,jd_).lt.3 ).or.( ptOXO(3,jd_).gt.LJp-2 ) ) then
       write(6,*) "[solveXOpt] ptOXO Error :: jO2 == ", ptOXO(3,jd_)
       ptOXO(3,jd_) = LJp / 2
    endif

    return
  end subroutine solveXOpt


  subroutine solveXpt_ij( psi, ptOXO, LIp, LJp, Flag__FixedXpt )
    implicit none
    integer         , intent(in)    :: LIp, LJp
    integer         , intent(inout) :: ptOXO(3,2)
    double precision, intent(in)    :: psi(LIp,LJp)
    logical         , intent(inout) :: Flag__FixedXpt
    integer                         :: iter
    integer                         :: iX, iO1, iO2, iX_
    integer                         :: jX, jO1, jO2, jX_
    integer         , parameter     :: iterMax = 30
    integer         , parameter     :: id_=1, jd_=2

    ! -- Initialize -- !
    iO1 = ptOXO(1,id_) ; jO1 = ptOXO(1,jd_)
    iX  = ptOXO(2,id_) ; jX  = ptOXO(2,jd_)
    iO2 = ptOXO(3,id_) ; jO2 = ptOXO(3,jd_)
    do iter=1, iterMax
       ! ---------- !
       ! -- LOOP -- !
       ! ---------- !
       iX_ = iX
       jX_ = jX
       ! -- i-direction -- !
       call solveXpt_i( psi, iX, jX, iO1, iO2, LIp, LJp )
       ! -- j-direction -- !
       call solveXpt_j( psi, iX, jX, jO1, jO2, LIp, LJp )
       if ( ( iX_.eq.iX ).and.( jX_.eq.jX ) ) exit
    enddo
    ! -- Check -- !
    call CheckXpt( psi, iX, jX, iO1, iO2, jO1, jO2, LIp, LJp, Flag__FixedXpt )

    ! -- Insert Back -- !
    ptOXO(1,id_) = iO1 ; ptOXO(1,jd_) = jO1
    ptOXO(2,id_) = iX  ; ptOXO(2,jd_) = jX
    ptOXO(3,id_) = iO2 ; ptOXO(3,jd_) = jO2
    return
  end subroutine solveXpt_ij
  
  
  subroutine solveXpt_i( psi, iX, jX, iO1, iO2, LIp, LJp )
    implicit none
    integer         , intent(in)    :: LIp, LJp, jX, iO1, iO2
    integer         , intent(inout) :: iX
    double precision, intent(in)    :: psi(LIp,LJp)
    integer                         :: i
    double precision                :: psiMin

    iX      = ( iO1+iO2 ) / 2
    psiMin  = +1.d8
    do i=iO1+2, iO2-2
       if ( psi(i,jX).lt.psiMin ) then
          psiMin = psi(i,jX)
          iX     = i
       endif
    enddo
    return
  end subroutine solveXpt_i


  subroutine solveXpt_j( psi, iX, jX, jO1, jO2, LIp, LJp )
    implicit none
    integer         , intent(in)    :: LIp, LJp, iX, jO1, jO2
    integer         , intent(inout) :: jX
    double precision, intent(in)    :: psi(LIp,LJp)
    integer                         :: j
    double precision                :: dpsi, dpsiMin

    jX      = ( jO1+jO2 )/2
    dpsiMin = +1.d8
    do j=4, LJp-3
       dpsi = abs( psi(iX,j+1) - psi(iX,j) )
       if ( dpsi.lt.dpsiMin ) then
          dpsiMin = dpsi
          jX      = j
       endif
    enddo
    return
  end subroutine solveXpt_j


  subroutine CheckXpt( psi, iX, jX, iO1, iO2, jO1, jO2, LIp, LJp, Flag__FixedXpt )
    implicit none
    integer         , intent(in)    :: LIp, LJp, jO1, jO2, iO1, iO2
    integer         , intent(inout) :: iX, jX
    double precision, intent(in)    :: psi(LIp,LJp)
    logical         , intent(out)   :: Flag__FixedXpt
    double precision                :: jrel, irel
    integer                         :: jXu, jXl, jXt

    ! -- i-Direction -- !
    irel   = dble( iX-iO1 ) / dble( iO2-iO1 )
    if ( ( irel.lt.0.1d0 ).or.( irel.gt.0.9d0 ) ) then
       iX = ( iO1+iO2 ) / 2
       call solveXpt_j( psi, iX, jX, jO1, jO2, LIp, LJp )
       Flag__FixedXpt = .true.
    endif

    ! -- j-Direction -- !
    jrel   = ( abs( dble(jX) - 0.5d0*( jO1+jO2 ) ) ) / dble( LJp )
    if ( jrel.gt.0.20d0 ) then
       jXl = max( ( jO1+jO2 )/2 - LJp/10,     3 )
       jXu = min( ( jO1+jO2 )/2 + LJp/10, LJp-1 )
       call solveXpt_j( psi(:,jXl:jXu), iX, jXt, jO1-jXl, jO2-jXl, LIp, jXu-jXl+1 )
       jX             = jXt + jXl - 1
       Flag__FixedXpt = .true.
    endif

    return
  end subroutine CheckXpt

  
end module getXOptMod
