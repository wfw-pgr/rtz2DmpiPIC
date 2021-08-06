module debugerMod
contains

  subroutine ampereconf
    use constants, only : LI, LJ, dt, dtdz, dtdr, valfe, myRank
    use constants, only : br_, bt_, bz_, er_, et_, ez_
    use variables, only : EBf, Jcr, rf, rhinv
    implicit none
    integer            :: i, j
    double precision   :: factor1, factor2, factor3, sum
    double precision   :: rotB(3,LI,LJ)

    ! --- [1] Initialization --- !
    factor1     = dtdr
    factor2     = dtdz
    factor3     = dt / vAlfe**2
    rotB(:,:,:) = 0.d0
    ! --- [2] rot(B)         --- !
    do j=2, LJ-1
       do i=2, LI-1
          rotB(1,i,j) = - factor2*( EBf(bt_,i+1,j) - EBf(bt_,i,j) ) &
               &                  - factor3*Jcr(1,i,j,0)
          rotB(2,i,j) = + factor2*( EBf(br_,i+1,j) - EBf(br_,i,j) ) &
               &        - factor1*( EBf(bz_,i,j+1) - EBf(bz_,i,j) ) &
               &                  - factor3*Jcr(2,i,j,0)
          rotB(3,i,j) = + factor1*( rf(j+1)*EBf(bt_,i,j+1) - rf(j)*EBf(bt_,i,j) )*rhinv(j) &
               &                  - factor3*Jcr(3,i,j,0)
       enddo
    enddo
    ! --- [3] Sum of Amplitude --- !
    sum = 0.d0
    do j=2, LJ-1
       do i=2, LI-1
          sum = sum + sqrt( rotB(1,i,j)**2 + rotB(2,i,j)**2 + rotB(3,i,j)**2 )
       enddo
    enddo

    if ( myRank.eq.0 ) then
       write(6,*)
       write(6,*) ' ampereconf :: sum =  ', sum
       write(6,*)
    endif

    return
  end subroutine ampereconf

  
  subroutine writeSingleEBMap( FileName )
    use variables, only : EBf
    implicit none
    character(*)       :: FileName

    open( 50, file=trim(FileName), form='unformatted' )
    write(50) EBf(1,:,:)
    close(50)
    write(6,'(2x,a)') '----------------------------------------------------------------------'
    write(6,'(2x,a10,a45,2x,a6)') "* SAVE :: ", trim(FileName), '[ OK ]'


    return
  end subroutine writeSingleEBMap
    

  
end module debugerMod
