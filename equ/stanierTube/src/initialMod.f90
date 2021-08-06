module InitialMod
contains

  ! =================================================================== !
  ! ===  InitVariables  ::  Initialization of variables             === !
  ! =================================================================== !
  subroutine InitVariables
    use constants, only : N1, N2, normType, vthcv, wpewce, fixedRange
    use constants, only : dx1_Debye, dx2_Debye, x1Min, x1Max, x2Min, x2Max, myRank
    use variables, only : dx1, dx2, x1, x2
    use variables, only : Bfd, Jcr, rhs, Avp, prs, rho
    implicit none
    integer            :: i, j
    double precision   :: lDebye
    
    ! ------------------------------------- !
    ! --- [1]  Initialize x1, x2  Grid  --- !
    ! ------------------------------------- !
    !  -- [1-1] dx Determination        --  !
    if ( normType.eq.'PIC' ) then
       lDebye =     vthcv / wpewce
       if ( fixedRange ) then
          dx1 = ( x1Max - x1Min ) / dble( N1-1 )
          dx2 = ( x2Max - x2Min ) / dble( N2-1 )
       else
          dx1 = dx1_Debye * lDebye
          dx2 = dx2_Debye * lDebye
       endif
       x1Min  = - dx1*dble(N1-1) * 0.5d0
       x1Max  = + dx1*dble(N1-1) * 0.5d0
       x2Max  = + dx2*dble(N2-1) + x2Min
    endif
    if ( normType.eq.'MHD' ) then
       dx1    = ( x1Max - x1Min ) / dble( N1-1 )
       dx2    = ( x2Max - x2Min ) / dble( N2-1 )
    endif
    !  -- [1-2] x1,x2 Axis              --  !
    do i=1, N1
       x1(i) = dx1*dble(i-1) + x1Min
    enddo
    do j=1, N2
       x2(j) = dx2*dble(j-1) + x2Min
    enddo
    ! ------------------------------------- !
    ! --- [2]  Initialize  Variables    --- !
    ! ------------------------------------- !
    Bfd(:,:,:) = 0.d0
    Jcr(:,:,:) = 0.d0
    rhs(:,:)   = 0.d0
    Avp(:,:)   = 0.d0
    prs(:,:)   = 0.d0
    rho(:,:)   = 0.d0
    if ( myRank.eq.0 ) then
       write(6,'(2x,a)') '[ InitVariables  @initialMod ]'
       write(6,'(5x,a)') '* variables are Initialized.'
       write(6,*)
    endif
    return
  end subroutine InitVariables


  ! =================================================================== !
  ! ===  setJBRHS  ::  set Initial Profile ( J, B, rhs )            === !
  ! =================================================================== !
  subroutine setJBRHS
    use constants, only : N1, N2, coordinate, normType, myRank
    use variables, only : x1, x2, dx1, dx2
    use variables, only : Jcr, Bfd, rhs
    Use constants, only : x1cnt, x2cnt, r_omega, x2min
    use constants, only : Bt0, valfe, Bv0, jCenter1, jCenter2
    implicit none
    integer            :: i, j, hN1
    double precision   :: rw, radius, rhat, inroot, x1cntr, x2cntr, x1len, x2len
    double precision   :: jCt, jmB, jmB1, jmB2, coef, Itotal, Bv0Sample
    double precision   :: Bt_coef(7), x2factor(N2)    

    ! ------------------------------------- !
    ! --- [1] Preparation               --- !
    ! ------------------------------------- !
    hN1     = int( N1 / 2 )
    x1len   = dx1*dble( N1-1 )
    x2len   = dx2*dble( N2-1 )
    rw      = r_omega * min( x1len, x2len )
    x1cntr  = x1cnt * x1len * 0.5d0
    x2cntr  = x2cnt * x2len + x2min
    if ( normType.eq.'MHD' ) then
       jmB1 = jCenter1
       jmB2 = jCenter2
    endif
    if ( normType.eq.'PIC' ) then
       jmB1 = jCenter1 / valfe**2
       jmB2 = jCenter2 / valfe**2
    endif

    ! ------------------------------------- !
    ! --- [2] Bt (Taylor Solution)      --- !
    ! ------------------------------------- !
    !  -- ( see A.Stanier 2013 PoP )    --  !
    !  -- [2-1] Series expansion Coeff. --  !
    Bt_coef(1) = + ( Bt0 / ( jmB1*rw ) )**2 
    Bt_coef(2) = + 47.0d0 / 360.d0
    Bt_coef(3) = -  0.5d0
    Bt_coef(4) = +  3.0d0 /  4.0d0
    Bt_coef(5) = -  5.0d0 /  9.0d0
    Bt_coef(6) = +  5.0d0 / 24.0d0
    Bt_coef(7) = -  1.0d0 / 30.0d0
    !  -- [2-2] 1/r Correction ( RZ )   --  !
    if ( coordinate.eq.'RZ' ) then
       if ( x2min.gt.0.d0 ) then   ! - r=0 exc. - !
          do j=1, N2
             x2factor(j) = ( x2cntr / x2(j) )**2
          enddo
       else                        ! - r=0 inc. - !
          x2factor(1) = 0.d0
          do j=2, N2
             x2factor(j) = ( x2cntr / x2(j) )**2
          enddo
       endif
    endif
    if ( coordinate.eq.'XZ' ) then
       do j=1, N2
          x2factor(j) = 1.d0
       enddo
    endif
    !  -- [2-3] J & B substitution      --  !
    do j=1, N2
       do i=1, N1

          radius = sqrt( ( abs(x1(i))-x1cntr )**2 + ( abs(x2(j))-x2cntr )**2 )
          if ( radius.lt.rw ) then  ! -( Inside FluxTube )- !
             if ( i.lt.hN1 ) then
                jCt        = jCenter1
                jmB        = jmB1
                Bt_coef(1) = + ( Bt0 / ( jmB1*rw ) )**2
             endif
             if ( i.gt.hN1 ) then
                jCt        = jCenter2
                jmB        = jmB2
                Bt_coef(1) = + ( Bt0 / ( jmB2*rw ) )**2
             endif
             rhat          = radius / rw
             Jcr(3,i,j)    = +jCt * ( 1.d0 - rhat**2 )**2
             inroot        =   Bt_coef(1)*x2factor(j) + Bt_coef(2)         &
                  &          + Bt_coef(3)*rhat**2     + Bt_coef(4)*rhat**4 &
                  &          + Bt_coef(5)*rhat**6     + Bt_coef(6)*rhat**8 &
                  &          + Bt_coef(7)*rhat**10
             if ( inroot.lt.0.d0 ) inroot = 0.d0  ! --: ( garbage Error prevention  ) !
             Bfd(3,i,j)    = -jmB * rw * sqrt( inroot )
             
          else                      ! -( Outside FluxTube )- !
             Jcr(3,i,j)    =  0.d0
             Bfd(3,i,j)    = - Bt0 * sqrt( x2factor(j) )
          endif
          
       enddo
    enddo
    
    ! ------------------------------------- !
    ! --- [3] Set R.H.S. of Poisson Eq. --- !
    ! ------------------------------------- !
    Itotal = 0.d0
    if ( normType.eq.'MHD' ) coef = - 1.d0
    if ( normType.eq.'PIC' ) coef = - 1.d0 / valfe**2
    do j=2, N2-1
       do i=2,N1-1
          rhs(i,j) =   coef * Jcr(3,i,j)
          Itotal   = Itotal + Jcr(3,i,j)*dx1*dx2
       enddo
    enddo
    
    ! ------------------------------------- !
    ! --- [4] Boundary Condition        --- !
    ! ------------------------------------- !
    if ( coordinate.eq.'XZ' ) then
       !  - @ Cartesian Geometry                   - !
       !  - No vertical Magnetic Field on the wall - !
       rhs( :, 1) = 0.d0
       rhs( :,N2) = 0.d0
       rhs( 1, :) = 0.d0
       rhs(N1, :) = 0.d0
    endif
    if ( coordinate.eq.'RZ' ) then
       !  - @ Cylindrical Geometry                 - !
       !  - EF :: Bv0 is imposed on the wall       - !
       do i=1, N1
          rhs(i,N2) = 0.5d0 * Bv0 * x2(N2)
       enddo
       do j=1, N2
          rhs( 1,j) = 0.5d0 * Bv0 * x2( j)
          rhs(N1,j) = 0.5d0 * Bv0 * x2( j)
       enddo
       if ( x2min.gt.0.d0 ) then ! -- Doubly Connected Case -- !
          do i=1, N1
             rhs(i, 1) = 0.5d0 * Bv0 * x2(1)
          enddo
       endif
       if ( myRank.eq.0 ) then
          Bv0Sample = - coef*Itotal / ( 16.d0*atan(1.d0)*x2cntr ) * 0.5d0
          write(6,'(2x,a)'         ) '[ setJBRHS       @initialMod ]'
          write(6,'(5x,a)'         ) '* Bv0 is set on Boundary'
          write(6,'(5x,a)'         ) '* Suggested value of Bv0 :: Bv0 = mu I / 4pi R0 * 0.5d0'
          write(6,'(5x,a)'         ) '* Suggested value of Bv0 :: Bv0 = mu I / 4pi R0 * 0.5d0'
          write(6,'(5x,a,1x,e12.5)') '                                = ', Bv0Sample
          write(6,*)
          write(6,'(2x,a)'         ) "----------------------------------------------------------------------"
          write(6,*)
          write(6,*)
          write(6,*)
       endif
    endif
    
    return
  end subroutine setJBRHS

  
end module InitialMod




    ! ! Bt ( Taylor solution ) :: see A.Stanier 2013 POP
    ! Bt_coef(1) = + ( Bt0 / jmB1 )**2 
    ! Bt_coef(2) = + 47.d0 * rw**2 / 3.6d2
    ! Bt_coef(3) = - 0.5d0
    ! Bt_coef(4) = + 3.d0 / ( 4.0d0 * rw**2 )
    ! Bt_coef(5) = - 5.d0 / ( 9.0d0 * rw**4 )
    ! Bt_coef(6) = + 5.d0 / ( 2.4d1 * rw**6 )
    ! Bt_coef(7) = - 1.d0 / ( 3.0d1 * rw**8 )

    ! if ( coordinate.eq.'RZ' ) then
    !    if ( x2min.gt.0.d0 ) then   ! - r=0 exc. - !
    !       do j=1, N2
    !          x2factor(j) = ( x2cntr / x2(j) )**2
    !       enddo
    !    else                        ! - r=0 inc. - !
    !       x2factor(1) = 0.d0
    !       do j=2, N2
    !          x2factor(j) = ( x2cntr / x2(j) )**2
    !       enddo
    !    endif
    ! endif
    ! if ( coordinate.eq.'XZ' ) then
    !    do j=1, N2
    !       x2factor(j) = 1.d0
    !    enddo
    ! endif
    
    ! ! --- [2] Define J & B on each grid --- !
    ! do j=1, N2
    !    do i=1, N1

    !       radius = sqrt( ( abs(x1(i))-x1cntr )**2 + ( abs(x2(j))-x2cntr )**2 )
    !       if ( radius.lt.rw ) then
    !          if ( i.lt.hN1 ) then
    !             jCt = jCenter1
    !             jmB = jmB1
    !             Bt_coef(1) = + ( Bt0 / jmB1 )**2 
    !          endif
    !          if ( i.gt.hN1 ) then
    !             jCt = jCenter2
    !             jmB = jmB2
    !             Bt_coef(1) = + ( Bt0 / jmB2 )**2 
    !          endif
    !          Jcr(3,i,j) = +jCt * ( 1.d0 - ( radius / rw )**2 )**2
    !          Bfd(3,i,j) = -jmB     * sqrt( &
    !               & + Bt_coef(1)*x2factor(j) + Bt_coef(2)           &
    !               & + Bt_coef(3)*radius**2   + Bt_coef(4)*radius**4 &
    !               & + Bt_coef(5)*radius**6   + Bt_coef(6)*radius**8 &
    !               & + Bt_coef(7)*radius**10  )
    !       else
    !          Jcr(3,i,j) =  0.d0
    !          Bfd(3,i,j) = - Bt0 * sqrt( x2factor(j) )
    !       endif
          
    !    enddo
    ! enddo
