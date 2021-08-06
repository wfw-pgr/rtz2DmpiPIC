module psifunc
contains
  
  subroutine pghat( psis, phat, ghat, dpdp, dgdp )

    use constants, only            : Nr, Nz
    use variables, only            : r, raxis
    use constants, only            : epsilon, alpha, Djh, beta_sep, Betap, lambda_
    use constants, only            : eta, c1, c2, psisw
    implicit none
    double precision, intent(in)  :: psis(Nz,Nr)
    double precision, intent(out) :: phat(Nz,Nr), ghat(Nz,Nr)
    double precision, intent(out) :: dpdp(Nz,Nr), dgdp(Nz,Nr)
    integer                       :: i, j
    double precision              :: etainv, K0, F0, beta, lamb
    
    if ( psisw.eq.1 ) then  ! -- Standard power law :: for ST or Spheromak -- !
       do j=1, Nr
          do i=1, Nz
             phat(i,j) =           psis(i,j)**( epsilon      )
             ghat(i,j) =           psis(i,j)**( alpha        )
             dpdp(i,j) = epsilon * psis(i,j)**( epsilon-1.d0 )
             dgdp(i,j) =   alpha * psis(i,j)**( alpha  -1.d0 )
          enddo
       enddo
    endif
    if ( psisw.eq.2 ) then  ! -- H-mode like profile :: for ST or tokamak -- !
       etainv = 1.d0 / eta
       do j=1, Nr
          do i=1, Nz
             phat(i,j) = ( c1 + c2*psis(i,j)**2 ) * ( 1.d0 - exp( - etainv * psis(i,j)**2 ) )
             dpdp(i,j) = ( 2.d0*c2*psis(i,j)    ) * ( 1.d0 - exp( - etainv * psis(i,j)**2 ) ) &
                  &      + 2.d0*etainv*( c1 + c2*psis(i,j)**2 ) * psis(i,j) * exp( - etainv * psis(i,j)**2 )
             ghat(i,j) =         psis(i,j)**( alpha        )
             dgdp(i,j) = alpha * psis(i,j)**( alpha  -1.d0 )
          enddo
       enddo
    endif
    if ( psisw.eq.3 ) then  ! -- FRC Hollowness distribution :: FRC -- !
       K0 = beta_sep * ( 1.d0 - 0.5d0 * Djh ) / ( 1.d0 - beta_sep )
       F0 = 1.d0 / K0
       do j=1, Nr
          do i=1, Nz
             if ( psis(i,j).ge.0.d0 ) then
                phat(i,j) =  K0 + psis(i,j) - 0.5d0*Djh*psis(i,j)**2
                dpdp(i,j) =       1.d0      -       Djh*psis(i,j)
             else
                phat(i,j) =  K0 * ( exp( F0*psis(i,j) ) )
                dpdp(i,j) =         exp( F0*psis(i,j) )
             endif
             ghat(i,j) =  0.d0
             dgdp(i,j) =  0.d0
          enddo
       enddo
    endif
    if ( psisw.eq.4 ) then  ! -- Spheromak -- !
       beta = Betap
       do j=1, Nr
          do i=1, Nz
             if ( psis(i,j).ge.0.d0 ) then
                phat(i,j) =  beta * psis(i,j)
                dpdp(i,j) =  beta
                ghat(i,j) =         lambda_ * psis(i,j)**2
                dgdp(i,j) =  2.d0 * lambda_ * psis(i,j)
             else
                phat(i,j) = 0.d0
                dpdp(i,j) = 0.d0
                ghat(i,j) = 0.d0
                dgdp(i,j) = 0.d0
             endif
          enddo
       enddo
    endif

    return
  end subroutine pghat

  
  subroutine pgWrite

    use constants, only : epsilon, alpha, jobdir
    use constants, only : eta, c1, c2, beta_sep, Djh, psisw, Betap, lambda_
    implicit none
    integer, parameter :: Npsi = 100
    integer            :: i
    double precision   :: dpsi, etainv, K0, F0, beta, lamb
    double precision   :: psi(Npsi), phat(Npsi), ghat(Npsi), dpdp(Npsi), dgdp(Npsi)
    
    !  --- [1] preparation  ---  !
    dpsi = 1.d0 / dble( Npsi-1 )
    do i=1, Npsi
       psi(i) = dble(i-1) * dpsi
    enddo

    !  --- [2] p(psi), g(psi) ---  !
    if ( psisw.eq.1 ) then
       do i=1, Npsi
          phat(i) =           psi(i)**( epsilon      )
          ghat(i) =           psi(i)**( alpha        )
          dpdp(i) = epsilon * psi(i)**( epsilon-1.d0 )
          dgdp(i) =   alpha * psi(i)**( alpha  -1.d0 )
       enddo
    endif
    if ( psisw.eq.2) then
       etainv = 1.d0 / eta
       do i=1, Npsi
          phat(i) = ( c1 + c2*psi(i)**2 ) * ( 1.d0 - exp( - etainv * psi(i)**2 ) )
          dpdp(i) = ( 2.d0*c2*psi(i)    ) * ( 1.d0 - exp( - etainv * psi(i)**2 ) ) &
               &    + 2.d0*etainv*( c1 + c2*psi(i)**2 ) * psi(i)*exp( - etainv * psi(i)**2 )
          ghat(i) =         psi(i)**( alpha        )
          dgdp(i) = alpha * psi(i)**( alpha  -1.d0 )
       enddo
    endif
    if ( psisw.eq.3 ) then
       K0 = beta_sep * ( 1.d0 - 0.5d0 * Djh ) / ( 1.d0 - beta_sep )
       F0 = 1.d0 / K0
       do i=1, Npsi
          phat(i) =  K0 + psi(i) - 0.5d0*Djh*psi(i)**2
          dpdp(i) =       1.d0    -       Djh*psi(i)
          ghat(i) =  0.d0
          dgdp(i) =  0.d0
       enddo
    endif
    if ( psisw.eq.4 ) then
       beta = Betap
       lamb = lambda_**2
       do i=1, Npsi
          phat(i) =  beta * psi(i)
          dpdp(i) =  beta
          ghat(i) =  lamb * psi(i)
          dgdp(i) =  lamb
       enddo
    endif
    
    !  --- [3] Write out p(psi), g(psi)  ---  !
    open(50,file=trim(jobdir)//'psifunc.dat',status='replace',form='formatted')
    do i=1, Npsi
       write(50,'(5(e14.7,1x))') psi(i), phat(i), ghat(i), dpdp(i), dgdp(i)
    enddo
    close(50)

    return
  end subroutine pgWrite
  
end module psifunc
