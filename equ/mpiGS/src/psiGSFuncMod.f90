module psiGSFuncMod
contains

  ! ====================================================== !
  ! === psi Function  P(psi), G(psi)                   === !
  ! ====================================================== !
  subroutine pghat( psih, phat, ghat, dpdp, dgdp )
    use variablesMod
    implicit none
    double precision, intent(in)  :: psih(LI,LJ)
    double precision, intent(out) :: phat(LI,LJ), ghat(LI,LJ)
    double precision, intent(out) :: dpdp(LI,LJ), dgdp(LI,LJ)
    integer                       :: i, j
    double precision              :: etainv, K0, F0, beta, lamb
    
    if ( psisw.eq.1 ) then  ! -- Standard power law :: for ST or Spheromak -- !
       do j=1, LJ
          do i=1, LI
             phat(i,j) =           psih(i,j)**( epsilon      )
             ghat(i,j) =           psih(i,j)**( alpha        )
             dpdp(i,j) = epsilon * psih(i,j)**( epsilon-1.d0 )
             dgdp(i,j) =   alpha * psih(i,j)**( alpha  -1.d0 )
          enddo
       enddo
    endif
    if ( psisw.eq.2 ) then  ! -- H-mode like profile :: for ST or tokamak -- !
       etainv = 1.d0 / eta
       do j=1, LJ
          do i=1, LI
             phat(i,j) = ( c1 + c2*psih(i,j)**2 ) * ( 1.d0 - exp( - etainv * psih(i,j)**2 ) )
             dpdp(i,j) = ( 2.d0*c2*psih(i,j)    ) * ( 1.d0 - exp( - etainv * psih(i,j)**2 ) ) &
                  &      + 2.d0*etainv*( c1 + c2*psih(i,j)**2 ) * psih(i,j) * exp( - etainv * psih(i,j)**2 )
             ghat(i,j) =         psih(i,j)**( alpha        )
             dgdp(i,j) = alpha * psih(i,j)**( alpha  -1.d0 )
          enddo
       enddo
    endif
    if ( psisw.eq.3 ) then  ! -- FRC Hollowness distribution :: FRC -- !
       K0 = beta_sep * ( 1.d0 - 0.5d0 * Djh ) / ( 1.d0 - beta_sep )
       F0 = 1.d0 / K0
       do j=1, LJ
          do i=1, LI
             if ( psih(i,j).ge.0.d0 ) then
                phat(i,j) =  K0 + psih(i,j) - 0.5d0*Djh*psih(i,j)**2
                dpdp(i,j) =       1.d0      -       Djh*psih(i,j)
             else
                phat(i,j) =  K0 * ( exp( F0*psih(i,j) ) )
                dpdp(i,j) =         exp( F0*psih(i,j) )
             endif
             ghat(i,j) =  0.d0
             dgdp(i,j) =  0.d0
          enddo
       enddo
    endif
    if ( psisw.eq.4 ) then  ! -- Spheromak -- !
       beta = beta_pol
       do j=1, LJ
          do i=1, LI
             if ( psih(i,j).ge.0.d0 ) then
                phat(i,j) =  beta * psih(i,j)
                dpdp(i,j) =  beta
                ghat(i,j) =         lambda * psih(i,j)**2
                dgdp(i,j) =  2.d0 * lambda * psih(i,j)
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
    use variablesMod
    implicit none
    integer, parameter :: Npsi = 100
    integer            :: i
    double precision   :: dpsi, etainv, K0, F0, beta, lamb
    double precision   :: psiF(Npsi), phat(Npsi), ghat(Npsi), dpdp(Npsi), dgdp(Npsi)
    
    !  --- [1] preparation  ---  !
    dpsi = 1.d0 / dble( Npsi-1 )
    do i=1, Npsi
       psiF(i) = dble(i-1) * dpsi
    enddo

    !  --- [2] p(psi), g(psi) ---  !
    if ( psisw.eq.1 ) then
       do i=1, Npsi
          phat(i) =           psiF(i)**( epsilon      )
          ghat(i) =           psiF(i)**( alpha        )
          dpdp(i) = epsilon * psiF(i)**( epsilon-1.d0 )
          dgdp(i) =   alpha * psiF(i)**( alpha  -1.d0 )
       enddo
    endif
    if ( psisw.eq.2) then
       etainv = 1.d0 / eta
       do i=1, Npsi
          phat(i) = ( c1 + c2*psiF(i)**2 ) * ( 1.d0 - exp( - etainv * psiF(i)**2 ) )
          dpdp(i) = ( 2.d0*c2*psiF(i)    ) * ( 1.d0 - exp( - etainv * psiF(i)**2 ) ) &
               &    + 2.d0*etainv*( c1 + c2*psiF(i)**2 ) * psiF(i)*exp( - etainv * psiF(i)**2 )
          ghat(i) =         psiF(i)**( alpha        )
          dgdp(i) = alpha * psiF(i)**( alpha  -1.d0 )
       enddo
    endif
    if ( psisw.eq.3 ) then
       K0 = beta_sep * ( 1.d0 - 0.5d0 * Djh ) / ( 1.d0 - beta_sep )
       F0 = 1.d0 / K0
       do i=1, Npsi
          phat(i) =  K0 + psiF(i) - 0.5d0*Djh*psiF(i)**2
          dpdp(i) =       1.d0    -       Djh*psiF(i)
          ghat(i) =  0.d0
          dgdp(i) =  0.d0
       enddo
    endif
    if ( psisw.eq.4 ) then
       beta = beta_pol
       lamb = lambda**2
       do i=1, Npsi
          phat(i) =  beta * psiF(i)
          dpdp(i) =  beta
          ghat(i) =  lamb * psiF(i)
          dgdp(i) =  lamb
       enddo
    endif
    
    !  --- [3] Write out p(psi), g(psi)  ---  !
    open(50,file=trim(jobdir)//'psifunc.dat',status='replace',form='formatted')
    do i=1, Npsi
       write(50,'(5(e14.7,1x))') psiF(i), phat(i), ghat(i), dpdp(i), dgdp(i)
    enddo
    close(50)

    return
  end subroutine pgWrite
  
end module psiGSFuncMod
