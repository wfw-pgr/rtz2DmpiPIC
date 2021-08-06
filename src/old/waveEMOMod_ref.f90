module waveEMOMod
contains
  
  subroutine setEMO
    use constants, only : LI, LJ, ns, npt, mode_num, wave_amplitude
    use variables, only : Ex, Ey, Ez, Bx, By, Bz
    use variables, only : x, y, vx, vy, vz, xpast, ypast
    use constants, only : dx, dt, wp, np, qm, cv, pch, vth_perp, vth_para, angle, vthe, vthi
    use variables, only : xleng, yleng, Bx0, By0, cvinv
    use variables, only : Exold, Eyold, Ezold
    use randomgen, only : unrndm, strndm, gaussdev
    use maxwell  , only : periodic_xy
    implicit none
    integer            :: i, j, k, n1, n2, seed
    double precision   :: gamma, sum
    double precision   :: thisBz, thisEy, B2inv, vdx, vdy, vdz
    double precision   :: thiskx, thisw, B1, E1, theta, rand(2)
    double precision, parameter :: pi = 4.d0 * atan( 1.d0 )
    
    ! (0) preparation : calc disp.rel.
    thiskx = 2.d0 * pi / ( xleng / mode_num )
    thisw  = sqrt( wp(1)**2 + wp(2)**2 + cv**2 * thiskx**2 )
    write(6,*) 'wp == ', wp(1)**2, wp(2)**2
    write(6,*) 'ck == ', cv * thiskx
    
    B1 = By0 * wave_amplitude
    E1 = thisw / ( cv * thiskx ) * B1
    write(6,*) 'EMW(O-mode) :: k == ', thiskx, ' w== ', thisw
    write(6,*) 'amp :: B0 = ', By0, 'B1 = ', B1, 'E1 = ', E1
    if ( angle.ne.90.d0 ) write(6,*) ' [ALERT]  angle is not 90 degree in EM-O mode wave [ALERT] '

    ! (1) set E
    do j=2, LJ-1
       do i=2, LI-1
          Ex(i,j) = 0.d0
          Ey(i,j) = E1 * cos( thiskx * dx * dble(i-2) )
          Ez(i,j) = 0.d0
       enddo
    enddo
    ! b.c.
    call periodic_xy(Ex)
    call periodic_xy(Ey)
    call periodic_xy(Ez)

    do j=1, LJ
       do i=1, LI
          Exold(i,j) = Ex(i,j)
          Eyold(i,j) = Ey(i,j)
          Ezold(i,j) = Ez(i,j)
       enddo
    enddo

    ! (2) set B
    do j=2, LJ-1
       do i=2, LI-1
          Bx(i,j) = Bx0
          By(i,j) = By0
          Bz(i,j) = B1 * cos( thiskx * dx * ( dble(i-2) + 0.5d0 ) )
       enddo
    enddo
    ! b.c.
    call periodic_xy(Bx)
    call periodic_xy(By)
    call periodic_xy(Bz)

    ! (3) set x  :: quiet start
    seed = 0
    do i=1, np(1)
       x(i) = xleng * unrndm( seed )
       x(i+np(1)) = x(i)
       y(i) = yleng * unrndm( seed )
       y(i+np(1)) = y(i)
    enddo

    do i=1, npt
       xpast(i) = x(i)
       ypast(i) = y(i)
    enddo

    ! (4) set v ( Lorentz eq. )    @ t' = 1.0
    n2= 0
    do k=1, ns
       n1 = n2
       n2 = n1 + np(k)
       do i=n1+1, n2

          vdx = 0.d0
          vdy = - qm(k) * E1 / thisw * sin( thiskx*x(i) - thisw*0.5d0*dt )
          vdz = 0.d0
          theta = 2.d0 * pi * unrndm( seed )
          ! vx(i) = vdx + vth_perp(k) * strndm(seed) * cos(theta)
          ! vy(i) = vdy + vth_para(k) * strndm(seed)
          vx(i) = vdx + vth_perp(k) * gaussdev( seed ) * cos(theta)
          vy(i) = vdy + vth_para(k) * gaussdev( seed )
          ! vz(i) = vdz + vth_perp(k) * strndm(seed) * sin(theta)
          vz(i) = vdz + vth_perp(k) * gaussdev( seed ) * sin(theta)

          gamma = 1.d0 / sqrt( 1.d0 - ( vx(i)**2 + vy(i)**2 + vz(i)**2 ) )
          vx(i) = gamma * vx(i)
          vy(i) = gamma * vy(i)
          vz(i) = gamma * vz(i)
          
       enddo
    enddo

    ! particle back step :: xypast calculation
    do i=1, npt
       gamma = dt / sqrt( 1.d0 + ( vx(i)**2 + vy(i)**2 + vz(i)**2 ) )
       xpast(i) = x(i) - gamma * vx(i)
       ypast(i) = y(i) - gamma * vy(i)
    enddo

    sum = 0.d0
    do i=1, np(1)
       sum = sum + vx(i)**2 + vy(i)**2 + vz(i)**2
    enddo
    sum = sqrt( sum / dble( np(1) ) )

    write(6,*) 'vthe check ', vthe, sum

    sum = 0.d0
    do i=np(1) + 1, np(1) + np(2)
       sum = sum + vx(i)**2 + vy(i)**2 + vz(i)**2
    enddo
    sum = sqrt( sum / dble( np(2) ) )
    
    write(6,*) 'vthi check ', vthi, sum

    return
  end subroutine setEMO

end module waveEMOMod
