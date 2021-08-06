module diagnosis
contains

  subroutine EquDiagnosis
    implicit none
    write(6,*)
    write(6,*)
    write(6,*) '    -----------------------------------------------------------  '
    write(6,*) '    -----------        Equilibrium Diagnosis        -----------  '
    write(6,*) '    -----------------------------------------------------------  '
    write(6,*)
    write(6,*)
    call ZerothIntegral
    call CrossSectionAnalysis
    call q_profile
    write(6,*)
    write(6,*) '    ------------------     END OF ANALYSIS   - ----------------  '
    write(6,*)
    return
  end subroutine EquDiagnosis


  subroutine ZerothIntegral
    
    use constants    , only : Nr, Nz, limposr, limposz, limidx, valfe, jobdir
    use constants    , only : alpha, epsilon, Beta_sep, Djh
    use variables    , only : r, z, dr, dz
    use variables    , only : prs, psi, psilim, psimin, Itotal
    use Field        , only : Br, Bt, Bz, Jt
    use NewtonRaphson, only : cSpline2D
    implicit none
    integer                :: i,j
    double precision       :: fx(Nr,Nz) , psis(Nz,Nr)
    double precision       :: Bpol1, Bpol2, Ppol, area, Betap1, Beta_Vol, pres, Itot, Phit, Phid, volm, absB
    double precision       :: separatrix, Bra, Bza, dBdx, dBdy, dS, dV, dVor, dpsiinv, pLength, Betap2, CrossCheck

    Bpol1    = 0.d0
    Bpol2    = 0.d0
    Ppol     = 0.d0
    area     = 0.d0    
    !  --- [1] calculate Area and Pressure  --- !
    fx(:,:)  = 1.d0
    area     = pol_integral( z, r, fx , Nz, Nr, psi, psilim )
    Ppol     = pol_integral( z, r, prs, Nz, Nr, psi, psilim )
    Ppol     = Ppol / area

    !  --- [2] calculate Poloidal Magnetic Pressure --- !
    !   -- [2-1] Define Bp as |B| @limiter  -- !
    call cSpline2D( z, r, Br, limposz(limidx), limposr(limidx), Bra, dBdx, dBdy )
    call cSpline2D( z, r, Bz, limposz(limidx), limposr(limidx), Bza, dBdx, dBdy )
    Bpol1    = 0.5d0 * ( Bra**2 + Bza**2 )
    Betap1   = valfe**2 * Ppol / Bpol1
    !  --- [2-2] Define Bp as Ip / l        -- !
    dpsiinv  = 1.d0 / ( psilim - psimin )
    do j=1, Nr
       do i=1, Nz
          psis(i,j) = ( psilim - psi(i,j) ) * dpsiinv
       enddo
    enddo
    pLength  = Perimeter( psis, Nz, Nr, dz, dr )
    Bpol2    = ( 0.5d0 * ( Itotal / pLength )**2 )
    Betap2   = valfe**2 * Ppol / Bpol2
    
    !  --- [3] Phi total ( Toroidal Flux ) --- !
    Phit     = pol_integral( z, r, Bt , Nz, Nr, psi, psilim )
    Itot     = pol_integral( z, r, Jt , Nz, Nr, psi, psilim )

    !  --- [4] Volume Integral  ---  !
    volm     = 0.d0
    pres     = 0.d0
    absB     = 0.d0
    dVor     = 2.d0 * ( 4.d0*atan(1.d0) ) * dr*dz
    do j=1, Nr
       do i=1, Nz
          if ( psi(i,j).le.psilim ) then
             dV   = r(j) * dVor
             volm = volm + dV
             pres = pres + prs(i,j)*dV
             absB = absB + 0.5d0 * ( Br(i,j)**2 + Bt(i,j)**2 + Bz(i,j)**2 ) * dV
          endif
       enddo
    enddo
    pres     = pres / volm * valfe**2
    absB     = absB / volm
    Beta_Vol = pres / absB
    
    write(6,*)
    write(6,*) '    -----------------------------------------------    '
    write(6,*) '    ----  In Separatrix Equilibrium Diagnosis  ----    '
    write(6,*) '    -----------------------------------------------    '
    write(6,*)
    write(6,'(6x,a32,2x,e13.6)') '                      Area :: ', area
    write(6,'(6x,a32,2x,e13.6)') '  1/2 Bp^2 (@Limiter     ) :: ', Bpol1
    write(6,'(6x,a32,2x,e13.6)') '  1/2 Bp^2 (@LineIntegral) :: ', Bpol2
    write(6,'(6x,a32,2x,e13.6)') '    Area Averaged Pressure :: ', Ppol
    write(6,'(6x,a32,2x,e13.6)') '                   Volmume :: ', Volm
    write(6,'(6x,a32,2x,e13.6)') '  Volume Averaged Pressure :: ', pres
    write(6,'(6x,a32,2x,e13.6)') '              1/2 B(vol)^2 :: ', absB
    write(6,*)
    write(6,'(6x,a32,2x,e13.6)') '    BetaP1 (@Limiter     ) :: ', Betap1
    write(6,'(6x,a32,2x,e13.6)') '    BetaP2 (@LineIntegral) :: ', Betap2
    write(6,'(6x,a32,2x,e13.6)') '    Beta_Volume            :: ', Beta_Vol
    write(6,*)
    write(6,'(6x,a32,2x,e13.6)') '             Toroidal Flux :: ', Phit
    write(6,'(6x,a32,2x,e13.6)') '                   I_total :: ', Itot
    write(6,*)

    open(20,file=trim(jobdir)//'equilibrium.dat',form='formatted',position='append')
    write(20,'(a20,1x,e13.6)') 'area'    , area
    write(20,'(a20,1x,e13.6)') 'Bpol'    , Bpol1
    write(20,'(a20,1x,e13.6)') 'Ppol'    , Ppol
    write(20,'(a20,1x,e13.6)') 'Betap1'  , Betap1
    write(20,'(a20,1x,e13.6)') 'Phit'    , Phit
    write(20,'(a20,1x,e13.6)') 'Itot'    , Itot
    write(20,'(a20,1x,e13.6)') 'alpha'   , alpha
    write(20,'(a20,1x,e13.6)') 'epsilon' , epsilon
    write(20,'(a20,1x,e13.6)') 'Beta_sep', Beta_sep
    write(20,'(a20,1x,e13.6)') 'Djh'     , Djh
    close(20)

    !  --- [4] Direct Summation ( cross-check ) --- !
    CrossCheck = 0
    if ( crossCheck.eq.1 ) then
       area = 0.d0
       pres = 0.d0
       Itot = 0.d0
       Phid = 0.d0
       dS   = dr*dz
       do j=1, Nr
          do i=1, Nz
             if ( psi(i,j).le.psilim ) then
                area = area + dS
                pres = pres + prs(i,j)*dS
                Itot = Itot +  Jt(i,j)*dS
                Phid = Phid +  Bt(i,j)*dS
             endif
          enddo
       enddo
       Betap1 = valfe**2 * pres / ( area * Bpol1 )

       write(6,*)
       write(6,*) '    -- Direct Summation ( w/o spline :: validation ) -- '
       write(6,'(6x,a32,2x,e13.6)') '             Poloidal Beta :: ', Betap1
       write(6,'(6x,a32,2x,e13.6)') '             Toroidal Flux :: ', Phid
       write(6,'(6x,a32,2x,e13.6)') '                   I_total :: ', Itot
       write(6,*)
    endif
    
    return
  end subroutine ZerothIntegral
        

  function gauss_quad_spline( x, x1, x2, fx )

    use NewtonRaphson, only       : cSpline1D
    implicit none
    double precision, intent(in) :: x(:), fx(:)
    double precision, intent(in) :: x1, x2
    double precision             :: gauss_quad_spline
    integer                      :: j, N
    double precision             :: xp(5), wp(5)
    double precision             :: dx, xm, xr, yplus, yminus, yderiv, sum
    ! abscissas and weights
    data xp/.1488743389,.4333953941,.6794095682,.8650633666,.9739065285/
    data wp/.2955242247,.2692667193,.2190863625,.1494513491,.0666713443/

    !  (0) preparation
    xm  = 0.5d0 * ( x2 + x1 )
    xr  = 0.5d0 * ( x2 - x1 )
    N   = size( x )
    sum = 0.d0
    
    !  (1) integration
    do j=1, 5
       dx = xr * xp(j)
       call cSpline1D( x, fx, N, xm+dx, yplus , yderiv )
       call cSpline1D( x, fx, N, xm-dx, yminus, yderiv )
       sum = sum + wp(j) * ( yplus + yminus )
    enddo
    gauss_quad_spline = xr * sum

    return
  end function gauss_quad_spline


  function pol_integral( z, r, var, Nz, Nr, psi, psilim )

    use NewtonRaphson, only       : cSplineNewtonRaphson1D
    implicit none
    integer         , intent(in) :: Nr, Nz
    double precision, intent(in) :: r(Nr), z(Nz), var(Nz,Nr)
    double precision, intent(in) :: psi(Nz,Nr), psilim
    integer                      :: i, j, inflag, flag, edger, edgel
    double precision             :: start, end, dfx
    double precision             :: rintegrated(Nz-1), fx(Nr)
    double precision             :: pol_integral

    inflag         = 0
    pol_integral   = 0.d0
    rintegrated(:) = 0.d0
    start          = 0.d0
    end            = 0.d0
    do i=1, Nz-1
       inflag = 0
       do j=1, Nr-1
          !if ( ( psi(i,j).le.psilim ).and.( inflag.eq.0 ) ) then
          !if ( ( ( psi(i,j)*psi(i+1,j) ).le.0.d0 ).and.( ( psi(i+1,j) - psi(i,j) ).le.0.d0 ) ) then
          if ( ( ( ( psi(i,j) - psilim ) * ( psi(i,j+1) - psilim ) ).le.0.d0 ) &
               .and.( ( psi(i,j+1) - psi(i,j) ).le.0.d0 ) ) then
             
             inflag  = 1
             start   = 0.5d0 * ( r(j) + r(j+1) )
             fx(:)   = psi(i,:) - psilim
             call cSplineNewtonRaphson1D( r, fx, start, dfx, flag )
             
          endif
          !if ( ( psi(i,j).gt.psilim ).and.( inflag.eq.1 ) ) then
          !if ( ( ( psi(i,j)*psi(i+1,j) ).le.0.d0 ).and.( ( psi(i+1,j) - psi(i,j) ).ge.0.d0 ) ) then

          if ( ( ( (psi(i,j)-psilim)*(psi(i,j+1)-psilim) ).le.0.d0 ) &
               .and.( ( psi(i,j+1) - psi(i,j) ).ge.0.d0 ) ) then
             inflag  = 0
             end     = 0.5d0 * ( r(j) + r(j+1) )
             fx(:)   = psi(i,:) - psilim
             call cSplineNewtonRaphson1D( r, fx, end, dfx, flag )
             
             rintegrated(i) = rintegrated(i) + gauss_quad_spline( r, start, end, var(i,:) )
             
          endif
       enddo
    enddo
    
    ! search edge
    inflag  = 0
    edgel   = 0
    do i=1, Nz-1
       do j=1, Nr-1
          if ( ( psi(i,j).lt.psilim ).and.( inflag.eq.0 ) ) then
             edgel  = i
             inflag = 1
             exit
          endif
       enddo
       if ( inflag .eq. 1 ) exit
    enddo
    
    inflag = 0
    edger  = 0
    do i=Nz-1, 1, -1
       do j=1, Nr-1
          if ( (psi(i,j).lt.psilim).and.(inflag.eq.0) ) then
             edger  = i+1
             inflag = 1
             exit
          endif
       enddo
       if ( inflag .eq. 1 ) exit
    enddo
    pol_integral = gauss_quad_spline( z(edgel:edger), z(edgel), z(edger), rintegrated(edgel:edger) )
 
    return
  end function pol_integral

  
  subroutine CrossSectionAnalysis
    !  (*) calculate ellipticity by kappa = ( Zmax - Zmin ) / 2a
    use constants    , only : Nr, Nz, mr, vthcv, TiTe, jobdir
    use variables    , only : r, z, psi
    use variables    , only : psimin, psilim, raxis
    use Field        , only : Br, Bt, Bz
    use NewtonRaphson, only : cSplineNewtonRaphson1D
    implicit none
    integer                :: i, j, j1, j2
    integer                :: flag, count
    double precision       :: rk, zk, psival, det, psirr, dpsidr, threshold, dl, s_
    double precision       :: LCFSrmin, LCFSrmax, LCFSzmin, LCFSzmax, LCFSrupr, LCFSrlwr
    double precision       :: Rgeo, triplus, triminus, aMinor, Triangularity, Ellipticity, AspectRatio
    double precision       :: psis(Nr), lambda_i(Nr)

    !  --- [0]  preparation  --- !
    threshold = 0.05 * ( psimin - psilim ) + psilim
    psis      = psi(Nz/2,:) - threshold
    
    !  --- [1]  LCFS Rmin, Rmax, aMinor --- !
    !    -- [1-1] Find Last Closed Flux Serface : Rmin
    LCFSrmin  = 0.3d0 * ( r(Nr) - r(1) )                     ! -- first estimation
    call cSplineNewtonRaphson1D( r, psis, LCFSrmin, dpsidr, flag )

    !    -- [1-2] Find Last Closed Flux Serface : RMax
    LCFSrmax  = 0.7d0 * ( r(Nr) - r(1) )                     ! -- first estimation
    call cSplineNewtonRaphson1D( r, psis, LCFSrmax, dpsidr, flag )
    
    aMinor    = 0.5d0 * ( LCFSrmax - LCFSrmin )

    !  --- [2] LCFS :: zmax, zmin, r_upper, r_lower
    do j=1, Nr
       do i=1, Nz
          if ( ( psi(i,j).le.threshold ).and.( z(i).le.LCFSzmin ) ) then
             LCFSzmin = z(i)
             LCFSrlwr = r(j)
          endif
          if ( ( psi(i,j).le.threshold ).and.( z(i).ge.LCFSzmax ) ) then
             LCFSzmax = z(i)
             LCFSrupr = r(j)
          endif
       enddo
    enddo

    !  --- [3] Ellipticity, AspectRatio, Triangularity --- !
    !   -- [3-1] Geometric Factor -- !
    Ellipticity   = ( LCFSzmax - LCFSzmin ) / ( 2.d0 * aMinor )
    AspectRatio   =   raxis / aMinor
    Rgeo          = ( LCFSrmax + LCFSrmin ) * 0.5d0
    triplus       = ( Rgeo - LCFSrupr ) / aMinor
    triminus      = ( Rgeo - LCFSrlwr ) / aMinor
    Triangularity = 0.5d0 * ( triplus + triminus )
    !   -- [3-2] s_ :: kinetic parameter -- !
    j1 = 1
    j2 = 1
    do j=1, Nr
       if ( abs( r(j)-raxis    ).lt.abs( r(j1)-raxis    ) ) j1 = j
       if ( abs( r(j)-LCFSrmax ).lt.abs( r(j2)-LCFSrmax ) ) j2 = j
    enddo ! -- Search rAxis & Separatix Grid :: (j1, j2) -- !
    do j=j1, j2
       lambda_i(j) = + sqrt( mr*TiTe ) * vthcv  &
            &        / sqrt( Br(Nz/2,j)**2 + Bt(Nz/2,j)**2 + Bz(Nz/2,j)**2 )
    enddo ! -- Larmor Radius    -- !
    dl = ( LCFSrmax - raxis ) / dble( j2-j1 )
    s_ = 0.d0
    do j=j1, j2
       s_ = s_ + r(j) * dl / ( LCFSrmax * lambda_i(j) )
    enddo ! -- Definition of s_ -- !

    !  --- [4] write out the results  --- !
    !    -- [4-1] to Standard Output  --- !
    write(6,*)
    write(6,*)
    write(6,*) '    -----------------------------------------------    '
    write(6,*) '    ----        Cross Section Analysis         ----    '
    write(6,*) '    -----------------------------------------------    '
    write(6,*)
    write(6,'(6x,a32,2x,e13.6)') ' LCFS Rmin      :: ', LCFSrmin
    write(6,'(6x,a32,2x,e13.6)') ' LCFS Rmax      :: ', LCFSrmax
    write(6,'(6x,a32,2x,e13.6)') ' LCFS Zmin      :: ', LCFSzmin
    write(6,'(6x,a32,2x,e13.6)') ' LCFS Zmax      :: ', LCFSzmax
    write(6,'(6x,a32,2x,e13.6)') ' LCFS Rupper    :: ', LCFSrupr
    write(6,'(6x,a32,2x,e13.6)') ' LCFS Rlower    :: ', LCFSrlwr
    write(6,*)
    write(6,'(6x,a32,2x,e13.6)') ' Rgeometry      :: ', Rgeo
    write(6,'(6x,a32,2x,e13.6)') ' triangular+    :: ', triplus
    write(6,'(6x,a32,2x,e13.6)') ' triangular-    :: ', triminus
    write(6,*)
    write(6,'(6x,a32,2x,e13.6)') ' aMinor         :: ', aMinor
    write(6,'(6x,a32,2x,e13.6)') ' Raxis          :: ', Raxis
    write(6,'(6x,a32,2x,e13.6)') ' Aspect Ratio   :: ', AspectRatio
    write(6,'(6x,a32,2x,e13.6)') ' Ellipticity    :: ', Ellipticity
    write(6,'(6x,a32,2x,e13.6)') ' Triangularity  :: ', Triangularity
    write(6,'(6x,a32,2x,e13.6)') ' s_ (kinetic)   :: ', s_
    write(6,*)
    write(6,*)
    !    -- [4-2] to File  --- !
    open(20,file=trim(jobdir)//'equilibrium.dat',form='formatted',position='append')
    write(20,'(a20,1x,e13.6)') 'aMinor'          , aMinor
    write(20,'(a20,1x,e13.6)') 'Raxis'           , Raxis
    write(20,'(a20,1x,e13.6)') 'AspectRatio'     , AspectRatio
    write(20,'(a20,1x,e13.6)') 'Ellipticity'     , Ellipticity
    write(20,'(a20,1x,e13.6)') 'Triangularity'   , Triangularity
    write(20,'(a20,1x,e13.6)') 'LCFS_Rmin'       , LCFSrmin
    write(20,'(a20,1x,e13.6)') 'LCFS_Rmax'       , LCFSrmax
    write(20,'(a20,1x,e13.6)') 'LCFS_Zmin'       , LCFSzmin
    write(20,'(a20,1x,e13.6)') 'LCFS_Zmax'       , LCFSzmax
    write(20,'(a20,1x,e13.6)') 'LCFS_Rupper'     , LCFSrupr
    write(20,'(a20,1x,e13.6)') 'LCFS_Rlower'     , LCFSrlwr
    write(20,'(a20,1x,e13.6)') 'psilim'          , psilim
    write(20,'(a20,1x,e13.6)') 'Raxis'           , Raxis
    write(20,'(a20,1x,e13.6)') 'psimin'          , psimin
    write(20,'(a20,1x,e13.6)') 'Rgeometry'       , Rgeo
    write(20,'(a20,1x,e13.6)') 'triangular+'     , triplus
    write(20,'(a20,1x,e13.6)') 'triangular-'     , triminus
    close(20)

    return
  end subroutine CrossSectionAnalysis

  
  subroutine q_profile

    use constants,          only : Nr, Nz, jobdir
    use variables,          only : r, z, dr, dz, psi, psilim, psimin
    use Field    ,          only : Bt
    implicit none
    integer         , parameter :: Nq = 30
    double precision, parameter :: pi = 4.d0 * atan(1.d0)
    integer                     :: i, j, k
    double precision            :: phi(Nq), psis(Nq), psil(Nq), dphi(Nq), psiaxis(Nq), qfactor(Nq)
    double precision            :: dpsi, q95, qaxis, dS, coef

    !  --- [1]  Preparation ---  !
    phi(:)     = 0.d0
    dphi(:)    = 0.d0
    psiaxis(:) = 0.d0
    dpsi       = ( psimin - psilim ) / dble(Nq-1)
    do k=1, Nq
       psiaxis(k) = dpsi * dble(k-1) + psilim
    enddo ! -- x-Axis -- !
    
    !  --- [2]  Toroidal Flux ---  !
    phi(:)     = 0.d0
    dS         = dr * dz
    do k=1, Nq
       do j=1, Nr
          do i=1, Nz
             if ( psi(i,j).le.psiaxis(k) ) then
                phi(k) = phi(k) + Bt(i,j)*dS
             endif
          enddo
       enddo
    enddo
    
    !  --- [3]  Safety Factor ---  !
    coef       = 1.d0 / ( dpsi * 2.d0 * pi )
    do k=1, Nq-1
       qfactor(k) = ( phi(k+1) - phi(k) ) * coef
    enddo    
    q95        = qfactor(1)
    qaxis      = qfactor(Nq-1)
    
    !  --- [3]  Write  ---  !
    write(6,*)
    write(6,*) '    ----        q-Profile  Analysis         ----    '
    write(6,*)
    write(6,'(6x,a32,2x,e13.6)') ' q95      :: ', q95
    write(6,'(6x,a32,2x,e13.6)') ' qaxis    :: ', qaxis
    write(6,'(6x,a32,2x,e13.6)') ' psilim   :: ', psilim
    write(6,'(6x,a32,2x,e13.6)') ' psimin   :: ', psimin
    write(6,*)
    write(6,*)
    open( 21,file=trim(jobdir)//'qprofile.dat',form='formatted',status='replace')
    do k=1, Nq-1
       psil(k) =   0.5d0*( psiaxis(k)+psiaxis(k+1) )
       psis(k) = ( 0.5d0*( psiaxis(k)+psiaxis(k+1) ) - psilim ) / ( psimin - psilim )
       write(21,'(4(e15.8,1x))') psil(k), psis(k), qfactor(k), phi(k)
       ! write( 6,'(5(e12.5,1x))') psiaxis(k), psil(k), psis(k), qfactor(k), phi(k)
    enddo
    close(21)

    open(20,file=trim(jobdir)//'equilibrium.dat',form='formatted',position='append')
    write(20,'(a20,1x,e13.6)'  ) ' q95'      , q95
    write(20,'(a20,1x,e13.6)'  ) ' qaxis'    , qaxis
    close(20)

    return
  end subroutine q_profile

  
  function Perimeter( map, N1, N2, dx1, dx2 )

    implicit none
    integer         , intent(in) :: N1, N2
    double precision, intent(in) :: map(N1,N2), dx1, dx2
    double precision             :: Perimeter
    integer                      :: i, j, Nisec(N1), iEdge(2)
    integer, parameter           :: jBtm = 3, jTop = 3
    integer, parameter           :: X    = 1,    Y = 2
    integer, parameter           :: Lwr  = 1,  Upr = 2
    double precision             :: lLeft, lRight, lTop, lBot
    double precision             :: Xisec(2,2,N1), Xnext(2,2)

    ! --- [1] InterSection Search ( Map = 0.0 ) --- !
    do i=1, N1
       Nisec(i) = 0
    enddo
    do i=1, N1
       do j=jBtm, N2-jTop
          if ( map(i,j)*map(i,j-1).le.0.d0 ) then
             Nisec(i) = Nisec(i) + 1
             if ( Nisec(i).le.2 ) then
                Xisec(X,Nisec(i),i) = dx1*dble(i-1)
                Xisec(Y,Nisec(i),i) = dx2*dble(j-1)  &
                     &  + abs(map(i,j-1)) / ( abs(map(i,j-1)) + abs(map(i,j)) ) * dx2
             endif
          endif
       enddo
    enddo

    ! --- [2] Calculation of Perimeter --- !
    !  -- [2-1] Edge Search -- !
    do i=1, N1, +1
       if ( Nisec(i).eq.2 ) then
          iEdge(1) = i
          exit
       endif
    enddo
    do i=N1, 1, -1
       if ( Nisec(i).eq.2 ) then
          iEdge(2) = i
          exit
       endif
    enddo
    !  -- [2-2] Left & Right  -- !
    lLeft  = ( Xisec(Y,Upr,iEdge(1)) - Xisec(Y,Lwr,iEdge(1)) )
    lRight = ( Xisec(Y,Upr,iEdge(2)) - Xisec(Y,Lwr,iEdge(2)) )
    !  -- [2-3] Top  & Bottom -- !
    Xnext(X,Lwr) = Xisec(X,Lwr,iEdge(1))
    Xnext(Y,Lwr) = Xisec(Y,Lwr,iEdge(1))
    Xnext(X,Upr) = Xisec(X,Upr,iEdge(1))
    Xnext(Y,Upr) = Xisec(Y,Upr,iEdge(1))
    lBot = 0.d0
    lTop = 0.d0
    do i=iEdge(1)+1, iEdge(2)
       if ( Nisec(i).eq.2 ) then
          lBot = lBot + sqrt( ( Xisec(X,Lwr,i) - Xnext(X,Lwr) )**2 &
               &            + ( Xisec(Y,Lwr,i) - Xnext(Y,Lwr) )**2 )
          lTop = lTop + sqrt( ( Xisec(X,Upr,i) - Xnext(X,Upr) )**2 &
               &            + ( Xisec(Y,Upr,i) - Xnext(Y,Upr) )**2 )
          Xnext(X,Lwr) = Xisec(X,Lwr,i)
          Xnext(Y,Lwr) = Xisec(Y,Lwr,i)
          Xnext(X,Upr) = Xisec(X,Upr,i)
          Xnext(Y,Upr) = Xisec(Y,Upr,i)
       endif
    enddo

    ! --- [3] Return Perimeter --- !
    Perimeter = lLeft + lRight + lBot + lTop

    return
  end function Perimeter

end module diagnosis
