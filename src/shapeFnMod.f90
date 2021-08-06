module shapeFnMod
contains
  
  function shapef_r( iposit, rposit, drinv )
    implicit none
    integer         , intent(in) :: iposit
    double precision, intent(in) :: rposit
    double precision, intent(in) :: drinv
    double precision             :: delta
    double precision             :: shapef_r(5)

    ! ----------------------------------- !
    ! --- [1]  CIC  (1st-order)       --- !
    ! ----------------------------------- !
    ! delta       = rposit*drinv - dble( iposit )
    ! shapef_z(1) = 0.d0
    ! shapef_z(2) = 0.d0
    ! shapef_z(3) = 1.d0 - abs(delta)
    ! shapef_z(4) = 0.d0
    ! shapef_z(5) = 0.d0
    ! shapef_z( 3 + int( sign( 1.d0, delta ) ) ) = abs(delta)
    
    ! ----------------------------------- !
    ! --- [2]  2nd-order b-Spline     --- !
    ! ----------------------------------- !
    delta       = rposit*drinv - dble( iposit )
    shapef_r(1) = 0.d0
    shapef_r(2) = 0.50d0*( 0.5d0-delta )**2
    shapef_r(3) = 0.75d0 - delta*delta
    shapef_r(4) = 0.50d0*( 0.5d0+delta )**2
    shapef_r(5) = 0.d0

    ! ----------------------------------- !
    ! --- [3]  Cylindrical (1st-order) -- !
    ! ----------------------------------- !
    ! rmndr       = rMin * drinv
    ! rptdr       = ( rposit + rMin ) * drinv
    ! nsign       = sign( 1.d0, rposit*drinv - dble( iposit ) )
    ! delta       =  0.5d0 * abs( ( rptdr + nsign*0.5d0 )**2 - ( dble( iposit ) + rmndr + 0.5d0*nsign )**2 ) / rptdr
    ! if ( iposit.eq.0 ) then
    !    delta_u = delta * rptdr * 2.d0
    !    delta_n = ( dble(iposit) + rmndr + 0.5d0*nsign  )**2 - rmndr**2
    !    delta_r = ( 2.d0*rmndr   - rptdr + 0.5d0        )**2 - rmndr**2
    !    delta   = delta_u / ( delta_u + delta_n + delta_r )
    ! endif
    ! if ( iposit.eq.LJ-2 ) then
    !    rmxdr   = iposit + rmndr
    !    delta_d = delta * rptdr * 2.d0
    !    delta_n = rmxdr**2 - ( rmxdr-0.5 )**2
    !    delta_r = rmxdr**2 - ( 2.d0 * rmxdr - rptdr - 0.5d0 )**2
    !    delta   = delta_d / ( delta_d + delta_n + delta_r )
    ! endif
    ! shapef_r(1) = 0.d0
    ! shapef_r(2) = 0.d0
    ! shapef_r(3) = 1.d0 - delta
    ! shapef_r(4) = 0.d0
    ! shapef_r(5) = 0.d0
    ! shapef_r( 3 + int( nsign ) ) = delta
    
    return
  end function shapef_r
  
  
  function shapef_z( iposit, rposit, dzinv )
    implicit none
    integer         , intent(in) :: iposit
    double precision, intent(in) :: rposit
    double precision, intent(in) :: dzinv
    double precision             :: shapef_z(5)
    double precision             :: delta

    ! ----------------------------------- !
    ! --- [1]  CIC  (1st-order)       --- !
    ! ----------------------------------- !
    ! delta       = rposit*dzinv - dble( iposit )
    ! shapef_z(1) = 0.d0
    ! shapef_z(2) = 0.d0
    ! shapef_z(3) = 1.d0 - abs(delta)
    ! shapef_z(4) = 0.d0
    ! shapef_z(5) = 0.d0
    ! shapef_z( 3 + int( sign( 1.d0, delta ) ) ) = abs(delta)
    
    ! ----------------------------------- !
    ! --- [2]  2nd-order b-Spline     --- !
    ! ----------------------------------- !
    delta       = rposit*dzinv - dble( iposit )
    shapef_z(1) = 0.d0
    shapef_z(2) = 0.50d0*( 0.5d0-delta )**2
    shapef_z(3) = 0.75d0 - delta*delta
    shapef_z(4) = 0.50d0*( 0.5d0+delta )**2
    shapef_z(5) = 0.d0
    
    return
  end function shapef_z

end module shapeFnMod
