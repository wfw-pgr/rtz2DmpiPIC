program main
  ! --- Constants & Variables --- !
  use constants , only : N1, N2, coordinate, x2min, RgridMode, BgridMode, myRank, PEtot, sftype
  use variables , only : Avp, rhs, dx1, dx2, x2
  ! --- Subroutines           --- !
  use cylPoisMod, only : BiCGSTABCyl2D_MPI
  use xyzPoisMod, only : BiCGSTABCrt2D_MPI
  use RgridRZMod, only : RgridField
  use BgridRZMod, only : BgridField
  use jxBchckMod, only : jxB_Force, checkNan
  use utilityMod, only : makeJobFile
  use InitialMod, only : InitVariables, setJBRHS
  use writingMod, only : WriteField, WriteParam, WriteF_MPI
  use lInterpMod, only : LinearInterpFWD, LinearInterpRET
  implicit none
  include "mpif.h"
  integer             :: ierr
  double precision    :: dx2Np1, x2Np1(0:N2), rhsNp1(N1,0:N2), AvpNp1(N1,0:N2)

  ! ------------------------------------- !
  ! --- [1] Preparation               --- !
  ! ------------------------------------- !
  !  -- [1-1] MPI Initialize          --  !
  call MPI_Init( ierr )
  call MPI_Comm_Rank( MPI_COMM_WORLD, myRank, ierr )
  call MPI_Comm_Size( MPI_COMM_WORLD, PEtot , ierr )
  !  -- [1-2] Opening Remark          --  !
  if ( myRank.eq.0 ) then
     write(6,*)
     write(6,'(2x,a)') "======================================================================"
     write(6,'(2x,a)') "=================     Flux-Tube Equilibrium Solver    ================"
     write(6,'(2x,a)') "======================================================================"
     write(6,*)
  endif
  call makeJobFile
  call InitVariables
  call setJBRHS

  ! ------------------------------------- !
  ! --- [2] Main Solver               --- !
  ! ------------------------------------- !
  !  -- [2-1]  (xyz) Cartesian Solver --  !
  if ( coordinate .eq. 'XZ' ) then
     call BiCGSTABCrt2D_MPI( Avp, rhs, dx1, dx2, N1, N2 )
  endif
  !  -- [2-2]  (rtz) Cartesian Solver --  !
  if ( coordinate .eq. 'RZ' ) then
     ! -- (1) Sinply Connected Case   --  !
     if ( x2min.eq.0.d0 ) then
        call LinearInterpFWD( x2, rhs, x2Np1, rhsNp1, N1, N2, sftype )
        dx2Np1         = x2Np1(2) - x2Np1(1)
        call BiCGSTABCyl2D_MPI( AvpNp1(1:N1,1:N2), rhsNp1(1:N1,1:N2), dx1, dx2Np1, N1, N2, x2min )
        AvpNp1(1:N1,0) = AvpNp1(1:N1,1)
        call LinearInterpRET( x2Np1, AvpNp1, x2, Avp, N1, N2 )
     endif
     ! -- (2) Doubly Connected Case   --  !
     if ( x2min.gt.0.d0 ) then
        call BiCGSTABCyl2D_MPI( Avp, rhs, dx1, dx2, N1, N2, x2min )
     endif
     
  end if
  !  -- [2-3] Check NaN & Loop End    --  !
  if ( myRank.eq.0 ) then
     call checkNan
     write(6,*)
     write(6,*)
     write(6,'(2x,a)') "----------------------------------------------------------------------"
     write(6,'(2x,a)') "--------------   Flux-Tube Equilibrium Solver - END -    -------------"
     write(6,'(2x,a)') "----------------------------------------------------------------------"
     write(6,*)
  endif

  ! ------------------------------------- !
  ! --- [3] Post Process              --- !
  ! ------------------------------------- !
  !  -- [3-1] Bgrid Case              --  !
  if ( ( myRank.eq.0 ).and.( BgridMode ) ) then
     write(6,*)
     write(6,'(2x,a)') '++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(6,'(2x,a)') '++   Grid System :: B-grid ( Staggered Grid )   ++'
     write(6,'(2x,a)') '++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(6,*)
     call BgridField( coordinate, 'B' )
     call WriteField( 'B' )
     call WriteParam( 'B' )
     call WriteF_MPI( 'B' )
     call jxB_Force ( 'B' )
  endif
  !  -- [3-2] Rgrid Case              --  !
  if ( ( myRank.eq.0 ).and.( RgridMode ) ) then
     write(6,*)
     write(6,'(2x,a)') '++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(6,'(2x,a)') '++   Grid System :: R-grid (   Regular Grid )   ++'
     write(6,'(2x,a)') '++++++++++++++++++++++++++++++++++++++++++++++++++'
     write(6,*)
     call RgridField( coordinate, 'R' )
     call WriteField( 'R' )
     call WriteParam( 'R' )
     call jxB_Force ( 'R' )
  endif
  call MPI_Barrier( MPI_COMM_WORLD, ierr )
  !  -- [3-3] Display Program End     --  !
  if ( myRank.eq.0 ) then
     write(6,*)
     write(6,'(2x,a)') "======================================================================"
     write(6,'(2x,a)') "=============     Flux-Tube Equilibrium Solver  - END -   ============"
     write(6,'(2x,a)') "======================================================================"
     write(6,*)
  endif
  call MPI_Finalize( ierr )
  
end program main
