  subroutine updatePointwiseEmissions (km, pBase, pTop, pEmis, nPts, pStart, &
                                       pEnd, hghte, area, &
                                       iPoint, jPoint, nhms, emissions_point, rc)
    implicit none

!   !ARGUMENTS:
    integer,                intent(in)  :: km     ! total model levels
    real, dimension(:),     intent(in)  :: pBase  ! base altitude (e.g., bottom of plume)
    real, dimension(:),     intent(in)  :: pTop   ! top altitude (e.g., top of plume)
    real, dimension(:),     intent(in)  :: pEmis  ! emission flux (e.g., kg/sec of species)
    integer,                intent(in)  :: nPts   ! number of events in file
    integer, dimension(:),  intent(in)  :: pStart ! HHMMSS to start emissions
    integer, dimension(:),  intent(in)  :: pEnd   ! HHMMSS to end emissions
    real, dimension(:,:,:), intent(in)  :: hghte  ! model level geopotential height [m]
    real, dimension(:,:),   intent(in)  :: area   ! grid cell area [m^2]
    integer, dimension(:),  intent(in)  :: iPoint ! i dimension location of emission on grid
    integer, dimension(:),  intent(in)  :: jPoint ! j dimension location of emission on grid
    integer,                intent(in)  :: nhms   ! model hour mintue second
    real, dimension(:,:,:), intent(inout)  :: emissions_point ![kg/kg]
    integer, optional,      intent(out)  :: rc  ! return code

!   !Local
    real, dimension(km)              :: point_column_emissions
    integer                          :: n, i, j
    real, dimension(:), allocatable  :: pEmis_

    integer :: status

!   Description: Returns 3D array of pointwise emissions.
!
!   Revision History:
!EOP
!-----------------------------------------------------------------------------
!    Begin...

     pEmis_ = pEmis

     do n = 1, nPts
        i = iPoint(n)
        j = jPoint(n)
        if( i<1 .OR. j<1 ) cycle    ! Point emission not in this sub-domain
        ! Emissions not occurring in current time step
        if(nhms < pStart(n) .or. nhms >= pEnd(n)) cycle

        call DistributePointEmission(km, hghte(i,j,:), pBase(n), &
                                     pTop(n), pEmis_(n), area(i,j), &
                                     point_column_emissions, __RC__)

        emissions_point(i,j,:) =  point_column_emissions
        end do

      __RETURN__(__SUCCESS__)
  end subroutine updatePointwiseEmissions
