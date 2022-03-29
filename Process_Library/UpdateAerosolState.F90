   subroutine UpdateAerosolState (emissions, emissions_surface, emissions_point, &
                                  sfrac, nPts, km, cdt, grav, nbins, delp, aero, rc)

! !USES:
  implicit NONE

! !INPUT PARAMETERS:
!   real, pointer, dimension(:,:)     :: emissions_surface
   real, dimension(:,:,:), intent(in)     :: emissions_surface
   real, dimension(:,:,:,:), intent(inout) :: emissions
   real, dimension(:,:,:), intent(in) :: emissions_point

   real, dimension(:), intent(in)  :: sfrac ! source fraction [1]
   integer, intent(in)             :: nPts  ! number of point emissions
   integer, intent(in)             :: km    ! total model levels
   real, intent(in)                :: cdt   ! chemistry model time-step [sec]
   real, intent(in)                :: grav  ! gravity [m/sec^2]
   integer, intent(in)             :: nbins ! number of aerosol size bins
   real, pointer, dimension(:,:,:), intent(in) :: delp  ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:,:), intent(inout)  :: aero ! aerosol [kg/kg]

! !OUTPUT PARAMETERS:
   integer, intent(out)             :: rc          ! Error return code:

! !DESCRIPTION: Updates internal state variables
!
! !REVISION HISTORY:
!
!  15May2020 - E.Sherman
!
! !Local Variables
   integer :: n, kmin


!EOP
!--------------------------------------------------------------------------------
!   Begin...

    rc = __SUCCESS__

    do n = 1, nbins
       emissions(:,:,km,n) = emissions_surface(:,:,n) * sfrac(n)
       if (nPts > 0) then
          kmin = 1
          emissions(:,:,:,n) = emissions(:,:,:,n) + emissions_point * sfrac(n)
       else
          kmin = km
       end if
       aero(:,:,kmin:km,n) = aero(:,:,kmin:km,n) + emissions(:,:,kmin:km,n) * &
                             cdt * grav / delp(:,:,kmin:km)
    end do

   end subroutine UpdateAerosolState
