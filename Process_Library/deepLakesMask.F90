   subroutine deepLakesMask (lons, lats, radToDeg, deep_lakes_mask, rc)

! !USES:
   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:,:), intent(inout) :: deep_lakes_mask

! !INPUT PARAMETERS:
   real, pointer, dimension(:,:), intent(in)   :: lats ! latitude [radians]
   real, pointer, dimension(:,:), intent(in)   :: lons ! longtitude [radians]
   real, intent(in)                       :: radToDeg

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc
!EOP

! !Local Variables
    real                :: dummylon
    integer             :: i, j
!EOP
!-------------------------------------------------------------------------
!  Begin...

   deep_lakes_mask = 1.0
   do j = 1, ubound(lons, 2)
      do i = 1, ubound(lons, 1)
                           dummylon = lons(i,j)*radToDeg
      if( dummylon < 0.0 ) dummylon = dummylon + 360.0
      ! The Great Lakes: lon = [91W,75W], lat = [40.5N, 50N]
      if ((dummylon > 267.0) .and. &
          (dummylon < 285.0) .and. &
          (lats(i,j)*radToDeg >  40.5) .and. &
          (lats(i,j)*radToDeg <  50.0)) deep_lakes_mask(i,j) = 0.0

       ! The Caspian Sea: lon = [45.0, 56], lat = 35, 48]
      if ((dummylon >  45.0) .and. &
          (dummylon <  56.0) .and. &
          (lats(i,j)*radToDeg >  35.0) .and. &
          (lats(i,j)*radToDeg <  48.0)) deep_lakes_mask(i,j) = 0.0
      end do
   end do

   __RETURN__(__SUCCESS__)
   end subroutine deepLakesMask
