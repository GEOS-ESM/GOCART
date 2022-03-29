   real function soilMoistureConvertVol2Grav(vsoil, sandfrac, rhop)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in) :: vsoil       ! volumetric soil moisture fraction [1]
   real, intent(in) :: sandfrac    ! fractional sand content [1]
   real, intent(in) :: rhop        ! dry dust density [kg m-3]

! !DESCRIPTION: Convert soil moisture fraction from volumetric to gravimetric.
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
   real :: vsat

!  !CONSTANTS:
   real, parameter :: rhow = 1000.    ! density of water [kg m-3]

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Saturated volumetric water content (sand-dependent) ! [m3 m-3]
   vsat = 0.489 - 0.00126 * ( 100. * sandfrac )

!  Gravimetric soil content
   soilMoistureConvertVol2Grav = vsoil * rhow / (rhop * (1. - vsat))

   end function soilMoistureConvertVol2Grav
