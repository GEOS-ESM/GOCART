   real function moistureCorrectionFecan(slc, sand, clay, rhop)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in) :: slc     ! liquid water content of top soil layer, volumetric fraction [1]
   real, intent(in) :: sand    ! fractional sand content [1]
   real, intent(in) :: clay    ! fractional clay content [1]
   real, intent(in) :: rhop    ! dry dust density [kg m-3]

! !DESCRIPTION: Compute correction factor to account for Fecal soil moisture
!
! !REVISION HISTORY:
!
!  02Apr2020, B.Baker/NOAA    - Original implementation
!  01Apr2020, R.Montuoro/NOAA - Adapted for GOCART process library

!  !Local Variables
   real :: grvsoilm
   real :: drylimit

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Convert soil moisture from volumetric to gravimetric
   grvsoilm = soilMoistureConvertVol2Grav(slc, sand, rhop)

!  Compute fecan dry limit
   drylimit = clay * (14.0 * clay + 17.0)

!  Compute soil moisture correction
   moistureCorrectionFecan = sqrt(1.0 + 1.21 * max(0., grvsoilm - drylimit)**0.68)

   end function moistureCorrectionFecan
