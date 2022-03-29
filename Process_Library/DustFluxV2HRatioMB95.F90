   real function DustFluxV2HRatioMB95(clay, kvhmax)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in) :: clay      ! fractional clay content [1]
   real, intent(in) :: kvhmax    ! maximum flux ratio [1]

!  !CONSTANTS:
   real, parameter :: clay_thresh = 0.2    ! clay fraction above which the maximum flux ratio is returned

! !DESCRIPTION: Computes the vertical-to-horizontal dust flux ratio according to
!               B.Marticorena, G.Bergametti, J.Geophys.Res., 100(D8), 16415–16430, 1995
!               doi:10.1029/95JD00690
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Adapted for GOCART process library
!
!EOP
!-------------------------------------------------------------------------
!  Begin...

   if (clay > clay_thresh) then
     DustFluxV2HRatioMB95 = kvhmax
   else
     DustFluxV2HRatioMB95 = 10.0**(13.4*clay-6.0)
   end if

   end function DustFluxV2HRatioMB95
