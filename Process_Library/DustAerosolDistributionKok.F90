   subroutine DustAerosolDistributionKok ( radius, rLow, rUp, distribution )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, dimension(:), intent(in)  :: radius      ! Dry particle bin effective radius [um]
   real, dimension(:), intent(in)  :: rLow, rUp   ! Dry particle bin edge radii [um]

! !OUTPUT PARAMETERS:
   real, dimension(:), intent(out) :: distribution    ! Normalized dust aerosol distribution [1]

! !DESCRIPTION: Computes lognormal aerosol size distribution for dust bins according to
!               J.F.Kok, PNAS, Jan 2011, 108 (3) 1016-1021; doi:10.1073/pnas.1014798108
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 01Apr2021 R.Montuoro/NOAA - Refactored for GOCART process library
!

! !Local Variables
   integer :: n, nbins
   real    :: diameter, dlam, dvol

! !CONSTANTS
   real, parameter    :: mmd    = 3.4          ! median mass diameter [um]
   real, parameter    :: stddev = 3.0          ! geometric standard deviation [1]
   real, parameter    :: lambda = 12.0         ! crack propagation length [um]
   real, parameter    :: factor = 1.e0 / (sqrt(2.e0) * log(stddev))  ! auxiliary constant

   character(len=*), parameter :: myname = 'DustAerosolDistributionKok'

!EOP
!-------------------------------------------------------------------------
!  Begin...

   distribution = 0.

!  Assume all arrays are dimensioned consistently
   nbins = size(radius)

   dvol = 0.
   do n = 1, nbins
     diameter = 2 * radius(n)
     dlam = diameter/lambda
     distribution(n) = diameter * (1. + erf(factor * log(diameter/mmd))) &
                     * exp(-dlam * dlam * dlam) * log(rUp(n)/rLow(n))
     dvol = dvol + distribution(n)
   end do

!  Normalize distribution
   do n = 1, nbins
     distribution(n) = distribution(n) / dvol
   end do

   end subroutine DustAerosolDistributionKok
