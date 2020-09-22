#define __SUCCESS__ 0
#define __VERIFY__(x) if(x/=0) then; if(present(rc)) rc=x; return; endif
#define __RC__ rc=status); __VERIFY__(status
#define __STAT__ stat=status); __VERIFY__(status
#define __IOSTAT__ iostat=status); __VERIFY__(status
#define __RETURN__(x) if (present(rc)) rc=x; return
#define __ASSERT__(expr) if(.not. (expr)) then; if (present(rc)) rc=-1; return; endif

module DustEmissionGOCART
   implicit none
   private

   public :: DustEmissionGOCART2G
   public :: u_threshold
   public :: bin_emissions
   public :: OCEAN, LAND

   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
contains
   real function u_threshold(radius, grav)
      real,    intent(in) :: radius ! particle radius [m]
      real,    intent(in) :: grav   ! gravity [m/sec^2]

      real, parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
      real, parameter ::  soil_density  = 2650.  ! km m-3
      real            ::  diameter          ! dust effective diameter [m]

      !  Calculate the threshold velocity of wind erosion [m/s] for each radius
      !  for a dry soil, as in Marticorena et al. [1997].
      !  The parameterization includes the air density which is assumed
      !  = 1.25 kg m-3 to speed the calculation.  The error in air density is
      !  small compared to errors in other parameters.

      diameter = 2. * radius

       u_threshold = 0.13 * sqrt( &
            & (soil_density*grav*diameter/air_dens) &
            & * (1.+6.e-7/(soil_density*grav*diameter**2.5)) &
            & / (1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.) &
            & )

   end function u_threshold

   function bin_emissions(u_thresh0, oro, u10m, v10m, gwettop, fraclake) result(emissions)
      real, allocatable :: emissions(:,:)
      real, intent(in) :: u_thresh0
      real, intent(in) :: oro(:,:)
      real, intent(in) :: u10m(:,:)
      real, intent(in) :: v10m(:,:)
      real, intent(in) :: gwettop(:,:)
      real, intent(in) :: fraclake(:,:)

      real    ::  w10m2, u_thresh
      integer :: i, j, ni, nj

      !  Get dimensions
      !  ---------------

      ni = size(u10m, 1)
      nj = size(u10m, 2)

      allocate(emissions(ni, nj))
      emissions = 0.0

      !     Spatially dependent part of calculation
      !     ---------------------------------------
      do j = 1, nj
         do i = 1, ni

            if ( oro(i,j) /= LAND ) cycle ! only over LAND gridpoints

            w10m2 = u10m(i,j)**2 + v10m(i,j)**2
            !           Modify the threshold depending on soil moisture as in Ginoux et al. [2001]

            if (gwettop(i,j) < 0.5) then
               u_thresh = u_thresh0*(1.2 + 0.2*log10(max(1.e-3,gwettop(i,j))))

               if (w10m2 > u_thresh**2) then
                  !     Emission of dust [kg m-2 s-1]
                  emissions(i,j) = (1.-fraclake(i,j)) * w10m2 * (sqrt(w10m2)-u_thresh)
               endif
            endif !(gwettop(i,j) < 0.5)
         end do ! i
      end do ! j
   end function bin_emissions

   subroutine DustEmissionGOCART2G(radius, fraclake, gwettop, oro, u10m, &
                                   v10m, Ch_DU, du_src, grav, &
                                   emissions, rc )
      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      real,          intent(in) :: radius(:)     ! particle radius [m]
      real, pointer, intent(in) :: fraclake(:,:) ! fraction of lake [1]
      real, pointer, intent(in) :: gwettop(:,:)  ! surface soil wetness [1]
      real, pointer, intent(in) :: oro(:,:)      ! land-ocean-ice mask [1]
      real, pointer, intent(in) :: u10m(:,:)     ! 10-meter eastward wind [m/sec]
      real, pointer, intent(in) :: v10m(:,:)     ! 10-meter northward wind [m/sec]
      real, pointer, intent(in) :: du_src(:,:)   ! dust emissions [(sec^2 m^5)/kg]
      real,          intent(in) :: Ch_DU         ! dust emission tuning coefficient [kg/(sec^2 m^5)]
      real,          intent(in) :: grav          ! gravity [m/sec^2]

      ! !OUTPUT PARAMETERS:
      real, pointer, intent(inout) :: emissions(:,:) ! Local emission [kg/(m^2 sec)]
      integer,        intent(out)  :: rc             ! Error return code:


      ! !DESCRIPTION: Computes the dust emissions for one time step
      !
      ! !REVISION HISTORY:
      !
      ! 11Feb2020 E.Sherman - First attempt at refactor
      !

      ! !Local Variables
      integer :: n, nbins
      real    :: u_thresh0

      !EOP
      !-------------------------------------------------------------------------
      !  Begin

      nbins = size(radius)
      emissions(:,:) = 0.

      do n = 1, nbins
         u_thresh0 = u_threshold(radius(n), grav)
         emissions = emissions + &
            Ch_DU * du_src * bin_emissions(u_thresh0, oro, u10m, v10m, gwettop, fraclake)
      end do

      rc = 0
   end subroutine DustEmissionGOCART2G
end module DustEmissionGOCART
