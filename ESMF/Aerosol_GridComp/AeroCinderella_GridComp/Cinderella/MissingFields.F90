#include "MAPL_Generic.h"

module MissingFields_mod
   use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

   use ESMF
   use MAPL

   use GEOS_UtilsMod

   implicit none
   private

   public get_AIRDENS
   public get_DELP
   public get_ZLE

   logical, parameter :: ice_flip   = .false.
   logical, parameter :: lake_flip  = .false.
   logical, parameter :: ocean_flip = .false.

   interface get_AIRDENS
      module procedure :: get_AIRDENS_real32
      module procedure :: get_AIRDENS_real64
   end interface get_AIRDENS

   interface get_DELP
      module procedure :: get_DELP_real32
      module procedure :: get_DELP_real64
   end interface get_DELP

   interface get_FRACI
      module procedure :: get_FRACI_real32
      module procedure :: get_FRACI_real64
   end interface get_FRACI

   interface get_FRLAKE
      module procedure :: get_FRLAKE_real32
      module procedure :: get_FRLAKE_real64
   end interface get_FRLAKE

   interface get_FROCEAN
      module procedure :: get_FROCEAN_real32
      module procedure :: get_FROCEAN_real64
   end interface get_FROCEAN

   interface get_ZLE
      module procedure :: get_ZLE_real32
      module procedure :: get_ZLE_real64
   end interface get_ZLE
contains
   subroutine get_AIRDENS_real32(PLE, T, AIRDENS)
      real(kind=REAL32), intent(in   ) :: PLE(:,:,:)
      real(kind=REAL32), intent(in   ) :: T(:,:,:)
      real(kind=REAL32), intent(  out) :: AIRDENS(:,:,:)

      AIRDENS = PLE/(MAPL_RDRY*T)
   end subroutine get_AIRDENS_real32

   subroutine get_AIRDENS_real64(PLE, T, AIRDENS)
      real(kind=REAL64), intent(in   ) :: PLE(:,:,:)
      real(kind=REAL64), intent(in   ) :: T(:,:,:)
      real(kind=REAL64), intent(  out) :: AIRDENS(:,:,:)

      AIRDENS = PLE/(MAPL_RDRY*T)
   end subroutine get_AIRDENS_real64

   subroutine get_AIRDENS_DRYP()
   end subroutine get_AIRDENS_DRYP

   subroutine get_DELP_real32(PRESS, DELP)
      real(kind=REAL32), intent(in   ) :: PRESS(:,:,:)
      real(kind=REAL32), intent(  out) :: DELP(:,:,:)

      integer :: i, n(3)

      n = shape(PRESS)

      do i=1, n(3)-1
         DELP(:,:,i) =  PRESS(:,:,i+1) - PRESS(:,:,i)
      end do
   end subroutine get_DELP_real32

   subroutine get_DELP_real64(PRESS, DELP)
      real(kind=REAL64), intent(in   ) :: PRESS(:,:,:)
      real(kind=REAL64), intent(  out) :: DELP(:,:,:)

      integer :: i, n(3)

      n = shape(PRESS)

      do i=1, n(3)-1
         DELP(:,:,i) =  PRESS(:,:,i+1) - PRESS(:,:,i)
      end do
   end subroutine get_DELP_real64

   subroutine get_FRACI_real32(sea_ice, FRACI, flip)
      real(kind=REAL32), intent(in   ) :: sea_ice
      real(kind=REAL32), intent(  out) :: FRACI
      logical, optional, intent(in   ) :: flip

      logical :: flip_values

      if (present(flip)) then
         flip_values = flip
      else
         flip_values = ice_flip
      end if

      if (flip_values) then
         FRACI = 1.0 - sea_ice
      else
         FRACI = sea_ice
      end if
   end subroutine get_FRACI_real32

   subroutine get_FRACI_real64(sea_ice, FRACI, flip)
      real(kind=REAL64), intent(in   ) :: sea_ice
      real(kind=REAL64), intent(  out) :: FRACI
      logical, optional, intent(in   ) :: flip

      logical :: flip_values

      if (present(flip)) then
         flip_values = flip
      else
         flip_values = ice_flip
      end if

      if (flip_values) then
         FRACI = 1.0 - sea_ice
      else
         FRACI = sea_ice
      end if
   end subroutine get_FRACI_real64

   subroutine get_FRLAKE_real32(lake_land, FRLAKE, flip)
      real(kind=REAL32), intent(in   ) :: lake_land
      real(kind=REAL32), intent(  out) :: FRLAKE
      logical, optional, intent(in   ) :: flip

      logical :: flip_values

      if (present(flip)) then
         flip_values = flip
      else
         flip_values = lake_flip
      end if

      if (flip_values) then
         FRLAKE = 1.0 - lake_land
      else
         FRLAKE = lake_land
      end if
   end subroutine get_FRLAKE_real32

   subroutine get_FRLAKE_real64(lake_land, FRLAKE, flip)
      real(kind=REAL64), intent(in   ) :: lake_land
      real(kind=REAL64), intent(  out) :: FRLAKE
      logical, optional, intent(in   ) :: flip

      logical :: flip_values

      if (present(flip)) then
         flip_values = flip
      else
         flip_values = lake_flip
      end if

      if (flip_values) then
         FRLAKE = 1.0 - lake_land
      else
         FRLAKE = lake_land
      end if
   end subroutine get_FRLAKE_real64

   subroutine get_FROCEAN_real32(land_sea, FROCEAN, flip)
      real(kind=REAL32), intent(in   ) :: land_sea
      real(kind=REAL32), intent(  out) :: FROCEAN
      logical, optional, intent(in   ) :: flip

      logical :: flip_values

      if (present(flip)) then
         flip_values = flip
      else
         flip_values = ocean_flip
      end if

      if (flip_values) then
         FROCEAN = 1.0 - land_sea
      else
         FROCEAN = land_sea
      end if
   end subroutine get_FROCEAN_real32

   subroutine get_FROCEAN_real64(land_sea, FROCEAN, flip)
      real(kind=REAL64), intent(in   ) :: land_sea
      real(kind=REAL64), intent(  out) :: FROCEAN
      logical, optional, intent(in   ) :: flip

      logical :: flip_values

      if (present(flip)) then
         flip_values = flip
      else
         flip_values = ocean_flip
      end if

      if (flip_values) then
         FROCEAN = 1.0 - land_sea
      else
         FROCEAN = land_sea
      end if
   end subroutine get_FROCEAN_real64

   subroutine get_PS()
   end subroutine get_PS

   subroutine get_RH2(T, PLE, Q, RH2)
      real(kind=REAL32), intent(in   ) :: T(:,:,:)
      real(kind=REAL32), intent(in   ) :: PLE(:,:,:)
      real(kind=REAL32), intent(in   ) :: Q(:,:,:)
      real(kind=REAL32), intent(  out) :: RH2(:,:,:)

      integer :: i, j, k, n(3)

      n = size(T)

      do i=1, n(1)
         do j=1, n(2)
            do k=1, n(3)
               RH2(i,j,k) = Q(i,j,k) / GEOS_Qsat(T(i,j,k), PLE(i,j,k))
            end do
         end do
      end do

   end subroutine

   subroutine get_DZ()
   end subroutine get_DZ

   subroutine get_ZLE_real32(inst_geop_interface, ZLE)
      real(kind=REAL32), intent(in   ) :: inst_geop_interface(:,:,:)
      real(kind=REAL32), intent(  out) :: ZLE(:,:,:)

      ZLE = inst_geop_interface/MAPL_GRAV
   end subroutine get_ZLE_real32

   subroutine get_ZLE_real64(inst_geop_interface, ZLE)
      real(kind=REAL64), intent(in   ) :: inst_geop_interface(:,:,:)
      real(kind=REAL64), intent(out) :: ZLE(:,:,:)

      ZLE = inst_geop_interface/MAPL_GRAV
   end subroutine get_ZLE_real64

   subroutine get_PFI_LSAN()
   end subroutine get_PFI_LSAN

   subroutine get_PFL_LSAN()
   end subroutine get_PFL_LSAN
end module MissingFields_mod