module PointEmissionGOCART
   implicit none
   private

   public DistributePointEmission
   public updatePointwiseEmissions

   public find_dz
   public find_level
   public find_wieghts
contains
   function find_dz(km, hghte) result(dz)
      real                :: dz(km)
      integer, intent(in) :: km
      real,    intent(in) :: hghte(:)

      integer :: k

      do k = km, 1, -1
         dz(k) = hghte(k)-hghte(k+1)
      end do
   end function find_dz

   integer function find_level(km, z, z_level)
      integer, intent(in) :: km
      real,    intent(in) :: z(:)
      real,    intent(in) :: z_level

      integer :: k

      do k = km, 1, -1
         if (z(k) >= z_level) then
            find_level = k
            exit
         end if
      end do
   end function find_level

   subroutine find_wieghts(k_bot, k_top, z_bot, z_top, z, dz, w)
      integer, intent(in) :: k_bot
      integer, intent(in) :: k_top
      real,    intent(in) :: z_bot
      real,    intent(in) :: z_top
      real,    intent(in) :: z(:)
      real,    intent(in) :: dz(:)

      real,    intent(inout) :: w(:)

      integer :: k

      do k = k_bot, k_top, -1
         if ((k < k_bot) .and. (k > k_top)) then
            w(k) = dz(k)
         else
            if (k == k_bot) then
               w(k) = (z(k) - z_bot)
            end if

            if (k == k_top) then
               w(k) = z_top - (z(k)-dz(k))
            end if
         end if
      end do
   end subroutine find_wieghts

   subroutine DistributePointEmission(km, hghte, z_bot, z_top, &
                                      emissions_point, area, &
                                      point_column_emissions, rc)
      ! !USES:
      implicit NONE

      ! !INPUT PARAMETERS:
      integer,            intent(in)  :: km       ! total model levels
      real, dimension(:), intent(in)  :: hghte    ! model level geopotential height [m]
      real,               intent(in)  :: z_bot, z_top ! base and top altitude respectively
      real,               intent(in)  :: area     ! grid cell area [m^2]
      real,               intent(in)  :: emissions_point ![kg/kg]


      ! !OUTPUT PARAMETERS:
      real, dimension(:), intent(out) ::  point_column_emissions ![kg/kg]
      integer, optional, intent(out) :: rc                       ! Error return code:


      ! !DESCRIPTION: Distributes piont emissions uniformily in the vertical in height coordinates.
      !
      ! !REVISION HISTORY:
      ! ??? A. Darmenov
      ! ??? P. Colarco
      !
      ! !Locals
      integer :: k_bot, k_top
      real, dimension(km) :: z, dz, w_ !dz units = meters

      !EOP
      !-------------------------------------------------------------------------
      ! Begin

      z(1:km) = hghte(1:km)

      k_bot = find_level(km, z, z_bot)
      k_top = find_level(km, z, z_top)

      !   find the weights
      w_ = 0

      if (k_bot == k_top) then
         if (z_top == z_bot) then ! for non-explosive volcanic emissions
            w_(k_bot) = tiny(0.)
         else
            w_(k_bot) = z_top - z_bot
         end if
      else
         dz = find_dz(km, hghte)
         call find_wieghts(k_bot, k_top, z_bot, z_top, z, dz, w_)
      end if

      point_column_emissions(:) = (w_ / sum(w_)) * emissions_point

      !  Logical consistency check
      if (k_top > k_bot) then
         rc = 1
      else
         rc = 0
      end if
   end subroutine DistributePointEmission

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
                                     point_column_emissions, rc)

        emissions_point(i,j,:) =  point_column_emissions
        end do
   rc = 0

  end subroutine updatePointwiseEmissions
end module PointEmissionGOCART
