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
! ??2020 E.Sherman - ported to process library
!
! !Locals
   integer :: k
   integer :: k_bot, k_top
   real    :: z_
   real, dimension(km) :: z, dz, w_ !dz units = meters

!EOP
!-------------------------------------------------------------------------
! Begin

!    z(1:km) = hghte(0:km-1)
    z(1:km) = hghte(1:km)

    do k = km, 1, -1
!       dz(k) = hghte(k-1)-hghte(k)
       dz(k) = hghte(k)-hghte(k+1)
    end do

!   find the bottom level
    do k = km, 1, -1
       if (z(k) >= z_bot) then
           k_bot = k
           exit
       end if
    end do

!   find the top level
    do k = k_bot, 1, -1
       if (z(k) >= z_top) then
           k_top = k
           exit
       end if
    end do

!   find the weights
    w_ = 0

!   if (k_top > k_bot) then
!       need to bail - something went wrong here
!   end if

    if (k_bot == k_top) then
       if (z_top == z_bot) then ! for non-explosive volcanic emissions
          w_(k_bot) = tiny(0.)
       else
          w_(k_bot) = z_top - z_bot
       end if
    else
     do k = k_bot, k_top, -1
        if ((k < k_bot) .and. (k > k_top)) then
             w_(k) = dz(k)
        else
             if (k == k_bot) then
                 w_(k) = (z(k) - z_bot)
             end if

             if (k == k_top) then
                 w_(k) = z_top - (z(k)-dz(k))
             end if
        end if
     end do
    end if

!   distribute emissions in the vertical
    point_column_emissions(:) = ((w_ / sum(w_)) * emissions_point) / area

      __RETURN__(__SUCCESS__)
    end subroutine DistributePointEmission
