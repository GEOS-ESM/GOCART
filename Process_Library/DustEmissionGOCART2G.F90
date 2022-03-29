   subroutine DustEmissionGOCART2G(radius, fraclake, gwettop, oro, u10m, &
                                   v10m, Ch_DU, du_src, grav, &
                                   emissions, rc )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in) :: radius(:)       ! particle radius [m]
   real, dimension(:,:), intent(in) :: fraclake ! fraction of lake [1]
   real, dimension(:,:), intent(in) :: gwettop  ! surface soil wetness [1]
   real, dimension(:,:), intent(in) :: oro      ! land-ocean-ice mask [1]
   real, dimension(:,:), intent(in) :: u10m     ! 10-meter eastward wind [m/sec]
   real, dimension(:,:), intent(in) :: v10m     ! 10-meter northward wind [m/sec]
   real, dimension(:,:), intent(in) :: du_src   ! dust emissions [(sec^2 m^5)/kg]
   real, intent(in) :: Ch_DU   ! dust emission tuning coefficient [kg/(sec^2 m^5)]
   real, intent(in) :: grav    ! gravity [m/sec^2]

! !OUTPUT PARAMETERS:
!   real, pointer, intent(inout)  :: emissions(:,:)    ! Local emission [kg/(m^2 sec)]
   real, intent(inout)  :: emissions(:,:,:)    ! Local emission [kg/(m^2 sec)]

   integer, intent(out) :: rc  ! Error return code:


! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
! 11Feb2020 E.Sherman - First attempt at refactor
!

! !Local Variables
   integer         ::  i, j, n
   real, parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
   real, parameter ::  soil_density  = 2650.  ! km m-3
   real            ::  diameter         ! dust effective diameter [m]
   real            ::  u_thresh0
   real            ::  u_thresh
   real            ::  w10m
   integer         ::  i1, i2, j1, j2, nbins
   integer         ::  dims(2)
!   real, allocatable ::  emissions_(:,:)

!EOP
!-------------------------------------------------------------------------
!  Begin

!  Initialize local variables
!  --------------------------
!   emissions(:,:,:) = 0.
!  Get dimensions
!  ---------------
   nbins = size(radius)
   dims = shape(u10m)
   i1 = 1; j1 = 1
   i2 = dims(1); j2 = dims(2)

!   allocate(emissions_(i2,j2))

!  Calculate the threshold velocity of wind erosion [m/s] for each radius
!  for a dry soil, as in Marticorena et al. [1997].
!  The parameterization includes the air density which is assumed
!  = 1.25 kg m-3 to speed the calculation.  The error in air density is
!  small compared to errors in other parameters.

!print*,'DustEmiss shape(emissions) = ',shape(emissions)

   do n = 1, nbins
      diameter = 2. * radius(n)

      u_thresh0 = 0.13 * sqrt(soil_density*grav*diameter/air_dens) &
                       * sqrt(1.+6.e-7/(soil_density*grav*diameter**2.5)) &
              / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)

!      emissions_(:,:) = 0.

!     Spatially dependent part of calculation
!     ---------------------------------------
      do j = j1, j2
         do i = i1, i2
            if ( oro(i,j) /= LAND ) cycle ! only over LAND gridpoints

            w10m = sqrt(u10m(i,j)**2.+v10m(i,j)**2.)
!           Modify the threshold depending on soil moisture as in Ginoux et al. [2001]
            if(gwettop(i,j) .lt. 0.5) then
               u_thresh = amax1(0.,u_thresh0* &
               (1.2+0.2*alog10(max(1.e-3,gwettop(i,j)))))

               if(w10m .gt. u_thresh) then
!                 Emission of dust [kg m-2 s-1]
                  emissions(i,j,n) = (1.-fraclake(i,j)) * w10m**2. * (w10m-u_thresh)
               endif
            endif !(gwettop(i,j) .lt. 0.5)
         end do ! i
      end do ! j
      emissions(:,:,n) = Ch_DU * du_src * emissions(:,:,n)
    end do ! n

   rc = __SUCCESS__

   end subroutine DustEmissionGOCART2G
