   subroutine DustEmissionK14( km, t_soil, w_top, rho_air,    &
                               z0, z, u_z, v_z, ustar,    &
                               f_land, f_snow,            &
                               f_src,                     &
                               f_sand, f_silt, f_clay,    &
                               texture, vegetation, gvf,  &
                               f_w, f_c, uts_gamma,       &
                               UNDEF, GRAV, VON_KARMAN,   &
                               opt_clay, Ch_DU,           &
                               emissions,                 &
                               u, u_t, u_ts,              &
                               R, H_w, f_erod,            &
                               rc )

! !USES:
   implicit none

! !INPUT PARAMETERS:
   integer, intent(in)              :: km         ! model levels
   real, dimension(:,:), intent(in) :: rho_air    ! air density
   real, dimension(:,:), intent(in) :: w_top      ! volumetric soil moisture in the top surface layer
   real, dimension(:,:), intent(in) :: t_soil     ! soil temperature
   real, dimension(:,:), intent(in) :: z0         ! aeolian aerodynamic roughness length
   real, dimension(:,:), intent(in) :: z, u_z, v_z! hight and wind at this height
   real, dimension(:,:), intent(in) :: ustar      ! friction velocity
   real, dimension(:,:), intent(in) :: f_land     ! land fraction
   real, dimension(:,:), intent(in) :: f_snow     ! snow fraction
   real, dimension(:,:), intent(in) :: f_src      ! dust source potential -- OBSOLETE
   real, dimension(:,:), intent(in) :: f_sand     ! sand fraction
   real, dimension(:,:), intent(in) :: f_silt     ! silt fraction
   real, dimension(:,:), intent(in) :: f_clay     ! clay fraction
   real, dimension(:,:), intent(in) :: texture    ! soil texture
   real, dimension(:,:), intent(in) :: vegetation ! vegetation categories (IGBP)
   real, dimension(:,:), intent(in) :: gvf        ! vegetation fraction

   integer, intent(in)              :: opt_clay   ! controls which clay&silt emissions term to use
   real, intent(in)                 :: Ch_DU      ! dust emission tuning coefficient [kg/(sec^2 m^5)]
   real,    intent(in)              :: f_w        ! factor to scale down soil moisture in the top 5cm to soil moisture in the top 1cm
   real,    intent(in)              :: f_c        ! scale down the wet sieving clay fraction to get it more in line with dry sieving measurements
   real,    intent(in)              :: uts_gamma  ! threshold friction velocity parameter 'gamma'
   real,    intent(in)              :: UNDEF      ! paramter for undefined varaible
   real,    intent(in)              :: GRAV       ! gravity
   real,    intent(in)              :: VON_KARMAN ! von karman constant

! !OUTPUT PARAMETERS:

   real, dimension(:,:,:), intent(out) :: emissions ! mass flux of emitted dust particles

   real, dimension(:,:), intent(out) :: u         ! aeolian friction velocity
   real, dimension(:,:), intent(out) :: u_t       ! threshold friction velocity
   real, dimension(:,:), intent(out) :: u_ts      ! threshold friction velocity over smooth surface

   real, dimension(:,:), intent(out) :: H_w       ! soil mosture correction
   real, dimension(:,:), intent(out) :: R         ! drag partition correction

   real, dimension(:,:), intent(out) :: f_erod    ! erodibility


   integer, intent(out) :: rc    ! Error return code:
                                 !  0 - all is well
                                 !  1 -

   character(len=*), parameter :: myname = 'DustEmissionK14'

!  !Local Variables

   real, dimension(:,:), allocatable :: w_g        ! gravimetric soil moisture
   real, dimension(:,:), allocatable :: w_gt       ! threshold gravimetric soil moisture

   real, dimension(:,:), allocatable :: f_veg      ! vegetation mask
   real, dimension(:,:), allocatable :: clay       ! 'corrected' clay fraction in '%'
   real, dimension(:,:), allocatable :: silt       ! 'corrected' silt fraction in '%'
   real, dimension(:,:), allocatable :: k_gamma    ! silt and clay term (gamma in K14 and I&K, 2017)
   real, dimension(:,:), allocatable :: z0s        ! smooth roughness length

   real, dimension(:,:), allocatable :: Dp_size    ! typical size of soil particles for optimal saltation
   real :: rho_p                              ! typical density of soil particles

   integer :: i, j, i1=1, i2, j1=1, j2, n

   real, parameter :: z0_valid = 0.08e-2      ! valid range of ARLEMS z0 is 0--0.08cm, z0 > 0.08cm is retreived but the data quality is low
   real, parameter :: z0_max = 6.25 * z0_valid! maximum roughness over arid surfaces
   real, parameter :: z0_    = 2.0e-4         ! representative aeolian aerodynamic roughness length z0 = 0.02cm

   real, parameter :: rho_water     = 1000.0  ! water density, 'kg m-3'
   real, parameter :: rho_soil_bulk = 1700.0  ! soil bulk density, 'kg m-3'
   real, parameter :: rho_soil      = 2500.0  ! soil particle density, 'kg m-3'

!  real, parameter :: f_w = 0.5               ! factor to scale down soil moisture in the top 5cm to soil moisture in the top 1cm
!  real, parameter :: f_c = 0.7               ! scale down the wet sieving clay fraction to get it more in line with dry sieving measurements

   ! Shao et al.
   real, parameter :: a_n = 0.0123
   real, parameter :: G   = 1.65e-4

   ! size of coarsest mode in the STATSGO/FAO soil type
   real, parameter :: Dc_soil(12) = (/ 710e-6, 710e-6, 125e-6, &
                                       125e-6, 125e-6, 160e-6, &
                                       710e-6, 125e-6, 125e-6, &
                                       160e-6, 125e-6,   2e-6 /)

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
!  15Aug2016, Darmenov - Initial implementation
!  15Dec2020, E.Sherman - Ported to GOCART2G process library

!EOP
!-------------------------------------------------------------------------
!  Begin...
   rc = 0
   i2 = ubound(t_soil,1)
   j2 = ubound(t_soil,2)

   allocate(w_g(i2,j2), w_gt(i2,j2), f_veg(i2,j2), clay(i2,j2), silt(i2,j2), k_gamma(i2,j2))
   allocate(z0s(i2,j2), Dp_size(i2,j2))

   ! typical size of soil particles for optimal saltation is about 75e-6m
   Dp_size = 75e-6

   ! typical density of soil particles, e.g. quartz grains
   rho_p = 2.65e3

   ! threshold friction velocity over smooth surface
   u_ts = UNDEF
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0))
       u_ts = sqrt(a_n * ( ((rho_p/rho_air) * GRAV * Dp_size) + uts_gamma / (rho_air * Dp_size)))
   end where

#if (0)
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  1) < 0.5) u_ts = u_ts * 1.176
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  2) < 0.5) u_ts = u_ts * 1.206
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  3) < 0.5) u_ts = u_ts * 1.234
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  4) < 0.5) u_ts = u_ts * 1.261
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  5) < 0.5) u_ts = u_ts * 1.272
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  6) < 0.5) u_ts = u_ts * 1.216
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  7) < 0.5) u_ts = u_ts * 1.211
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  8) < 0.5) u_ts = u_ts * 1.266
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture -  9) < 0.5) u_ts = u_ts * 1.222
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture - 10) < 0.5) u_ts = u_ts * 1.146
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture - 11) < 0.5) u_ts = u_ts * 1.271
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0) .and. abs(texture - 12) < 0.5) u_ts = u_ts * 1.216
#endif

   ! gravimetric soil moisture : scaled down to represent the values in the top 1cm and converted to '%'
   w_g = UNDEF
   where (f_land > 0.0)
#if (1)
       ! following Zender
       ! Q_s_ = 0.489 - 0.126*f_sand
       ! rho_soil_bulk = rho_soil*(1 - Q_s_)
       ! w_g = 100 * f_w * (rho_water / rho_soil_bulk) * w_top
       ! ...the equivalent one-liner
       w_g = 100 * f_w * rho_water / rho_soil / (1.0 - (0.489 - 0.126*f_sand)) * w_top
#else
       w_g = 100 * f_w * (rho_water / rho_soil_bulk) * w_top
#endif
   end where

   ! soil moisture correction following Fecan
   clay = UNDEF
   silt = UNDEF
   w_gt = UNDEF
   where ((f_land > 0.0) .and. (f_clay <= 1.0) .and. (f_clay >= 0.0))
       clay = f_c * f_clay
       silt = f_silt + (1.0-f_c)*f_clay   ! move the excess clay to the silt fraction

       w_gt = 14.0*clay*clay + 17.0*clay  ! w_gt in '%'
   end where

   H_w  = 1.0
#if (1)
   ! Fecan, 1999
   where ((f_land > 0.0) .and. (w_g > w_gt))
       H_w = sqrt(1.0 + 1.21*(w_g - w_gt)**0.68)
   end where
#else
   ! Shao, 1996
   where ((f_land > 0.0) .and. (w_top <= 1.0) .and. (w_top >= 0.0))
       H_w = exp(22.7*f_w *w_top)
   end where
#endif


   select case (opt_clay)
   case (1)
       ! following Ito and Kok, 2017
       k_gamma = 0.05

       where ((f_land > 0.0) .and. (clay < 0.2) .and. (clay >= 0.05))
           k_gamma = clay
       end where

       where ((f_land > 0.0) .and. (clay >= 0.2) .and. (clay <= 1.0))
           k_gamma = 0.2
       end where
   case (2)
       ! following Ito and Kok, 2017
       k_gamma = 1.0/1.4

       where ((f_land > 0.0) .and. (clay < 0.2) .and. (clay >= 0.0))
           k_gamma = 1.0 / (1.4 - clay - silt)
       end where

       where ((f_land > 0.0) .and. (clay >= 0.2) .and. (clay <= 1.0))
           k_gamma = 1.0 / (1.0 + clay - silt)
       end where
   case default
       ! following Kok et al, 2014
       k_gamma = 0.0

       where ((f_land > 0.0) .and. (clay <= 1.0) .and. (clay >= 0.0))
           k_gamma = clay
       end where
   end select


   ! roughness over smooth surface
   z0s = 125e-6
   do j = j1, j2
       do i = i1, i2
           if (texture(i,j) > 0 .and. texture(i,j) < 13) then
               z0s(i,j) = Dc_soil(nint(texture(i,j)))
           end if
       end do
   end do

   z0s = z0s / 30.0    ! z0s = MMD/x, where typically x is in the range 24--30; x=10 was recommended
                       ! as a more appropriate value for this parameter in a recent article

   ! drag partition correction
   R = 1.0
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > z0s))
#if (1)
       ! MacKinnon et al, 2004
       R = 1.0 - log(z0/z0s)/log(0.7 * (122.55/z0s)**0.8)
#else
       ! King et al, 2005, Darmenova et al, 2009, and K14-S1 use the corrected MB expression
       R = 1.0 - log(z0/z0s)/log(0.7 * (0.1/z0s)**0.8)
#endif
   end where


   ! *soil* friction velocity, see Equations 5, S.10, S11 in Kok et al, 2014 and the supplement paper
   u = UNDEF
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0))
#if (1)
       u = ustar
#else
       u = VON_KARMAN / log(z/z0) * sqrt(u_z*u_z + v_z*v_z)
#endif
       u = R * u           ! correction for roughness elements
   end where


   ! *soil* threshold friction velocity, Section 2.2 in Kok et al, 2014
   u_t = UNDEF
   where ((f_land > 0.0) .and. (z0 < z0_max) .and. (z0 > 0.0))
       u_t = u_ts * H_w    ! apply moisture correction
   end where


   ! erodibility
   f_erod = UNDEF
   where (f_land > 0.0)
       f_erod = 1.0
   end where

   ! erodibility parameterization - Laurent et al., 2008
   where ((f_land > 0.0) .and. (z0 > 3.0e-5) .and. (z0 < z0_max))
       f_erod = 0.7304 - 0.0804*log10(100*z0)
   end where

   ! bedrock
   where (abs(texture - 15) < 0.5) f_erod = 0.0


   ! vegetation mask
   f_veg = 0.0
   where ((f_land > 0.0) .and. abs(vegetation -  7) < 0.1) f_veg = 1.0  ! open shrublands
!  where ((f_land > 0.0) .and. abs(vegetation -  9) < 0.1) f_veg = 1.0  ! savannas
!  where ((f_land > 0.0) .and. abs(vegetation - 10) < 0.1) f_veg = 1.0  ! grasslands
!  where ((f_land > 0.0) .and. abs(vegetation - 12) < 0.1) f_veg = 1.0  ! croplands
   where ((f_land > 0.0) .and. abs(vegetation - 16) < 0.1) f_veg = 1.0  ! barren or sparsely vegetated

   ! vegetation mask: modulate with vegetation fraction
   where (f_land > 0.0 .and. gvf >= 0.0 .and. gvf < 0.8) f_veg = f_veg * (1 - gvf)


   ! final erodibility
   f_erod = f_erod * f_veg * f_land * (1.0 - f_snow)

   ! ...kludge to deal with high emissions in Australia
   where (f_src >= 0.0) f_erod = f_src * f_erod

   call VerticalDustFluxK14( i1, i2, j1, j2, km, &
                             u, u_t, rho_air,    &
                             f_erod, k_gamma,    &
                             emissions(:,:,1) )

   ! duplicate dust emissions across the 3rd dimension for use in call to UpdateAerosolState
   ! UpdateAerosolState expects surface dust emissions array of 3 dimensions(x, y, bin).
   emissions(:,:,1) = emissions(:,:,1) * Ch_DU
   do n = 2, size(emissions, dim=3)
      emissions(:,:,n) = emissions(:,:,1)
   end do

   end subroutine DustEmissionK14
