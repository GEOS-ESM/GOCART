#define __SUCCESS__ 0
#define __FAIL__ 1
#define __VERIFY__(x) if(x/=0) then; if(present(rc)) rc=x; return; endif
#define __VERIFY_NO_OPT__(x) if(x/=0) then; rc=x; return; endif
#define __RC__ rc=status); __VERIFY__(status
#define __RC_NO_OPT__ rc=status); __VERIFY_NO_OPT__(status
#define __STAT__ stat=status); __VERIFY__(status
#define __IOSTAT__ iostat=status); __VERIFY__(status
#define __RETURN__(x) if (present(rc)) rc=x; return
#define __ASSERT__(expr) if(.not. (expr)) then; if (present(rc)) rc=-1; return; endif
!-------------------------------------------------------------------------
!
! !MODULE: GOCART2G_Process -- GOCART2G process library
!
! !INTERFACE:
   module  GOCART2G_Process

! !USES:
!  Only instrinsic fortran types and functions are allowed.
   use GOCART2G_MieMod
   use, intrinsic :: iso_fortran_env, only: IOSTAT_END

   implicit none
   private

!
! !PUBLIC MEMBER FUNCTIONS:
!
   public DustAerosolDistributionKok
   public DustEmissionFENGSHA
   public DustEmissionGOCART2G
   public DustEmissionK14
   public DustFluxV2HRatioMB95
   public moistureCorrectionFecan
   public soilMoistureConvertVol2Grav
   public DistributePointEmission
   public updatePointwiseEmissions
   public Chem_Settling
   public Chem_SettlingSimple
   public Chem_Settling2Gorig
   public Chem_SettlingSimpleOrig
   public DryDeposition
   public WetRemovalGOCART2G
   public UpdateAerosolState
   public Aero_Compute_Diags
   public jeagleSSTcorrection
   public deepLakesMask
   public weibullDistribution
   public SeasaltEmission
   public wetRadius
   public hoppelCorrection
   public CAEmission
   public phobicTophilic
   public NIheterogenousChem
   public SulfateDistributeEmissions
   public DMSemission
   public SUvolcanicEmissions
   public SulfateUpdateOxidants
   public SU_Wet_Removal
   public SU_Compute_Diags
   public SulfateChemDriver
   public get_HenrysLawCts
   public NIthermo
   public Chem_UtilResVal
   public Chem_UtilIdow
   public Chem_UtilCdow
   public Chem_BiomassDiurnal
   public ReadPointEmissions
   public EmissionReader


   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   integer, parameter     :: DP = kind(1.0d0)

   type :: EmissionReader
      private
      integer, allocatable :: unit
   contains
      procedure :: open
      procedure :: close
      procedure :: rewind => rewind_reader
      procedure :: is_end_marker
      procedure :: read_table
      procedure :: next_line
      procedure :: count_words
      procedure :: scan_to_label
      procedure :: get_dims
   end type EmissionReader

   type KeywordEnforcer
   end type KeywordEnforcer

!
! !DESCRIPTION:
!
!  This module contains and implements all necessary process calculations for GOCART.
!
! !REVISION HISTORY:
!
!  11Feb2020  E.Sherman, A.da Silva, T.Clune, A.Darmenov - Ported/consolidated/refactored GOCART
!                   physics and chemistry code into a single process library that only uses
!                   intrinsic Fortran functions.
!
!  01Apr2021  R.Montuoro/NOAA - Added FENGSHA dust scheme and related methods.
!
!
!EOP
!-------------------------------------------------------------------------
CONTAINS

!=====================================================================================
!BOP
!
! !IROUTINE:  DustAerosolDistributionKok - Compute Kok's dust size aerosol distribution
!
! !INTERFACE:
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

!===============================================================================
!BOP
!
! !IROUTINE: soilMoistureConvertVol2Grav - volumetric to gravimetric soil moisture
!
! !INTERFACE:
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

!===============================================================================
!BOP
!
! !IROUTINE: moistureCorrectionFecan - Correction factor for Fecan soil moisture
!
! !INTERFACE:
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

!===============================================================================
!BOP
!
! !IROUTINE: DustFluxV2HRatioMB95 - vertical-to-horizontal dust flux ratio (MB95)
!
! !INTERFACE:
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

!==================================================================================
!BOP
!
! !IROUTINE: DustEmissionFENGSHA - Compute dust emissions using NOAA/ARL FENGSHA model
!
! !INTERFACE:
   subroutine DustEmissionFENGSHA(fraclake, fracsnow, oro, slc, clay, sand, silt,  &
                                  ssm, rdrag, airdens, ustar, uthrs, alpha, gamma, &
                                  kvhmax, grav, rhop, distribution, emissions, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, dimension(:,:), intent(in) :: fraclake ! fraction of lake [1]
   real, dimension(:,:), intent(in) :: fracsnow ! surface snow area fraction [1]
   real, dimension(:,:), intent(in) :: slc      ! liquid water content of soil layer, volumetric fraction [1]
   real, dimension(:,:), intent(in) :: oro      ! land-ocean-ice mask [1]
   real, dimension(:,:), intent(in) :: clay     ! fractional clay content [1]
   real, dimension(:,:), intent(in) :: sand     ! fractional sand content [1]
   real, dimension(:,:), intent(in) :: silt     ! fractional silt content [1]
   real, dimension(:,:), intent(in) :: ssm      ! erosion map [1]
   real, dimension(:,:), intent(in) :: rdrag    ! drag partition [1/m]
   real, dimension(:,:), intent(in) :: airdens  ! air density at lowest level [kg/m^3]
   real, dimension(:,:), intent(in) :: ustar    ! friction velocity [m/sec]
   real, dimension(:,:), intent(in) :: uthrs    ! threshold velocity [m/2]
   real,                 intent(in) :: alpha    ! scaling factor [1]
   real,                 intent(in) :: gamma    ! scaling factor [1]
   real,                 intent(in) :: kvhmax   ! max. vertical to horizontal mass flux ratio [1]
   real,                 intent(in) :: grav     ! gravity [m/sec^2]
   real, dimension(:),   intent(in) :: rhop            ! soil class density [kg/m^3]
   real, dimension(:),   intent(in) :: distribution    ! normalized dust binned distribution [1]

! !OUTPUT PARAMETERS:
   real,    intent(out) :: emissions(:,:,:)     ! binned surface emissions [kg/(m^2 sec)]
   integer, intent(out) :: rc                   ! Error return code: __SUCCESS__ or __FAIL__

! !DESCRIPTION: Compute dust emissions using NOAA/ARL FENGSHA model
!
! !REVISION HISTORY:
!
! 22Feb2020 B.Baker/NOAA    - Original implementation
! 29Mar2021 R.Montuoro/NOAA - Refactored for process library
!

! !Local Variables
   logical               :: skip
   integer               :: i, j, n, nbins
   integer, dimension(2) :: ilb, iub
   real                  :: alpha_grav
   real                  :: fracland
   real                  :: h
   real                  :: kvh
   real                  :: q
   real                  :: rustar
   real                  :: total_emissions
   real                  :: u_sum, u_thresh

! !CONSTANTS:
   real, parameter       :: ssm_thresh = 1.e-02    ! emit above this erodibility threshold [1]

!EOP
!-------------------------------------------------------------------------
!  Begin

   rc = __SUCCESS__

!  Get dimensions and index bounds
!  -------------------------------
   nbins = size(emissions, dim=3)
   ilb = lbound(ustar)
   iub = ubound(ustar)

!  Initialize emissions
!  --------------------
   emissions = 0.

!  Prepare scaling factor
!  ----------------------
   alpha_grav = alpha / grav

!  Compute size-independent factors for emission flux
!  ---------------------------
   do j = ilb(2), iub(2)
     do i = ilb(1), iub(1)
       ! skip if we are not on land
       ! --------------------------
       skip = (oro(i,j) /= LAND)
       ! threshold and sanity check for surface input
       ! --------------------------------------------
       if (.not.skip) skip = (ssm(i,j) < ssm_thresh) &
         .or. (clay(i,j) < 0.) .or. (sand(i,j) < 0.) &
         .or. (rdrag(i,j) < 0.)

       if (.not.skip) then
         fracland = max(0., min(1., 1.-fraclake(i,j))) &
                  * max(0., min(1., 1.-fracsnow(i,j)))

         ! Compute vertical-to-horizontal mass flux ratio
         ! ----------------------------------------------
         kvh = DustFluxV2HRatioMB95(clay(i,j), kvhmax)

         ! Compute total emissions
         ! -----------------------
         total_emissions = alpha_grav * fracland * (ssm(i,j) ** gamma) &
                         * airdens(i,j) * kvh

         !  Compute threshold wind friction velocity using drag partition
         !  -------------------------------------------------------------
         rustar = rdrag(i,j) * ustar(i,j)

         !  Now compute size-dependent total emission flux
         !  ----------------------------------------------
         do n = 1, nbins
           ! Fecan moisture correction
           ! -------------------------
           h = moistureCorrectionFecan(slc(i,j), sand(i,j), clay(i,j), rhop(n))

           ! Adjust threshold
           ! ----------------
           u_thresh = uthrs(i,j) * h

           u_sum = rustar + u_thresh

           ! Compute Horizontal Saltation Flux according to Eq (9) in Webb et al. (2020)
           ! ---------------------------------------------------------------------------
           q = max(0., rustar - u_thresh) * u_sum * u_sum

           ! Distribute emissions to bins and convert to mass flux (kg s-1)
           ! --------------------------------------------------------------
           emissions(i,j,n) = distribution(n) * total_emissions * q
         end do

       end if

     end do
   end do

   end subroutine DustEmissionFENGSHA

!==================================================================================
!BOP
! !IROUTINE: DustEmissionGOCART2G

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

!==================================================================================
!BOP
! !IROUTINE: DustEmissionK14

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

   allocate(w_g(i2,j2), w_gt(i2,j2), f_veg(i2,j2), clay(i2,j2), silt(i2,j2), k_gamma(i2,j2), source=0.0)
   allocate(z0s(i2,j2), Dp_size(i2,j2), source=0.0)

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

!==================================================================================
!BOP
! !IROUTINE: VerticalDustFluxK14

   subroutine VerticalDustFluxK14( i1, i2, j1, j2, km, &
                                   u, u_t, rho_air, &
                                   f_erod, k_gamma,     &
                                   emissions )

! !USES:
   implicit none
! !INPUT PARAMETERS:
   integer, intent(in) ::  i1, i2, j1, j2, km

   real, dimension(:,:), intent(in) :: u           ! friction velocity, 'm s-1'
   real, dimension(:,:), intent(in) :: u_t         ! threshold friction velocity, 'm s-1'
   real, dimension(:,:), intent(in) :: rho_air     ! air density, 'kg m-3'
   real, dimension(:,:), intent(in) :: f_erod      ! erodibility
   real, dimension(:,:), intent(in) :: k_gamma     ! clay and silt dependent term that modulates the emissions

! !OUTPUT PARAMETERS:
   real, intent(out)    :: emissions(:,:)          ! total vertical dust mass flux, 'kg m-2 s-1'

   character(len=*), parameter :: myname = 'VerticalDustFluxK14'

! !Local Variables
   integer :: i, j
   real    :: u_st                     ! standardized threshold friction velocity
   real    :: C_d                      ! dust emission coefficient
   real    :: f_ust                    ! numerical term

  ! parameters from Kok et al. (2012, 2014)
   real, parameter :: rho_a0 = 1.225   ! standard atmospheric density at sea level, 'kg m-3'
   real, parameter :: u_st0  = 0.16    ! the minimal value of u* for an optimally erodible soil, 'm s-1'
   real, parameter :: C_d0   = 4.4e-5  ! C_d0 = (4.4 +/- 0.5)*1e-5
   real, parameter :: C_e    = 2.0     ! C_e  = 2.0 +/- 0.3
   real, parameter :: C_a    = 2.7     ! C_a  = 2.7 +/- 1.0

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
!  11Oct2011, Darmenov - For now use the GOCART emission scheme to
!                        calculate the total emission
!
!EOP
!-------------------------------------------------------------------------

  emissions = 0.0 ! total emission

  !  Vertical dust flux
  !  ------------------
  do j = j1, j2
      do i = i1, i2

          if ((f_erod(i,j) > 0.0) .and. (u(i,j) > u_t(i,j))) then
              u_st  = u_t(i,j) * sqrt(rho_air(i,j) / rho_a0)
              u_st  = max(u_st, u_st0)

              f_ust = (u_st - u_st0)/u_st0
              C_d = C_d0 * exp(-C_e * f_ust)

              emissions(i,j) = C_d * f_erod(i,j) * k_gamma(i,j) * rho_air(i,j)  &
                                   * ((u(i,j)*u(i,j) - u_t(i,j)*u_t(i,j)) / u_st) &
                                   * (u(i,j) / u_t(i,j))**(C_a * f_ust)
          end if

      end do
  end do

  ! all done
  end subroutine VerticalDustFluxK14

!==================================================================================
!BOP

! !IROUTINE: updatePointwiseEmissions

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

    integer :: status

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
                                     point_column_emissions, __RC__)

        emissions_point(i,j,:) =  point_column_emissions
        end do

      __RETURN__(__SUCCESS__)
  end subroutine updatePointwiseEmissions

!==================================================================================
!BOP
! !IROUTINE: DistributePointEmissions

! !INTERFACE:

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


!==================================================================================
!BOP
! !IROUTINE: Chem_SettlingSimple

   subroutine Chem_SettlingSimple ( km, klid, flag, cdt, grav, &
                                    radiusInp, rhopInp, int_qa, tmpu, &
                                    rhoa, rh, hghte, delp, fluxout,  &
                                    vsettleOut, correctionMaring, rc)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)    :: km     ! total model levels
   integer, intent(in)    :: klid   ! index for pressure lid
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt
   real, intent(in)    :: grav   ! gravity [m/sec^2]
   real, intent(in)  :: radiusInp  ! particle radius [microns]
   real, intent(in)  :: rhopInp    ! soil class density [kg/m^3]
   real, dimension(:,:,:), intent(inout) :: int_qa  ! aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu   ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa   ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in)  :: hghte  ! geopotential height [m]
   real, pointer, dimension(:,:,:), intent(in)  :: delp   ! pressure level thickness [Pa]

! !OUTPUT PARAMETERS:

   real, pointer, dimension(:,:), intent(inout)  :: fluxout ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, optional, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -
!  Optionally output the settling velocity calculated
   real, pointer, optional, dimension(:,:,:)  :: vsettleOut

!  Optionally correct the settling velocity following Maring et al, 2003
   logical, optional, intent(in)    :: correctionMaring

   character(len=*), parameter :: myname = 'SettlingSimple'

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop)
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, n, dk

   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

   integer         :: i1=1, i2, j1=1, j2
   integer         :: dims(3)

   real, pointer, dimension(:,:)   :: hsurf
   real(kind=DP), dimension(:,:), allocatable  :: cmass_before, cmass_after
   real, allocatable    :: dz(:,:,:)
   real, dimension(:,:,:), allocatable  :: radius, rhop, qa
   real, dimension(:,:,:), allocatable  :: vsettle   ! fall speed [m s-1]
   real ::  ONE_OVER_G
   integer :: status

!EOP
!-------------------------------------------------------------------------

!  Get dimensions
!  ---------------
   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

   ONE_OVER_G = 1.0/grav

   hsurf => hghte(i1:i2,j1:j2,km)

   allocate(dz(i2,j2,km), radius(i2,j2,km), rhop(i2,j2,km), vsettle(i2,j2,km), qa(i2,j2,km), source=0.0)
   allocate(cmass_before(i2,j2), cmass_after(i2,j2), source=0.0_DP)

   qa = int_qa

   if(associated(fluxout)) fluxout(:,:) = 0.0

!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1

!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo

!  If radius le 0 then get out
   if(radiusInp .le. 0.) then
      status = 100
      __RETURN__(STATUS)
   end if

!   Find the column dry mass before sedimentation
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k) * delp(i,j,k) * ONE_OVER_G
          enddo
       enddo
    enddo

!   Particle swelling
    call ParticleSwelling(i1, i2, j1, j2, km, rh, radiusInp, rhopInp, radius, rhop, flag)

!   Settling velocity of the wet particle
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             call Chem_CalcVsettle(radius(i,j,k), rhop(i,j,k), rhoa(i,j,k), &
                                   tmpu(i,j,k), vsettle(i,j,k), grav)
          end do
       end do
    end do

    if(present(correctionMaring)) then
       if (correctionMaring) then
          vsettle = max(1.0e-9, vsettle - v_upwardMaring)
       endif
    endif

    if(present(vsettleOut)) then
       vsettleOut = vsettle
    endif

!   Time integration
    call SettlingSolver(i1, i2, j1, j2, km, cdt, delp, dz, vsettle, qa)

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k) * delp(i,j,k) * ONE_OVER_G
          enddo
       enddo
    enddo

    if( associated(fluxout) ) then
       fluxout(:,:) = (cmass_before - cmass_after)/cdt
    endif

    int_qa = qa

   __RETURN__(__SUCCESS__)

   end subroutine Chem_SettlingSimple

!==================================================================================
!BOP
! !IROUTINE: Chem_Settling

   subroutine Chem_Settling ( km, klid, bin, flag, cdt, grav, &
                              radiusInp, rhopInp, int_qa, tmpu, &
                              rhoa, rh, hghte, delp, fluxout,  &
                              vsettleOut, correctionMaring, rc)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)    :: km     ! total model levels
   integer, intent(in)    :: klid   ! index for pressure lid
   integer, intent(in)    :: bin    ! aerosol bin index
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt
   real, intent(in)    :: grav   ! gravity [m/sec^2]
   real, intent(in)  :: radiusInp  ! particle radius [microns]
   real, intent(in)  :: rhopInp    ! soil class density [kg/m^3]
   real, dimension(:,:,:), intent(inout) :: int_qa  ! aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu   ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa   ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in)  :: hghte  ! geopotential height [m]
   real, pointer, dimension(:,:,:), intent(in)  :: delp   ! pressure level thickness [Pa]

! !OUTPUT PARAMETERS:

   real, pointer, dimension(:,:,:), intent(inout)  :: fluxout ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, optional, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -
!  Optionally output the settling velocity calculated
   real, pointer, optional, dimension(:,:,:)  :: vsettleOut

!  Optionally correct the settling velocity following Maring et al, 2003
   logical, optional, intent(in)    :: correctionMaring

   character(len=*), parameter :: myname = 'Settling'

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop)
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
!  21Apr2021  E.Sherman Ported to process library
!  15May2019  Darmenov  Refactor and speed up code
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, n, dk

   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

   integer         :: i1=1, i2, j1=1, j2
   integer         :: dims(3)

   real, pointer, dimension(:,:)   :: hsurf
   real(kind=DP), dimension(:,:), allocatable  :: cmass_before, cmass_after
   real, allocatable    :: dz(:,:,:)
   real, dimension(:,:,:), allocatable  :: radius, rhop, qa
   real, dimension(:,:,:), allocatable  :: vsettle   ! fall speed [m s-1]
   real ::  ONE_OVER_G
   integer :: status

!EOP
!-------------------------------------------------------------------------

!  Get dimensions
!  ---------------
   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

   ONE_OVER_G = 1.0/grav

   hsurf => hghte(i1:i2,j1:j2,km)

   allocate(dz(i2,j2,km), radius(i2,j2,km), rhop(i2,j2,km), vsettle(i2,j2,km), qa(i2,j2,km), source=0.0)
   allocate(cmass_before(i2,j2), cmass_after(i2,j2), source=0.0_DP)

   qa = int_qa

   if(associated(fluxout)) fluxout(:,:,bin) = 0.0

!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1

!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo

!  If radius le 0 then get out
   if(radiusInp .le. 0.) then
      status = 100
      __RETURN__(STATUS)
   end if

!   Find the column dry mass before sedimentation
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k) * delp(i,j,k) * ONE_OVER_G
          enddo
       enddo
    enddo

!   Particle swelling
    call ParticleSwelling(i1, i2, j1, j2, km, rh, radiusInp, rhopInp, radius, rhop, flag)

!   Settling velocity of the wet particle
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             call Chem_CalcVsettle(radius(i,j,k), rhop(i,j,k), rhoa(i,j,k), &
                                   tmpu(i,j,k), vsettle(i,j,k), grav)
          end do
       end do
    end do

    if(present(correctionMaring)) then
       if (correctionMaring) then
          vsettle = max(1.0e-9, vsettle - v_upwardMaring)
       endif
    endif

    if(present(vsettleOut)) then
       vsettleOut = vsettle
    endif

!   Time integration
    call SettlingSolver(i1, i2, j1, j2, km, cdt, delp, dz, vsettle, qa)

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k) * delp(i,j,k) * ONE_OVER_G
          enddo
       enddo
    enddo

    if( associated(fluxout) ) then
       fluxout(:,:,bin) = (cmass_before - cmass_after)/cdt
    endif

    int_qa = qa

   __RETURN__(__SUCCESS__)

   end subroutine Chem_Settling


!==================================================================================
!BOP
! !IROUTINE: Chem_CalcVsetle

   subroutine Chem_CalcVsettle ( radius, rhop, rhoa, tmpu, &
                                 vsettle, grav )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   real, intent(in)    :: radius              ! Particle radius [m]
   real, intent(in)    :: rhop                ! Particle density [kg m-3]
   real, intent(in)    :: rhoa                ! Layer air density [kg m-3]
   real, intent(in)    :: tmpu                ! Layer temperature [K]
   real, intent(in)    :: grav                ! Gravity [m s-2]

! !OUTPUT PARAMETERS:

   real, intent(out)   :: vsettle                 ! Layer fall speed [m s-1]

   character(len=*), parameter :: myname = 'Vsettle'

! !DESCRIPTION: Calculates the aerosol settling velocity and Brownian diffusion
!               coefficient
!               Follows discussions in Seinfeld and Pandis, Pruppacher and
!               Klett, and the coding in CARMA (Toon et al., 1988)
!               Should work satisfactorily for al reasonable sized aerosols
!               (up to Reynolds number 300)
!
! !REVISION HISTORY:
!
!  06Nov2003  Colarco   Initial version.
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   real :: rmu                       ! Dynamic viscosity [kg m-1 s-1]
   real :: vt                        ! Thermal velocity of air molecule [m s-1]
   real :: rmfp                      ! Air molecule mean free path [m]
   real :: bpm                       ! Cunningham slip correction factor
   real :: rkn                       ! Knudsen number
   real :: re, x, y                  ! reynold's number and parameters
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]
   real, parameter :: pi = 3.141529265

   real, parameter :: f_vt = 8*kb/pi/m_air
   real, parameter :: f_diff_coef = kb/(6*pi)
   real, parameter :: two_over_nine = 2./9.

   real, parameter :: a0 = -3.18657
   real, parameter :: a1 =  0.992696
   real, parameter :: a2 = -1.53193e-3
   real, parameter :: a3 = -9.870593e-4
   real, parameter :: a4 = -5.78878e-4
   real, parameter :: a5 =  8.55176e-5
   real, parameter :: a6 = -3.27815e-6


!  Dynamic viscosity from corrected Sutherland's Equation
   rmu = 1.8325e-5*(416.16/(tmpu+120.))*(tmpu/296.16)**1.5

!  Thermal velocity of air molecule
   vt = sqrt(tmpu * f_vt)

!  Air molecule mean free path
   rmfp = 2*rmu/(rhoa*vt)

!  Knudsen number
   rkn = rmfp/radius

!  Cunningham slip correction factor
!  bpm = 1. + 1.246*rkn + 0.42*rkn*exp(-0.87/rkn)
!  linearized form, Binkowski and Shankar (equation A27, 1995)
   bpm = 1 + 1.246*rkn

!  Brownian diffusion coefficient
!  diff_coef = tmpu*bpm/(rmu*radius) * f_diff_coef

!  Fall speed (assumes Reynolds # < 0.01)
   vsettle = two_over_nine*rhop*radius*radius*grav*bpm/rmu

!  Check the Reynold's number to see if we need a drag correction
!  First guess at Reynold's number using Stoke's calculation
   re = 2.*rhoa*radius*vsettle/rmu

!  If Re > 0.01 then apply drag correction following Pruppacher and
!  Klett regime 2 (eq. 10-142).  Assuming reasonable aerosols we
!  do not consider that particle Re may exceed 300.
   if(re .gt. 0.01) then
      x  = log(24.*re/bpm)
      y  = a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*(a5 + a6*x)))))
      re = exp(y)*bpm
      vsettle = 0.5*rmu*re/(rhoa*radius)
   endif

   end subroutine Chem_CalcVsettle

!==================================================================================
!BOP
! !IROUTINE: SettlingSolver

  subroutine SettlingSolver(i1, i2, j1, j2, km, cdt, delp, dz, vs, qa)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real,    intent(in) :: cdt

    real, dimension(i1:i2,j1:j2,km), intent(in) :: delp
    real, dimension(i1:i2,j1:j2,km), intent(in) :: dz
    real, dimension(i1:i2,j1:j2,km), intent(in) :: vs

    real, dimension(i1:i2,j1:j2,km), intent(inout) :: qa


    ! local
    integer :: i, j, iit
    integer :: nSubSteps

    real, dimension(i1:i2, j1:j2, km) :: tau

    real, dimension(km) :: dp_
    real, dimension(km) :: tau_

    real :: dt, dt_cfl


    tau = vs/dz

    do j = j1, j2
      do i = i1, i2

          dp_  = delp(i,j,:)
          tau_ = tau(i,j,:)

          dt_cfl  = 1 / maxval(tau_)

          if (dt_cfl > cdt) then
              ! no need for time sub-splitting
              nSubSteps = 1
              dt = cdt
          else
              nSubSteps = ceiling(cdt / dt_cfl)
              dt = cdt/nSubSteps
          end if

          do iit = 1, nSubSteps
              qa(i,j,   1) = qa(i,j,   1) * (1 - dt*tau_(1))
              qa(i,j,2:km) = qa(i,j,2:km) + ( (dp_(1:km-1)/dp_(2:km))*(dt*tau_(1:km-1))*qa(i,j,1:km-1) ) &
                                          - dt*tau_(2:km)*qa(i,j,2:km)
          end do

      enddo
    enddo

   end subroutine SettlingSolver

!==================================================================================
!BOP
! !IROUTINE: ParticleSwelling

   subroutine ParticleSwelling (i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop, flag)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    integer, intent(in) :: flag

    real, intent(in) :: radius_dry
    real, intent(in) :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    select case (flag)
        case (0)
            radius = radius_dry
            rhop   = rhop_dry

        case (1)
            call ParticleSwelling_Fitzgerald(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

        case (2)
            call ParticleSwelling_Gerber(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

        case (3)
            call ParticleSwelling_Gerber_AmmoniumSulfate(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

        case (4)
            call ParticleSwelling_PK2007(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

        case default
            radius = radius_dry
            rhop   = rhop_dry
    end select

   end subroutine ParticleSwelling

   subroutine ParticleSwelling_Fitzgerald(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in) :: radius_dry
    real, intent(in) :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    ! the following parameters relate to the swelling of seasalt like particles
    ! following Fitzgerald, Journal of Applied Meteorology, 1975.
    real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
    real, parameter :: alphaNaCl = 1.35

    real :: alpha, alpha1, alpharat, beta, theta, f1, f2
    real :: sat, rrat

    integer :: i, j, k

!   Adjust particle size for relative humidity effects,
!   based on Fitzgerald, Journal of Applied Meteorology, 1975

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

       radius(i,j,k) = radius_dry
       rhop(i,j,k) = rhop_dry

       sat = rh(i,j,k)

       if (sat > 0.80) then
!       parameterization blows up for RH > 0.995, so set that as max
!       rh needs to be scaled 0 - 1
        sat = min(0.995, sat)

!       Calculate the alpha and beta parameters for the wet particle
!       relative to amonium sulfate
        beta = exp( (0.00077*sat) / (1.009-sat) )

        if(sat .le. 0.97) then
         theta = 1.058
        else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
        endif

        alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )

! no need of this calculations, because epsilon == 1
!       f1 = 10.2 - 23.7*sat + 14.5*sat*sat
!       f2 = -6.7 + 15.5*sat -  9.2*sat*sat
!       alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
!       alpha = alphaNaCl * (alpha1*alpharat)
! instead, it is faster to do
        alpha = alphaNaCl * alpha1

        radius(i,j,k) = alpha * radius_dry**beta

        rrat = radius_dry/radius(i,j,k)
        rrat = rrat*rrat*rrat

        rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow
       endif

     end do
    end do
   end do

   end subroutine ParticleSwelling_Fitzgerald

   subroutine ParticleSwelling_Gerber(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in)  :: radius_dry
    real, intent(in)  :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    ! parameters from Gerber 1985 (units require radius in cm, see the variable rcm)
    real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424

    real :: sat, rrat, rcm

    integer :: i, j, k

    ! Adjust the particle size for relative humidity effects,
    ! based on Gerber 1985
    do k = 1, km
     do j = j1, j2
      do i = i1, i2

       sat = max(rh(i,j,k), tiny(1.0)) ! to avoid zero FPE
       sat = min(0.995, sat)

       rcm = radius_dry*100. ! radius in 'cm'

       radius(i,j,k) = 0.01 * ( c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                                + rcm*rcm*rcm )**(1./3.)

       rrat = radius_dry/radius(i,j,k)
       rrat = rrat*rrat*rrat

       rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow

      end do
     end do
    end do

   end subroutine ParticleSwelling_Gerber

   subroutine ParticleSwelling_Gerber_AmmoniumSulfate(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in) :: radius_dry
    real, intent(in) :: rhop_dry
    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    ! parameters for ammonium sulfate from Gerber 1985 (units require radius in cm, see the variable rcm)
    real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428

    real :: sat, rrat, rcm

    integer :: i, j, k


    ! Adjust the particle size for relative humidity effects,
    ! based on Gerber parameterization for Ammonium Sulfate

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

       sat = max(rh(i,j,k), tiny(1.0)) ! to avoid zero FPE
       sat = min(0.995, sat)

       rcm = radius_dry*100. ! radius in 'cm'
       radius(i,j,k) = 0.01 * ( SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                                + rcm*rcm*rcm )**(1./3.)

       rrat = radius_dry/radius(i,j,k)
       rrat = rrat*rrat*rrat

       rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow

      end do
     end do
    end do

   end subroutine ParticleSwelling_Gerber_AmmoniumSulfate

   subroutine ParticleSwelling_PK2007(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop)

    implicit none

    integer, intent(in) :: i1, i2
    integer, intent(in) :: j1, j2
    integer, intent(in) :: km

    real, dimension(i1:i2,j1:j2,km), intent(in)  :: rh

    real, intent(in) :: radius_dry
    real, intent(in) :: rhop_dry

    real, dimension(i1:i2,j1:j2,km), intent(out) :: radius  ! radius  of the wet particle
    real, dimension(i1:i2,j1:j2,km), intent(out) :: rhop    ! density of the wet particle


    ! local
    real, parameter :: rhow = 1000.  ! density of water [kg m-3]

    real :: sat, rrat

    integer :: i, j, k


    ! Adjust the particle size for relative humidity effects,
    ! based on Petters and Kreidenweis (ACP2007) parameterization

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

       sat = rh(i,j,k)
       sat = min(0.99, sat)

       radius(i,j,k) = radius_dry * (1+1.19*sat/(1-sat))**(1./3.)

       rrat = radius_dry/radius(i,j,k)
       rrat = rrat*rrat*rrat

       rhop(i,j,k) = rrat*rhop_dry + (1 - rrat)*rhow

     end do
    end do
   end do

   end subroutine ParticleSwelling_PK2007


!==================================================================================
!BOP
! !IROUTINE: Chem_Settling2Gorig

   subroutine Chem_Settling2Gorig (km, klid, flag, bin, int_qa, grav, delp, &
                                   radiusInp, rhopInp, cdt, tmpu, rhoa, &
                                   rh, hghte, fluxout, vsettleOut, correctionMaring, rc)

! !USES:
   implicit none

! !INPUT PARAMETERS:
   integer,    intent(in)    :: km     ! total model levels
   integer,    intent(in)    :: klid   ! index for pressure lid
   integer,    intent(in)    :: flag   ! flag to control particle swelling (see note)
   integer,    intent(in)    :: bin    ! aerosol bin index
   real,       intent(in)    :: grav   ! gravity [m/sec^2]
   real,       intent(in)    :: cdt    ! chemistry model time-step
   real, intent(in)  :: radiusInp  ! particle radius [microns]
   real, intent(in)  :: rhopInp    ! soil class density [kg/m^3]
   real, dimension(:,:,:), intent(inout) :: int_qa  ! aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu   ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa   ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in)  :: hghte  ! geopotential height [m]
   real, pointer, dimension(:,:,:), intent(in)  :: delp   ! pressure level thickness [Pa]
   logical, optional, intent(in)   :: correctionMaring


! !OUTPUT PARAMETERS:
   real, pointer, dimension(:,:,:), intent(inout) :: fluxout ! Mass lost by settling to surface [kg/(m^2 sec)]

!  Optionally output the settling velocity calculated
   real, dimension(:,:,:), optional, intent(out)  :: vsettleOut !Layer fall speed [m/sec]

   integer, optional, intent(out)   :: rc

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop)
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
! 11Feb2020 E.Sherman - First attempt at port/refactor
!

!  !Local
   real,    parameter     :: rhow = 1000.  ! Density of water [kg m-3]
!  parameter from Gerber 1985 (units require radius in cm, see rcm)
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424
!  parameters for ammonium sulfate
   real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

   integer         :: i1=1, i2, j1=1, j2
   integer         :: dims(3)
   integer         :: nSubSteps, ijl, dk
   integer         :: i, j, k, iit, n

   real, allocatable    :: dz(:,:,:)
   real :: radius, rhop   ! particle radius and density passed to
                          ! fall velocity calculation
   real :: minTime, qmin
   real :: sat, rrat
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2
   real :: rcm
   real :: diff_coef                 ! Brownian diffusion coefficient [m2 s-1]
   real(kind=DP)  :: dt_settle, gravDP
   real, pointer, dimension(:,:)   :: hsurf

   real, dimension(:,:,:), allocatable  :: vsettle   ! fall speed [m s-1]
   real(kind=DP), dimension(:,:,:), allocatable   :: dzd, vsd, qa, qa_temp
   real(kind=DP), dimension(:,:), allocatable  :: cmass_before, cmass_after
   real(kind=DP) :: qdel, qsrc, d_p, d_pm1

   integer :: status

!EOP
!-------------------------------------------------------------------------

!  Get dimensions
!  ---------------
   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)
   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   gravDP = grav

   hsurf => hghte(i1:i2,j1:j2,km)

!  Allocate arrays
!  ---------------
   allocate(dz, mold=rhoa)
   allocate(vsettle(i2,j2,km), source=0.0)
   allocate(dzd(i2,j2,km), vsd(i2,j2,km), qa(i2,j2,km), qa_temp(i2,j2,km), source=0.0_DP)
   allocate(cmass_before(i2,j2), cmass_after(i2,j2), source=0.0_DP)

!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1

!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
      dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo
   dzd = dz

   qa(:,:,:) = int_qa(:,:,:)

   radius = radiusInp
   rhop = rhopInp

   if(associated(fluxout)) fluxout(:,:,bin) = 0.0
   cmass_before(:,:) = 0.d0
   cmass_after(:,:) = 0.d0

!  If radius le 0 then get out
   if(radius .le. 0.) then
      status = 100
      __RETURN__(STATUS)
   end if

   do k = klid, km
      do j = j1, j2
         do i = i1, i2
!           Find the column dry mass before sedimentation
            cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k)/gravDP * delp(i,j,k)

!           Adjust the particle size for relative humidity effects
            sat = max(rh(i,j,k),tiny(1.0)) ! to avoid zero FPE

!           Fitzgerald
            select case(flag)
            case(1)
               if (sat >= 0.80) then
!              if(flag .eq. 1 .and. sat .ge. 0.80) then
!              parameterization blows up for RH > 0.995, so set that as max
!              rh needs to be scaled 0 - 1
                  sat = min(0.995,sat)
!                 Calculate the alpha and beta parameters for the wet particle
!                 relative to amonium sulfate
                  beta = exp( (0.00077*sat) / (1.009-sat) )
                  if(sat .le. 0.97) then
                     theta = 1.058
                  else
                     theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
                  endif
                  alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
                  f1 = 10.2 - 23.7*sat + 14.5*sat**2.
                  f2 = -6.7 + 15.5*sat - 9.2*sat**2.
                  alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
                  alpha = alphaNaCl * (alpha1*alpharat)
!                 radius is the radius of the wet particle
                  radius = alpha * radiusInp**beta
                  rrat = (radiusInp/radius)**3.
                  rhop = rrat*rhopInp + (1.-rrat)*rhow
               end if
            case(2)
               sat = min(0.995,sat)
               rcm = radiusInp*100.
               radius = 0.01 * (c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                                + rcm**3.)**(1./3.)
               rrat = (radiusInp/radius)**3.
               rhop = rrat*rhopInp + (1.-rrat)*rhow
            case(3)
!              Gerber parameterization for Ammonium Sulfate
               sat = min(0.995,sat)
               rcm = radiusInp*100.
               radius = 0.01 * (SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                                + rcm**3.)**(1./3.)
               rrat = (radiusInp/radius)**3.
               rhop = rrat*rhopInp + (1.-rrat)*rhow
            case(4)
!              Petters and Kreidenweis (ACP2007) parameterization
               sat = min(0.99,sat)
               radius = (radiusInp**3 * (1+1.19*sat/(1-sat)))**(1./3.)
               rrat = (radiusInp/radius)**3
               rhop = rrat*rhopInp + (1.-rrat)*rhow
            end select

!            Calculate the settling velocity
             call Chem_CalcVsettle2Gorig(radius, rhop, rhoa(i,j,k), tmpu(i,j,k), &
                                     grav, diff_coef, vsettle(i,j,k))
         end do !do i
      end do !do j
   end do !do k

   if(present(correctionMaring)) then
      if ((correctionMaring) .and. (radiusInp .le. (0.5*diameterMaring))) then
         vsettle = max(1.0e-9, vsettle - v_upwardMaring)
      endif
   endif

   vsd = vsettle

   if(present(vsettleOut)) then
      vsettleOut = vsettle
   endif

!   Loop over sub-timestep
    do j = j1, j2
       do i = i1, i2
      !   Determine global min time to cross grid cell
          qmin = minval(dz(i,j,:)/vsettle(i,j,:))
          minTime = min(cdt,qmin)
      !   Now, how many iterations do we need to do?
          if ( minTime < 0 ) then
             nSubSteps = 0
          else if(minTime .ge. cdt) then
             nSubSteps = 1
             dt_settle = cdt
          else
             nSubSteps = cdt/minTime+1
             dt_settle = cdt/nSubSteps
          endif

          do iit = 1, nSubSteps
!          Try a simple forward Euler scheme
           qdel = qa(i,j,klid)*dt_settle*vsd(i,j,klid)/dzd(i,j,klid)
           qa(i,j,klid) = qa(i,j,klid) - qdel

!             do k = 2, km
             do k = klid+1, km
               d_p  = delp(i,j,k)
               d_pm1 = delp(i,j,k-1)
               qsrc = qdel * d_pm1 / d_p
               qdel = qa(i,j,k)*dt_settle*vsd(i,j,k)/dzd(i,j,k)
               qa(i,j,k) = qa(i,j,k) - qdel + qsrc
            end do
         end do  !itt
      end do
    end do

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = klid, km
     do j = j1, j2
      do i = i1, i2
       cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k)/ gravDP * delp(i,j,k)
      enddo
     enddo
    enddo

    if( associated(fluxout) ) then
       fluxout(:,:,bin) = (cmass_before - cmass_after)/cdt
    endif

    int_qa = qa

   __RETURN__(__SUCCESS__)

   end subroutine Chem_Settling2Gorig

!=========================================================================================

!BOP
!
! !IROUTINE:  Chem_CalcVsettle2Gorig - Calculate the aerosol settling velocity
!
! !INTERFACE:
   subroutine Chem_CalcVsettle2Gorig ( radius, rhop, rhoa, tmpu, grav, &
                                      diff_coef, vsettle )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   real, intent(in) :: radius   ! Particle radius [m]
   real, intent(in) :: rhop     ! Particle density [kg/m^3]
   real, intent(in) :: rhoa     ! Layer air density [kg/m^3]
   real, intent(in) :: tmpu     ! Layer temperature [K]
   real, intent(in) :: grav     ! gravity [m/sec^2]

! !OUTPUT PARAMETERS:
   real, intent(out)   :: diff_coef  ! Brownian diffusion
                                     ! coefficient [m2/sec]
   real, intent(out)   :: vsettle    ! Layer fall speed [m s-1]

   character(len=*), parameter :: myname = 'Vsettle'

! !DESCRIPTION: Calculates the aerosol settling velocity and Brownian diffusion
!               coefficient
!               Follows discussions in Seinfeld and Pandis, Pruppacher and
!               Klett, and the coding in CARMA (Toon et al., 1988)
!               Should work satisfactorily for al reasonable sized aerosols
!               (up to Reynolds number 300)
!
! !REVISION HISTORY:
!
!  06Nov2003  Colarco   Initial version.
!  23Jan2003  da Silva  Standardization
!

! !Local Variables
   real(kind=DP) :: rmu                   ! Dynamic viscosity [kg m-1 s-1]
   real(kind=DP) :: vt                    ! Thermal velocity of air molecule [m s-1]
   real(kind=DP) :: rmfp                  ! Air molecule mean free path [m]
   real(kind=DP) :: bpm                   ! Cunningham slip correction factor
   real(kind=DP) :: rkn                   ! Knudsen number
   real(kind=DP) :: re, x, y              ! reynold's number and parameters
   real, parameter :: kb = 1.3807e-23     ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26  ! Mass of <avg> air molecule [kg]
   real, parameter :: pi = 3.141529265

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Dynamic viscosity from corrected Sutherland's Equation
   rmu = 1.8325e-5*(416.16/(tmpu+120.))*(tmpu/296.16)**1.5

!  Thermal velocity of air molecule
   vt = sqrt(8.*kb*tmpu/pi/m_air)

!  Air molecule mean free path
   rmfp = 2.*rmu/rhoa/vt

!  Knudsen number
   rkn = rmfp/radius

!  Cunningham slip correction factor
   bpm = 1. + 1.246*rkn + 0.42*rkn*exp(-0.87/rkn)

!  Brownian diffusion coefficient
   diff_coef = kb*tmpu*bpm/3./pi/rmu/(2.*radius)

!  Fall speed (assumes Reynolds # < 0.01)
   vsettle = 2./9.*rhop*radius**2.*grav*bpm/rmu

!  Check the Reynold's number to see if we need a drag correction
!  First guess at Reynold's number using Stoke's calculation
   re = 2.*rhoa*radius*vsettle/rmu

!  If Re > 0.01 then apply drag correction following Pruppacher and
!  Klett regime 2 (eq. 10-142).  Assuming reasonable aerosols we
!  do not consider that particle Re may exceed 300.
   if(re .gt. 0.01) then
    x = log(24.*re/bpm)
    y = -3.18657 + 0.992696   *x     - .00153193   *x**2. &
                 - 0.000987059*x**3. - .000578878  *x**4. &
                 + 8.55176E-05*x**5. -  3.27815E-06*x**6.
    re = exp(y)*bpm
    vsettle = rmu*re/2./rhoa/radius
   endif


   end subroutine Chem_CalcVsettle2Gorig

!==================================================================================
!BOP
! !IROUTINE: Chem_SettlingSimpleOrig

   subroutine Chem_SettlingSimpleOrig ( km, klid, flag, grav, cdt, radiusInp, rhopInp, &
                                        int_qa, tmpu, rhoa, rh, delp, hghte, &
                                        fluxout, rc, vsettleOut, correctionMaring )


! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer,    intent(in)  :: km        ! total model levels
   integer,    intent(in)  :: klid      ! index for pressure level lid
   integer,    intent(in)  :: flag      ! flag to control particle swelling (see note)
   real,       intent(in)  :: grav      ! gravity [m/sec^2]
   real,       intent(in)  :: cdt       ! chemistry model time-step [sec]
   real, intent(in)        :: radiusInp ! aerosol radius  microns
   real, intent(in)        :: rhopInp   ! aerosol density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(inout) :: int_qa ! aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(in) :: tmpu ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa  ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: rh    ! relative humidyt [1]
   real, pointer, dimension(:,:,:), intent(in) :: hghte ! geopotential height [m]
   real, pointer, dimension(:,:,:), intent(in) :: delp  ! air pressure thickness [Pa]
   logical, optional, intent(in)    :: correctionMaring

! !OUTPUT PARAMETERS:
   real, pointer, dimension(:,:), intent(inout)  :: fluxout
   real, dimension(:,:,:), optional, intent(out)  :: vsettleOut
   integer, optional, intent(out)   :: rc

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
! 17Aug2020 E.Sherman - Ported and modified from Chem_SettlingMod.F90

! !Local Variables
   real,    parameter     :: rhow = 1000.  ! Density of water [kg m-3]
!  parameter from Gerber 1985 (units require radius in cm, see rcm)
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424
!  parameters for ammonium sulfate
   real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

   integer         :: i1=1, i2, j1=1, j2
   integer         :: nSubSteps, ijl, dk
   integer         :: i, j, k, iit

   real, allocatable    :: dz(:,:,:)
   real :: radius, rhop   ! particle radius and density passed to
                          ! fall velocity calculation
   real :: minTime, qmin
   real :: sat, rrat
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2
   real :: rcm
   real :: diff_coef                 ! Brownian diffusion coefficient [m2 s-1]
   real(kind=DP)  :: dt_settle, gravDP
!   real, allocatable, dimension(:,:)   :: hsurf
   real, pointer, dimension(:,:)   :: hsurf

   real, dimension(:,:,:), allocatable  :: vsettle   ! fall speed [m s-1]
   real(kind=DP), dimension(:,:,:), allocatable   :: dzd, vsd, qa, qa_temp
   real(kind=DP), dimension(:,:), allocatable  :: cmass_before, cmass_after, qdel, &
        d_p, d_pm1, qsrc

!EOP
!-------------------------------------------------------------------------
!  Begin

   gravDP = grav

   i2 = ubound(hghte,1)
   j2 = ubound(hghte,2)

   hsurf => hghte(i1:i2,j1:j2,km)

!  Allocate arrays
!  ---------------
   allocate(dz, mold=rhoa)
   allocate(vsettle(i2,j2,km), source=0.0)
   allocate(dzd(i2,j2,km), vsd(i2,j2,km), qa(i2,j2,km), qa_temp(i2,j2,km), source=0.0_DP)
   allocate(cmass_before(i2,j2), cmass_after(i2,j2), qdel(i2,j2), d_p(i2,j2), &
            d_pm1(i2,j2), qsrc(i2,j2), source=0.0_DP)

!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1

!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo
   dzd = dz

!   qa = w_c%qa(ibin)%data3d
   qa = int_qa

    radius = radiusInp
    rhop = rhopInp

!   Reset a (large) minimum time to cross a grid cell in settling
    minTime = cdt

    if( associated(fluxout) ) fluxout = 0.0
    cmass_before(:,:) = 0.d0
    cmass_after(:,:) = 0.d0

!   If radius le 0 then get out
    if(radius .le. 0.) return

    do k = klid, km
     do j = j1, j2
      do i = i1, i2

!      Find the column dry mass before sedimentation
       cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k)/gravDP * delp(i,j,k)

!      Adjust the particle size for relative humidity effects
       sat = max(rh(i,j,k),tiny(1.0)) ! to avoid zero FPE

!      Fitzgerald
       if(flag .eq. 1 .and. sat .ge. 0.80) then
!       parameterization blows up for RH > 0.995, so set that as max
!       rh needs to be scaled 0 - 1
        sat = min(0.995,sat)
!       Calculate the alpha and beta parameters for the wet particle
!       relative to amonium sulfate
        beta = exp( (0.00077*sat) / (1.009-sat) )
        if(sat .le. 0.97) then
         theta = 1.058
        else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
        endif
        alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
        f1 = 10.2 - 23.7*sat + 14.5*sat**2.
        f2 = -6.7 + 15.5*sat - 9.2*sat**2.
        alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
        alpha = alphaNaCl * (alpha1*alpharat)
!       radius is the radius of the wet particle
        radius = alpha * radiusInp**beta
        rrat = (radiusInp/radius)**3.
        rhop = rrat*rhopInp + (1.-rrat)*rhow
       elseif(flag .eq. 2) then   ! Gerber
        sat = min(0.995,sat)
        rcm = radiusInp*100.
        radius = 0.01 * (   c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                          + rcm**3.)**(1./3.)
        rrat = (radiusInp/radius)**3.
        rhop = rrat*rhopInp + (1.-rrat)*rhow
       elseif(flag .eq. 3) then
!       Gerber parameterization for Ammonium Sulfate
        sat = min(0.995,sat)
        rcm = radiusInp*100.
        radius = 0.01 * (   SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                      + rcm**3.)**(1./3.)
        rrat = (radiusInp/radius)**3.
        rhop = rrat*rhopInp + (1.-rrat)*rhow
       elseif(flag .eq. 4) then
!       Petters and Kreidenweis (ACP2007) parameterization
        sat = min(0.99,sat)
        radius = (radiusInp**3 * (1+1.19*sat/(1-sat)))**(1./3.)
        rrat = (radiusInp/radius)**3
        rhop = rrat*rhopInp + (1.-rrat)*rhow
       endif

!      Calculate the settling velocity
       call Chem_CalcVsettle2Gorig(radius, rhop, rhoa(i,j,k), tmpu(i,j,k), &
                                   grav, diff_coef, vsettle(i,j,k))
      end do
     end do
    end do

    if(present(correctionMaring)) then
     if ((correctionMaring) .and. (radiusInp .le. (0.5*diameterMaring))) then
       vsettle = max(1.0e-9, vsettle - v_upwardMaring)
     endif
    endif

    vsd = vsettle

    if(present(vsettleOut)) then
       vsettleOut = vsettle
    endif

!   Determine global min time to cross grid cell
    qmin = minval(dz/vsettle)
    minTime = min(minTime,qmin)

!   Now, how many iterations do we need to do?
    if ( minTime < 0 ) then
         nSubSteps = 0
!         call mpout_log(myname,'no Settling because minTime = ', minTime )
    else if(minTime .ge. cdt) then
     nSubSteps = 1
     dt_settle = cdt
    else
     nSubSteps = cdt/minTime+1
     dt_settle = cdt/nSubSteps
    endif

!   Loop over sub-timestep
    do iit = 1, nSubSteps

!     Try a simple forward Euler scheme
     qdel = qa(i1:i2,j1:j2,klid)*dt_settle*vsd(i1:i2,j1:j2,klid)/dzd(i1:i2,j1:j2,klid)
     qa(i1:i2,j1:j2,klid) = qa(i1:i2,j1:j2,klid) - qdel

     do k = klid+1, km
      d_p   = delp(i1:i2,j1:j2,k)
      d_pm1 = delp(i1:i2,j1:j2,k-1)
      qsrc = qdel * d_pm1 / d_p
      qdel = qa(i1:i2,j1:j2,k)*dt_settle*vsd(i1:i2,j1:j2,k)/dzd(i1:i2,j1:j2,k)
      qa(i1:i2,j1:j2,k) = qa(i1:i2,j1:j2,k) - qdel + qsrc
     enddo

    end do  ! iit

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = klid, km
     do j = j1, j2
      do i = i1, i2
       cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k)/ gravDP * delp(i,j,k)
      enddo
     enddo
    enddo

    if( associated(fluxout) ) then
       fluxout(i1:i2,j1:j2) &
        = (cmass_before(i1:i2,j1:j2) - cmass_after(i1:i2,j1:j2))/cdt
    endif

    int_qa = qa

      __RETURN__(__SUCCESS__)
   end subroutine Chem_SettlingSimpleOrig


!============================================================================
!BOP
!
! !IROUTINE: DryDeposition - Calculate aerosol dry deposition for lowest layer
!
! !INTERFACE:
!
   subroutine DryDeposition ( km, tmpu, rhoa, hghte, oro, ustar, pblh, shflux, &
                              von_karman, cpd, grav, z0h, drydepf, rc, &
                              radius, rhop, u10m, v10m, fraclake, gwettop )

! !USES:
  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km       ! total model levels
   real, pointer, dimension(:,:,:), intent(in) :: tmpu    ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa    ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: hghte   ! top of layer geopotential height [m]
   real, pointer, dimension(:,:), intent(in)   :: oro     ! orography flag
   real, pointer, dimension(:,:), intent(in)   :: ustar   ! friction speed [m/sec]
   real, pointer, dimension(:,:), intent(in)   :: pblh    ! PBL height [m]
   real, pointer, dimension(:,:), intent(in)   :: shflux  ! sfc. sens. heat flux [W m-2]
   real, intent(in)                :: von_karman ! Von Karman constant [unitless]
   real, intent(in)                :: cpd       ! thermodynamic constant, specific heat of something?
   real, intent(in)                :: grav      ! gravity [m/sec^2]
   real, pointer, dimension(:,:)   :: z0h       ! rough height, sens. heat [m]

! !OUTPUT PARAMETERS:
   real, intent(inout)        :: drydepf(:,:)     ! Deposition frequency [1/sec]
   integer, intent(out)          :: rc          ! Error return code:

! !OPTIONAL PARAMETERS:
!  If these parameters are provided we compute a "resuspension" term as
!  if the particles are lifted like dust
   real, optional                            :: radius    ! particle radius [m]
   real, optional                            :: rhop      ! particle density [kg/m^3]
   real, pointer, dimension(:,:), optional   :: u10m      ! 10-m u-wind component [m/sec]
   real, pointer, dimension(:,:), optional   :: v10m      ! 10-m v-wind component [m/sec]
   real, pointer, dimension(:,:), optional   :: fraclake  ! fraction covered by water [1]
   real, pointer, dimension(:,:), optional   :: gwettop   ! fraction soil moisture [1]


! !DESCRIPTION: Calculates the deposition velocity for aerosols in the lowest
!               model layer.
!
! !REVISION HISTORY:
!
!  08Jan2010 - Colarco, based on GOCART implementation, does not
!                       include any size dependent deposition term

! !Local Variables
   character(len=*), parameter :: myname = 'DryDeposition'
   integer, parameter :: DP=kind(1.0d0)
   integer :: i, j
   integer :: dims(3)
   integer :: i1=1, i2, j1=1, j2
   real, parameter :: rhow = 1000.      ! density of water [kg m-3]
   real, parameter :: coll_size = 0.002 ! collector size [m]
   real, allocatable :: dz(:,:)     ! lowest layer thickness
   real, allocatable :: rmu(:,:)    ! dynamic viscosity [kg m-1 s-1]
   real, allocatable :: Ra(:,:)     ! aerodynamic resistance
   real, allocatable :: Rs(:,:)     ! surface resistance
   real, allocatable :: vdep(:,:)   ! Deposition speed [m s-1]
   real, allocatable :: obk(:,:)    ! Obukhov Length [m]

   real(kind=DP) :: Rttl        ! total surface resistance

   real(kind=DP) :: R2, w10m, u_thresh0
   real(kind=DP) :: vds, vdsmax, czh
   real(kind=DP) :: frac, cz, psi_h, eps, logmfrac, z0h_min, z0h_
   real(kind=DP) :: one = 1.0, zero = 0.0
!
!EOP
!-------------------------------------------------------------------------
!  Begin...

   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

   allocate(dz(i2,j2),rmu(i2,j2),Ra(i2,j2),Rs(i2,j2),vdep(i2,j2), &
            obk(i2,j2), source=0.0)

!  Calculate the viscosity and thickness of the surface level
   dz = hghte(:,:,km-1) - hghte(:,:,km)
   rmu = 1.8325e-5*(416.16/(tmpu(i1:i2,j1:j2,km)+120.)) &
                  *(tmpu(i1:i2,j1:j2,km)/296.16)**1.5

   z0h_min = 100. * tiny(1.0)  ! because sometimes we may get z0h=0.

!  Calculate the Obukhov length scale
!  -----------------------------------
   call ObukhovLength2G( i1, i2, j1, j2, von_karman, cpd, grav, &
                        tmpu(:,:,km), rhoa(:,:,km), shflux, ustar, &
                        obk )

!  Aerodynamic Resistance
!  psi_h and Ra are equations 2, 4-5 of Walcek et al. 1986 Atmospheric Environment
!  ----------------------------
   do j = j1, j2
    do i = i1, i2

      cz = dz(i,j) / 2.
      frac = cz / obk(i,j)
      if(frac .gt. 1.) frac = 1.
      if(frac .gt. 0. .and. frac .le. 1.) then
       psi_h = -5.0*frac
      else if (frac .lt. 0.) then
       eps = min(one,-frac)
       logmfrac = log(eps)
       psi_h = exp(0.598 + 0.39*logmfrac - 0.09*(logmfrac)**2.)
      endif

      z0h_ = max ( z0h(i,j), z0h_min )

      Ra(i,j) = (log(cz/z0h_) - psi_h) / (von_karman*ustar(i,j))

    enddo
   enddo

!  Surface Resistance term for aerosols
!  Rs formulation from eqn. 15 - 18 of Walcek et al. 1986 Atmospheric Environment
!  Loop over space
!  -------------------------
   do j = j1, j2
    do i = i1, i2

!     Calculate the surface resistance term
      vds = 0.002*ustar(i,j)
!     Set to small value of vds if ustar too small
      vds = max(vds, 0.002 * 0.00001)
      if(obk(i,j) .lt. 0.) vds = vds*(1.+(-300./obk(i,j))**0.6667)
      czh = pblh(i,j)/obk(i,j)
      if(czh .lt. -30.) vds = 0.0009*ustar(i,j)*(-czh)**0.6667
!     vdsMax is from Table 2 of Walcek et al. 1986
!     There are actually seasonal and regionally varying values,
!     but for most of the world a value of 1.0 cm s-1 is used.
      vdsMax = 0.01

      Rs(i,j) = 1./min(vds,vdsmax)

      if(Rs(i,j) .gt. 9999.) Rs(i,j) = 9999.
      if(Rs(i,j) .lt. 1.)    Rs(i,j) = 1.

!     If doing dust over land, possibly re-emit
!     Logic is to check on optional provided parameter and modify R2
      R2 = 1.
      if(present(fraclake) .and. present(u10m) .and. present(v10m) .and. &
         present(radius) .and. present(rhop) .and. present(gwettop)) then

!      Calculate the threshold velocity for dust emissions
       u_thresh0 = 0.13 * sqrt(rhop*grav*2.*radius/rhoa(i,j,km)) &
                        * sqrt(1.+6.e-7/(rhop*grav*(2.*radius)**2.5)) &
              / sqrt(1.928*(1331.*(100.*2.*radius)**1.56+0.38)**0.092 - 1.)
       w10m = sqrt(u10m(i,j)**2. + v10m(i,j)**2.)

!      Calculate the coefficient for resuspension
       if(oro(i,j) .eq. OCEAN) then
        R2 = 1.
       else
        R2 = fraclake(i,j)+(1.-fraclake(i,j)) &
                           *( gwettop(i,j)+(1.-gwettop(i,j)) &
                             *exp(-max(zero,(w10m-u_thresh0))))
       endif
      endif
!     Now what is the deposition velocity
      Rttl = Ra(i,j) + Rs(i,j)

      vdep(i,j) = 1./Rttl*R2

!     Set a minimum value of deposition velocity
      vdep(i,j) = max(vdep(i,j),1.e-4)

!     Save the dry deposition frequency for the chemical removal terms
!     in units of s-1
      drydepf(i,j) = max(0.,vdep(i,j) / dz(i,j))

    end do  ! i
   end do   ! j

   rc = __SUCCESS__

   end subroutine DryDeposition

!====================================================================================
! !IROUTINE: ObukhovLength2G - Calculate the Obukhov length scale stability parameter
!
! !INTERFACE:
!
!  =========================================================================
!  Calculate the Obukhov length scale
!  Wesely and Hicks (1977) Journal of Air Pollution Control Association
!  Equation 9.  Note: we have adopted this from GOCART, which neglected
!  the latent heat of evaporation term.  Also, we are using surface
!  mid-layer values of air density and temperature, not absolute surface
!  values (and not the potential temperature, either).

   subroutine ObukhovLength2G ( i1, i2, j1, j2, von_karman, cpd, grav, &
                              t, rhoa, shflux, ustar, &
                              obk )

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2
   real, intent(in) :: von_karman ! Von Karman constant [unitless]
   real, intent(in) :: cpd        ! thermodynamic constant, specific heat of something?
   real, intent(in) :: grav       ! gravity [m/sec^2]
   real, dimension(i1:i2,j1:j2)  :: t         ! temperature [K]
   real, dimension(i1:i2,j1:j2)  :: rhoa      ! air density [kg/m^3]
   real, pointer, dimension(:,:) :: ustar     ! friction speed [m/sec]
   real, pointer, dimension(:,:) :: shflux    ! sfc. sens. heat flux [W/m^2]

! !OUTPUT PARAMETERS
   real, dimension(i1:i2,j1:j2)  :: obk       ! Obukhov length [m]

!  Local

!  Calculate the Monin-Obhukov length:
!          -Air density * Cp * T(surface) * Ustar^3
!   OBK = -------------------------------------------
!               vK * g * Sensible heat flux
!  vK = 0.4               von Karman constant
!  Cp = 1000 J kg-1 K-1   specific heat of air at constant pressure
!  If OBK < 0 the air is unstable; if OBK > 0 the air is stable
!  For sensible heat flux of zero OBK goes to infinity (set to 1.e5)

   obk = 1.e5
   where(abs(shflux) > 1.e-32) &
       obk =   - rhoa * cpd * t * ustar**3. &
             / (von_karman * grav * shflux)

   return
   end subroutine ObukhovLength2G

!==================================================================================
!BOP

! !IROUTINE: WetRemovalGOCART2G
   subroutine WetRemovalGOCART2G ( km, klid, n1, n2, bin_ind, cdt, aero_type, kin, grav, fwet, &
                                   aerosol, ple, tmpu, rhoa, pfllsan, pfilsan, &
                                   precc, precl, fluxout, rc )

! !USES:
  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km  ! total model levels
   integer, intent(in) :: klid ! index for pressure lid
   integer, intent(in) :: n1  ! total number of bins (probably can be removed)
   integer, intent(in) :: n2  ! total number of bins (probably can be removed)
   integer, intent(in) :: bin_ind ! bin index (usually the loop iteration)
   real, intent(in)    :: cdt     ! chemistry model time-step [sec]
   character(len=*)    :: aero_type
   logical, intent(in) :: KIN ! true for aerosol
   real, intent(in)    :: grav    ! gravity [m/sec^2]
   real, intent(in)    :: fwet
   real, dimension(:,:,:), intent(inout) :: aerosol  ! internal state aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(in)  :: ple     ! pressure level thickness [Pa]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu    ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa    ! moist air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in)  :: pfllsan ! 3D flux of liquid nonconvective precipitation [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:), intent(in)  :: pfilsan ! 3D flux of ice nonconvective precipitation [kg/(m^2 sec)]
   real, pointer, dimension(:,:), intent(in)    :: precc   ! surface convective rain flux [kg/(m^2 sec)]
   real, pointer, dimension(:,:), intent(in)    :: precl   ! Non-convective precipitation [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:)  :: fluxout ! tracer loss flux [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, intent(out)             :: rc          ! Error return code:

! !DESCRIPTION: Calculates the updated species concentration due to wet
!               removal.  As written, intended to function for large
!               scale (not convective) wet removal processes
!
! !REVISION HISTORY:
!
!  08Jan2010 - Colarco, based on GOCART implementation, does not
!                       include any size dependent term
!
! !Local Variables
   character(len=*), parameter :: myname = 'WetRemovalGOCART2G'
   integer, parameter :: DP=kind(1.0d0)
   integer  ::  i, j, k, n, LH, kk, ios, nbins
   integer  :: i1=1, i2, j1=1, j2, dims(3)
   real, allocatable, dimension(:,:,:) :: pdog   ! air mass factor dp/g [kg m-2]
   real, allocatable, dimension(:,:,:) :: delz   ! box height  dp/g/rhoa [m]
   real :: pls, pcv, pac             ! ls, cv, tot precip [mm day-1]
   real :: qls(km), qcv(km)          ! ls, cv portion of moisture tendency[kg m-3 s-1]
   real :: qmx, qd, A                ! temporary variables on moisture
   real :: F, B, BT                  ! temporary variables on cloud, freq.
   real :: WASHFRAC, WASHFRAC_F_14
   real, allocatable :: fd(:,:)      ! flux across layers [kg m-2]
   real, allocatable :: dpfli(:,:,:)  ! vertical gradient of LS ice+rain precip flux
   real, allocatable :: DC(:)        ! scavenge change in mass mixing ratio
   real, allocatable, dimension(:,:,:) :: c_h2o, cldliq, cldice

!  Rain parameters from Liu et al.
   real, parameter :: B0_ls = 1.0e-4
   real, parameter :: F0_ls = 1.0
   real, parameter :: XL_ls = 5.0e-4
   real, parameter :: B0_cv = 1.5e-3
   real, parameter :: F0_cv = 0.3
   real, parameter :: XL_cv = 2.0e-3
   real, parameter :: k_wash = 1.d0  ! first order washout rate, constant, [cm^-1]
!  Duration of rain: ls = model timestep, cv = 1800 s (<= cdt)
   real            :: Td_ls
   real, parameter :: Td_cv = 1800.
   real(kind=DP), PARAMETER   :: R = 8.2057d-2  ! universal gas constant [L*atm/moles/K]
   real(kind=DP), PARAMETER   :: INV_T0 = 1d0 / 298d0
   real(kind=DP), PARAMETER   :: conv_NH3 = 5.69209978831d-1 ! 0.6*SQRT(0.9) for ice to gas ratio
   real(kind=DP)  :: k_rain, Kstar298, H298_R, I2G, L2G, C_TOT, F_L, F_I
   real(kind=DP)  :: PP, LP

   logical :: snow_scavenging

!  Efficiency of dust wet removal (since dust is really not too hygroscopic)
!  Applied only to in-cloud scavenging
   real :: effRemoval

   rc = __SUCCESS__

!EOP
!-----------------------------------------------------------------------------
!  Begin...

   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

!  Allocate arrays
   allocate(c_h2o(i2,j2,km), cldliq(i2,j2,km), cldice(i2,j2,km), pdog(i2,j2,km), &
            delz(i2,j2,km), dpfli(i2,j2,km), source=0.0)

!  Initialize local variables
!  --------------------------
!  c_h2o, cldliq, and cldice are respectively intended to be the
!  water mixing ratio (liquid or vapor?, in or out of cloud?)
!  cloud liquid water mixing ratio
!  cloud ice water mixing ratio
   c_h2o  = (10d0**(-2663.5d0/tmpu(:,:,:) + 12.537d0 ) ) /  &
                   (ple(:,:,0:km-1)+ple(:,:,1:km)) /2d0
   cldliq = 0.d0
   where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.d0) / 20.d0 )
   where(tmpu >= 268.) cldliq = 1.d-6
   cldice = 1.d-6 - cldliq

   Td_ls = cdt
   nbins = n2-n1+1

   allocate(fd(km,nbins),source=0.0,stat=ios)
   allocate(dc(nbins),source=0.0,stat=ios)

   if( associated(fluxout) ) fluxout(i1:i2,j1:j2,bin_ind) = 0.0

!  Accumulate the 3-dimensional arrays of rhoa and pdog
   pdog = (ple(:,:,1:km)-ple(:,:,0:km-1)) / grav
   delz = pdog / rhoa
   dpfli = pfllsan(:,:,1:km)-pfllsan(:,:,0:km-1)+pfilsan(:,:,1:km)-pfilsan(:,:,0:km-1)
   if (.not. KIN) then              ! Gases
      if (aero_type == 'NH3') then  ! Only for NH3 at present
      ! values adopted in Umich/IMPACT and GMI, effective Henry's law coefficient at pH=5
        Kstar298 = 1.05d6
        H298_R = -4.2d3
      else
        !$omp critical (G2G_proc_1)
        print *, 'stop in WetRemoval, need Kstar298 and H298_R'
        !$omp end critical (G2G_proc_1)
        rc = __FAIL__
        return
      endif
   endif

!  Snow scavenging flag
   snow_scavenging = .true.

   if ( (aero_type == 'OC'      ) .or. &
        (aero_type == 'sea_salt') .or. &
        (aero_type == 'sulfur'  ) .or. &
        (aero_type == 'seasalt' ) .or. &
        (aero_type == 'sulfate' ) .or. &
        (aero_type == 'NH3'     ) .or. &
        (aero_type == 'NH4a'    ) .or. &
        (aero_type == 'nitrate' ) .or. &
        (aero_type == 'bromine' ) .or. &
        (aero_type == 'dust'    ) ) then
     snow_scavenging = .false.
   end if

!  Loop over spatial indices
   do j = j1, j2
    do i = i1, i2

!    Check for total precipitation amount
!    Assume no precip in column if precl+precc = 0
     pac = precl(i,j) + precc(i,j)
     if(pac .le. 0.) cycle
     pls = precl(i,j)
     pcv = precc(i,j)

!    Initialize the precipitation fields
     qls(:)  = 0.
     qcv(:)  = 0.
     fd(:,:) = 0.

!    Find the highest model layer experiencing rainout.  Assumes no
!    scavenging if T < 258 K
     LH = 0
     do k = klid, km
      if(dpfli(i,j,k) .gt. 0. ) then
        LH = k
        exit
      endif
     end do

     if(LH .lt. 1) cycle

     do k = LH, km
      qls(k) = dpfli(i,j,k)/pdog(i,j,k)*rhoa(i,j,k)
     end do

!    Loop over vertical to do the scavenging!
     do k = LH, km

!-----------------------------------------------------------------------------
!   (1) LARGE-SCALE RAINOUT:
!       Tracer loss by rainout = TC0 * F * exp(-B*dt)
!         where B = precipitation frequency,
!               F = fraction of grid box covered by precipitating clouds.
!       We assume that tracer scavenged by rain is falling down to the
!       next level, where a fraction could be re-evaporated to gas phase
!       if Qls is less then 0 in that level.
!-----------------------------------------------------------------------------
      if (qls(k) .gt. tiny(0.)) then
       F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
       k_rain  = B0_ls/F0_ls +1./(F0_ls*XL_ls/qls(k))
       if ( kin ) then     ! Aerosols
          B = k_rain
       else                ! Gases
        ! ice to gas ratio
          if ( c_h2o(i,j,k) > 0.d0) then
             I2G = (cldice(i,j,k) / c_h2o(i,j,k)) * conv_NH3
          else
             I2G = 0.d0
          endif
          L2G = cldliq(i,j,k) * R * tmpu(i,j,k) * &
                  Kstar298 * EXP( -H298_R * ( ( 1d0 / tmpu(i,j,k) ) - INV_T0 ) )
        ! fraction of NH3 in liquid & ice phases
          C_TOT = 1d0 + L2G + I2G
          F_L = L2G / C_TOT
          F_I = I2G / C_TOT
        ! compute kg, the retention factor for liquid NH3 is 0 at T < 248K and
        ! 0.05 at 248K < T < 268K
          if (tmpu(i,j,k) >=268d0) then
             B = k_rain * ( F_L+F_I )
          elseif ( (248d0 < tmpu(i,j,k)) .and. (tmpu(i,j,k) < 268d0) ) then
             B = k_rain * ( (0.05*F_L)+F_I )
          else
             B = k_rain * F_I
          endif
       endif ! kin
       BT = B * Td_ls
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >
!      Adjust du level:
       do n = 1, nbins
! supress scavenging at cold T except for HNO3
        if (tmpu(i,j,k) < 258d0 .and. .not.snow_scavenging) then
            F = 0.d0
        endif

        effRemoval = fwet
        DC(n) = aerosol(i,j,k) * F * effRemoval *(1.-exp(-BT))
        if (DC(n).lt.0.) DC(n) = 0.
        aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
        if (aerosol(i,j,k) .lt. 1.0E-32) aerosol(i,j,k) = 1.0E-32
       end do
!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = DC(n)*pdog(i,j,k)
       end do

      end if                                    ! if Qls > 0  >>>

!-----------------------------------------------------------------------------
! * (2) LARGE-SCALE WASHOUT:
! *     Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------
      if(k .gt. LH .and. qls(k) .ge. 0.) then
       if(qls(k) .lt. qls(k-1)) then
!       Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1,LH,-1
         if (Qls(kk).gt.0.) then
          Qmx = max(Qmx,Qls(kk))
         else
          exit
         end if
        end do

        F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(Qmx*cdt/Td_ls))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

!       Aerosols
        Qd = Qmx /rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Gases
        if ( .not. KIN ) then
           IF ( tmpu(i,j,k) >= 268d0 ) THEN
            !------------------------
            ! T >= 268K: Do washout
            !------------------------
            ! Rainwater content in the grid box (Eq. 17, Jacob et al, 2000)
            PP = (PFLLSAN(i,j,k)/1000d0 + PFILSAN(i,j,k)/917d0 )*100d0 ! from kg H2O/m2/s to cm3 H2O/cm2/s
            LP = ( PP * cdt ) / ( F * delz(i,j,k)*100.d0 )  ! DZ*100.d0 in cm
            ! Compute liquid to gas ratio for H2O2, using the appropriate
            ! parameters for Henry's law -- also use rainwater content Lp
            ! (Eqs. 7, 8, and Table 1, Jacob et al, 2000)
            !CALL COMPUTE_L2G( Kstar298, H298_R, tmpu(i,j,k), LP, L2G )
            L2G = Kstar298 * EXP( -H298_R*((1d0/tmpu(i,j,k))-INV_T0) ) &
                  * LP * R * tmpu(i,j,k)
            ! Washout fraction from Henry's law (Eq. 16, Jacob et al, 2000)
            WASHFRAC = L2G / ( 1d0 + L2G )
            ! Washout fraction / F from Eq. 14, Jacob et al, 2000
            ! Note: WASHFRAC_F_14 should match what's used for HNO3 (hma, 13aug2011)
            WASHFRAC_F_14 = 1d0 - EXP( -K_WASH * ( PP / F ) * cdt )
            ! Do not let the Henry's law washout fraction exceed
            IF ( WASHFRAC > WASHFRAC_F_14 ) THEN
              WASHFRAC = WASHFRAC_F_14
            ENDIF
           ELSE
            !------------------------
            ! T < 268K: No washout
            !------------------------
            WASHFRAC = 0d0
           ENDIF
        endif

!       Adjust du level:
        do n = 1, nbins
         if ( KIN ) then
            DC(n) = aerosol(i,j,k) * F * (1.-exp(-BT))
         else
            DC(n) = aerosol(i,j,k) * F * WASHFRAC
         endif
         if (DC(n).lt.0.) DC(n) = 0.
         aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
         if (aerosol(i,j,k) .lt. 1.0E-32) &
          aerosol(i,j,k) = 1.0E-32
         if( associated(fluxout)) then
          fluxout(i,j,bin_ind) = fluxout(i,j,bin_ind)+DC(n)*pdog(i,j,k)/cdt

         endif
        end do

       end if
      end if                                    ! if ls washout  >>>

!-----------------------------------------------------------------------------
!  (3) CONVECTIVE RAINOUT:
!      Tracer loss by rainout = dd0 * F * exp(-B*dt)
!        where B = precipitation frequency,
!              F = fraction of grid box covered by precipitating clouds.
!-----------------------------------------------------------------------------

      if (qcv(k) .gt. 0.) then
       F  = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qcv(k)*cdt/Td_cv))
       B  = B0_cv
       BT = B * Td_cv
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >

!      Adjust du level:
       do n = 1, nbins
        effRemoval = fwet
        DC(n) = aerosol(i,j,k) * F * effRemoval * (1.-exp(-BT))
        if (DC(n).lt.0.) DC(n) = 0.
        aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
        if (aerosol(i,j,k) .lt. 1.0E-32) aerosol(i,j,k) = 1.0E-32
       end do

!------  Flux down:  unit is kg. Including both ls and cv.
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + DC(n)*pdog(i,j,k)
       end do

      end if                                  ! if Qcv > 0   >>>

!-----------------------------------------------------------------------------
!  (4) CONVECTIVE WASHOUT:
!      Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------

      if (k.gt.LH .and. Qcv(k).ge.0.) then
       if (Qcv(k).lt.Qcv(k-1)) then
!-----  Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1, LH, -1
         if (Qcv(kk).gt.0.) then
          Qmx = max(Qmx,Qcv(kk))
         else
          exit
         end if
        end do

        F = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qmx*cdt/Td_cv))
        if (F.lt.0.01) F = 0.01

!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx / rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust du level:
        do n = 1, nbins
         DC(n) = aerosol(i,j,k) * F * (1.-exp(-BT))
         if (DC(n).lt.0.) DC(n) = 0.
         aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
         if (aerosol(i,j,k) .lt. 1.0E-32) &
          aerosol(i,j,k) = 1.0E-32
         if( associated(fluxout)) then
          fluxout(i,j,bin_ind) = fluxout(i,j,bin_ind)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if cv washout  >>>

!-----------------------------------------------------------------------------
!  (5) RE-EVAPORATION.  Assume that SO2 is re-evaporated as SO4 since it
!      has been oxidized by H2O2 at the level above.
!-----------------------------------------------------------------------------
!     Add in the flux from above, which will be subtracted if reevaporation occurs
      if(k .gt. LH) then
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + Fd(k-1,n)
       end do

!      Is there evaporation in the currect layer?
       if (dpfli(i,j,k) .lt. 0.) then
!       Fraction evaporated = H2O(k)evap / H2O(next condensation level).
        if (dpfli(i,j,k-1) .gt. 0.) then

          A =  abs(  dpfli(i,j,k) /  dpfli(i,j,k-1)  )
          if (A .gt. 1.) A = 1.

!         Adjust tracer in the level
          do n = 1, nbins
           DC(n) =  Fd(k-1,n) / pdog(i,j,k) * A
           aerosol(i,j,k) = aerosol(i,j,k) + DC(n)
           aerosol(i,j,k) = max(aerosol(i,j,k),1.e-32)
           Fd(k,n)  = Fd(k,n) - DC(n)*pdog(i,j,k)
          end do

        endif
       endif                                   ! if -moistq < 0
      endif
     end do  ! k


     do n = 1, nbins
      if( associated(fluxout)) then
       fluxout(i,j,bin_ind) = fluxout(i,j,bin_ind)+Fd(km,n)/cdt
      endif
     end do

    end do   ! i
   end do    ! j

   deallocate(fd,DC,dpfli,stat=ios)

   end subroutine WetRemovalGOCART2G

!=============================================================================
!BOP

! !IROUTINE: UpdateAerosolState
   subroutine UpdateAerosolState (emissions, emissions_surface, emissions_point, &
                                  sfrac, nPts, km, cdt, grav, nbins, delp, aero, rc)

! !USES:
  implicit NONE

! !INPUT PARAMETERS:
!   real, pointer, dimension(:,:)     :: emissions_surface
   real, dimension(:,:,:), intent(in)     :: emissions_surface
   real, dimension(:,:,:,:), intent(inout) :: emissions
   real, dimension(:,:,:), intent(in) :: emissions_point

   real, dimension(:), intent(in)  :: sfrac ! source fraction [1]
   integer, intent(in)             :: nPts  ! number of point emissions
   integer, intent(in)             :: km    ! total model levels
   real, intent(in)                :: cdt   ! chemistry model time-step [sec]
   real, intent(in)                :: grav  ! gravity [m/sec^2]
   integer, intent(in)             :: nbins ! number of aerosol size bins
   real, pointer, dimension(:,:,:), intent(in) :: delp  ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:,:), intent(inout)  :: aero ! aerosol [kg/kg]

! !OUTPUT PARAMETERS:
   integer, intent(out)             :: rc          ! Error return code:

! !DESCRIPTION: Updates internal state variables
!
! !REVISION HISTORY:
!
!  15May2020 - E.Sherman
!
! !Local Variables
   integer :: n, kmin


!EOP
!--------------------------------------------------------------------------------
!   Begin...

    rc = __SUCCESS__

    do n = 1, nbins
       emissions(:,:,km,n) = emissions_surface(:,:,n) * sfrac(n)
       if (nPts > 0) then
          kmin = 1
          emissions(:,:,:,n) = emissions(:,:,:,n) + emissions_point * sfrac(n)
       else
          kmin = km
       end if
       aero(:,:,kmin:km,n) = aero(:,:,kmin:km,n) + emissions(:,:,kmin:km,n) * &
                             cdt * grav / delp(:,:,kmin:km)
    end do

   end subroutine UpdateAerosolState

!==============================================================================

!BOP
!
! !IROUTINE:  Aero_Compute_Diags - Calculate aerosol diagnostics
!
! !INTERFACE:
   subroutine Aero_Compute_Diags (mie, km, klid, nbegin, nbins, rlow, rup, &
                                  wavelengths_profile, wavelengths_vertint, aerosol, &
                                  grav, tmpu, rhoa, rh, u, v, delp, ple,tropp, &
                                  sfcmass, colmass, mass, exttau, stexttau, scatau, stscatau,&
                                  sfcmass25, colmass25, mass25, exttau25, scatau25, &
                                  fluxu, fluxv, conc, extcoef, scacoef, bckcoef,&
                                  exttaufm, scataufm, angstrom, aerindx, NO3nFlag, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   type(GOCART2G_Mie),  intent(in) :: mie        ! mie table
   integer, intent(in) :: km, nbegin, nbins
   integer,    intent(in)    :: klid   ! index for pressure lid
   real, optional, dimension(:), intent(in)    :: rlow   ! bin radii - low bounds
   real, optional, dimension(:), intent(in)    :: rup    ! bin radii - upper bounds
   real, dimension(:), intent(in)    :: wavelengths_profile
   real, dimension(:), intent(in)    :: wavelengths_vertint
   real, dimension(:,:,:,:), intent(in) :: aerosol     !
   real, intent(in) :: grav
   real, pointer, dimension(:,:,:), intent(in) :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa  ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: delp  ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: rh    ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in) :: u     ! east-west wind [m/s]
   real, pointer, dimension(:,:,:), intent(in) :: v     ! north-south wind [m/s]
   real, pointer, dimension(:,:,:), intent(in) :: ple   ! level edge air pressure [Pa]
   real, pointer, dimension(:,:), intent(in)   :: tropp ! tropopause pressure [Pa]
   logical, optional, intent(in)               :: NO3nFlag

! !OUTPUT PARAMETERS:
!  Total mass
   real, optional, dimension(:,:), intent(inout)   :: sfcmass   ! sfc mass concentration kg/m3
   real, optional, dimension(:,:), intent(inout)   :: colmass   ! col mass density kg/m2
   real, optional, dimension(:,:,:), intent(inout) :: mass      ! 3d mass mixing ratio kg/kg
   real, optional, dimension(:,:,:), intent(inout) :: conc      ! 3d mass concentration, kg/m3
!  Total optical properties
   real, optional, dimension(:,:,:), intent(inout)   :: exttau    ! ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: stexttau  ! stratospheric ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: scatau    ! sct. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: stscatau  ! stratospheric sct. AOT at 550 nm
   real, optional, dimension(:,:), intent(inout)   :: sfcmass25 ! sfc mass concentration kg/m3 (pm2.5)
   real, optional, dimension(:,:), intent(inout)   :: colmass25 ! col mass density kg/m2 (pm2.5)
   real, optional, dimension(:,:,:), intent(inout) :: mass25    ! 3d mass mixing ratio kg/kg (pm2.5)
   real, optional, dimension(:,:,:), intent(inout)   :: exttau25  ! ext. AOT at 550 nm (pm2.5)
   real, optional, dimension(:,:,:), intent(inout)   :: scatau25  ! sct. AOT at 550 nm (pm2.5)
   real, optional, dimension(:,:),  intent(inout)  :: aerindx   ! TOMS UV AI
   real, optional, dimension(:,:), intent(inout)   :: fluxu     ! Column mass flux in x direction
   real, optional, dimension(:,:), intent(inout)   :: fluxv     ! Column mass flux in y direction
   real, optional, dimension(:,:,:,:), intent(inout) :: extcoef   ! 3d ext. coefficient, 1/m
   real, optional, dimension(:,:,:,:), intent(inout) :: scacoef   ! 3d scat.coefficient, 1/m
   real, optional, dimension(:,:,:,:), intent(inout) :: bckcoef   ! 3d backscatter coefficient, m-1 sr-1
   real, optional, dimension(:,:,:), intent(inout)   :: exttaufm  ! fine mode (sub-micron) ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)   :: scataufm  ! fine mode (sub-micron) sct. AOT at 550 nm
   real, optional, dimension(:,:), intent(inout)   :: angstrom  ! 470-870 nm Angstrom parameter
   integer, optional, intent(out)   :: rc        ! Error return code:
                                                 !  0 - all is well
                                                 !  1 -

! !DESCRIPTION: Calculates some simple 2d diagnostics from the dust fields
!
! !REVISION HISTORY:
!
!  16APR2004, Colarco
!  11MAR2010, Nowottnick
!  11AUG2020, E.Sherman - refactored to work for multiple aerosols

! !Local Variables
   character(len=*), parameter :: myname = 'Aero_Compute_Diags'
   integer :: i, j, k, n, w, ios, status
   integer :: i1 =1, i2, j1=1, j2
   integer :: ilam470, ilam870
   real, allocatable, dimension(:,:,:) :: tau, ssa, bck
!   real :: fPMfm(nbins)  ! fraction of bin with particles diameter < 1.0 um
!   real :: fPM25(nbins)  ! fraction of bin with particles diameter < 2.5 um
   real, dimension(:), allocatable :: fPMfm  ! fraction of bin with particles diameter < 1.0 um
   real, dimension(:), allocatable :: fPM25  ! fraction of bin with particles diameter < 2.5 um
   logical :: do_angstrom
   real, dimension(:,:), allocatable :: tau470, tau870
   logical   :: NO3nFlag_ !local version of the input

!EOP
!-------------------------------------------------------------------------
!  Begin...

   if( present(NO3nFlag) ) then
      NO3nFlag_ = NO3nFlag
   else
      NO3nFlag_ = .false.
   end if

!  Initialize local variables
!  --------------------------
   i2 = size(rhoa,1)
   j2 = size(rhoa,2)
   allocate(fPMfm(nbins),source=0.0)
   allocate(fPM25(nbins),source=0.0)

!  Get the wavelength indices
!  --------------------------

   ilam470 = mie%getChannel(4.70e-7)
   if(ilam470 <= 0) ilam470 = 0

   ilam870 = mie%getChannel(8.70e-7)
   if(ilam870 <= 0) ilam870 = 0

!  Determine if going to do Angstrom parameter calculation
!  -------------------------------------------------------
   do_angstrom = .false.
!  If both 470 and 870 channels provided (and not the same) then
!  possibly will do Angstrom parameter calculation
   if(ilam470 .ne. 0 .and. &
      ilam870 .ne. 0 .and. &
      ilam470 .ne. ilam870) do_angstrom = .true.

   if( present(angstrom) )  then
      if (do_angstrom ) then
         allocate(tau470(i1:i2,j1:j2), tau870(i1:i2,j1:j2), source=0.0)
      end if
   end if

!  Compute the fine mode (sub-micron) and PM2.5 bin-wise fractions
!  ------------------------------------
   if (present(rlow) .and. present(rup)) then
      call Aero_Binwise_PM_Fractions(fPMfm, 0.50, rlow, rup, nbins)   ! 2*r < 1.0 um
      call Aero_Binwise_PM_Fractions(fPM25, 1.25, rlow, rup, nbins)   ! 2*r < 2.5 um
   end if

   if (present(aerindx))  aerindx = 0.0  ! for now

!  Calculate the diagnostic variables if requested
!  -----------------------------------------------
!  Calculate the surface mass concentration
   if( present(sfcmass) ) then
      sfcmass(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         sfcmass(i1:i2,j1:j2) &
              =   sfcmass(i1:i2,j1:j2) &
              + aerosol(i1:i2,j1:j2,km,n)*rhoa(i1:i2,j1:j2,km)
      end do
   endif
   if( present(sfcmass25) ) then
      sfcmass25(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         sfcmass25(i1:i2,j1:j2) &
              =   sfcmass25(i1:i2,j1:j2) &
              + aerosol(i1:i2,j1:j2,km,n)*rhoa(i1:i2,j1:j2,km)*fPM25(n)
      end do
   endif

!  Calculate the aerosol column loading
   if( present(colmass) ) then
      colmass(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
       do k = klid, km
        colmass(i1:i2,j1:j2) &
         =   colmass(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav
       end do
      end do
   endif
   if( present(colmass25) ) then
      colmass25(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         do k = klid, km
            colmass25(i1:i2,j1:j2) &
             = colmass25(i1:i2,j1:j2) &
             + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*fPM25(n)
       end do
      end do
   endif

!  Calculate the total mass concentration
   if( present(conc) ) then
      conc(i1:i2,j1:j2,1:km) = 0.
      do n = nbegin, nbins
         conc(i1:i2,j1:j2,1:km) &
             = conc(i1:i2,j1:j2,1:km) &
             + aerosol(i1:i2,j1:j2,1:km,n)*rhoa(i1:i2,j1:j2,1:km)
      end do
   endif

!  Calculate the total mass mixing ratio
   if( present(mass) ) then
      mass(i1:i2,j1:j2,1:km) = 0.
      do n = nbegin, nbins
       mass(i1:i2,j1:j2,1:km) &
         =   mass(i1:i2,j1:j2,1:km) &
           + aerosol(i1:i2,j1:j2,1:km,n)
      end do
   endif
   if( present(mass25) ) then
      mass25(i1:i2,j1:j2,1:km) = 0.
      do n = nbegin, nbins
       mass25(i1:i2,j1:j2,1:km) &
         =   mass25(i1:i2,j1:j2,1:km) &
           + aerosol(i1:i2,j1:j2,1:km,n)*fPM25(n)
      end do
   endif

!  Calculate the column mass flux in x direction
   if( present(fluxu) ) then
      fluxu(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         do k = klid, km
           fluxu(i1:i2,j1:j2) &
            = fluxu(i1:i2,j1:j2) &
            + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
         end do
      end do
   endif

!  Calculate the column mass flux in y direction
   if( present(fluxv) ) then
      fluxv(i1:i2,j1:j2) = 0.
      do n = nbegin, nbins
         do k = klid, km
           fluxv(i1:i2,j1:j2) &
           = fluxv(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
         end do
      end do
   endif

   allocate(tau(i1:i2,j1:j2,km),source = 0.)
   allocate(ssa(i1:i2,j1:j2,km),source = 0.)
   allocate(bck(i1:i2,j1:j2,km),source = 0.)
!  Calculate the extinction and/or scattering AOD
   if( present(extcoef)  .or. &
       present(scacoef)  .or. &
       present(bckcoef))   then

      if( present(extcoef) ) extcoef = 0.
      if( present(scacoef) ) scacoef = 0.
      if( present(bckcoef) ) bckcoef = 0.

      do n = nbegin, nbins
        do w = 1, size(wavelengths_profile)
          call mie%Query(wavelengths_profile(w),n,   &
                         aerosol(:,:,:,n)*delp/grav, &
                         rh, tau=tau, ssa=ssa, bbck=bck,__RC__)
!         Calculate the total ext. and scat. coefficients
          if ( present(extcoef) ) then
             extcoef(:,:,:,w) = extcoef(:,:,:,w) + &
                               tau * (grav * rhoa / delp)
          endif
          if ( present(scacoef) ) then
             scacoef(:,:,:,w) = scacoef(:,:,:,w) + &
                               ssa * tau * (grav * rhoa / delp)
          endif
          !calculate the backscatter coefficient
          if ( present(bckcoef) ) then
             bckcoef(:,:,:,w) = bckcoef(:,:,:,w) + &
                               bck * aerosol(:,:,:,n)*rhoa
          endif
       enddo !wavelengths_profile
      enddo !nbins
    end if !present(extcoef)...

   if( present(exttau) .or. &
       present(stexttau) .or. &
       present(scatau) .or. &
       present(stscatau) .or. &
       present(exttau25) .or. &
       present(exttaufm) .or. &
       present(scatau25) .or. &
       present(scataufm) ) then

      if( present(exttau) ) exttau = 0.
      if( present(stexttau) ) stexttau = 0.
      if( present(scatau) ) scatau = 0.
      if( present(stscatau) ) stscatau = 0.

      if( present(exttau25) ) exttau25 = 0.
      if( present(scatau25) ) scatau25 = 0.

      if( present(exttaufm) ) exttaufm = 0.
      if( present(scataufm) ) scataufm = 0.

      do w = 1, size(wavelengths_vertint)
        do n = nbegin, nbins
           call mie%Query(wavelengths_vertint(w), n,  &
                          aerosol(:,:,:,n)*delp/grav, &
                          rh, tau=tau, ssa=ssa, __RC__)
           do k = klid, km
!             Integrate in the vertical
              if( present(exttau) ) exttau(:,:,w) = exttau(:,:,w) + tau(:,:,k)
              if( present(stexttau) ) then
                 where (ple(:,:,k) .le. tropp) 
                    stexttau(:,:,w) = stexttau(:,:,w) + tau(:,:,k)
                 elsewhere(ple(:,:,k) .gt. tropp .and. ple(:,:,k-1) .lt. tropp) 
                    stexttau(:,:,w) = stexttau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)
                 endwhere
              endif

              if( present(exttaufm) ) then
                 if( NO3nFlag_ ) then
                    exttaufm(:,:,w) = exttaufm(:,:,w) + tau(:,:,k)
                 else
                    exttaufm(:,:,w) = exttaufm(:,:,w) + tau(:,:,k) * fPMfm(n)
                 end if
              end if

              if( present(exttau25) ) then
                 if( NO3nFlag_ ) then
                    exttau25(:,:,w) = exttau25(:,:,w) + tau(:,:,k)
                 else
                    exttau25(:,:,w) = exttau25(:,:,w) + tau(:,:,k) * fPM25(n)
                 end if
              end if

              if( present(scatau) ) scatau(:,:,w) = scatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
              if( present(stscatau) ) then
                 where (ple(:,:,k) .le. tropp) 
                    stscatau(:,:,w) = stscatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
                 elsewhere(ple(:,:,k) .gt. tropp .and. ple(:,:,k-1) .lt. tropp) 
                    stscatau(:,:,w) = stscatau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)*ssa(:,:,k)
                 endwhere
              endif
              if( present(scataufm) ) then
                 if( NO3nFlag_ ) then
                    scataufm(:,:,w) = scataufm(:,:,w) + tau(:,:,k) * ssa(:,:,k)
                 else
                    scataufm(:,:,w) = scataufm(:,:,w) + tau(:,:,k) * ssa(:,:,k) * fPMfm(n)
                 end if
              end if

              if( present(scatau25) ) then
                 if( NO3nFlag_ ) then
                    scatau25(:,:,w) = scatau25(:,:,w) + tau(:,:,k) * ssa(:,:,k)
                 else
                    scatau25(:,:,w) = scatau25(:,:,w) + tau(:,:,k) * ssa(:,:,k) * fPM25(n)
                 end if
              end if
           enddo !k
        enddo !wavelengths_vertint
      enddo !nbins
   endif !present(exttau)...

!  Calculate the 470-870 Angstrom parameter
   if( present(angstrom) .and. do_angstrom ) then

      angstrom(i1:i2,j1:j2) = 0.
!     Set tau to small number by default
      tau470(i1:i2,j1:j2) = tiny(1.0)
      tau870(i1:i2,j1:j2) = tiny(1.0)

      do n = nbegin, nbins

!       Select the name for species
        call mie%Query(4.70E-7, n,                 &
                       aerosol(:,:,:,n)*delp/grav, &
                       rh, tau=tau, __RC__)
        do k = klid, km
          tau470 = tau470 + tau(:,:,k)
        enddo

        call mie%Query(8.70E-7, n,                 &
                       aerosol(:,:,:,n)*delp/grav, &
                       rh, tau=tau, __RC__)
        do k = klid, km
          tau870 = tau870 + tau(:,:,k)
        enddo
      enddo  ! nbins

      angstrom(i1:i2,j1:j2) = &
        -log(tau470(i1:i2,j1:j2)/tau870(i1:i2,j1:j2)) / &
         log(470./870.)
   endif
   deallocate(tau,ssa)
   __RETURN__(__SUCCESS__)
   end subroutine Aero_Compute_Diags
!====================================================================

!BOP
!
! !IROUTINE:  Aero_Binwise_PM_Fractions - Calculate bin-wise PM fractions
!
! !INTERFACE:
   subroutine Aero_Binwise_PM_Fractions(fPM, rPM, r_low, r_up, nbins)

! !USES:
  implicit NONE

! !INPUT/OUTPUT PARAMETERS:
  real, dimension(:), intent(inout) :: fPM     ! bin-wise PM fraction (r < rPM)

! !INPUT PARAMETERS:
   real,    intent(in)              :: rPM     ! PM radius
   integer, intent(in)              :: nbins   ! number of bins
   real, dimension(:), intent(in)   :: r_low   ! bin radii - low bounds
   real, dimension(:), intent(in)   :: r_up    ! bin radii - upper bounds

! !Local Variables

   integer :: n

   character(len=*), parameter :: myname = 'Aero_Binwise_PM_Fractions'
!EOP
!-------------------------------------------------------------------------
!  Begin...

   do n = 1, nbins
     if(r_up(n) < rPM) then
       fPM(n) = 1.0
     else
       if(r_low(n) < rPM) then
!        Assume dm/dlnr = constant, i.e., dm/dr ~ 1/r
         fPM(n) = log(rPM/r_low(n)) / log(r_up(n)/r_low(n))
       else
         fPM(n) = 0.0
       endif
     endif
   enddo

   end subroutine Aero_Binwise_PM_Fractions

!======================================================================================

!BOP
!
! !IROUTINE:  deepLakesMask
!
! !INTERFACE:
   subroutine deepLakesMask (lons, lats, radToDeg, deep_lakes_mask, rc)

! !USES:
   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real, dimension(:,:), intent(inout) :: deep_lakes_mask

! !INPUT PARAMETERS:
   real, pointer, dimension(:,:), intent(in)   :: lats ! latitude [radians]
   real, pointer, dimension(:,:), intent(in)   :: lons ! longtitude [radians]
   real, intent(in)                       :: radToDeg

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc
!EOP

! !Local Variables
    real                :: dummylon
    integer             :: i, j
!EOP
!-------------------------------------------------------------------------
!  Begin...

   deep_lakes_mask = 1.0
   do j = 1, ubound(lons, 2)
      do i = 1, ubound(lons, 1)
                           dummylon = lons(i,j)*radToDeg
      if( dummylon < 0.0 ) dummylon = dummylon + 360.0
      ! The Great Lakes: lon = [91W,75W], lat = [40.5N, 50N]
      if ((dummylon > 267.0) .and. &
          (dummylon < 285.0) .and. &
          (lats(i,j)*radToDeg >  40.5) .and. &
          (lats(i,j)*radToDeg <  50.0)) deep_lakes_mask(i,j) = 0.0

       ! The Caspian Sea: lon = [45.0, 56], lat = 35, 48]
      if ((dummylon >  45.0) .and. &
          (dummylon <  56.0) .and. &
          (lats(i,j)*radToDeg >  35.0) .and. &
          (lats(i,j)*radToDeg <  48.0)) deep_lakes_mask(i,j) = 0.0
      end do
   end do

   __RETURN__(__SUCCESS__)
   end subroutine deepLakesMask

!========================================================================================
!BOP
!
! !IROUTINE:  jeagleSSTcorrection - Apply SST correction following Jaegle et al. 2011
!
! !INTERFACE:
   subroutine jeagleSSTcorrection(sstEmisFlag, fsstemis, ts, rc)

! !USES:
  implicit NONE

! !INPUT/OUTPUT PARAMETERS:
  real, dimension(:,:), intent(inout) :: fsstemis     !

! !INPUT PARAMETERS:
   integer, intent(in)                       :: sstEmisFlag  ! 1 or 2
   real, dimension(:,:), intent(in)          :: ts  ! surface temperature (K)

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc
!EOP

! !Local Variables
   real, allocatable, dimension(:,:) :: tskin_c
!EOP
!-------------------------------------------------------------------------
!  Begin...

   fsstemis = 1.0

   if (sstemisFlag == 1) then          ! SST correction folowing Jaegle et al. 2011
      fsstemis = 0.0

      allocate( tskin_c, mold=fsstemis )
      tskin_c  = ts - 273.15
      fsstemis = (0.3 + 0.1*tskin_c - 0.0076*tskin_c**2 + 0.00021*tskin_c**3)

      where(fsstemis < 0.0) fsstemis = 0.0

      deallocate( tskin_c )
   else if (sstemisFlag == 2) then     ! GEOS5 SST correction
      fsstemis = 0.0

      allocate( tskin_c, mold=fsstemis )
      tskin_c  = ts - 273.15

      where(tskin_c < -0.1) tskin_c = -0.1    ! temperature range (0, 36) C
      where(tskin_c > 36.0) tskin_c = 36.0    !

      fsstemis = (-1.107211 -0.010681*tskin_c -0.002276*tskin_c**2 + 60.288927*1.0/(40.0 - tskin_c))
      where(fsstemis < 0.0) fsstemis = 0.0
      where(fsstemis > 7.0) fsstemis = 7.0

      deallocate( tskin_c )
   end if

   __RETURN__(__SUCCESS__)
   end subroutine jeagleSSTcorrection

!=====================================================================================
!BOP
!
! !IROUTINE: weibullDistribution - Apply a Weibull distribution to emissions wind speeds
!
! !INTERFACE:
   subroutine weibullDistribution (gweibull, weibullFlag, u10m, v10m, rc)

! !USES:
   implicit NONE

! !INPUT/OUTPUT PARAMETERS:
   real(kind=DP), dimension(:,:), intent(inout)    :: gweibull

! !INPUT PARAMETERS:
   logical, intent(in)                    :: weibullFlag
   real, dimension(:,:), intent(in)       :: u10m
   real, dimension(:,:), intent(in)       :: v10m


! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc

! !Descrption: The Weibull distribution correction ends up being a multiplicative constant
!  (g) times our present source function (see Eq. 12 in Fan & Toon, 2011 and notes for
!  (9/22/11). This constant is derived from the incomplete and complete forms of the gamma
!  function, hence the utilities pasted below.  The Weibull function and shape
!  parameters (k, c) assumed are from Justus 1978.

!EOP

! !Local Variables
   real(kind=DP)                 :: a, c, k, wt, x
   real(kind=DP), dimension(:,:), allocatable :: wm
   integer     :: i, j

   integer :: status

!EOP
!-------------------------------------------------------------------------
!  Begin...

   gweibull = 1.0

   allocate(wm(ubound(u10m, 1),ubound(u10m, 2)))
   wm = sqrt(u10m**2 + v10m**2)   ! mean wind speed
   wt = 4.d0                      ! a threshold (Fan & Toon, 2011)

   if (weibullFlag) then
       gweibull = 0.0

   do j = 1, ubound(u10m, 2)
      do i = 1, ubound(u10m, 1)
         if (wm(i,j) > 0.01) then
            k = 0.94d0 * sqrt(wm(i,j))         ! Weibull shape parameter
            c = wm(i,j) / gamma(1.d0 + 1.d0/k) ! Weibull shape parameter
            x = (wt / c) ** k
            a = 3.41d0 / k + 1.d0
            gweibull(i,j)  = (c / wm(i,j))**3.41d0 * igamma(a,x,__RC__)
         end if
      end do ! i
   end do ! j
   endif

   deallocate(wm)

   __RETURN__(__SUCCESS__)
   end subroutine weibullDistribution

!=====================================================================================

 DOUBLE PRECISION function igamma(A, X, rc)
!-----------------------------------------------------------------------
! incomplete (upper) Gamma function
! \int_x^\infty t^{A-1}\exp(-t) dt
!-----------------------------------------------------------------------
 IMPLICIT NONE
 double precision, intent(in) ::        A
 DOUBLE PRECISION, INTENT(IN) ::      X

 integer, intent(out) :: rc
! LOCAL VARIABLE
 DOUBLE PRECISION :: XAM, GIN,  S, R, T0
 INTEGER K
      rc = __SUCCESS__

        XAM=-X+A*LOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'IGAMMA: a and/or x too large, X = ', X
           WRITE(*,*) 'A = ', A
           rc = __FAIL__
           return
        ENDIF

        IF (X.EQ.0.0) THEN
           IGAMMA=GAMMA(A)

        ELSE IF (X.LE.1.0+A) THEN
           S=1.0/A
           R=S
           DO  K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (ABS(R/S).LT.1.0e-15) EXIT
           END DO
           GIN=EXP(XAM)*S
           IGAMMA=GAMMA(A)-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0
           DO K=60,1,-1
              T0=(K-A)/(1.0+K/(X+T0))
           end do

           IGAMMA=EXP(XAM)/(X+T0)

        ENDIF

 end function igamma

!=====================================================================================

! !IROUTINE:  SeasaltEmission - Master driver to compute the sea salt emissions
!
! !INTERFACE:
!
   subroutine SeasaltEmission ( rLow, rUp, method, u10m, v10m, ustar, pi, &
                                memissions, nemissions, rc )

! !DESCRIPTION: Calculates the seasalt mass emission flux every timestep.
!  The particular method (algorithm) used for the calculation is based
!  on the value of "method" passed on input.  Mostly these algorithms are
!  a function of wind speed and particle size (nominally at 80% RH).
!  Routine is called once for each size bin, passing in the edge radii
!  "rLow" and "rUp" (in dry radius, units of um).  Returned in the emission
!  mass flux [kg m-2 s-1].  A sub-bin assumption is made to break (possibly)
!  large size bins into a smaller space.
!
! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real, intent(in)             :: rLow, rUp   ! Dry particle bin edge radii [um]
   real, intent(in)             :: u10m(:,:)   ! 10-meter eastward wind [m s-1]
   real, intent(in)             :: v10m(:,:)   ! 10-m northward wind [m s-1]
   real, target, intent(in)     :: ustar(:,:)  ! friction velocity [m s-1]
   integer, intent(in)          :: method      ! Algorithm to use
   real, intent(in)             :: pi          ! pi constant

! !INOUTPUT PARAMETERS:
   real, dimension(:,:), intent(inout) :: memissions      ! Mass Emissions Flux [kg m-2 s-1]
   real, dimension(:,:), intent(inout) :: nemissions      ! Number Emissions Flux [# m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, intent(out)          :: rc              ! Error return code:
                                                    !  0 - all is well
                                                    !  1 -
! !Local Variables
   integer       :: ir
   real, pointer :: w(:,:)                          ! Intermediary wind speed [m s-1]
   real          :: r, dr                           ! sub-bin radius spacing (dry, um)
   real          :: rwet, drwet                     ! sub-bin radius spacing (rh=80%, um)
   real          :: aFac, bFac, scalefac, rpow, exppow, wpow
   real, allocatable, dimension(:,:), target  :: w10m  ! 10-m wind speed [m s-1]

! !CONSTANTS
   real, parameter    :: r80fac = 1.65     ! ratio of radius(RH=0.8)/radius(RH=0.) [Gerber]
   real, parameter    :: rhop = 2200.      ! dry seasalt density [kg m-3]
!   real, parameter    :: pi = 3.1415       ! ratio of circumference to diameter of circle
   integer, parameter :: nr = 10                    ! Number of (linear) sub-size bins

   character(len=*), parameter :: myname = 'SeasaltEmission'

!EOP
!-------------------------------------------------------------------------
!  Begin...

   rc = __SUCCESS__

!  Define 10-m wind speed
   allocate(w10m, mold=u10m)
   w10m = sqrt(u10m*u10m + v10m*v10m)
!  Define the sub-bins (still in dry radius)
   dr = (rUp - rLow)/nr
   r  = rLow + 0.5*dr

!  Loop over size bins
   nemissions = 0.
   memissions = 0.

   do ir = 1, nr

    rwet  = r80fac * r
    drwet = r80fac * dr

    select case(method)

     case(1)  ! Gong 2003
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 1.
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41
      w        => w10m

     case(2)  ! Gong 1997
      aFac     = 3.
      bFac     = (0.380-log10(rwet))/0.650
      scalefac = 1.
      rpow     = 1.05
      exppow   = 1.19
      wpow     = 3.41
      w        => w10m

     case(3)  ! GEOS5 2012
      aFac     = 4.7*(1.+30.*rwet)**(-0.017*rwet**(-1.44))
      bFac     = (0.433-log10(rwet))/0.433
      scalefac = 33.0e3
      rpow     = 3.45
      exppow   = 1.607
      wpow     = 3.41 - 1.
      w        => ustar

     case default
      !$omp critical (G2G_proc_4)
      print *, 'GOCART2G_Process.F90 - SeasaltEmission - missing algorithm method'
      !$omp end critical (G2G_proc_4)
      rc = __FAIL__
      return

    end select

!   Number emissions flux (# m-2 s-1)
    nemissions = nemissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )

!   Mass emissions flux (kg m-2 s-1)
    scalefac = scalefac * 4./3.*pi*rhop*r**3.*1.e-18
    memissions = memissions + SeasaltEmissionGong( rwet, drwet, w, scalefac, aFac, bFac, rpow, exppow, wpow )

    r = r + dr

   end do

   deallocate(w10m)
  end subroutine SeasaltEmission


! Function to compute sea salt emissions following the Gong style
! parameterization.  Functional form is from Gong 2003:
!  dN/dr = scalefac * 1.373 * (w^wpow) * (r^-aFac) * (1+0.057*r^rpow) * 10^(exppow*exp(-bFac^2))
! where r is the particle radius at 80% RH, dr is the size bin width at 80% RH, and w is the wind speed

  function SeasaltEmissionGong ( r, dr, w, scalefac, aFac, bFac, rpow, exppow, wpow )

   real, intent(in)    :: r, dr     ! Wet particle radius, bin width [um]
   real, pointer, intent(in)    :: w(:,:)    ! Grid box mean wind speed [m s-1] (10-m or ustar wind)
   real, intent(in)    :: scalefac, aFac, bFac, rpow, exppow, wpow
   real                :: SeasaltEmissionGong(size(w,1),size(w,2))

!  Initialize
   SeasaltEmissionGong = 0.

!  Particle size distribution function
   SeasaltEmissionGong = scalefac * 1.373*r**(-aFac)*(1.+0.057*r**rpow) &
                         *10**(exppow*exp(-bFac**2.))*dr
!  Apply wind speed function
   SeasaltEmissionGong = w**wpow * SeasaltEmissionGong

  end function SeasaltEmissionGong

!============================================================================================

!BOP
!
! !IROUTINE: wetRadius - Compute the wet radius of sea salt particle
!
! !INTERFACE:
   subroutine wetRadius (radius, rhop, rh, flag, radius_wet, rhop_wet, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)  :: radius    ! dry radius [m]
   real, intent(in)  :: rhop      ! dry density [kg m-3]
   real, intent(in)  :: rh        ! relative humidity [0-1]
   integer           :: flag      ! 1 (Fitzgerald, 1975)
                                  ! 2 (Gerber, 1985)

! !OUTPUT PARAMETERS:
   real, intent(out) :: radius_wet ! humidified radius [m]
   real, intent(out) :: rhop_wet   ! wet density [kg m-3]
   integer, intent(out) :: rc

! !Local Variables
   real :: sat, rcm, rrat
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

!  The following parameters relate to the swelling of seasalt like particles
!  following Fitzgerald, Journal of Applied Meteorology, 1975.
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2

!  parameter from Gerber 1985 (units require radius in cm, see rcm)
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424

!EOP
!------------------------------------------------------------------------------------
!  Begin...

   rc = __SUCCESS__

!  Default is to return radius as radius_wet, rhop as rhop_wet
   radius_wet = radius
   rhop_wet   = rhop

!  Make sure saturation ratio (RH) is sensible
   sat = max(rh,tiny(1.0)) ! to avoid zero FPE

!  Fitzgerald Scheme
   if(flag .eq. 1 .and. sat .ge. 0.80) then
!     parameterization blows up for RH > 0.995, so set that as max
!     rh needs to be scaled 0 - 1
      sat = min(0.995,sat)
!     Calculate the alpha and beta parameters for the wet particle
!     relative to amonium sulfate
      beta = exp( (0.00077*sat) / (1.009-sat) )
      if(sat .le. 0.97) then
         theta = 1.058
      else
         theta = 1.058 - (0.0155*(sat-0.97)) /(1.02-sat**1.4)
      endif
      alpha1 = 1.2*exp( (0.066*sat) / (theta-sat) )
      f1 = 10.2 - 23.7*sat + 14.5*sat**2.
      f2 = -6.7 + 15.5*sat - 9.2*sat**2.
      alpharat = 1. - f1*(1.-epsilon) - f2*(1.-epsilon**2.)
      alpha = alphaNaCl * (alpha1*alpharat)
!     radius_wet is the radius of the wet particle
      radius_wet = alpha * radius**beta
      rrat       = (radius/radius_wet)**3.
      rhop_wet   = rrat*rhop + (1.-rrat)*rhow
   elseif(flag .eq. 2) then   ! Gerber
      sat = min(0.995,sat)
      rcm = radius*100.
      radius_wet = 0.01 * (c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                           + rcm**3.)**(1./3.)
      rrat       = (radius/radius_wet)**3.
      rhop_wet   = rrat*rhop + (1.-rrat)*rhow
   endif

 end subroutine wetRadius

!===============================================================================

!BOP
!
! !IROUTINE: hoppelCorrection
!
! !INTERFACE:
   subroutine hoppelCorrection (radius, rhop, rh, dz, ustar, rhFlag, &
                                airdens, t, grav, karman, fhoppel, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)     :: radius    ! dry radius [m]
   real, intent(in)     :: rhop      ! dry density [kg m-3]
   integer, intent(in)  :: rhFlag    ! 1 (Fitzgerald, 1975)
                                     ! 2 (Gerber, 1985)
   real, dimension(:,:), intent(in)  :: rh    ! relative humidity [0-1]
   real, dimension(:,:), intent(in)  :: dz    ! surface layer height [m]
   real, dimension(:,:), intent(in)  :: ustar ! surface velocity scale [m s-1]
   real, dimension(:,:), intent(in)  :: airdens ! air density [kg/m^3]s
   real, dimension(:,:), intent(in)  :: t  ! temperature [k]
   real, intent(in)  :: grav    ! gravity [m/sec^2]
   real, intent(in)  :: karman  ! Von Karman constant [unitless]


! !INOUTPUT PARAMETERS:
   real, dimension(:,:), intent(inout) :: fhoppel

! !OUTPUT PARAMETERS:
   integer, intent(out) :: rc

! !Local Variables
   real    :: radius_wet ! humidified radius [m]
   real    :: rhop_wet   ! wet density [kg m-3]
   real    :: diff_coef
   real, allocatable, dimension(:,:) ::  vsettle
   integer :: i, j

   integer :: status


!EOP
!------------------------------------------------------------------------------------
!  Begin..

   rc = __SUCCESS__
   fhoppel = 1.0
   allocate(vsettle, mold=rh)
   vsettle=0.0

   do j = 1, ubound(rh,2)
      do i = 1, ubound(rh,1)
         call wetRadius (radius, rhop, rh(i,j), rhFlag, &
                         radius_wet, rhop_wet, __RC_NO_OPT__)
         call Chem_CalcVsettle2Gorig (radius_wet, rhop_wet, airdens(i,j), t(i,j), &
                                      GRAV, diff_coef, vsettle(i,j))
         fhoppel(i,j) = (10./dz(i,j)) ** (vsettle(i,j)/KARMAN/ustar(i,j))
      end do
   end do


   deallocate(vsettle)

   end subroutine hoppelCorrection

!===============================================================================
!BOP
!
! !IROUTINE:  CAEmission - Adds Carbonaceous Aerosol emission for one timestep
!             We have emissions from 6 sources, which are distributed
!             differently in the vertical
!             1) biomass burning - uniformly mixed in PBL
!             2) biofuel sources - emitted into lowest 100 m
!             3) anthropogenic l1 - emitted into lowest 100 m
!             4) anthropogenic l2 - emitted into 100 - 500 m levels
!             5) terpene          - emitted to surface (hydrophilic only)
!             6) point sources    - emitted in altitudes specified in input
!
! !INTERFACE:
!

   subroutine CAEmission (mie, km, nbins, cdt, grav, prefix, ratPOM, eAircraftfuel, aircraft_fuel_src, &
                           aviation_lto_src, aviation_cds_src, aviation_crs_src, &
                           fHydrophobic, pblh, tmpu, rhoa, rh, aerosolPhilic, aerosolPhobic, &
                           delp, aviation_layers, &
                            biomass_src, terpene_src, eocant1_src, eocant2_src, oc_ship_src, biofuel_src, &
                           OC_emis, OC_emisAN, OC_emisBB, OC_emisBF, OC_emisBG, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   type(GOCART2G_Mie),  intent(in) :: mie        ! mie table
   integer, intent(in) :: km     ! total model levels
   integer, intent(in) :: nbins  ! number of aerosol size bins
   real, intent(in)    :: cdt    ! chemistry model time-step [sec]
   real, intent(in)    :: grav   ! gravity [m/sec^2]
   character(len=2), intent(in)  :: prefix ! varaible name prefix
   real, intent(in)    :: ratPOM
   real, intent(in)    :: eAircraftFuel ! Aircraft emission factor: go from kg fuel to kg C
   real, dimension(:), intent(in)  :: aviation_layers ! Heights [m] of LTO, CDS and CRS aviation emissions layers
   real, pointer, dimension(:,:), intent(in)    :: pblh  ! PBL height [m]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa  ! air density [kg m-3]
   real, pointer, dimension(:,:,:), intent(in)  :: rh    ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in)  :: delp  ! pressure level thickness [Pa]
   real, dimension(:,:,:), intent(in) :: aircraft_fuel_src ! aircraft fuel source [1]
   real, dimension(:,:), intent(in) :: aviation_cds_src ! Climb/Descent aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_crs_src ! Cruise aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_lto_src ! Landing/Take-off aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: biomass_src
   real, dimension(:,:), intent(in) :: terpene_src
   real, dimension(:,:), intent(in) :: eocant1_src  ! anthropogenic emissions
   real, dimension(:,:), intent(in) :: eocant2_src  ! anthropogenic emissions
   real, dimension(:,:), intent(in) :: oc_ship_src  ! ship emissions
   real, dimension(:,:), intent(in) :: biofuel_src  ! biofuel emissions
   real, intent(in) :: fHydrophobic

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: aerosolPhobic
   real, dimension(:,:,:), intent(inout) :: aerosolPhilic
   real, pointer, dimension(:,:,:)  :: OC_emis  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisAN  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisBB  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisBF  ! OC emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: OC_emisBG  ! OC emissions, kg/m2/s
   integer, optional, intent(out) :: rc         ! Error return code:
                                                !  0 - all is well
                                                !  1 -
   character(len=*), parameter :: myname = 'CAEmission'

! !DESCRIPTION: Updates the OC concentration with emissions every timestep
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!    June2020 E.Sherman - moved to process library
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, n, ios, ijl
   integer  ::  n1, n2
   integer  :: i1=1, i2, j1=1, j2
!  pressure at 100m, 500m, & PBLH
   real, dimension(:,:), allocatable :: p100, p500, pPBL
   real, dimension(:,:), allocatable :: p0, z0, ps
   real :: p1, z1, dz, delz, delp_, f100, f500, fPBL, fBot
   real :: qmax, qmin, eBiofuel, eBiomass, eTerpene, eAnthro

   real, dimension(:,:), allocatable :: factor, srcHydrophobic, srcHydrophilic
   real, dimension(:,:), allocatable :: srcBiofuel, srcBiomass, srcAnthro, srcBiogenic
   real                         :: srcTmp, zpbl, maxAll

   real, dimension(:,:,:), allocatable :: emis_aviation
   real, dimension(:,:,:), allocatable :: srcAviation
   real                            :: z_lto_bot, z_lto_top
   real                            :: z_cds_bot, z_cds_top
   real                            :: z_crs_bot, z_crs_top

   real, dimension(:,:), allocatable          :: f_bb_        ! scaling factor for BB emissions based on maximum allowed exttau
   real, dimension(:,:), allocatable          :: exttau_bb_   ! increment of exttau due to BB during the current time step
   real, allocatable, dimension(:,:,:,:) :: qa_bb_       ! increment of qa due to BB during the current time step (nbins,i1:i2,j1:j2:km) 
                                                         ! W.Jiang note, changed to (i1:i2,j1:j2,km,nbins) for efficiency
   real                                  :: cutoff_bb_exttau
   integer                               :: idx
   integer                               :: ilam550, status
   real                                  :: wavelength550
   real, dimension(:,:,:), allocatable   :: tau
   character(len=255)                    :: qname
   real, parameter                       :: max_bb_exttau = 30.0

!  Indices for point emissions
   real, dimension(km)          :: point_column_emissions

!  Source function terms for SOA from Anthropogenic VOCs
   real :: srcSOAanthro = 0.0
!  Initialize local variables
!  --------------------------
   i2 = size(rhoa,1)
   j2 = size(rhoa,2)
   ijl = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

   allocate(factor(i2,j2), srcHydrophobic(i2,j2), srcHydrophilic(i2,j2), srcBiofuel(i2,j2), &
            srcBiomass(i2,j2), srcAnthro(i2,j2), srcBiogenic(i2,j2), f_bb_(i2,j2), exttau_bb_(i2,j2), source=0.0)

!  Emission factors scaling from source files to desired mass quantity
   eBiomass = ratPOM
   eBiofuel = ratPOM
   eTerpene = ratPOM
   eAnthro  = ratPOM

!  Zero diagnostic accumulators
     if(associated(OC_emis)) OC_emis = 0.0
     if(associated(OC_emisAN)) OC_emisAN = 0.0
     if(associated(OC_emisBF)) OC_emisBF = 0.0
     if(associated(OC_emisBB)) OC_emisBB = 0.0
     if(associated(OC_emisBG)) OC_emisBG = 0.0

!  Distribute aircraft emissions from LTO, CDS and CRS layers
!  ----------------------------------------------------------
   z_lto_bot = max(1e-3, aviation_layers(1))
   z_lto_top = max(2e-3, aviation_layers(2))

   z_cds_bot = max(2e-3, aviation_layers(2))
   z_cds_top = max(3e-3, aviation_layers(3))

   z_crs_bot = max(3e-3, aviation_layers(3))
   z_crs_top = max(4e-3, aviation_layers(4))

   allocate(emis_aviation, mold=tmpu)
   allocate(srcAviation, mold=tmpu)
   emis_aviation = 0.0
   srcAviation   = 0.0

   call distribute_aviation_emissions(delp, rhoa, z_lto_bot, z_lto_top, aviation_lto_src, emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_cds_bot, z_cds_top, aviation_cds_src, emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_crs_bot, z_crs_top, aviation_crs_src, emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation


!  Determine surface pressure
!  AMS Note: pass this in
!  --------------------------
   allocate(ps, mold=pblh)
   allocate(p0, mold=pblh)
   allocate(z0, mold=pblh)
   allocate(p100, mold=pblh)
   allocate(p500, mold=pblh)
   allocate(pPBL, mold=pblh)
!AOO initialization
   p0=0.;z0=0.;p100=0.;p500=0.;pPBL=0.
!AOO end initialization
   ps = 0.0
   p0 = 0.0
   z0 = 0.0
   p100 = 0.0
   p500 = 0.0
   pPBL = 0.0
   do k = 1, km
    ps(i1:i2,j1:j2) = ps(i1:i2,j1:j2) + delp(i1:i2,j1:j2,k)
   end do

!  Find the pressure of the 100m, 500m, and PBLH altitudes
!  AMS Note: this could be greatly simplified by using ze/zm and having a
!      generic routine from the bottom up with an early exit condition
!  -----------------------------------------------------------------------
   p0 = ps
   z0(i1:i2,j1:j2) = 0.
   do k = km, 1, -1
    do j = j1, j2
     do i = i1, i2
      p1 = p0(i,j) - delp(i,j,k)
      dz = delp(i,j,k)/rhoa(i,j,k)/grav
      z1 = z0(i,j)+dz
      if(z0(i,j) .lt. 100 .and. z1 .ge. 100.) then
       delz = z1-100.
       delp_ = delz*rhoa(i,j,k)*grav
       p100(i,j) = p1+delp_
      endif
      if(z0(i,j) .lt. 500 .and. z1 .ge. 500.) then
       delz = z1-500.
       delp_ = delz*rhoa(i,j,k)*grav
       p500(i,j) = p1+delp_
      endif
      zpbl = max ( pblh(i,j), 100. )
      if(z0(i,j) .lt. zpbl .and. z1 .ge. zpbl ) then
       delz = z1-zpbl
       delp_ = delz*rhoa(i,j,k)*grav
       pPBL(i,j) = p1+delp_
      endif
      p0(i,j) = p1
      z0(i,j) = z1
     end do
    end do
   end do

!   Limit biomass burning emissions
!   -------------------------------
    allocate(qa_bb_(i1:i2,j1:j2,km,nbins))
    qa_bb_ = 0.0

    p0 = ps
K_LOOP_BB: do k = km, 1, -1

!   First determine emissions for this layer
!   ----------------------------------------
    maxAll = 0.0
    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - delp(i,j,k)

!     Pressure @ PBL height
!     ---------------------
      fPBL = 0.
      if(p1 .ge. pPBL(i,j)) fPBL = delp(i,j,k)/(ps(i,j)-pPBL(i,j))
      if(p1 .lt. pPBL(i,j) .and. p0(i,j) .ge. pPBL(i,j)) fPBL = (p0(i,j)-pPBL(i,j))/(ps(i,j)-pPBL(i,j))

!     Sources by class in kg m-2 s-1
!     ------------------------------
      srcBiomass(i,j)  = fPBL * eBiomass * biomass_src(i,j)

      srcHydrophobic(i,j) =     fHydrophobic  * srcBiomass(i,j)
      srcHydrophilic(i,j) = (1.-fHydrophobic) * srcBiomass(i,j)

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

     end do ! i
    end do  ! j


!   Update concentrations at this layer
!   The "1" element is hydrophobic
!   The "2" element is hydrophilic
!   -----------------------------------
    factor = cdt * grav / delp(:,:,k)

    qa_bb_(:,:,k,1) = factor * srcHydrophobic
    qa_bb_(:,:,k,2) = factor * srcHydrophilic

   end do K_LOOP_BB

!   Get the wavelength indices
!   --------------------------
!   Must provide ilam550 for AOT calculation
    ilam550 = mie%getChannel(5.50e-7)
    if (ilam550 <=0) ilam550 = 1
    wavelength550 = mie%getWavelength(ilam550, __RC__)

!  Calculate the extinction and/or scattering AOD

   exttau_bb_(i1:i2,j1:j2) = 0.0
   allocate(tau(i1:i2,j1:j2,km), source = 0.)
   do n = 1, nbins
!     Select the name for species and the index
     call mie%Query(wavelength550, n,           &
              qa_bb_(:,:,:,n)*delp(:,:,:)/grav, &
              rh, tau=tau, __RC__)
     do k = 1, km
!        Integrate in the vertical
        exttau_bb_(:,:) = exttau_bb_(:,:) + tau(:,:,k)
     enddo
   enddo  ! nbins

   f_bb_ = 1.0
   cutoff_bb_exttau = (cdt / (24 * 3600.0)) * max_bb_exttau

   do j = j1, j2
    do i = i1, i2
     if (exttau_bb_(i,j) > cutoff_bb_exttau) then
      f_bb_(i,j) = cutoff_bb_exttau / exttau_bb_(i,j)
     end if
    enddo
   enddo

   deallocate(qa_bb_, tau)

!  Now update the tracer mixing ratios with the aerosol sources
!  ------------------------------------------------------------
   p0 = ps
K_LOOP: do k = km, 1, -1

!   First determine emissions for this layer
!   ----------------------------------------
    maxAll = 0.0
    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - delp(i,j,k)

!     Pressure @ 100m
!     ---------------
      f100 = 0.
      if(p1 .ge. p100(i,j)) f100 = delp(i,j,k)/(ps(i,j)-p100(i,j))
      if(p1 .lt. p100(i,j) .and. p0(i,j) .ge. p100(i,j)) &
       f100 = (p0(i,j)-p100(i,j))/(ps(i,j)-p100(i,j))

!     Pressure @ 500m
!     ---------------
      f500 = 0.
      if ( p0(i,j) .ge. p100(i,j) .and. p1 .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = (p100(i,j)-p1)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = delp(i,j,k)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .ge. p500(i,j) .and. p1 .lt. p500(i,j)) &
       f500 = (p0(i,j)-p500(i,j))/(p100(i,j)-p500(i,j))

!     Pressure @ PBL height
!     ---------------------
      fPBL = 0.
      if(p1 .ge. pPBL(i,j)) fPBL = delp(i,j,k)/(ps(i,j)-pPBL(i,j))
      if(p1 .lt. pPBL(i,j) .and. p0(i,j) .ge. pPBL(i,j)) &
       fPBL = (p0(i,j)-pPBL(i,j))/(ps(i,j)-pPBL(i,j))

!     Terpene is tree-top emission; only add in bottom layer
!     ------------------------------------------------------
      if ( k .eq. km ) then
         fBot = 1.0
      else
         fBot = 0.0
      end if

!     Sources by class in kg m-2 s-1
!     ------------------------------
      srcBiofuel(i,j)  = f100 * eBiofuel * biofuel_src(i,j)
      srcAnthro(i,j)   = f100 * eAnthro  * eocant1_src(i,j) &
                       + f500 * eAnthro  * eocant2_src(i,j) &
                       + f100 * eAnthro  * oc_ship_src(i,j) &
                       +        eAnthro  * srcAviation(i,j,k) &
                       +        eAnthro  * eAircraftFuel * aircraft_fuel_src(i,j,k)
      if ((prefix == 'OC') .or. (prefix == 'BR')) then
         srcBiomass(i,j)  = fPBL * eBiomass * biomass_src(i,j) * f_bb_(i,j)
      else
         srcBiomass(i,j)  = fPBL * eBiomass * biomass_src(i,j)
      end if

      srcBiogenic(i,j) = fBot * eTerpene * terpene_src(i,j) !Black carbon has no biogenic source. Should be zeros.

      srcTmp = srcBiofuel(i,j) + srcAnthro(i,j) + srcBiomass(i,j)

      srcHydrophobic(i,j) =     fHydrophobic  * srcTmp
      srcHydrophilic(i,j) = (1.-fHydrophobic) * srcTmp + srcBiogenic(i,j)

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

     end do ! i
    end do  ! j

!   Update concentrations at this layer
!   The "1" element is hydrophobic
!   The "2" element is hydrophilic
!   -----------------------------------
    factor = cdt * grav / delp(:,:,k)

    aerosolPhobic(:,:,k) = aerosolPhobic(:,:,k) &
                             + factor * srcHydrophobic

    aerosolPhilic(:,:,k) = aerosolPhilic(:,:,k) &
                             + factor * srcHydrophilic

!   Fill in diagnostics if requested
!   --------------------------------
    if ( associated(OC_emis)) &
                    OC_emis(:,:,1) = OC_emis(:,:,1) + srcHydrophobic

    if ( associated(OC_emis)) &
                    OC_emis(:,:,2) = OC_emis(:,:,2) + srcHydrophilic

    if ( associated(OC_emisBF)) &
                    OC_emisBF  = OC_emisBF  + srcBiofuel

    if ( associated(OC_emisBB)) &
                    OC_emisBB  = OC_emisBB  + srcBiomass

    if ( associated(OC_emisAN)) &
                    OC_emisAN  = OC_emisAN  + srcAnthro

    if ( associated(OC_emisBG)) &
                    OC_emisBG  = OC_emisBG + srcBiogenic
   end do K_LOOP

   __RETURN__(__SUCCESS__)
   end subroutine CAEmission

   subroutine distribute_aviation_emissions(delp, rhoa, z_bot, z_top, emissions_layer, emissions, i1, i2, j1, j2, km, grav)

    implicit none

    integer, intent(in) :: i1, i2, j1, j2, km

    real, dimension(:,:,:), intent(in) :: delp
    real, dimension(:,:,:), intent(in) :: rhoa
    real, dimension(:,:),   intent(in) :: emissions_layer
    real, intent(in)                   :: z_bot
    real, intent(in)                   :: z_top
    real, dimension(:,:,:), intent(out):: emissions
    real, intent(in)                   :: grav
!   local
    integer :: i, j, k
    integer :: k_bot, k_top
    real    :: z_
    real, dimension(km) :: z, dz, w_

    do j = j1, j2
        do i = i1, i2
            ! find level height
            z = 0.0
            z_= 0.0

            do k = km, 1, -1
                dz(k) = delp(i,j,k)/rhoa(i,j,k)/grav
                z_    = z_ + dz(k)
                z(k)  = z_
            end do

            ! find the bottom level
            do k = km, 1, -1
                if (z(k) >= z_bot) then
                    k_bot = k
                    exit
                end if
            end do

            ! find the top level
            do k = k_bot, 1, -1
                if (z(k) >= z_top) then
                    k_top = k
                    exit
                end if
            end do

            ! find the weights
            w_ = 0

!           if (k_top > k_bot) then
!               need to bail - something went wrong here
!           end if

            if (k_bot .eq. k_top) then
                w_(k_bot) = z_top - z_bot
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
            ! distribute emissions in the vertical
            emissions(i,j,:) = (w_ / sum(w_)) * emissions_layer(i,j)
        end do
    end do

    end subroutine distribute_aviation_emissions

!============================================================================

!BOP
!
! !IROUTINE: phobicTophilic
!
! !INTERFACE:
   subroutine phobicTophilic (aerosol_phobic, aerosol_philic, aerosol_toHydrophilic, &
                              km, cdt, grav, delp, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)   :: km   ! total model level
   real, intent(in)      :: cdt  ! chemistry model time-step [sec]
   real, intent(in)      :: grav ! [m/sec^2]
   real, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: aerosol_phobic   ! OCphobic [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: aerosol_philic   ! OCphilic [kg kg-1]
   real, dimension(:,:), pointer   :: aerosol_toHydrophilic ! OCHYPHIL [kg m-2 s-1]
! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc

! !Local Variables
   integer :: i, j, k
   real :: qUpdate, delq

!EOP
!------------------------------------------------------------------------------------
!  Begin...

   if(associated(aerosol_toHydrophilic)) aerosol_toHydrophilic = 0.0

   do k = 1, km
    do j = 1, ubound(delp, 2)
     do i = 1, ubound(delp, 1)
      qUpdate = aerosol_phobic(i,j,k)*exp(-4.63e-6*cdt)
      qUpdate = max(qUpdate,1.e-32)
      delq = max(0.,aerosol_phobic(i,j,k)-qUpdate)
      aerosol_phobic(i,j,k) = qUpdate
      aerosol_philic(i,j,k) = aerosol_philic(i,j,k)+delq
      if(associated(aerosol_toHydrophilic)) &
       aerosol_toHydrophilic(i,j) = aerosol_toHydrophilic(i,j) &
        + delq*delp(i,j,k)/grav/cdt
     end do
    end do
   end do

   __RETURN__(__SUCCESS__)
  end subroutine phobicTophilic


!============================================================================
!BOP
!
! !IROUTINE: NIheterogenousChemOpt
!
! !INTERFACE:
   subroutine NIheterogenousChem (NI_phet, xhno3, UNDEF, AVOGAD, AIRMW, PI, RUNIV, rhoa, tmpu, relhum, delp, &
                                  DU, SS, rmedDU, rmedSS, fnumDU, fnumSS,                                    &
                                  km, klid, cdt, grav, fMassHNO3, fMassNO3, nNO3an1, nNO3an2,                &
                                  nNO3an3, HNO3_conc, HNO3_sfcmass, HNO3_colmass, rc)


! !DESCRIPTION: Nitrogen heterogeneous chemistry - Optimized by A. Darmenov

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)                    :: UNDEF          ! set an undefined value (MAPL_UNDEF)
   real, intent(in)                    :: AVOGAD         ! Avogadro's number [1/kmol]
   real, intent(in)                    :: AIRMW          ! molecular weight of air [kg/kmol]
   real, intent(in)                    :: PI             ! pi constant
   real, intent(in)                    :: RUNIV          ! ideal gas constant [J/(Kmole*K)]
   real, dimension(:,:,:), intent(in)  :: rhoa           ! Layer air density [kg/m^3]
   real, dimension(:,:,:), intent(in)  :: tmpu           ! Layer temperature [K]
   real, dimension(:,:,:), intent(in)  :: relhum         ! relative humidity [1]
   real, dimension(:,:,:), intent(in)  :: delp           ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:,:), intent(in) :: DU   ! dust aerosol [kg/kg]
   real, pointer, dimension(:,:,:,:), intent(in) :: SS   ! sea salt aerosol [kg/kg]
   real, dimension(:) ,intent(in)      :: rmedDU         ! dust aerosol radius [um]
   real, dimension(:) ,intent(in)      :: rmedSS         ! sea salt aerosol radius [um]
   real, dimension(:) ,intent(in)      :: fnumDU         ! number of dust particles per kg mass
   real, dimension(:) ,intent(in)      :: fnumSS         ! number of sea salt particles per kg mass
   integer, intent(in)                 :: km             ! number of model levels
   integer, intent(in)                 :: klid           ! index for pressure lid
   real, intent(in)                    :: cdt            ! chemistry model timestep (sec)
   real, intent(in)                    :: grav           ! gravity (m/sec)
   real, intent(in)                    :: fMassHNO3      ! gram molecular weight
   real, intent(in)                    :: fMassNO3       ! gram molecular weight

! !INOUTPUT PARAMETERS:
   real, pointer, dimension(:,:,:), intent(inout)  :: NI_phet   ! Nitrate Production from Het Chem [kg/(m^2 sec)]
   real, dimension(:,:,:), intent(inout)  :: xhno3     ! buffer for NITRATE_HNO3 [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:), intent(inout)  :: HNO3_conc ! Nitric Acid Mass Concentration [kg/m^3]
   real, pointer, dimension(:,:), intent(inout)    :: HNO3_sfcmass ! Nitric Acid Surface Mass Concentration [kg/m^3]
   real, pointer, dimension(:,:), intent(inout)    :: HNO3_colmass ! Nitric Acid Column Mass Density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(inout)  :: nNO3an1 ! Nitrate bin 1 [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout)  :: nNO3an2 ! Nitrate bin 2 [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout)  :: nNO3an3 ! Nitrate bin 3 [kg/kg]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc

! !Local Variables
   real, dimension(:,:,:), allocatable :: kan1, kan2, kan3, kan
   real, dimension(:,:,:), allocatable :: deltahno3

   integer :: i1, j1, i2, j2, n, i, j, k
   integer :: nbinsDU        ! number of dust bins
   integer :: nbinsSS        ! number of sea salt bins


!
! !REVISION HISTORY:
!
! ???? Optimized NI Het Chem - A. Darmenov
! 15Dec2020 - Ported to process library - E. Sherman

!EOP
!------------------------------------------------------------------------------------
!  Begin..

   nbinsDU = size(DU,4)
   nbinsSS = size(SS,4)

!  Heterogeneous chemistry
!  -----------------------
!  Heterogeneous chemistry wants to know about GOCART dust and sea
!  salt tracers.  This code is not at the moment generalized as it
!  seems very wedded to the traditional GOCART arrangement (5 dust,
!  5 sea salt) and the particulars of the nitrate aerosol arrangement.

   j1 = lbound(tmpu, 2)
   j2 = ubound(tmpu, 2)
   i1 = lbound(tmpu, 1)
   i2 = ubound(tmpu, 1)

   allocate(kan1, mold=tmpu)
   allocate(kan2, mold=tmpu)
   allocate(kan3, mold=tmpu)
   allocate(kan, mold=tmpu)

   kan1 = 0.0
   kan2 = 0.0
   kan3 = 0.0
   kan  = UNDEF

   DUST_HETEROGENOUS_CHEM: if (associated(DU)) then
      DUST_REACTION_RATES: do n = 1, nbinsDU
         kan = 0.0
         call HNO3_reaction_rate(i1, i2, j1, j2, km, klid, &
                                 rmedDU(n), fnumDU(n), &
                                 rhoa, tmpu, relhum, DU(:,:,:,n), kan, &
                                 AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)

         select case(n)
            case (1)
               kan1 = kan1 + kan
            case (2)
               kan2 = kan2 + kan
            case (3)
               kan2 = kan2 + kan
            case (4)
               kan3 = kan3 + kan
            case (5)
               kan3 = kan3 + kan
         end select

      end do DUST_REACTION_RATES
   end if DUST_HETEROGENOUS_CHEM


   SALT_HETEROGENOUS_CHEM: if (associated(SS)) then
      SALT_REACTION_RATES: do n = 1, nbinsSS
         kan = 0.0
         call SSLT_reaction_rate(i1, i2, j1, j2, km, klid, &
                                 rmedSS(n), fnumSS(n), &
                                 rhoa, tmpu, relhum, SS(:,:,:,n), kan, &
                                 AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)

         select case(n)
            case (1)
               kan1 = kan1 + kan
            case (2)
               kan1 = kan1 + kan
            case (3)
               kan2 = kan2 + kan
            case (4)
               kan2 = kan2+ kan
            case (5)
               kan3 = kan3 + kan
         end select

      end do SALT_REACTION_RATES
   end if SALT_HETEROGENOUS_CHEM

!  Compute the nitric acid loss (but don't actually update)
   kan = max(0.0, (kan1 + kan2 + kan3))

   call apportion_reaction_rate(i1, i2, j1, j2, km, kan1, kan)
   call apportion_reaction_rate(i1, i2, j1, j2, km, kan2, kan)
   call apportion_reaction_rate(i1, i2, j1, j2, km, kan3, kan)

   allocate(deltahno3, mold=kan)
   deltahno3 = xhno3 * fMassHNO3 / AIRMW * (1.0 - exp(-kan*cdt))
   deltahno3 = max(0.0, deltahno3)

   xhno3 = xhno3 - deltahno3 * AIRMW / fMassHNO3

   nNO3an1 = nNO3an1 + kan1 * deltahno3 * fMassNO3 / fMassHNO3
   nNO3an2 = nNO3an2 + kan2 * deltahno3 * fMassNO3 / fMassHNO3
   nNO3an3 = nNO3an3 + kan3 * deltahno3 * fMassNO3 / fMassHNO3

   if(associated(NI_phet)) then
      NI_phet(:,:,1) = (1.0 / (grav*cdt)) * sum(kan1*deltahno3*delp, dim=3)
      NI_phet(:,:,2) = (1.0 / (grav*cdt)) * sum(kan2*deltahno3*delp, dim=3)
      NI_phet(:,:,3) = (1.0 / (grav*cdt)) * sum(kan3*deltahno3*delp, dim=3)
   end if

!  Output diagnostic HNO3
!  ----------------------
!  Calculate the HNO3 mass concentration
   if( associated(HNO3_conc) ) then
      HNO3_conc = xhno3 * fMassHNO3 / AIRMW * rhoa
   endif
!  Calculate the HNO3 surface mass concentration
   if( associated(HNO3_sfcmass) ) then
      HNO3_sfcmass(i1:i2,j1:j2) = xhno3(i1:i2,j1:j2,km) * fMassHNO3 / AIRMW * rhoa(i1:i2,j1:j2,km)
   endif
!  Calculate the HNO3 column loading
   if( associated(HNO3_colmass) ) then
      HNO3_colmass(i1:i2,j1:j2) = 0.
      do k = klid, km
        HNO3_colmass(i1:i2,j1:j2) &
         =   HNO3_colmass(i1:i2,j1:j2) + xhno3(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      end do
   endif

   __RETURN__(__SUCCESS__)
   end subroutine NIheterogenousChem

!============================================================================
!BOP
!
! !IROUTINE: HNO3_reaction_rate
!
! !INTERFACE:
   subroutine HNO3_reaction_rate(i1, i2, j1, j2, km, klid, rmed, fnum, rhoa, temp, rh, q, kan, &
                                 AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)

! !DESCRIPTION:

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)  :: i1, i2, j1, j2        ! grid dimension
   integer, intent(in)                 :: km     ! model levels
   integer, intent(in)                 :: klid   ! index for pressure lid
   real, intent(in)                    :: rmed   ! aerosol radius [um]
   real, intent(in)                    :: fnum   ! number of aerosol particles per kg mass
   real, dimension(:,:,:), intent(in)  :: rhoa   ! Layer air density [kg/m^3]
   real, dimension(:,:,:), intent(in)  :: temp   ! Layer temperature [K]
   real, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
   real, dimension(:,:,:), intent(in)  :: q      ! aerosol
   real, intent(in)                    :: AVOGAD ! Avogadro's number [1/kmol]
   real, intent(in)                    :: AIRMW  ! molecular weight of air [kg/kmol]
   real, intent(in)                    :: PI     ! pi constant
   real, intent(in)                    :: RUNIV  ! ideal gas constant [J/(Kmole*K)]
   real, intent(in)                    :: fMassHNO3      ! gram molecular weight

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(out) :: kan

! !Local Variables
   integer :: i, j, k

   real :: f_sad
   real :: f_ad
   real :: radius
   real :: ad
   real :: sad

!EOP
!------------------------------------------------------------------------------------
!  Begin..

   f_ad = 1.e-6 * AVOGAD / AIRMW  ! air number density # cm-3 per unit air density

   ! surface area per unit air density and unit aerosol mass mixing ratio
   f_sad  = 0.01 * 4 * PI * rmed**2 * fnum

   ! radius in 'cm'
   radius = 100 * rmed

   do k = klid, km
      do j = j1, j2
         do i = i1, i2
            ad   = f_ad  * rhoa(i,j,k)             ! air number density # cm-3
            sad  = f_sad * rhoa(i,j,k) * q(i,j,k)  ! surface area density cm2 cm-3

            kan(i,j,k) = sktrs_hno3(temp(i,j,k), rh(i,j,k), sad, ad, radius, PI, &
                                    RUNIV, fMassHNO3)
         end do
      end do
   end do

   end subroutine HNO3_reaction_rate

!============================================================================
!BOP
!
! !IROUTINE: SSLT_reaction_rate
!
! !INTERFACE:
   subroutine SSLT_reaction_rate(i1, i2, j1, j2, km, klid, rmed, fnum, rhoa, temp, rh, q, kan, &
                                 AVOGAD, AIRMW, PI, RUNIV, fMassHNO3)

! !DESCRIPTION:

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)  :: i1, i2, j1, j2        ! grid dimension
   integer, intent(in)                 :: km     ! model levels
   integer, intent(in)                 :: klid   ! index for pressure lid
   real, intent(in)                    :: rmed   ! aerosol radius [um]
   real, intent(in)                    :: fnum   ! number of aerosol particles per kg mass
   real, dimension(:,:,:), intent(in)  :: rhoa   ! Layer air density [kg/m^3]
   real, dimension(:,:,:), intent(in)  :: temp   ! Layer temperature [K]
   real, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
   real, dimension(:,:,:), intent(in)  :: q      ! aerosol
   real, intent(in)                    :: AVOGAD ! Avogadro's number [1/kmol]
   real, intent(in)                    :: AIRMW  ! molecular weight of air [kg/kmol]
   real, intent(in)                    :: PI     ! pi constant
   real, intent(in)                    :: RUNIV  ! ideal gas constant [J/(Kmole*K)]
   real, intent(in)                    :: fMassHNO3      ! gram molecular weight

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(out) :: kan

! !Local Variables
   integer :: i, j, k

   real :: f_sad
   real :: f_ad
   real :: radius
   real :: ad
   real :: sad

!EOP
!------------------------------------------------------------------------------------
!  Begin..

      f_ad = 1.e-6 * AVOGAD / AIRMW  ! air number density # cm-3 per unit air density

      ! surface area per unit air density and unit aerosol mass mixing ratio
      f_sad  = 0.01 * 4 * PI * rmed**2 * fnum

      ! radius in 'cm'
      radius = 100 * rmed

      do k = klid, km
       do j = j1, j2
         do i = i1, i2
          ad   = f_ad  * rhoa(i,j,k)             ! air number density # cm-3
          sad  = f_sad * rhoa(i,j,k) * q(i,j,k)  ! surface area density cm2 cm-3

          kan(i,j,k) = sktrs_sslt(temp(i,j,k), sad, ad, radius, PI, RUNIV, fMassHNO3)
         end do
       end do
      end do

   end subroutine SSLT_reaction_rate

!============================================================================
!BOP
!
! !IROUTINE: apportion_reaction_rate
!
! !INTERFACE:
   subroutine apportion_reaction_rate (i1, i2, j1, j2, km, kan, kan_total)

! !DESCRIPTION:

! !USES:
   implicit NONE


   integer, intent(in) :: i1, i2, j1, j2, km

   real, dimension(i1:i2,j1:j2,km), intent(inout) :: kan
   real, dimension(i1:i2,j1:j2,km), intent(in)    :: kan_total
!EOP
!------------------------------------------------------------------------------------
!  Begin..

   where (kan_total > tiny(kan_total))
       kan = kan / kan_total
   else where
       kan = 0.0
   end where

   end subroutine apportion_reaction_rate

!============================================================================
!BOP
!
! !IROUTINE: sktrs_hno3
!
! !INTERFACE:
   function sktrs_hno3 ( tk, rh, sad, ad, radA, pi, rgas, fMassHNO3 )

! !DESCRIPTION:
! Below are the series of heterogeneous reactions
! The reactions sktrs_hno3n1, sktrs_hno3n2, and sktrs_hno3n3 are provided
! as given by Huisheng Bian.  As written they depend on knowing the GOCART
! structure and operate per column but the functions themselves are
! repetitive.  I cook up a single sktrs_hno3 function which is called per
! grid box per species with an optional parameter gamma being passed.
! Following is objective:
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface

   implicit none

! !INPUT PARAMETERS:
   real, intent(in) ::  tk   ! temperature [K]
   real, intent(in) ::  rh   ! fractional relative humidity [0 - 1]
   real, intent(in) ::  sad  ! aerosol surface area density [cm2 cm-3]
   real, intent(in) ::  ad   ! air number concentration [# cm-3]
   real, intent(in) ::  radA ! aerosol radius [cm]

   real  :: pi   ! pi constant
   real  :: rgas ! ideal gas constant [J/(K*mol)]
   real  :: fMassHNO3 ! gram molecular weight of HNO3
   real :: sktrs_hno3

! !Local Variables
   real, parameter :: fmassHNO3_hno3 = 63.013

!   REAL,  PARAMETER :: GAMMA_HNO3 = 0.1
   REAL,  PARAMETER :: GAMMA_HNO3 = 1.0e-3
!   REAL,  PARAMETER :: GAMMA_HNO3 = 5.0e-4

   real :: dfkg
   real :: avgvel
   real :: gamma
   real :: f_rh
   real :: sqrt_tk
!   real(kind=DP) :: pi_dp = pi
!   real(kind=DP) :: rgas_dp = rgas

!   real, parameter :: p_dfkg   = sqrt(3.472e-2 + 1.0/fmassHNO3)
!   real, parameter :: p_avgvel = sqrt(8.0 * rgas_dp * 1000.0 / (pi_dp * fmassHNO3))

   real(kind=DP) :: pi_dp
   real(kind=DP) :: rgas_dp

   real :: p_dfkg
   real :: p_avgvel

!EOP
!------------------------------------------------------------------------------------
!  Begin..

   pi_dp = pi
   rgas_dp = rgas
   p_dfkg   = sqrt(3.472e-2 + 1.0/fmassHNO3)
   p_avgvel = sqrt(8.0 * rgas_dp * 1000.0 / (pi_dp * fmassHNO3))

      ! RH factor - Figure 1 in Duncan et al. (2010)
      f_rh = 0.03

      if (rh >= 0.1 .and. rh < 0.3)       then
         f_rh = 0.03 + 0.8  * (rh - 0.1)
      else if (rh >= 0.3 .and. rh < 0.5 ) then
         f_rh = 0.19 + 2.55 * (rh - 0.3)
      else if (rh >= 0.5 .and. rh < 0.6)  then
         f_rh = 0.7  + 3.0  * (rh - 0.5)
      else if (rh >= 0.6 .and. rh < 0.7)  then
         f_rh = 1.0  + 3.0  * (rh - 0.6)
      else if (rh >= 0.7 .and. rh < 0.8)  then
         f_rh = 1.3  + 7.0  * (rh - 0.7)
      else if (rh >= 0.8 )                then
         f_rh = 2.0
      end if

!     Following uptake coefficients of Liu et al.(2007)
      gamma = gamma_hno3 * f_rh

      sqrt_tk = sqrt(tk)

!     calculate gas phase diffusion coefficient (cm2/s)
      dfkg = 9.45e17 / ad * sqrt_tk * p_dfkg

!     calculate mean molecular speed (cm/s)
      avgvel = 100.0 * sqrt_tk * p_avgvel

!     calculate rate coefficient
      sktrs_hno3 = sad / ( 4.0 / (gamma * avgvel) + radA / dfkg )

      END FUNCTION sktrs_hno3

!============================================================================
!BOP
!
! !IROUTINE: sktrs_sslt
!
! !INTERFACE:
   function sktrs_sslt ( tk, sad, ad, radA, pi, rgas, fMassHNO3 )

! !DESCRIPTION:
! Below are the series of heterogeneous reactions
! The reactions sktrs_hno3n1, sktrs_hno3n2, and sktrs_hno3n3 are provided
! as given by Huisheng Bian.  As written they depend on knowing the GOCART
! structure and operate per column but the functions themselves are
! repetitive.  I cook up a single sktrs_hno3 function which is called per
! grid box per species with an optional parameter gamma being passed.
! Following is objective:
! loss rate (k = 1/s) of species on aerosol surfaces
!
! k = sad * [ radA/Dg +4/(vL) ]^(-1)
!
! where
! Dg = gas phase diffusion coefficient (cm2/s)
! L = sticking coefficient (unitless)  = gamma
! v = mean molecular speed (cm/s) = [ 8RT / (pi*M) ]^1/2
!
! radA/Dg = uptake by gas-phase diffusion to the particle surface
! 4/(vL) = uptake by free molecular collisions of gas molecules with the surface

   implicit none

! !INPUT PARAMETERS:
   real  :: tk   ! temperature [K]
   real  :: sad  ! aerosol surface area density [cm2 cm-3]
   real  :: ad   ! air number concentration [# cm-3]
   real  :: radA ! aerosol radius [cm]
   real  :: sktrs_sslt
   real  :: pi   ! pi constant
   real  :: rgas ! ideal gas constant [J/(K*mol)]
   real  :: fMassHNO3 ! gram molecular weight of HNO3
!   real(kind=DP), optional  :: gammaInp ! optional uptake coefficient (e.g., 0.2 for SS, else calculated)

!  Locals
   REAL,  PARAMETER :: GAMMA_SSLT = 0.1e0

   real :: dfkg
   real :: avgvel
   real :: sqrt_tk

   real(kind=DP) :: pi_dp
   real(kind=DP) :: rgas_dp

   real :: p_dfkg
   real :: p_avgvel

!EOP
!------------------------------------------------------------------------------------
!  Begin..
   pi_dp = pi
   rgas_dp = rgas

   p_dfkg   = sqrt(3.472e-2 + 1.0/fmassHNO3)
   p_avgvel = sqrt(8.0 * rgas_dp * 1000.0 / (pi_dp * fmassHNO3))

!  Initialize
   sqrt_tk = sqrt(tk)

!     calculate gas phase diffusion coefficient (cm2/s)
      dfkg = 9.45e17 / ad * sqrt_tk * p_dfkg

!     calculate mean molecular speed (cm/s)
      avgvel = 100.0 * sqrt_tk * p_avgvel

!     calculate rate coefficient
      sktrs_sslt = sad / ( 4.0 / (gamma_sslt * avgvel) + radA / dfkg )

   end function sktrs_sslt

!==================================================================================
!BOP
! !IROUTINE: SulfateDistributeEmissions

  subroutine SulfateDistributeEmissions ( km, nbins, cdt, grav, nymd, nhms, &
                                          fMassSO4, fMassSO2, fSO4ant, eAircraftFuel, &
                                          nSO2, nSO4, &
                                          so2anthro_l1_src, so2anthro_l2_src, &
                                          so2biomass_src, dmso_conc, &
                                          so2ship_src, so4ship_src, &
                                          aircraft_fuel_src, &
                                          so2, so4, &
                                          oro, u10m, v10m, hghte, pblh, &
                                          tmpu, rhoa, delp, nVolc, &
                                          SU_emis, SU_SO4eman, SU_SO2eman, SU_SO2embb, &
!                                          maskString, gridMask, &
                                          aviation_layers,   &
                                          aviation_lto_src, &
                                          aviation_cds_src, &
                                          aviation_crs_src, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km, nbins ! number model layers, and number of species respectively
   real, intent(in)    :: cdt, grav ! model time, and gravity respectively
   integer, intent(in) :: nymd, nhms
   real, intent(in)    :: fMassSO4  ! gram molecular weight of SO4
   real, intent(in)    :: fMassSO2  ! gram molecular weight of SO2
   real, intent(in)    :: fSO4ant ! Fraction of anthropogenic emissions that are SO4
   integer, intent(in) :: nSO2  ! index of SO2 relative to other sulfate tracers
   integer, intent(in) :: nSO4  ! index of SO2 relative to other sulfate tracers
   real, intent(in)    :: eAircraftFuel ! Aircraft emission factor: go from kg fuel to kg SO2
   real, dimension(:,:), intent(in) :: so2anthro_l1_src ! anthropogenic source surface[1]
   real, dimension(:,:), intent(in) :: so2anthro_l2_src ! anthropogenic source [1]
   real, dimension(:,:), intent(in) :: so2biomass_src ! biomass burning source [1]
   real, dimension(:,:), intent(in) :: dmso_conc ! DMS source [1]
   real, dimension(:,:), intent(in) :: so2ship_src ! SO2 ship emissions [1]
   real, dimension(:,:), intent(in) :: so4ship_src ! SO4 ship emissions [1]
   real, dimension(:,:,:), intent(in) :: aircraft_fuel_src ! aircraft fuel source [1]

   real, pointer, dimension(:,:), intent(in)    :: oro   ! orography flag
   real, pointer, dimension(:,:), intent(in)    :: u10m  ! 10-m u-wind component [m s-1]
   real, pointer, dimension(:,:), intent(in)    :: v10m  ! 10-m v-wind component [m s-1]
   real, pointer, dimension(:,:,:), intent(in)  :: hghte ! top of layer geopotential height [m]
   real, pointer, dimension(:,:), intent(in)    :: pblh
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa  ! Layer air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]
   integer, intent(in) :: nVolc     ! number of volcanic emissions
   real, dimension(:), intent(in)  :: aviation_layers ! Heights [m] of LTO, CDS and CRS aviation emissions layers
   real, dimension(:,:), intent(in) :: aviation_cds_src ! Climb/Descent aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_crs_src ! Cruise aircraft fuel emission [1]
   real, dimension(:,:), intent(in) :: aviation_lto_src ! Landing/Take-off aircraft fuel emission [1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: so2, so4 ! Sulfate  internal state varaibles [kg/kg]
   real, pointer, dimension(:,:,:)  :: SU_emis      ! SU emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: SU_SO4eman  ! SO4 anthro emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: SU_SO2eman  ! SO2 anthro emissions, kg/m2/s
   real, pointer, dimension(:,:)  :: SU_SO2embb  ! SO2 bioburn emissions, kg/m2/s

!  OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc    ! Error return code:
                                             !  0 - all is well

! !DESCRIPTION: SulfateDistributeEmissions - Adds sulfate source emission for one timestep
!               We have emissions from 4 sources, which are distributed
!               differently in the vertical
!               1) biomass burning - uniformly mixed in PBL (SO2)
!               2) anthropogenic l1 - emitted into lowest 100 m (SO2,SO4)
!               3) anthropogenic l2 - emitted into 100 - 500 m levels (SO2,SO4)
!               4) volcanic emissions
!               Additionally have a source of DMS from transfer from seawater
!               into lowest model layer
!               Consider factors in conversion: we estimate that 5% of sulfur
!               from anthropogenic sources (by mass) goes directly to SO4.
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco -  Based on Ginoux
!  17July2020, Sherman - Refactored for GOCART2G. Only uses intrinsic Fortran

! !Local Variables
   integer  ::  i, j, k
   integer  :: i1=1, j1=1, i2, j2

   real, dimension(:,:), allocatable :: srcSO2, srcSO4, srcDMS, srcSO4anthro, &
                                        srcSO2anthro, srcSO2bioburn
   real, allocatable, dimension(:,:)    :: hsurf

   real :: p1, z1, dz, deltaz, deltap, f100, f500, fPblh
   real :: zpbl
                          ! pressure at 100m, 500m, & PBLH
   real, dimension(:,:), allocatable :: p100, p500, pPblh, p0, z0, ps

   real, dimension(:,:,:), allocatable :: emis_aviation
   real, dimension(:,:,:), allocatable :: srcAviation
   real  :: z_lto_bot, z_lto_top
   real  :: z_cds_bot, z_cds_top
   real  :: z_crs_bot, z_crs_top

!EOP
!-------------------------------------------------------------------------
!  Begin

   i2 = size(rhoa,1)
   j2 = size(rhoa,2)
   allocate(hsurf(i1:i2,j1:j2))
   hsurf = hghte(i1:i2,j1:j2,km)

   allocate(srcSO2(i2,j2), srcSO4(i2,j2), srcDMS(i2,j2), srcSO4anthro(i2,j2), &
            srcSO2anthro(i2,j2), srcSO2bioburn(i2,j2), source=0.0)

!  Initialize local variables
!  --------------------------
   srcSO2 = 0.0
   srcSO4 = 0.0
   srcDMS = 0.0
!AOO initialization
   srcSO4anthro=0.;srcSO2anthro=0.;srcSO2bioburn=0.
!AOO end initialization

   if ((nVolc <= 0) .and. associated(SU_emis)) SU_emis = 0.0 !SU_emis is usually set to zero in SUvolcanicEmissions.
!                                               !If there are no volcanic emissions, we need to set it to zero here.
   if (associated(SU_SO4eman)) SU_SO4eman = 0.0
   if (associated(SU_SO2eman)) SU_SO2eman = 0.0
   if (associated(SU_SO2embb)) SU_SO2embb = 0.0

!  Distribute aircraft emissions from LTO, CDS and CRS layers
!  ----------------------------------------------------------
   z_lto_bot = max(1e-3, aviation_layers(1))
   z_lto_top = max(2e-3, aviation_layers(2))

   z_cds_bot = max(2e-3, aviation_layers(2))
   z_cds_top = max(3e-3, aviation_layers(3))

   z_crs_bot = max(3e-3, aviation_layers(3))
   z_crs_top = max(4e-3, aviation_layers(4))

   allocate(emis_aviation, mold=tmpu)
   allocate(srcAviation, mold=tmpu)
   emis_aviation = 0.0
   srcAviation   = 0.0

   call distribute_aviation_emissions(delp, rhoa, z_lto_bot, z_lto_top, aviation_lto_src, &
                                      emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_cds_bot, z_cds_top, aviation_cds_src, &
                                      emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

   call distribute_aviation_emissions(delp, rhoa, z_crs_bot, z_crs_top, aviation_crs_src, &
                                      emis_aviation, i1, i2, j1, j2, km, grav)
   srcAviation = srcAviation + emis_aviation

!  Find the pressure of the 100m, 500m, and PBLH altitudes
   allocate(ps, mold=pblh)
   allocate(p0, mold=pblh)
   allocate(z0, mold=pblh)
   allocate(p100, mold=pblh)
   allocate(p500, mold=pblh)
   allocate(pPblh, mold=pblh)
!AOO initialization
   p0=0.;z0=0.;p100=0.;p500=0.;pPblh=0.
!AOO end initialization

   ps = 0.0
   p0 = 0.0
   z0 = 0.0
   p100 = 0.0
   p500 = 0.0
   pPblh = 0.0
   do k = 1, km
    ps(i1:i2,j1:j2) = ps(i1:i2,j1:j2) + delp(i1:i2,j1:j2,k)
   end do
   p0 = ps
   z0(i1:i2,j1:j2) = 0.
   do k = km, 1, -1
    do j = j1, j2
     do i = i1, i2
      p1 = p0(i,j) - delp(i,j,k)
      dz = delp(i,j,k)/rhoa(i,j,k)/grav
      z1 = z0(i,j)+dz
      if(z0(i,j) .lt. 100 .and. z1 .ge. 100.) then
       deltaz = z1-100.
       deltap = deltaz*rhoa(i,j,k)*grav
       p100(i,j) = p1+deltap
      endif
      if(z0(i,j) .lt. 500 .and. z1 .ge. 500.) then
       deltaz = z1-500.
       deltap = deltaz*rhoa(i,j,k)*grav
       p500(i,j) = p1+deltap
      endif
      zpbl = max ( pblh(i,j), 100. )
      if(z0(i,j) .lt. zpbl .and. z1 .ge. zpbl) then
       deltaz = z1-zpbl
       deltap = deltaz*rhoa(i,j,k)*grav
       pPblh(i,j) = p1+deltap
      endif
      p0(i,j) = p1
      z0(i,j) = z1
     end do
    end do
   end do

!  Now update the tracer mixing ratios with the aerosol sources
   p0 = ps
   z0 = hsurf
   do k = km, 1, -1

    do j = j1, j2
     do i = i1, i2

      p1 = p0(i,j) - delp(i,j,k)
      f100 = 0.
      if(p1 .ge. p100(i,j)) f100 = delp(i,j,k)/(ps(i,j)-p100(i,j))
      if(p1 .lt. p100(i,j) .and. p0(i,j) .ge. p100(i,j)) &
       f100 = (p0(i,j)-p100(i,j))/(ps(i,j)-p100(i,j))

      f500 = 0.
      if(p0(i,j) .ge. p100(i,j) .and. p1 .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = (p100(i,j)-p1)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .lt. p100(i,j) .and. p1 .ge. p500(i,j)) &
       f500 = delp(i,j,k)/(p100(i,j)-p500(i,j))
      if(p0(i,j) .ge. p500(i,j) .and. p1 .lt. p500(i,j)) &
       f500 = (p0(i,j)-p500(i,j))/(p100(i,j)-p500(i,j))

      fPblh = 0.
      if(p1 .ge. pPblh(i,j)) fPblh = delp(i,j,k)/(ps(i,j)-pPblh(i,j))
      if(p1 .lt. pPblh(i,j) .and. p0(i,j) .ge. pPblh(i,j)) &
       fPblh = (p0(i,j)-pPblh(i,j))/(ps(i,j)-pPblh(i,j))

!     All source from files specified in kg SO2 m-2 s-1 (unless filename
!     indicates otherwise!).
      srcSO4anthro(i,j) = fSO4ant * fMassSO4/fMassSO2 * &
                (   f100 * so2anthro_l1_src(i,j) &
                  + f500 * so2anthro_l2_src(i,j)  )
      srcSO2anthro(i,j) = (1.-fSO4ant) * &
                (   f100 * so2anthro_l1_src(i,j) &
                  + f500 * so2anthro_l2_src(i,j)  )

      srcSO2bioburn(i,j) = fPblh*so2biomass_src(i,j)

!     Add the ship emissions to anthro
      srcSO2anthro(i,j) = srcSO2anthro(i,j) + f100*so2ship_src(i,j)
      srcSO4anthro(i,j) = srcSO4anthro(i,j) + f100*so4ship_src(i,j)

!     Add the aircraft fuel emissions to anthro SO2
      srcSO2anthro(i,j) = srcSO2anthro(i,j) + &
       eAircraftFuel * aircraft_fuel_src(i,j,k)

      srcSO2anthro(i,j) = srcSO2anthro(i,j) + srcAviation(i,j,k)

      srcSO4(i,j) = srcSO4anthro(i,j)
      srcSO2(i,j) = srcSO2anthro(i,j)+srcSO2bioburn(i,j)

      so2(i,j,k)  =   so2(i,j,k) + srcSO2(i,j)*cdt*grav/delp(i,j,k)
      so4(i,j,k)  =   so4(i,j,k) + srcSO4(i,j)*cdt*grav/delp(i,j,k)

      p0(i,j) = p1

     end do ! i
    end do  ! j

    if (associated(SU_emis)) SU_emis(:,:,nSO2) = SU_emis(:,:,nSO2) + srcSO2
    if (associated(SU_emis)) SU_emis(:,:,nSO4) = SU_emis(:,:,nSO4) + srcSO4
    if (associated(SU_SO4eman)) SU_SO4eman = SU_SO4eman + srcSO4anthro
    if (associated(SU_SO2eman)) SU_SO2eman = SU_SO2eman + srcSO2anthro
    if (associated(SU_SO2embb)) SU_SO2embb = SU_SO2embb + srcSO2bioburn

   end do ! k

   __RETURN__(__SUCCESS__)
  end subroutine SulfateDistributeEmissions

!==================================================================================
!BOP
! !IROUTINE: DMSemission

   subroutine DMSemission (km, cdt, grav, tmpu, u10m, v10m, oro, delp, &
                           fMassDMS, dmso_conc, dms, SU_emis, ndms, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km  ! number model layers, and number of species respectively
   real, intent(in)    :: cdt ! model time step [seconds]
   real, intent(in)    :: grav ! gravity [m sec-1]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu  ! temperature [K]
   real, pointer, dimension(:,:), intent(in)    :: u10m  ! 10-m u-wind component [m s-1]
   real, pointer, dimension(:,:), intent(in)    :: v10m  ! 10-m v-wind component [m s-1]
   real, pointer, dimension(:,:), intent(in)    :: oro   ! orography flag
   real, pointer, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]
   real, dimension(:,:), intent(in) :: dmso_conc ! DMS source [1]
   integer, intent(in) :: ndms      ! index of DMS relative to other sulfate tracers
   real, intent(in)    :: fMassDMS  ! gram molecular weight of DMS


! !INOUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: dms ! dms [kg kg-1]
   real, pointer, dimension(:,:,:), intent(inout)  :: SU_emis   ! SU emissions, kg/m2/s

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc    ! Error return code:
                                             !  0 - all is well

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
! 11Feb2020 E.Sherman - First attempt at refactor
!

! !Local Variables
   integer         :: i, j, k
   integer         :: i1=1, j1=1, i2, j2

   real, dimension(:,:), allocatable :: srcDMS
   real :: sCO2, schmidt, w10m, akw, sst

!EOP
!-------------------------------------------------------------------------
!  Begin

!  Add in the DMS source
!  ---------------------
!  DMS emissions go into the lowest model layer only
!  The transfer of DMS from the ocean surface to the atmosphere is
!  a function of surface temperature and wind speed.
!  For now we use the lowest atmospheric temperature (really want SST)
!  and the 10-m wind speed.
!  This code follows from GOCART with the following notes:
!  :the Schmidt number for CO2 is assumed to be 600
!  :the Schmidt number of DMSo follows Saltzman et al., 1993
!  :the Schmidt number dependence breaks for high SST
!  :following www.knmi.nl/~velthove/TM/input we introduce a maximum
!   temperature of 28 C for the calculation
!  :the w10m dependence is from Liss and Merlivat (1986)
!  All this needs some thorough checking!

    i2 = size(tmpu,1)
    j2 = size(tmpu,2)

    allocate(srcDMS(i2,j2))
    srcDMS = 0.

    k = km
    sCO2 = 600.
    do j = j1, j2
       do i = i1, i2
          sst = tmpu(i,j,k)-273.15
          if(sst .gt. 28.) sst = 28.
!         only valid for ocean and warm enough temperatures
          if( (oro(i,j) /= OCEAN) .or. (sst .lt. -20.)) cycle
            schmidt = 2764.0 - 147.12*sst + 3.726*(sst**2.) - 0.038*(sst**3.)
!           w10m is the 10-m wind speed in m s-1
            w10m = sqrt(u10m(i,j)**2. + v10m(i,j)**2.)
          if(w10m .le. 3.6) then
             akw = 0.17*w10m*((sCO2/schmidt)**0.667)
          else if (w10m .le. 13.) then
             akw = (2.85*w10m - 9.65)*sqrt(sCO2/schmidt)
          else
             akw = (5.90*w10m - 49.3)*sqrt(sCO2/schmidt)
          endif
!         This parameterization has put akw in units cm hr-1 -> goto m s-1
          akw = akw/100./3600.
!         DMSo concentration is nMol/L
!         Want to put the source into units of kg m-2 s-1
          srcDMS(i,j) = akw * (fmassDMS/1000.)*(dmso_conc(i,j)*1.e-9/1.e-3)
          dms(i,j,k) =  dms(i,j,k) + srcDMS(i,j)*cdt*grav/delp(i,j,k)
       end do
    end do

    if( associated(SU_emis )) SU_emis(:,:,ndms) = srcDMS


      __RETURN__(__SUCCESS__)
   end subroutine DMSemission


!==================================================================================
!BOP
! !IROUTINE: SUvolcanicEmissions

   subroutine SUvolcanicEmissions (nVolc, vStart, vEnd, vSO2, vElev, vCloud, iPoint, &
                                   jPoint, nhms, SO2EMVN, SO2EMVE, SO2, nSO2, SU_emis, km, cdt, grav,&
                                   hghte, delp, area, vLat, vLon, rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: nVolc     ! number of emissions
   integer, dimension(:), intent(in) :: vStart ! emission start time [sec]
   integer, dimension(:), intent(in) :: vEnd   ! emission end time [sec]
   real, dimension(:), intent(in)    :: vSO2   ! volcanic emission from file [kg]
   real, dimension(:), intent(in)    :: vCloud ! top elevation of emissions [m]
   integer, dimension(:), intent(in) :: iPoint, jPoint ! sub-domain locations of volcanos
   integer, intent(in) :: nhms ! current model time [sec]
   integer, intent(in) :: nSO2   ! index of SO2 relative to other sulfate tracers
   integer, intent(in) :: km   ! number of model levels
   real, intent(in)    :: cdt  ! model time step [sec]
   real, pointer, dimension(:,:,:) :: hghte     ! top of layer geopotential height [m]
   real, intent(in)    :: grav ! gravity [m sec-1]
!   real, dimension(:,:,:), intent(in) :: airdens ! layer air density [kg/m^3]
   real, dimension(:,:,:), intent(in) :: delp  ! pressure thickness [Pa]
   real, dimension(:,:), intent(in)   :: area  ! area of grid cell [m^2]
   real, dimension(:), intent(in)     :: vLat  ! latitude specified in file [degree]
   real, dimension(:), intent(in)     :: vLon  ! longitude specified in file [degree]
! !INOUT PARAMETERS:
  real, pointer, dimension(:,:), intent(inout) :: SO2EMVN ! non-explosive volcanic emissions [kg m-2 s-1]
  real, pointer, dimension(:,:), intent(inout) :: SO2EMVE ! explosive volcanic emissions [kg m-2 s-1]
  real, pointer, dimension(:,:,:), intent(inout) :: SO2 ! SO2 [kg kg-1]
  real, pointer, dimension(:,:,:), intent(inout) :: SU_emis      ! SU emissions, kg/m2/s
  real, dimension(:), intent(inout) ::  vElev ! bottom elevation of emissions [m]

! !OUTPUT PARAMETERS:
  integer, optional, intent(out)   :: rc    ! Error return code:
                                            !  0 - all is well

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! 22July2020 E.Sherman
!
! !Local Variables
   integer  ::  i, j, it
   real, dimension(:,:,:), allocatable  :: emissions_point
   real :: so2volcano

   real :: hup, hlow, dzvolc, dz, z1, k
   real :: deltaSO2v
   real, dimension(:,:), allocatable :: z0
   real, allocatable, dimension(:,:) :: srcSO2volc
   real, allocatable, dimension(:,:) :: srcSO2volce

!EOP
!-------------------------------------------------------------------------
!  Begin

   if (nVolc > 0) then

   allocate(srcSO2volc, mold=area)
   allocate(srcSO2volce, mold=area)
   srcSO2volc = 0.
   srcSO2volce = 0.

   if (associated(SU_emis)) SU_emis = 0.0
   if (associated(SO2EMVN)) SO2EMVN = 0.
   if (associated(SO2EMVE)) SO2EMVE = 0.

   allocate(z0, mold=area)
   z0 = hghte(:,:,km)
!AOO initialization
   z0=0.
!AOO end initialization

    do it = 1, nVolc
       so2volcano = 0.
       i = iPoint(it)
       j = jPoint(it)

!      Skip this volcano?
       if (i<1 .or. j<1) cycle ! volcano not in sub-domain

!      Check time against time range of eruption
       if(nhms < vStart(it) .or. nhms >= vEnd(it)) cycle

!      Emissions per volcano
       if(area(i,j) > 1.) then
          so2volcano = vSO2(it) / area(i,j)     ! to kg SO2/sec/m2
          so2volcano = max(so2volcano,tiny(so2volcano))
       endif

!        Distribute in the vertical
!        Database provides altitude of top of volcano cone (vElev) and altitude
!        of plume top (vCloud).  If vCloud != vElev then distribute emissions
!        in top 1/3 of column extending from vElev to vCloud (case of explosive
!        eruption), else put emissions in grid cell containing vElev (degassing)
!        --------------------------
         hup  = vCloud(it)
         hlow = vElev(it)
         if (hup .ne. hlow) then
            hlow = hup - (hup-hlow)/3.
         endif

!        Diagnostic - sum of volcanos
!        ----------------------------
         if (hup .eq. hlow) then
            srcSO2volc(i,j) = srcSO2volc(i,j) + so2volcano
         else
            srcSO2volce(i,j) = srcSO2volce(i,j) + so2volcano
         endif

         dzvolc = hup-hlow
         do k = km, 1, -1
            z1 = hghte(i,j,k-1) ! geopotential altitude at gridbox top
            dz = z1-z0(i,j)     ! thickness of gridbox
            deltaSO2v = 0.

!           Volcano is above this level
!           ---------------------------
            if(z1 .lt. hlow) then
               z0(i,j) = z1
               cycle
            end if

!           Volcano is below this level (except at surface)
!           -----------------------------------------------
            if(z0(i,j) .gt. hup .and. k .ne. km) then
               z0(i,j) = z1
               cycle
            end if

!           Volcano is in this level
!           ------------------------
            if( (k .eq. km .and. z0(i,j) .gt. hup) .or. &     ! below surface
                 (z0(i,j) .le. hlow .and. z1 .ge. hup) ) then ! in level
               deltaSO2v = so2volcano

!           Volcano only partly in level                       ! Cell:
!           ----------------------------
            else if (z0(i,j) .lt. hlow .and. z1 .lt. hup) then ! has bottom of cloud
               deltaSO2v = (z1-hlow)/dzvolc*so2volcano

            else if (z0(i,j) .gt. hlow .and. z1 .gt. hup) then ! has top of cloud
               deltaSO2v = (hup-z0(i,j))/dzvolc*so2volcano

            else                                               ! is filled with cloud
               deltaSO2v = dz/dzvolc*so2volcano
            end if

            z0(i,j) = z1
            so2(i,j,k) = so2(i,j,k) + deltaSO2v*cdt*grav/delp(i,j,k)

      end do ! k
   enddo     ! it
  end if ! nVolc > 0

  if (associated(SO2EMVN)) SO2EMVN = SO2EMVN + srcSO2volc
  if (associated(SO2EMVE)) SO2EMVE = SO2EMVE + srcSO2volce
  if (associated(SU_emis)) SU_emis(:,:,nSO2) = SU_emis(:,:,nSO2) + srcSO2volc + srcSO2volce

  __RETURN__(__SUCCESS__)
  end subroutine SUvolcanicEmissions

!==================================================================================
!BOP
! !IROUTINE: SulfateUpdateOxidants

   subroutine SulfateUpdateOxidants (nymd_current, nhms_current, lonRad, latRad, &
                                     rhoa, km, cdt, nymd_last, &
                                     undefval, radToDeg, nAvogadro, pi, airMolWght, &
                                     oh_clim, no3_clim, h2o2_clim, &
                                     xoh, xno3, xh2o2, recycle_h2o2, rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)    :: nymd_current, &   ! current model NYMD
                             nhms_current      ! current model NHMS
   real, dimension(:,:), intent(in)   :: lonRad, latRad ! model grid lon and lat
   real, dimension(:,:,:), intent(in) :: rhoa           ! layer air density [kg/m^3]
   integer, intent(in)    :: km         ! number of model levels
   real, intent(in)       :: cdt        ! chemistry model time-step
   integer, intent(inout) :: nymd_last  ! NYMD of last emission update
   real, intent(in)       :: undefval   ! value for undefined values
   real, intent(in)       :: radToDeg   ! radian to degrees conversion
   real, intent(in)       :: nAvogadro  ! Avogadro's number [molecules per mole of air]
   real, intent(in)       :: pi         ! pi constant
   real, intent(in)       :: airMolWght ! air molecular weight [kg/Kmole]
   real, pointer, dimension(:,:,:) :: oh_clim, &   ! climatological OH
                                      no3_clim, &  ! climatological NO3
                                      h2o2_clim  ! climatological H2O2
   real, dimension(:,:,:), intent(inout) :: xoh, xno3, xh2o2 ! returned oxidant values
   logical, intent(inout) :: recycle_h2o2

! !OUTPUT PARAMETERS:
  integer, optional, intent(out)   :: rc    ! Error return code:
                                            !  0 - all is well

! !DESCRIPTION: Update Oxidant Fields for Sulfate
!               We have 3 oxidant fields (OH, NO3, H2O2) which may come
!               from either a climatological file or from interactive GMI.
!               IF from climatology, update (reset) values from climatology
!               if necessary (e.g., for a new day) and set to current values
!               needed by chemistry.
!               IF from GMI read as is

!
! !REVISION HISTORY:
! ???        ???       - Legacy code
! 23July2020 E.Sherman - ported/refactored for use in process library.
!

! !Local Variables
   integer :: i, j, k, jday
   real    :: qmax, xhour, xhouruse
   real, dimension(:,:), allocatable  :: cossza, sza
   real, dimension(:,:), allocatable  :: tcosz, tday, tnight
   integer :: n, ndystep
   integer :: i1=1, j1=1, i2, j2

!EOP
!-------------------------------------------------------------------------
!  Begin...

    i2 = size(rhoa,1)
    j2 = size(rhoa,2)

    allocate(cossza(i1:i2,j1:j2), sza(i1:i2,j1:j2), tcosz(i1:i2,j1:j2), &
             tday(i1:i2,j1:j2), tnight(i1:i2,j1:j2), source=0.0)

! Update emissions/production if necessary (daily)
!  -----------------------------------------------
!   Oxidant fields
!   The expectation here is that OH is being read in the form
!   volume mixing ratio from a file (so, like GMI would provide).
!   Below, in the scaling by solar zenith angle, we convert from
!   VMR to # cm-3 expected by the chemistry.
    where(1.01*oh_clim(i1:i2,j1:j2,1:km) > undefval) oh_clim(i1:i2,j1:j2,1:km) = 0.
    where(     oh_clim(i1:i2,j1:j2,1:km) < 0       ) oh_clim(i1:i2,j1:j2,1:km) = 0.

    where(1.01*no3_clim(i1:i2,j1:j2,1:km) > undefval) no3_clim(i1:i2,j1:j2,1:km) = 0.
    where(     no3_clim(i1:i2,j1:j2,1:km) < 0       ) no3_clim(i1:i2,j1:j2,1:km) = 0.

    where(1.01*h2o2_clim(i1:i2,j1:j2,1:km) > undefval) h2o2_clim(i1:i2,j1:j2,1:km) = 0.
    where(     h2o2_clim(i1:i2,j1:j2,1:km) < 0       ) h2o2_clim(i1:i2,j1:j2,1:km) = 0.

!   The first time through the reads we will save the h2o2 monthly
!   average in the instantaneous field
!   ---------------------------------
    if (nymd_last == nymd_current) then
       xh2o2 = h2o2_clim
       nymd_last = nymd_current
    end if

!   Find the day number of the year and hour (needed for later doing sza)
!   ----------------------------------
    jday = idaynum(nymd_current)
    xhour = (  real(nhms_current/10000)*3600. &
             + real(mod(nhms_current,10000)/100)*60. &
             + real(mod(nhms_current,100)) &
             ) / 3600.

!   Recycle H2O2 to input on 3 hour boundaries if not coupled to GMI
!   ----------------------------------
    if (recycle_h2o2) then
       xh2o2 = h2o2_clim
       recycle_h2o2 = .false.
    end if

!   If not getting instantaneous values from GMI, update for time of day.
!   ---------------------------------------------------------------------
!   OH
    xoh = oh_clim
    cossza(:,:) = 0.

!   Want to find the sum of the cos(sza) for use in scaling OH diurnal variation
!   tcosz is the sum of cossza over the whole day
!   tday is the time of day spent in light
!   Requires integrating over future times, so cannot use w_c%cosz
    xHourUse = xHour
    ndystep = 86400. / cdt
    tcosz(:,:) = 0.
    tday(:,:) = 0.
    do n = 1, ndystep
       call szangle(jday,xHourUse,lonRad,latRad,PI,radToDeg,sza,cossza, i2, j2)
       tcosz = tcosz + cossza
       xHourUse = xHourUse + cdt/3600.
       if(xHourUse .gt. 24.) xHourUse = xHourUse - 24.
!      Find the daylight portion of the day
       do j = j1, j2
          do i = i1, i2
             if(cossza(i,j) .gt. 0.) tday(i,j) = tday(i,j) + cdt
          end do
       end do
    end do

!   Find the cos(sza) now for use in scaling OH and NO3
    call szangle(jday,xHour,lonRad,latRad,PI,radToDeg,sza,cossza, i2, j2)

    tnight(i1:i2,j1:j2) = (86400.-tday(i1:i2,j1:j2))

    do k = 1, km
       where (tcosz(i1:i2,j1:j2) > 0)
          xoh(i1:i2,j1:j2,k) = oh_clim(i1:i2,j1:j2,k)*(86400./cdt)*cossza(i1:i2,j1:j2) / tcosz(i1:i2,j1:j2)
       elsewhere
          xoh(i1:i2,j1:j2,k) = 0.00
       end where
    end do
    where(xoh(i1:i2,j1:j2,1:km) < 0.00) xoh(i1:i2,j1:j2,1:km) = 0.00

!   To go from volume mixing ratio to # cm-3 (expected in chemistry)
!   include the following line
    xoh = xoh * 1000.*rhoa / airMolWght * nAvogadro * 1.e-6

!   NO3
    xno3 = no3_clim
    cossza(:,:) = 0.
    call szangle(jday,xHour,lonRad,latRad,PI,radToDeg,sza,cossza, i2, j2)

!   If there is daylight then no3 is small (assume zero) and the
!   average is distributed only over the night time portion

    do k=1,km
       where(cossza(i1:i2,j1:j2) > 0 .OR. tnight(i1:i2,j1:j2) < tiny(1.0))
          xno3(i1:i2,j1:j2,k) = 0.00
       elsewhere
          xno3(i1:i2,j1:j2,k) = no3_clim(i1:i2,j1:j2,k) * 86400./ tnight(i1:i2,j1:j2)
       end where
    end do

    __RETURN__(__SUCCESS__)
   end subroutine SulfateUpdateOxidants

!==================================================================================
!BOP
! !IROUTINE: szangle

   subroutine szangle (jday, xhour, lonRad, latRad, PI, radToDeg, sza, cossza, i2, j2)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
  integer, intent(in) :: jday ! day # of the year
  real, intent(in)  :: xhour
  real, dimension(:,:), intent(in) :: lonRad, latRad   ! model grid lon and lat
  real, intent(in) :: PI, radToDeg
  real, dimension(:,:), intent(inout) :: cossza, sza
  integer, intent(in) :: i2, j2 ! size of i and j grid dimensions

! !OUTPUT PARAMETERS:

! !DESCRIPTION: given locations and hour find the sza
!               from legacy GOCART (source?)

!
! !REVISION HISTORY:
! 29July2004 P.Colarco - legacy code
! 23July2020 E.Sherman - ported to process library.

! !Local Variables
   integer :: i, j, i1=1, j1=1
   real :: a0, a1, a2, a3, b1, b2, b3, r, dec
   real :: timloc, ahr, xlon, rlat

   a0 = 0.006918
   a1 = 0.399912
   a2 = 0.006758
   a3 = 0.002697
   b1 = 0.070257
   b2 = 0.000907
   b3 = 0.000148
   r  = 2.*pi*float(jday-1)/365. ! where jday is day # of the year


!EOP
!-------------------------------------------------------------------------
!  Begin

!  dec is the solar declination in radians
   dec = a0 - a1*cos(   r) + b1*sin(   r) &
            - a2*cos(2.*r) + b2*sin(2.*r) &
            - a3*cos(3.*r) + b3*sin(3.*r)

   do j = j1, j2
     do i = i1, i2
!    timloc is the local time in hours
     xlon = lonRad(i,j)*radToDeg
     timloc = xhour + xlon/15.
     if(timloc .lt. 0.)  timloc = timloc+24.
     if(timloc .gt. 24.) timloc = timloc-24.
!    ahr is the hour angle in radians
     ahr = abs(timloc - 12.)*15.*pi/180.
      rlat = latRad(i,j)
      cossza(i,j) =   sin(rlat)*sin(dec) &
                    + cos(rlat)*cos(dec)*cos(ahr)

      cossza(i,j)    = min(max(cossza(i,j),-1.0),1.0) !ALT make sure cos stays between -1.0 and 1.0
      sza(i,j)    = acos(cossza(i,j)) * radToDeg
      if(cossza(i,j) .lt. 0.) cossza(i,j) = 0.
     end do
   end do

   end subroutine szangle

!==================================================================================
!BOP
! !IROUTINE: idaynum

   integer function idaynum (nymd)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
  integer :: nymd

! !OUTPUT PARAMETERS:

! !DESCRIPTION: Given nymd compute the day number of the year.

!
! !REVISION HISTORY:
! 29July2004 P.Colarco - Legacy code
! 23July2020 E.Sherman - moved from SulfateChemDriverMod.F90 for use in process library.

! !Local Variables

   integer :: yyyy, mm, dd, imon, isleapyr
   integer :: ndays(12)

   data ndays /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

   yyyy = nymd / 10000
   mm = mod(nymd,10000) / 100
   dd = mod(nymd,100)

!EOP
!-------------------------------------------------------------------------
!  Begin...

!  Is it a leap year?
   isleapyr = 0
   if(mod(yyyy,4) .eq. 0) then
    isleapyr = 1
    if(mod(yyyy,100) .eq. 0) then
     isleapyr = 0
     if(mod(yyyy,400) .eq. 0) then
      isleapyr = 1
     endif
    endif
   endif

!  What day number
   idaynum = 0
   if(mm .eq. 1) then
    idaynum = dd
   else
    do imon = 1, mm-1
     if(imon .eq. 2 .and. isleapyr .eq. 1) then
      idaynum = idaynum+29
     else
      idaynum = idaynum + ndays(imon)
     endif
    enddo
    idaynum = idaynum + dd
   endif

   return
   end function idaynum

!==================================================================================
!BOP
! !IROUTINE: SU_Wet_Removal

   subroutine SU_Wet_Removal ( km, nbins, klid, cdt, kin, grav, airMolWght, delp, fMassSO4, fMassSO2, &
                               h2o2_int, ple, rhoa, precc, precl, pfllsan, pfilsan, tmpu, &
                               nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, MSA, &
                               fluxout, pSO4_colflux, pSO4wet_colflux, &
                               pso4, pso4wet, rc )


! !USES:
   implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: km, nbins  ! number of model levels and number of species respectively
   integer, intent(in) :: klid  ! index for pressure lid
   real, intent(in)    :: cdt   ! chemisty model timestep
   logical, intent(in) :: KIN   ! true for aerosol
   real, intent(in)    :: grav  ! gravity [m/sec]
   real, intent(in)    :: airMolWght ! air molecular weight [kg]
   real, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]
   real, intent(in) :: fMassSO4, fMassSO2
   real, dimension(:,:,:) :: h2o2_int
   real, pointer, dimension(:,:,:), intent(in) :: ple     ! level edge air pressure
   real, pointer, dimension(:,:,:), intent(in) :: rhoa    ! air density, [kg m-3]
   real, pointer, dimension(:,:), intent(in)   :: precc   ! total convective precip, [mm day-1]
   real, pointer, dimension(:,:), intent(in)   :: precl   ! total large-scale prec,  [mm day-1]
   real, pointer, dimension(:,:,:), intent(in) :: pfllsan ! 3D flux of liquid nonconvective precipitation [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:), intent(in) :: pfilsan ! 3D flux of ice nonconvective precipitation [kg/(m^2 sec)]
   real, pointer, dimension(:,:,:), intent(in) :: tmpu    ! temperature, [K]
   integer, intent(in) :: nDMS, nSO2, nSO4, nMSA !index position of sulfates
   real, dimension(:,:,:), intent(inout) :: DMS ! [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO2 ! [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO4 ! [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: MSA ! [kg/kg]

! !OUTPUT PARAMETERS:
   real, pointer, dimension(:,:,:),intent(inout) :: fluxout
   real, pointer, dimension(:,:),intent(inout)  :: pSO4_colflux
   real, pointer, dimension(:,:),intent(inout)  :: pSO4wet_colflux
   real, pointer, dimension(:,:,:),intent(inout) :: pso4
   real, pointer, dimension(:,:,:),intent(inout) :: pso4wet
   integer, optional, intent(out)   :: rc    ! Error return code:
                                            !  0 - all is well


! !DESCRIPTION: Updates the SU concentration due to chemistry
!  The SU grid component is currently established with 4 different
!  species (bins) following this convection:
!   1) DMS
!   2) SO2
!   3) SO4
!   4) MSA
!  Accordingly we have 4 chemical cycles to follow through, which are
!  sub-subroutines under this one.
!  The chemistry is a function of OH, NO3, and H2O2 concentrations
!  as well as DMS, SO2, SO4, MSA concentrations.  It is also a function
!  of solar zenith angle and temperature.  We pass in temperature.  SZA
!  will be a function of time of day and lat/lon.  For now we simply add
!  this to the grid component before calculating it.  I bet this is
!  somewhere else in the model.

!
! !REVISION HISTORY:
!
!
! !Local Variables
   integer :: status

   integer :: i1=1, j1=1, i2, j2
   integer :: dims(3)

   integer  ::  i, j, k, iit, n, LH, kk, ios
!   real :: pdog(i1:i2,j1:j2,km)      ! air mass factor dp/g [kg m-2]
   real, dimension(:,:,:), allocatable :: pdog           ! air mass factor dp/g [kg m-2]
   real*8 :: Td_ls, Td_cv              ! ls and cv timescales [s]
   real*8 :: pls, pcv, pac             ! ls, cv, tot precip [mm day-1]
   real*8 :: qls(km), qcv(km)          ! ls, cv portion of moisture tendency [kg m-3 s-1]
   real*8 :: qmx, qd, A                ! temporary variables on moisture
   real*8 :: F, B, BT                  ! temporary variables on cloud, freq.
   real*8, allocatable :: fd(:,:)      ! flux across layers [kg m-2]
   real*8, allocatable :: dpfli(:,:,:) !
   real*8, allocatable :: DC(:)        ! scavenge change in mass mixing ratio
!   real :: c_h2o(i1:i2,j1:j2,km), cldliq(i1:i2,j1:j2,km), cldice(i1:i2,j1:j2,km)
   real, dimension(:,:,:), allocatable :: c_h2o, cldliq, cldice
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]

!  Rain parameters (from where?)
   real, parameter :: B0_ls = 1.0e-4
   real, parameter :: F0_ls = 1.0
   real, parameter :: XL_ls = 5.0e-4
   real, parameter :: B0_cv = 1.5e-3
   real, parameter :: F0_cv = 0.3
   real, parameter :: XL_cv = 2.0e-3
   real, parameter :: one = 1.0, zero = 0.0

   integer :: nbins_ = 0 ! nbins needs to be redefined in case MSA is not being computed

!  Conversion of SO2 mmr to SO2 vmr (since H2O2 is carried around like
!  a volume mixing ratio)
   real*8 :: fmr, SO2Soluble
   fMR = airMolWght / fMassSO2

!EOP
!-------------------------------------------------------------------------
!  Begin

   rc = __SUCCESS__

   allocate(c_h2o, mold=rhoa)
   allocate(cldliq, mold=rhoa)
   allocate(cldice, mold=rhoa)
!  Initialize local variables
!  --------------------------
!  c_h2o, cldliq, and cldice are respectively intended to be the
!  water mixing ratio (liquid or vapor?, in or out of cloud?)
!  cloud liquid water mixing ratio
!  cloud ice water mixing ratio
   c_h2o  = (10d0**(-2663.5d0/tmpu(:,:,:) + 12.537d0 ) ) /  &
                   (ple(:,:,0:km-1)+ple(:,:,1:km)) /2d0
   cldliq = 0.d0
   where(tmpu >  248.) cldliq = 1.d-6 * ( ( tmpu - 248.d0) / 20.d0 )
   where(tmpu >= 268.) cldliq = 1.d-6
   cldice = 1.d-6 - cldliq

   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

!  check if doing MSA and define nbins_ accordingly
!   if (associated(MSA)) then
!      nbins_ = nbins
!   else
!      nbins_ = nbins - 1
!   end if

   do n = 1, nbins
    if( associated(fluxout)) fluxout(:,:,n) = 0.0
   end do
   if( associated(pso4wet_colflux)) pso4wet_colflux(i1:i2,j1:j2) = 0.
   if( associated(pso4wet)) pso4wet(i1:i2,j1:j2,1:km) = 0.

!  Allocate the dynamic arrays
   allocate(fd(km,nbins),__STAT__)
   allocate(dc(nbins),__STAT__)
   allocate(dpfli(i1:i2, j1:j2, km),__STAT__)
!AOO initialization
    fd=0.d0;dc=0.d0;dpfli=0.d0
!AOO end initialization

!  Duration of rain: ls = model timestep, cv = 1800 s (<= cdt)
   Td_ls = cdt
   Td_cv = 1800.

!  Accumulate the 3-dimensional arrays of rhoa and pdog
   pdog = delp/grav

   dpfli = pfllsan(:,:,1:km)-pfllsan(:,:,0:km-1)+pfilsan(:,:,1:km)-pfilsan(:,:,0:km-1)

!  Loop over spatial indices
   do j = j1, j2
    do i = i1, i2

!    Check for total precipitation amount
!    Assume no precip in column if precl+precc = 0
     pac = precl(i,j) + precc(i,j)
     if(pac .le. 0.) cycle
     pls = precl(i,j)
     pcv = precc(i,j)

!    Initialize the precipitation fields
     qls(:)  = 0.
     qcv(:)  = 0.
     fd(:,:) = 0.
     Dc(:)   = 0.

!    Find the highest model layer experiencing rainout.  Assumes no
!    scavenging if T < 258 K
     LH = 0
     do k = klid, km
      if(dpfli(i,j,k) .gt. 0. .and. tmpu(i,j,k) .gt. 258.) then
       LH = k
       exit
      endif
     end do

     if(LH .lt. 1) cycle

     do k = LH, km
      qls(k) = dpfli(i,j,k)/pdog(i,j,k)*rhoa(i,j,k)
     end do

!    Loop over vertical to do the scavenging!
     do k = LH, km

!-----------------------------------------------------------------------------
!   (1) LARGE-SCALE RAINOUT:
!       Tracer loss by rainout = TC0 * F * exp(-B*dt)
!         where B = precipitation frequency,
!               F = fraction of grid box covered by precipitating clouds.
!       We assume that tracer scavenged by rain is falling down to the
!       next level, where a fraction could be re-evaporated to gas phase
!       if Qls is less then 0 in that level.
!-----------------------------------------------------------------------------
      if (qls(k) .gt. tiny(0.)) then
       F  = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(qls(k)*cdt/Td_ls))
       B  = B0_ls/F0_ls +1./(F0_ls*XL_ls/qls(k))
       BT = B * Td_ls
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >
!      What is the soluble amount of SO2?
       SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = SO4(i,j,k) * F * (1.-exp(-BT))
       if (associated(MSA)) then
          DC(nMSA) = MSA(i,j,k) * F * (1.-exp(-BT))
       endif

!      Adjust H2O2 concentration in cloudy portion of cell
       if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
!        gcSU%h2o2_int(i,j,k) = max(zero,(1.-F)*gcSU%h2o2_int(i,j,k))
! GOCART removes all
        h2o2_int(i,j,k) = 0.
       else
        h2o2_int(i,j,k) &
          = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
       endif

       do n = 1, nbins
        if (DC(n).lt.0.) DC(n) = 0.
       end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = DC(n) * pdog(i,j,k)
       end do

      end if                                    ! if Qls > 0  >>>

!-----------------------------------------------------------------------------
! * (2) LARGE-SCALE WASHOUT:
! *     Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------
      if(k .gt. LH .and. qls(k) .ge. 0.) then
       if(qls(k) .lt. qls(k-1)) then
!       Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1,LH,-1
         if (Qls(kk).gt.0.) then
          Qmx = max(Qmx,Qls(kk))
         else
          exit
         end if
        end do

        F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(Qmx*cdt/Td_ls))
        if (F.lt.0.01) F = 0.01

!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx /rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!      What is the soluble amount of SO2?
       SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = SO4(i,j,k) * F * (1.-exp(-BT))
       if (associated(MSA)) then
          DC(nMSA) = MSA(i,j,k) * F * (1.-exp(-BT))
       end if

!       Adjust H2O2 concentration in cloudy portion of cell
        if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
         h2o2_int(i,j,k) = max(zero,(one-F)*h2o2_int(i,j,k))
!  GOCART removes all
!         gcSU%h2o2_int(i,j,k) = 0.
        else
         h2o2_int(i,j,k) &
           = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
        endif

        do n = 1, nbins
         if (DC(n).lt.0.) DC(n) = 0.
        end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

        do n = 1, nbins
         if( associated(fluxout) ) then
          fluxout(i,j,n) = fluxout(i,j,n)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if ls washout  >>>

!-----------------------------------------------------------------------------
!  (3) CONVECTIVE RAINOUT:
!      Tracer loss by rainout = dd0 * F * exp(-B*dt)
!        where B = precipitation frequency,
!              F = fraction of grid box covered by precipitating clouds.
!-----------------------------------------------------------------------------
      if (qcv(k) .gt. 0.) then
       F  = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qcv(k)*cdt/Td_cv))
       B  = B0_cv
       BT = B * Td_cv
       if (BT.gt.10.) BT = 10.               !< Avoid overflow >

!      Adjust SO2 for H2O2 oxidation
       SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
       if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!      Adjust SU amounts
       DC(nDMS) = 0.
       DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
       DC(nSO4) = SO4(i,j,k) * F * (1.-exp(-BT))
       if (associated(MSA)) then
          DC(nMSA) = MSA(i,j,k) * F * (1.-exp(-BT))
       end if
       DC(nSO4) = 0.
       if (associated(MSA)) then
          DC(nMSA) = 0.
       end if

!      Adjust H2O2 concentration in cloudy portion of cell
       if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
        h2o2_int(i,j,k) = max(zero,(one-F)*h2o2_int(i,j,k))
       else
        h2o2_int(i,j,k) &
          = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
       endif

       do n = 1, nbins
        if (DC(n).lt.0.) DC(n) = 0.
       end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

!      Flux down:  unit is kg m-2
!      Formulated in terms of production in the layer.  In the revaporation step
!      we consider possibly adding flux from above...
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + DC(n)*pdog(i,j,k)
       end do

      end if                                  ! if Qcv > 0   >>>

!-----------------------------------------------------------------------------
!  (4) CONVECTIVE WASHOUT:
!      Occurs when rain at this level is less than above.
!-----------------------------------------------------------------------------

      if (k.gt.LH .and. Qcv(k).ge.0.) then
       if (Qcv(k).lt.Qcv(k-1)) then
!-----  Find a maximum F overhead until the level where Qls<0.
        Qmx   = 0.
        do kk = k-1, LH, -1
         if (Qcv(kk).gt.0.) then
          Qmx = max(Qmx,Qcv(kk))
         else
          exit
         end if
        end do

        F = F0_cv / (1. + F0_cv*B0_cv*XL_cv/(Qmx*cdt/Td_cv))
        if (F.lt.0.01) F = 0.01
!-----------------------------------------------------------------------------
!  The following is to convert Q(k) from kgH2O/m3/sec to mm/sec in order
!  to use the Harvard formula.  Convert back to mixing ratio by multiplying
!  by rhoa.  Multiply by pdog gives kg/m2/s of precip.  Divide by density
!  of water (=1000 kg/m3) gives m/s of precip and multiply by 1000 gives
!  units of mm/s (omit the multiply and divide by 1000).
!-----------------------------------------------------------------------------

        Qd = Qmx / rhoa(i,j,k)*pdog(i,j,k)
        if (Qd.ge.50.) then
         B = 0.
        else
         B = Qd * 0.1
        end if
        BT = B * cdt
        if (BT.gt.10.) BT = 10.

!       Adjust SO2 for H2O2 oxidation
        SO2Soluble = min(fmr*SO2(i,j,k),h2o2_int(i,j,k)*one)/fmr
        if(SO2Soluble .lt. 0.) SO2Soluble = 0.

!       Adjust SU amounts
        DC(nDMS) = 0.
        DC(nSO2) = SO2Soluble * F * (1.-exp(-BT))
! Sulfate scavenged in moist
!        DC(nSO4) = w_c%qa(n1+nSO4-1)%data3d(i,j,k) * F * (1.-exp(-BT))
!        DC(nMSA) = w_c%qa(n1+nMSA-1)%data3d(i,j,k) * F * (1.-exp(-BT))
        DC(nSO4) = 0.
        if (associated(MSA)) then
           DC(nMSA) = 0.
        end if

!       Adjust H2O2 concentration in cloudy portion of cell
        if(fmr*SO2(i,j,k) .gt. h2o2_int(i,j,k)) then
         h2o2_int(i,j,k) = max(zero,(one-F)*h2o2_int(i,j,k))
        else
         h2o2_int(i,j,k) &
           = h2o2_int(i,j,k) - F*fmr*SO2(i,j,k)
        endif

        do n = 1, nbins
         if (DC(n).lt.0.) DC(n) = 0.
        end do

       call updateAerosol(DMS(i,j,k), DC(nDMS))
       call updateAerosol(SO2(i,j,k), DC(nSO2))
       call updateAerosol(SO4(i,j,k), DC(nSO4))
       if (associated(MSA)) then
          call updateAerosol(MSA(i,j,k), DC(nMSA))
       end if

        do n = 1, nbins
         if( associated(fluxout) ) then
          fluxout(i,j,n) = fluxout(i,j,n)+DC(n)*pdog(i,j,k)/cdt
         endif
        end do

       end if
      end if                                    ! if cv washout  >>>

!-----------------------------------------------------------------------------
!  (5) RE-EVAPORATION.  Assume that SO2 is re-evaporated as SO4 since it
!      has been oxidized by H2O2 at the level above.
!-----------------------------------------------------------------------------
!     Add in the flux from above, which will be subtracted if reevaporation occurs
      if(k .gt. LH) then
       do n = 1, nbins
        Fd(k,n) = Fd(k,n) + Fd(k-1,n)
       end do

!      Is there evaporation in the currect layer?
       if (dpfli(i,j,k) .lt. 0.) then
!       Fraction evaporated = H2O(k)evap / H2O(next condensation level).
        if (dpfli(i,j,k-1) .gt. 0.) then

          A =  abs(  dpfli(i,j,k) /  dpfli(i,j,k-1)  )
          if (A .gt. 1.) A = 1.

!         Adjust tracer in the level
!         For the SO2 tracer we do not allow re-evaporation.
!         We compute DC(nSO2) solely to add this to DC(nSO4) and to remove
!         from Fd(k,nSO2)
!         Instead, the SO2 gets re-evaporated to the SO4 bin because of
!         previous H2O2 oxidation

          DC(nDMS) = 0.
          DC(nSO2) = Fd(k-1,nSO2) / pdog(i,j,k) * A
          DC(nSO4) = Fd(k-1,nSO4) / pdog(i,j,k) * A
          if (associated(MSA)) then
             DC(nMSA) = Fd(k-1,nMSA) / pdog(i,j,k) * A
          end if

          do n = 1, nbins
           if (DC(n).lt.0.) DC(n) = 0.
          end do

          if (associated(MSA)) then
             MSA(i,j,k) = MSA(i,j,k) + DC(nMSA)
          end if
!         SO2 gets added to SO4, but remember to remove the SO2 from FD!
          SO4(i,j,k) = SO4(i,j,k) + DC(nSO4) + DC(nSO2)*fMassSO4/fMassSO2
          if( associated(pso4wet_colflux)) &
             pso4wet_colflux(i,j) = pso4wet_colflux(i,j) &
              + DC(nSO2)*fMassSO4/fMassSO2 / cdt * delp(i,j,k)/grav
          if( associated(pso4wet) ) &
             pso4wet(i,j,k) = DC(nSO2)*fMassSO4/fMassSO2 / cdt

          if( associated(pso4_colflux)) &
             pso4_colflux(i,j) = pso4_colflux(i,j) &
              + DC(nSO2)*fMassSO4/fMassSO2 / cdt * delp(i,j,k)/grav
          if( associated(pso4) ) &
             pso4(i,j,k) = pso4(i,j,k) + DC(nSO2)*fMassSO4/fMassSO2 / cdt


!         Adjust the flux out of the bottom of the layer--remove SO2 here!
          DMS(i,j,k) = max(DMS(i,j,k),tiny(1.0))
          Fd(k,nDMS) = Fd(k,nDMS) - DC(nDMS)*pdog(i,j,k)
          SO2(i,j,k) = max(SO2(i,j,k),tiny(1.0))
          Fd(k,nSO2) = Fd(k,nSO2) - DC(nSO2)*pdog(i,j,k)
          SO4(i,j,k) = max(SO4(i,j,k),tiny(1.0))
          Fd(k,nSO4) = Fd(k,nSO4) - DC(nSO4)*pdog(i,j,k)
          if (associated(MSA)) then
             MSA(i,j,k) = max(MSA(i,j,k),tiny(1.0))
             Fd(k,nMSA) = Fd(k,nMSA) - DC(nMSA)*pdog(i,j,k)
          end if
        endif
       endif                                   ! if -moistq < 0
      endif
     end do  ! k

     do n = 1, nbins
      if( associated(fluxout) ) then
       fluxout(i,j,n) = fluxout(i,j,n)+Fd(km,n)/cdt
      endif
     end do

    end do   ! i
   end do    ! j

   deallocate(fd,DC,dpfli,stat=ios)


   __RETURN__(__SUCCESS__)
   contains
     subroutine updateAerosol (aerosol, DC)

     ! !USES:
     implicit NONE
     ! !INPUT PARAMETERS:
      real, intent(inout) :: aerosol
      real*8, intent(in)  :: DC

        aerosol = aerosol - DC
        if (aerosol .lt. 1.0E-32) aerosol = 1.0E-32

      end subroutine updateAerosol

   end subroutine SU_Wet_Removal

!==================================================================================
!BOP
! !IROUTINE: SU_Compute_Diags

   subroutine SU_Compute_Diags ( km, klid, rmed, sigma, rhop, grav, pi, nSO4, mie, &
                                 wavelengths_profile, wavelengths_vertint, &
                                 tmpu, rhoa, delp, ple, tropp,rh, u, v, &
                                 DMS, SO2, SO4, MSA, &
                                 dmssfcmass, dmscolmass, &
                                 msasfcmass, msacolmass, &
                                 so2sfcmass, so2colmass, &
                                 so4sfcmass, so4colmass, &
                                 exttau, stexttau,scatau, stscatau,so4mass, so4conc, extcoef, &
                                 scacoef, bckcoef, angstrom, fluxu, fluxv, sarea, snum, rc )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km    ! number of model levels
   integer,    intent(in)    :: klid   ! index for pressure lid
   real, intent(in)    :: rmed  ! mean radius [um]
   real, intent(in)    :: sigma ! Sigma of lognormal number distribution
   real, intent(in)    :: rhop  ! dry particle density [kg m-3]
   real, intent(in)    :: grav  ! gravity [m/sec]
   real, intent(in)    :: pi    ! pi constant
   integer, intent(in) :: nSO4  ! index of SO4 relative to other internal variables
   type(GOCART2G_Mie), intent(in) :: mie   ! mie table
   real, dimension(:), intent(in)  :: wavelengths_profile
   real, dimension(:), intent(in)  :: wavelengths_vertint
   real, pointer, dimension(:,:,:), intent(in) :: tmpu    ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa    ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: delp    ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: ple   ! level edge air pressure [Pa]
   real, pointer, dimension(:,:), intent(in)   :: tropp ! tropopause pressure [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: rh      ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in) :: u       ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:), intent(in) :: v       ! north-south wind [m s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: DMS  ! dimethyl sulfide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO2  ! sulfer dioxide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO4  ! sulfate aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout)  :: MSA  ! methanesulphonic acid [kg/kg]
   real, optional, dimension(:,:),   intent(inout)  :: dmssfcmass ! sfc mass concentration [kg/m3]
   real, optional, dimension(:,:),   intent(inout)  :: dmscolmass ! col mass density [kg/m2]
   real, optional, dimension(:,:),   intent(inout)  :: msasfcmass ! sfc mass concentration [kg/m3]
   real, optional, dimension(:,:),   intent(inout)  :: msacolmass ! col mass density [kg/m2]
   real, optional, dimension(:,:),   intent(inout)  :: so2sfcmass ! sfc mass concentration [kg/m3]
   real, optional, dimension(:,:),   intent(inout)  :: so2colmass ! col mass density [kg/m2]
   real, optional, dimension(:,:),   intent(inout)  :: so4sfcmass ! sfc mass concentration [kg/m3]
   real, optional, dimension(:,:),   intent(inout)  :: so4colmass ! col mass density [kg/m2]
   real, optional, dimension(:,:,:), intent(inout)  :: exttau     ! ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)  :: stexttau   ! Stratosphere ext. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)  :: scatau     ! sct. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)  :: stscatau   ! Stratosphere sct. AOT at 550 nm
   real, optional, dimension(:,:,:), intent(inout)  :: so4mass    ! 3D sulfate mass mr
   real, optional, dimension(:,:,:), intent(inout)  :: so4conc    ! 3D mass concentration, [kg/m3]
   real, optional, dimension(:,:,:,:), intent(inout)  :: extcoef    ! 3D ext. coefficient, [1/m]
   real, optional, dimension(:,:,:,:), intent(inout)  :: scacoef    ! 3D scat.coefficient, [1/m]
   real, optional, dimension(:,:,:,:), intent(inout)  :: bckcoef    ! 3D backscatter coefficient, [m-1 sr-1]
   real, optional, dimension(:,:),   intent(inout)  :: angstrom   ! 470-870 nm Angstrom parameter
   real, optional, dimension(:,:),   intent(inout)  :: fluxu      ! Column mass flux in x direction
   real, optional, dimension(:,:),   intent(inout)  :: fluxv      ! Column mass flux in y direction
   real, optional, dimension(:,:,:), intent(inout)  :: sarea      ! Sulfate surface area density [m2 m-3]
   real, optional, dimension(:,:,:), intent(inout)  :: snum       ! Sulfate number density [# m-2]
   integer, optional, intent(out)   :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -


! !DESCRIPTION: Calculates some simple 2d diagnostics from the SU fields
!
! !REVISION HISTORY:
!
!  16APR2004, Colarco
!  29july2020, E.Sherman - refactored for process library

! !Local Variables
   integer :: i, j, k, w, i1=1, j1=1, i2, j2, status
   real, dimension(:,:,:), allocatable :: tau, ssa, bck
   real, dimension(:,:), allocatable :: tau470, tau870
   integer    :: ilam470, ilam870
   logical :: do_angstrom
   real :: rh_, gf, rwet, svol


!EOP
!-------------------------------------------------------------------------
!  Begin
   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

   allocate(tau470(i1:i2,j1:j2), tau870(i1:i2,j1:j2), source=0.0)

!  Get the wavelength indices
!  --------------------------

   ilam470 = mie%getChannel(4.70e-7)
   if(ilam470 <= 0) ilam470 = 0

   ilam870 = mie%getChannel(8.70e-7)
   if(ilam870 <= 0) ilam870 = 0

!  Determine if going to do Angstrom parameter calculation
!  -------------------------------------------------------
   do_angstrom = .false.
!  If both 470 and 870 channels provided (and not the same) then
!  possibly will do Angstrom parameter calculation
   if(ilam470 .ne. 0 .and. &
      ilam870 .ne. 0 .and. &
      ilam470 .ne. ilam870) do_angstrom = .true.


!  Calculate the diagnostic variables if requested
!  -----------------------------------------------

!  Calculate the surface mass concentration
   if( present(so4sfcmass) ) then
      so4sfcmass(i1:i2,j1:j2) = 0.
      so4sfcmass(i1:i2,j1:j2) &
       =  SO4(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( present(so2sfcmass) ) then
      so2sfcmass(i1:i2,j1:j2) = 0.
      so2sfcmass(i1:i2,j1:j2) &
       =   SO2(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( present(dmssfcmass) ) then
      dmssfcmass(i1:i2,j1:j2) = 0.
      dmssfcmass(i1:i2,j1:j2) &
       =   DMS(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( present(msasfcmass) .and. associated(MSA)) then
      msasfcmass(i1:i2,j1:j2) = 0.
      msasfcmass(i1:i2,j1:j2) &
       =   MSA(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif


!  Initialize the diagnostic variables
!  -----------------------------------

!  Calculate the column loading
   if( present(so4colmass) ) then
      so4colmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       so4colmass(i1:i2,j1:j2) &
        =   so4colmass(i1:i2,j1:j2) &
          + SO4(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( present(so2colmass) ) then
      so2colmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       so2colmass(i1:i2,j1:j2) &
        =   so2colmass(i1:i2,j1:j2) &
          + SO2(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( present(dmscolmass) ) then
      dmscolmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       dmscolmass(i1:i2,j1:j2) &
        =   dmscolmass(i1:i2,j1:j2) &
          + DMS(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( present(msacolmass) .and. associated(MSA)) then
      msacolmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       msacolmass(i1:i2,j1:j2) &
        =   msacolmass(i1:i2,j1:j2) &
          + MSA(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif


!  Calculate the mass concentration of sulfate
   if( present(so4conc) ) then
      so4conc(i1:i2,j1:j2,1:km) = 0.
      so4conc(i1:i2,j1:j2,1:km) = SO4(i1:i2,j1:j2,1:km)*rhoa(i1:i2,j1:j2,1:km)
   endif

!  Mass mixing ratio of sulfate
   if( present(so4mass) ) then
      so4mass(i1:i2,j1:j2,1:km) = 0.
      so4mass(i1:i2,j1:j2,1:km) = SO4(i1:i2,j1:j2,1:km)
   endif

!  Calculate the column mass flux in x direction
   if( present(fluxu) ) then
      fluxu(i1:i2,j1:j2) = 0.
       do k = klid, km
        fluxu(i1:i2,j1:j2) &
         =   fluxu(i1:i2,j1:j2) &
           + SO4(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
       end do
   endif

!  Calculate the column mass flux in y direction
   if( present(fluxv) ) then
      fluxv(i1:i2,j1:j2) = 0.
       do k = klid, km
        fluxv(i1:i2,j1:j2) &
         =   fluxv(i1:i2,j1:j2) &
           + SO4(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
       end do
   endif

!  Calculate the extinction and/or scattering AOD
   allocate(tau(i1:i2,j1:j2,km), source = 0.)
   allocate(ssa(i1:i2,j1:j2,km), source = 0.)
   allocate(bck(i1:i2,j1:j2,km), source = 0.)
   if( present(extcoef) .or. present(scacoef) .or. &
       present(bckcoef)) then

      if (present(extcoef)) extcoef = 0.
      if (present(scacoef)) scacoef = 0.
      if (present(bckcoef)) bckcoef = 0.

      do w = 1, size(wavelengths_profile)
         call mie%Query(wavelengths_profile(w), 1, & ! Only SO4 exists in the MieTable, so its index is 1
                        SO4*delp/grav, rh,         &
                        tau=tau, ssa=ssa, bbck=bck,__RC__)

!         Calculate the total ext. and scat. coefficients
         if( present(extcoef) ) then
              extcoef(:,:,:,w) = extcoef(:,:,:,w) + &
                              tau * (grav * rhoa / delp)
         endif
         if( present(scacoef) ) then
              scacoef(:,:,:,w) = scacoef(:,:,:,w) + &
                              ssa * tau * (grav * rhoa / delp)
         endif
         if( present(bckcoef) ) then
              bckcoef(:,:,:,w) = bckcoef(:,:,:,w) + &
                              bck * SO4 * rhoa 
         endif
      enddo
   endif

   if( present(exttau) .or. present(stexttau) .or. &
       present(scatau) .or. present(stscatau)) then

      if (present(exttau)) exttau = 0.
      if (present(stexttau)) stexttau = 0.
      if (present(scatau)) scatau = 0.
      if (present(stscatau)) stscatau = 0.

      do w = 1, size(wavelengths_vertint)
         call mie%Query(wavelengths_vertint(w), 1,  & ! Only SO4 exists in the MieTable, so its index is 1
                        SO4*delp/grav, rh,          &
                        tau=tau, ssa=ssa, __RC__)

         do k = klid, km
!           Integrate in the vertical
            if ( present(exttau) ) then
               exttau(:,:,w) = exttau(:,:,w) + tau(:,:,k)
            endif

            if (present(stexttau) ) then
               where (ple(:,:,k) .le. tropp) 
                  stexttau(:,:,w) = stexttau(:,:,w) + tau(:,:,k)
               elsewhere(ple(:,:,k-1) .lt. tropp) 
                 stexttau(:,:,w)  = stexttau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)
               endwhere
            endif

            if ( present(scatau) ) then
               scatau(:,:,w) = scatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
            endif

            if ( present(stscatau) ) then
               where (ple(:,:,k) .le. tropp) 
                  stscatau(:,:,w) = stscatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
               elsewhere(ple(:,:,k-1) .lt. tropp) 
                  stscatau(:,:,w) = stscatau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)*ssa(:,:,k)
               endwhere
           endif
         enddo
      enddo
   endif

!  Calculate the 470-870 Angstrom parameter
   if( present(angstrom) .and. do_angstrom ) then

      angstrom(i1:i2,j1:j2) = 0.
!     Set tau to small number by default
      tau470(i1:i2,j1:j2) = tiny(1.0)
      tau870(i1:i2,j1:j2) = tiny(1.0)

      call mie%Query(4.70E-7,  1,       & ! Only SO4 exists in the MieTable, so its index is 1
                     SO4*delp/grav, rh, &
                     tau=tau, __RC__)
      do k = klid, km
         tau470 = tau470 + tau(:,:,k)
      enddo

      call mie%Query(8.70E-7, 1,       &
                     SO4*delp/grav,rh, &
                     tau=tau, __RC__)
      do k = klid, km
         tau870 = tau870 + tau(:,:,k)
      enddo

!      enddo  ! nbins
      angstrom(i1:i2,j1:j2) = &
        -log(tau470(i1:i2,j1:j2)/tau870(i1:i2,j1:j2)) / &
         log(470./870.)
   endif

!  Calculate the sulfate surface area density [m2 m-3], possibly for use in
!  StratChem or other component.  Assumption here is a specified effective
!  radius (gcSU%radius for sulfate) and standard deviation of lognormal
!  distribution.  Hydration is by grid box provided RH and is follows Petters
!  and Kreeidenweis (ACP2007)
   if(present(sarea) .or. present(snum)) then
!        rmed   = w_c%reg%rmed(n1+nSO4-1)                    ! median radius, m
        if(rmed > 0.) then
!         sigma  = w_c%reg%sigma(n1+nSO4-1)                  ! width of lognormal distribution
         do k = klid, km
         do j = j1, j2
          do i = i1, i2
           rh_ = min(0.95,rh(i,j,k))
           gf = (1. + 1.19*rh_/(1.-rh_) )                   ! ratio of wet/dry volume, eq. 5
           rwet = rmed * gf**(1./3.)                      ! wet effective radius, m
!          Wet particle volume m3 m-3
           svol = SO4(i,j,k) * rhoa(i,j,k) / rhop * gf
!          Integral of lognormal surface area m2 m-3
           if(present(sarea)) sarea(i,j,k) = 3./rwet*svol*exp(-5./2.*alog(sigma)**2.)
!          Integral of lognormal number density # m-3
           if(present(snum)) snum(i,j,k) = svol / (rwet**3) * exp(-9/2.*alog(sigma)**2.) * 3./4./pi
          enddo
         enddo
        enddo
       endif
   endif
   deallocate(tau,ssa)
   __RETURN__(__SUCCESS__)
   end subroutine SU_Compute_Diags

!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver

   subroutine SulfateChemDriver (km, klid, cdt, PI, radToDeg, von_karman, &
                                 airMolWght, nAvogadro, cpd, grav, &
                                 fMassMSA, fMassDMS, fMassSO2, fMassSO4, &
                                 nymd, nhms, lonRad, latRad, &
                                 dms, so2, so4, msa, &
                                 nDMS, nSO2, nSO4, nMSA, &
                                 xoh, xno3, xh2o2, h2o2_init, &
                                 delp, tmpu, cloud, rhoa, hghte, &
                                 ustar, shflux, oro, pblh, z0h, &
                                 SU_dep, SU_PSO2, SU_PMSA, &
                                 SU_PSO4, SU_PSO4g, SU_PSO4aq, &     ! 2d diagnostics
                                 pso2, pmsa, pso4, pso4g, pso4aq, drydepositionfrequency, & ! 3d diagnostics
                                 rc)


! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: PI     ! pi constnat
   real, intent(in)    :: radToDeg ! radians to degree conversion
   real, intent(in)    :: von_karman ! Von Karman constant [unitless]
   real, intent(in)    :: nAvogadro  ! Avogadro's number [1/kmol]
   real, intent(in)    :: airMolWght ! molecular weight of air [kg/kmol]
   real, intent(in)    :: cpd
   real, intent(in)    :: grav   ! gravity [m/sec]
   real, intent(in)    :: fMassMSA, fMassDMS, fMassSO2, fMassSO4 ! gram molecular weights of species
   integer, intent(in) :: nymd   ! model year month day
   integer, intent(in) :: nhms   ! model hour mintue second
   real, dimension(:,:), intent(in) :: lonRad   ! model grid lon [radians]
   real, dimension(:,:), intent(in) :: latRad   ! model grid lat [radians]
   real, dimension(:,:,:), intent(inout) :: dms  ! dimethyl sulfide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: so2  ! sulfer dioxide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: so4  ! sulfate aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: msa  ! methanesulphonic acid [kg/kg]
   integer, intent(in) :: nDMS, nSO2, nSO4, nMSA ! index position of sulfates
   real, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: tmpu   ! temperature [K]
   real, dimension(:,:,:), intent(in) :: cloud  ! cloud fraction for radiation [1]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa   ! layer air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: hghte  ! top of layer geopotential height [m]
   real, pointer, dimension(:,:), intent(in)   :: ustar  ! surface velocity scale [m/sec]
   real, pointer, dimension(:,:), intent(in)   :: shflux ! sensible heat flux from turbulence [w/m^2]
   real, pointer, dimension(:,:), intent(in)   :: oro    ! land-ocean-ice mask
   real, pointer, dimension(:,:), intent(in)   :: pblh   ! planetary boundary layer height [m]
   real, pointer, dimension(:,:), intent(in)   :: z0h    ! surface roughness for heat [m]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: xoh, xno3, xh2o2 ! OH, NO3, H2O2 respectievly [kg/kg]
   real, dimension(:,:,:) :: h2o2_init ! private H2O2 that is saved and used to initialize [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO2 ! SO2 Prod from DMS Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PMSA ! MSA Prod from DMS Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO4 ! SO4 Prod from All SO2 Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO4g ! SO4 Prod from Gaseous SO2 Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO4aq ! SO4 Prod from Aqueous SO2 Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso2 ! SO2 Prod from DMS oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pmsa ! MSA Prod from DMS oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso4 ! SO4 Prod from all SO2 oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso4g ! SO4 Prod from gaseous SO2 oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso4aq ! SO4 Prod from aqueous SO2 oxidation [kg m-2 s-1]
   real, dimension(:,:), allocatable, intent(out) :: drydepositionfrequency

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -

! !DESCRIPTION: Updates the SU concentration due to chemistry
!  The SU grid component is currently established with 4 different
!  species (bins) following this convection:
!   1) DMS
!   2) SO2
!   3) SO4
!   4) MSA
!  Accordingly we have 4 chemical cycles to follow through, which are
!  sub-subroutines under this one.
!  The chemistry is a function of OH, NO3, and H2O2 concentrations
!  as well as DMS, SO2, SO4, MSA concentrations.  It is also a function
!  of solar zenith angle and temperature.  We pass in temperature.  SZA
!  will be a function of time of day and lat/lon.  For now we simply add
!  this to the grid component before calculating it.  I bet this is
!  somewhere else in the model.

!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!  30july2020 E.Sherman - ported to process library
!

! !Local Variables
   real, dimension(:,:), allocatable :: cossza, sza
   integer :: k, jday, i2, j2
   real, dimension(:,:,:), allocatable :: pSO2_DMS, pMSA_DMS, pSO4g_SO2, pSO4aq_SO2
   real    :: xhour
   integer :: status


!EOP
!-------------------------------------------------------------------------
!  Begin

   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

   allocate(drydepositionfrequency, mold=oro)
   allocate(cossza, mold=oro)
   allocate(sza, mold=oro)
!AOO initialization
   drydepositionfrequency=0.;cossza=0.;sza=0.
!AOO end initialization

   drydepositionfrequency = 0.0
   cossza = 0.0
   sza = 0.0

!  Reset the production terms
   allocate(pSO2_DMS, mold=tmpu)
   allocate(pMSA_DMS, mold=tmpu)
   allocate(pSO4g_SO2, mold=tmpu)
   allocate(pSO4aq_SO2, mold=tmpu)
   pSO2_DMS = 0.
   pMSA_DMS = 0.
   pSO4g_SO2 = 0.
   pSO4aq_SO2 = 0.

   if( associated(su_pSO2) )  su_pSO2   = 0.
   if( associated(su_pMSA) )  su_pMSA   = 0.
   if( associated(su_pSO4) )  su_pSO4   = 0.
   if( associated(su_pSO4g) )  su_pSO4g  = 0.
   if( associated(su_pSO4aq) )  su_pSO4aq = 0.
   if( associated(pSO2) )     pSO2   = 0.
   if( associated(pMSA) )     pMSA   = 0.
   if( associated(pSO4) )     pSO4   = 0.
   if( associated(pSO4g) )    pSO4g  = 0.
   if( associated(pSO4aq) )   pSO4aq = 0.


!  Find the cossza
!  ----------------------------------
   jday = idaynum(nymd)
   xhour = (  real(nhms/10000)*3600. &
            + real(mod(nhms,10000)/100)*60. &
            + real(mod(nhms,100)) &
           ) / 3600.

   call szangle (jday, xhour, lonRad, latRad, PI, radToDeg, sza, cossza, i2, j2)
!  Reset the dry deposition fluxes & frequencies
   if( associated(su_dep) ) su_dep = 0.0

   call DryDeposition ( km, tmpu, rhoa, hghte, oro, ustar, pblh, shflux, &
                        von_karman, cpd, grav, z0h, drydepositionfrequency, __RC__)


!  Now call the chemistry packages...
!  ----------------------------------
!  DMS source and oxidation to SO2 and MSA
   call SulfateChemDriver_DMS (km, klid, cdt, airMolWght, nAvogadro, cpd,&
                               fMassMSA, fMassDMS, fMassSO2, &
                               dms, nDMS, xoh, xno3, &
                               cossza, tmpu, rhoa, &
                               pSO2_DMS, pMSA_DMS, SU_dep, &
                               __RC__)

   if( associated(pSO2) )  pSO2 = pSO2_DMS
   if( associated(su_pSO2)) then
     do k = klid, km
      su_pSO2(:,:) = su_pSO2(:,:) + pSO2_DMS(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

   if( associated(pMSA) )  pMSA = pMSA_DMS
   if( associated(su_pMSA)) then
     do k = klid, km
      su_pMSA(:,:) = su_pMSA(:,:) + pMSA_DMS(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

!  SO2 source and oxidation to SO4
   call SulfateChemDriver_SO2 (km, klid, cdt, airMolWght, nAvogadro, cpd, grav, &
                               fMassSO4, fMassSO2, &
                               so2, nSO2, xoh, xh2o2, &
                               tmpu, rhoa, delp, oro, cloud, drydepositionfrequency, &
                               pSO2_DMS, pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                               __RC__)

   if( associated(pSO4g) )  pSO4g = pSO4g_SO2
   if( associated(su_pSO4g)) then
     do k = klid, km
      su_pSO4g(:,:) = su_pSO4g(:,:) + pSO4g_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

   if( associated(pSO4aq) )  pSO4aq = pSO4aq_SO2
   if( associated(su_pSO4aq)) then
     do k = klid, km
      su_pSO4aq(:,:) = su_pSO4aq(:,:) + pSO4aq_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

   if( associated(pSO4) ) pSO4 = pSO4g_SO2 + pSO4aq_SO2
   if( associated(su_pSO4)) then
     do k = klid, km
      su_pSO4(:,:) = su_pSO4(:,:) + pSO4g_SO2(:,:,k)*delp(:,:,k)/grav &
                     + pSO4aq_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

!  SO4 source and loss
   call SulfateChemDriver_SO4 (km, klid, cdt, grav, so4, nSO4, delp, &
                               drydepositionfrequency, pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                               __RC__)

!  MSA source and loss
   if( associated(msa)) then
      call SulfateChemDriver_MSA (km, klid, cdt, grav, msa, nMSA, delp, &
                                  drydepositionfrequency, pMSA_DMS, SU_dep, &
                                  __RC__)
   end if

!  Save the h2o2 value after chemistry
   h2o2_init = xh2o2

   __RETURN__(__SUCCESS__)
   end subroutine SulfateChemDriver

!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_DMS

   subroutine SulfateChemDriver_DMS (km, klid, cdt, airMolWght, nAvogadro, cpd, &
                                     fMassMSA, fMassDMS, fMassSO2, &
                                     qa, nDMS, xoh, xno3, &
                                     cossza, tmpu, rhoa, &
                                     pSO2_DMS, pMSA_DMS, SU_dep, &
                                     rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: nAvogadro  ! Avogadro's number [1/kmol]
   real, intent(in)    :: airMolWght ! molecular weight of air [kg/kmol]
   real, intent(in)    :: cpd
   real, intent(in)    :: fMassMSA, fMassDMS, fMassSO2 ! gram molecular weights of species
   integer, intent(in) :: nDMS       !index position of sulfates
   real, dimension(:,:,:), intent(in) :: xoh, xno3  ! OH, NO3 respectievly [kg/kg]
   real, dimension(:,:), intent(in)   :: cossza
   real, dimension(:,:,:), intent(in) :: tmpu   ! temperature [K]
   real, dimension(:,:,:), intent(in) :: rhoa   ! layer air density [kg/m^3]


! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), allocatable,  intent(out) :: pSO2_DMS ! SO2 production from DMS oxidation [kg kg-1 s-1]
   real, dimension(:,:,:), allocatable,  intent(out) :: pMSA_DMS ! MSA production from DMS oxidation [kg kg-1 s-1]
   integer, optional, intent(out)   :: rc

! !DESCRIPTION: Computes the production of SO2 and MSA due to DMS oxidation
!
!   R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel
!          k1 = { 1.7d-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }
!          a = 0.75, b = 0.25
!
!   R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel
!          k2 = 1.2e-11*exp(-260/T)
!
!      DMS_OH = DMS0 * exp(-(r1+r2)*NDT1)
!          where DMS0 is the DMS concentration at the beginning,
!          r1 = k1*[OH], r2 = k2*[OH]
!
!   R3:    DMS + NO3 ->   SO2 + ...
!          k3 = 1.9e-13*exp(500/T)
!
!      DMS = DMS_OH * exp(-r3*NDT1)
!          where r3 = k3*[NO3]
!
!   R4:    DMS + X   ->   SO2 + ...
!          assume to be at the rate of DMS+OH and DMS+NO3 combined.
!
!   The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved
!   for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They
!   are in unit of MixingRatio/second.
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!
!  03Aug2020 E.Sherman - ported to process library

! !Local Variables
   integer :: i, j, k, i1=1, j1=1, i2, j2
   real*8  :: Fx, b, eff
   real*8  :: rk1, rk2, rk3, rk4
   real*8  :: tk, o2, oh, no3, air
   real*8  :: dms, dms0, dms_oh

   data Fx  / 1.0 /
   data b   / 0.25 /
   data eff / 1. /

!EOP
!-------------------------------------------------------------------------
! Begin...

   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

   allocate(pSO2_DMS, mold=tmpu)
   allocate(pMSA_DMS, mold=tmpu)
   pSO2_DMS = 0.0
   pMSA_DMS = 0.0

!  spatial loop
   do k = klid, km
    do j = j1, j2
     do i = i1, i2

      rk1 = 0.
      rk2 = 0.
      rk3 = 0.
      rk4 = 0.

      tk  = tmpu(i,j,k)
      oh  = xoh(i,j,k)
!     air molecules in # cm-3
      air = 1000.*rhoa(i,j,k) / airMolWght * nAvogadro * 1.e-6
!     oxygen molecules in # cm-3
      o2 = 0.21 * air
!     no3 -> go from volume mixing ratio to # cm-3
      no3 = xno3(i,j,k) * air

!     initial DMS concentration (kg kg-1)
      dms0 = qa(i,j,k)
      dms0 = max(dms0,tiny(dms0))

!     1 & 2) DMS + OH: RK1 = addition, RK2 = abstraction
      if(oh .gt. 0.) then
       rk1 = (1.7d-42 * exp(7810./tk) * o2) / &
             (1. + 5.5e-31 * exp(7460./tk) * o2) * oh
       rk2 = (1.2e-11 * exp(-260./tk)) * oh
      endif

!     3) DMS + NO3: only happens at night
      if(cossza(i,j) .le. 0.) then
       rk3 = (1.9e-13 * exp(500./tk)) * no3
      endif

!     Now do the DMS loss
      dms_oh = dms0   * exp( -(rk1+rk2)* Fx * cdt)
      dms    = dms_oh * exp( -(rk3)    * Fx * cdt)

!     SO2 and MSA production terms
!     MSA is formed from the DMS+OH addition step
!     Production should go as mass mixing ratio change in MSA
      if( (rk1+rk2) .eq. 0.) then
       pMSA_DMS(i,j,k) = 0.
      else
       pMSA_DMS(i,j,k) =  (dms0 - dms_oh) * b*rk1/((rk1+rk2)*Fx) * eff &
                         * (fMassMSA/fMassDMS) / cdt
      endif

!     Everything else goes into SO2 formation step
      pSO2_DMS(i,j,k) = ( dms0 - dms - &
                          pMSA_DMS(i,j,k)*cdt*(fMassDMS/fMassMSA) &
                        ) * (fMassSO2/fMassDMS) / cdt


!     4) Dry deposition of DMS (not in GOCART?)
!      if(k .eq. km) rk4 = drydepf(i,j)
!      dms0 = dms
!      dms  = dms0 * exp(-rk4*cdt)
!      dms    = max(dms,1.e-32)

!     Update the mass mixing ratio and the dry deposition flux out of DMS
      dms    = max(dms,tiny(dms))
      qa(i,j,k) = dms

     end do ! i
    end do  ! j
    if(k .eq. km .and. associated(SU_dep) ) SU_dep(:,:,nDMS) = 0.
   end do   ! k


   __RETURN__(__SUCCESS__)
   end subroutine SulfateChemDriver_DMS


!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_SO2

   subroutine SulfateChemDriver_SO2 (km, klid, cdt, airMolWght, nAvogadro, cpd, grav, &
                                     fMassSO4, fMassSO2, &
                                     qa, nSO2, xoh, xh2o2, &
                                     tmpu, rhoa, delp, oro, cloud, drydepf, &
                                     pSO2_DMS, pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                                     rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: nAvogadro  ! Avogadro's number [1/kmol]
   real, intent(in)    :: airMolWght ! molecular weight of air [kg/kmol]
   real, intent(in)    :: cpd
   real, intent(in)    :: grav       ! gravity [m/sec]
   real, intent(in)    :: fMassSO4, fMassSO2 ! gram molecular weights of species
   integer, intent(in) :: nSO2       !index position of sulfates
   real, dimension(:,:,:), intent(in) :: tmpu   ! temperature [K]
   real, dimension(:,:,:), intent(in) :: rhoa   ! layer air density [kg/m^3]
   real, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]
   real, dimension(:,:,:), intent(in) :: cloud  ! cloud fraction for radiation [1]
   real, dimension(:,:), intent(in)   :: drydepf  ! dry deposition frequency [s-1]
   real, pointer, dimension(:,:), intent(in) :: oro  ! land-ocean-ice mask
   real, dimension(:,:,:), intent(in) :: pSO2_DMS ! SO2 production from DMS oxidation [kg m-2 s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: xoh, xh2o2  ! OH, H2O2 respectievly [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), allocatable, intent(out) :: pSO4g_SO2 ! SO4 production - gas phase [kg kg-1 s-1]
   real, dimension(:,:,:), allocatable, intent(out) :: pSO4aq_SO2 ! SO4 production - aqueous [kg kg-1 s-1]
   integer, optional, intent(out)   :: rc

! !DESCRIPTION: Computes the concentration of SO2 and production of SO4
!
!  SO2 production:
!    DMS + OH, DMS + NO3 (saved in SU_ChemDrv_DMS)
!
!  SO2 loss:
!    SO2 + OH  -> SO4
!    SO2       -> drydep
!    SO2 + H2O2 or O3 (aq) -> SO4
!
!  SO2 = SO2_0 * exp(-bt)
!      + PSO2_DMS*dt/bt * [1-exp(-bt)]
!    where b is the sum of the reaction rate of SO2 + OH and the dry
!    deposition rate of SO2, PSO2_DMS is SO2 production from DMS in
!    MixingRatio/timestep.
!
!  If there is cloud in the gridbox (fraction = fc), then the aqueous
!  phase chemistry also takes place in cloud. The amount of SO2 oxidized
!  by H2O2 in cloud is limited by the available H2O2; the rest may be
!  oxidized due to additional chemistry, e.g, reaction with O3 or O2
!  (catalyzed by trace metal).
!
! !REVISION HISTORY:
!  06Nov2003, Colarco - Based on Ginoux!
!  15Jul2010, Colarco - modularized
!  03Aug2020 E.Sherman - ported to process library


! !Local Variables
   integer :: i, j, k, j2, i2
   real*8  :: rk1, rk2, rk, rkt, f1
   real*8  :: L1, L2, Ld, SO2, SO2_cd, fc, fMR
   real*8  :: oh, h2o2, SO20, tk, air, k0, ki, kk
   real, dimension(:,:), allocatable :: fout

   data ki / 1.5e-12 /

!EOP
!-------------------------------------------------------------------------
! Begin...

   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

   allocate(pSO4g_SO2, mold=tmpu)
   allocate(pSO4aq_SO2, mold=tmpu)
   allocate(fout(i2,j2))
   pSO4g_SO2 = 0.0
   pSO4aq_SO2 = 0.0
   fout = 0.0

!  Conversion of SO2 mmr to SO2 vmr
   fMR = airMolWght / fMassSO2

!  Initialize flux variable
   fout = 0.

!  spatial loop
   do k = klid, km
    do j = 1, j2
     do i = 1, i2

      rk1 = 0.
      rk2 = 0.
      L1  = 0.
      L2  = 0.
      Ld  = 0.

      tk   = tmpu(i,j,k)
      oh   = xoh(i,j,k)
      h2o2 = max(xh2o2(i,j,k),tiny(xh2o2(i,j,k)))

!     air molecules in # cm-3
      air  = 1000.*rhoa(i,j,k) / airMolWght * nAvogadro * 1.e-6
!     1) SO2 + OH(g) in s-1
      k0 = 3.0e-31 * (300./tk)**4.3
      kk = k0 * air / ki
      f1 = (1. + (log10(kk))**2.)**(-1.)
      rk1 = ( (k0*air/(1.+kk)) * 0.6**f1) * oh

!     2) rk2 is the loss of SO2 due to dry deposition.
      if(k .eq. km) then
!      drydepf calculated for aerosol
!      follow Walcek: ocean drydepf_so2 = 10*drydepf_aer
!      or if land drydepf_so2 = 3*drydepf_aer
       if(oro(i,j) .eq. OCEAN) then
        rk2 = 10.*drydepf(i,j)
       else
        rk2 = 3.*drydepf(i,j)
       endif
!       rk2 = drydepf(i,j)
      else
       rk2 = 0.
      endif

      rk = (rk1 + rk2)
      rkt = rk*cdt

!     Update the SO2 concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.

!     initial SO2 concentration (kg kg-1) after adding source
      SO20 = qa(i,j,k) + pSO2_DMS(i,j,k)*cdt
      SO20 = max(SO20,tiny(SO20))

      if(rk .gt. 0.) then
       SO2_cd =  SO20 * exp(-rkt)
       L1     = (SO20 - SO2_cd) * rk1/rk
       if(k .eq. km) then
        Ld    = (SO20 - SO2_cd) * rk2/rk
        fout(i,j) = Ld * delp(i,j,km)/grav/cdt
       else
        Ld    = 0.
       endif
      else
       SO2_cd = SO20
       L1     = 0.
      endif

!     Update SO2 concentration after cloud chemistry, if it occurs
      fc = cloud(i,j,k)
      if(fc .gt. 0. .and. SO2_cd .gt. 0. .and. tk .gt. 258.) then
!      Check on H2O2 vmr -> is SO2 vmr greater?
       if(fMr * SO2_cd .gt. h2o2) then
        fc = fc*(h2o2/(fMR*SO2_cd))
        h2o2 = h2o2*(1.-cloud(i,j,k))
       else
        h2o2 = h2o2*(1. - cloud(i,j,k)*(fMR*SO2_cd)/h2o2)
       endif
       SO2 = SO2_cd*(1.-fc)
!      aqueous loss rate (mixing ratio/timestep)
       L2 = SO2_cd * fc
      else
       SO2 = SO2_cd
       L2 = 0.
      endif

!     Ideally you would update the H2O2 mixing ratio at this point
!     and then reset it periodically
      xh2o2(i,j,k) = max(h2o2,tiny(h2o2))

      SO2 = max(SO2,tiny(SO2))
      qa(i,j,k) = SO2
      pSO4g_SO2(i,j,k) = L1 * (fMassSO4/fMassSO2) / cdt
      pSO4aq_SO2(i,j,k) = L2 * (fMassSO4/fMassSO2) / cdt

     end do
    end do
   end do

   if( associated(SU_dep) ) SU_dep(:,:,nSO2) = fout

   __RETURN__(__SUCCESS__)
   end subroutine SulfateChemDriver_SO2

!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_SO4

   subroutine SulfateChemDriver_SO4 (km, klid, cdt, grav, qa, nSO4, delp, drydepf, &
                                     pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                                     rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: grav   ! gravity [m/sec]
   integer, intent(in) :: nSO4   ! index position of sulfate
   real, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]
   real, dimension(:,:), intent(in)   :: drydepf    ! dry deposition frequency [s-1]
   real, dimension(:,:,:), intent(in) :: pSO4g_SO2  ! SO4 production - gas phase [kg kg-1 s-1]
   real, dimension(:,:,:), intent(in) :: pSO4aq_SO2 ! SO4 production - aqueous [kg kg-1 s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc

! !DESCRIPTION:
!
!  SO4 production:
!    The only production term is due to SO2 oxidation.
!    SO4 = SO4_0 * exp(-kt) + pSO4_SO2/kt * (1.-exp(-kt))
!     where k is the dry deposition
!
! !REVISION HISTORY:
!
!  15Jul2010, Colarco - Modularized
!  06Nov2003, Colarco
!  Based on Ginoux
!
!  03Aug2020 E.Sherman - ported to process library
!
! !Local Variables
   integer :: i, j, k, i2, j2
   real*8  :: rk, rkt, Ld
   real*8  :: SO4, SO40, pSO4
   real, dimension(:,:), allocatable :: fout

!EOP
!-------------------------------------------------------------------------

! Begin...

   j2 = ubound(qa, 2)
   i2 = ubound(qa, 1)

   allocate(fout(i2,j2))
   fout=0.   !AOO initialization

!  Initialize flux variable
   fout = 0.

!  spatial loop
   do k = klid, km
    do j = 1, j2
     do i = 1, i2

      pSO4 = pSO4g_SO2(i,j,k)+pSO4aq_SO2(i,j,k)

!     initial SO4 concentration (kg kg-1)
      SO40 = qa(i,j,k)
      SO40 = max(SO40,tiny(SO40))

!     Update the SO4 concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.
!     RK is the dry deposition frequency
      if(k .eq. km) then
       RK = drydepf(i,j)
       RKT = RK*cdt
       SO4 = (SO40 + pSO4*cdt) * exp(-rkt)
       Ld  = (SO40 - SO4 + pSO4*cdt)
       fout(i,j) = Ld * delp(i,j,km)/grav/cdt
      else
       SO4 = SO40 + pSO4*cdt
       Ld = 0.
      endif

      SO4 = max(SO4,tiny(SO4))
      qa(i,j,k) = SO4

     end do
    end do
   end do

   if( associated(SU_dep) ) SU_dep(:,:,nSO4) = fout

   __RETURN__(__SUCCESS__)
   end subroutine SulfateChemDriver_SO4

!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_MSA

   subroutine SulfateChemDriver_MSA (km, klid, cdt, grav, qa, nMSA, delp, drydepf, &
                                     pMSA_DMS, SU_dep, &
                                     rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: grav   ! gravity [m/sec]
   integer, intent(in) :: nMSA   ! index position of sulfate
   real, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]
   real, dimension(:,:), intent(in)   :: drydepf   ! dry deposition frequency [s-1]
   real, dimension(:,:,:), intent(in) :: pMSA_DMS  ! MSA production - gas phase [kg kg-1 s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc

! !DESCRIPTION:
!
!  MSA production:
!    The only production term is due to DMS oxidation.
!    MSA = MSA_0 * exp(-kt) + pSO4_SO2/kt * (1.-exp(-kt))
!     where k is the dry deposition
!
! !REVISION HISTORY:
!
!  15Jul2010, Colarco -- modularized
!  06Nov2003, Colarco
!  Based on Ginoux
!
!  03Aug2020 E.Sherman - ported to process library

! !Local Variables
   integer :: i, j, k, i2, j2
   real*8  :: rk, rkt, Ld
   real*8  :: MSA, MSA0
   real, dimension(:,:), allocatable :: fout

!EOP
!-------------------------------------------------------------------------
! Begin...

   j2 = ubound(qa, 2)
   i2 = ubound(qa, 1)

   allocate(fout(i2,j2))
   fout = 0.0

!  spatial loop
   do k = klid, km
    do j = 1, j2
     do i = 1, i2

!     initial MSA concentration (kg kg-1)
      MSA0 = qa(i,j,k)
      MSA0 = max(MSA0,tiny(MSA0))

!     Update the MSA concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.
!     RK is the dry deposition frequency
      if(k .eq. km) then
       RK = drydepf(i,j)
       RKT = RK*cdt
       MSA = (MSA0 + pMSA_DMS(i,j,k)*cdt) * exp(-rkt)
       Ld  = (MSA0 + pMSA_DMS(i,j,k)*cdt - MSA)
       fout(i,j) = Ld * delp(i,j,km)/grav/cdt
      else
       MSA = MSA0 + pMSA_DMS(i,j,k)*cdt
       Ld = 0.
      endif

      MSA = max(MSA,tiny(MSA))
      qa(i,j,k) = MSA

     end do
    end do
   end do

   if( associated(SU_dep) ) SU_dep(:,:,nMSA) = fout

   __RETURN__(__SUCCESS__)
   end subroutine SulfateChemDriver_MSA

!==================================================================================
!BOP
! !IROUTINE: get_HenrysLawCts

   subroutine get_HenrysLawCts(name,c1,c2,c3,c4,rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   character(len=*), intent(in) :: name

! !OUTPUT PARAMETERS:
   real, intent(out):: c1,c2,c3,c4
   integer, optional, intent(out) :: rc                   ! Error return code:


! !DESCRIPTION: Provides Henry's Law coefficients for species.
!
! !REVISION HISTORY:
!
! 05Aug2020 E.Sherman - Moved over from Henrys_Law_cts.F90
!

! !Local Variables
   integer :: l,found

  INTEGER,PARAMETER :: nspecies_HL=051
  REAL   ,PARAMETER :: notfound = -1.

  !Name of species
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies_HL) :: spc_name=(/ &
      'O3  ' & !001
     ,'H2O2' & !002
     ,'NO  ' & !003
     ,'NO2 ' & !004
     ,'NO3 ' & !005
     ,'N2O5' & !006
     ,'HONO' & !007
     ,'HNO3' & !008
     ,'HNO4' & !009
     ,'SO2 ' & !010
     ,'SULF' & !011
     ,'CO  ' & !012
     ,'CO2 ' & !013
     ,'N2  ' & !014
     ,'O2  ' & !015
     ,'H2O ' & !016
     ,'H2  ' & !017
     ,'O3P ' & !018
     ,'O1D ' & !019
     ,'HO  ' & !020
     ,'HO2 ' & !021
     ,'CH4 ' & !022
     ,'ETH ' & !023
     ,'ALKA' & !024
     ,'ALKE' & !025
     ,'BIO ' & !026
     ,'ARO ' & !027
     ,'HCHO' & !028
     ,'ALD ' & !029
     ,'KET ' & !030
     ,'CRBO' & !031
     ,'ONIT' & !032
     ,'PAN ' & !033
     ,'OP1 ' & !034
     ,'OP2 ' & !035
     ,'ORA1' & !036
     ,'ORA2' & !037
     ,'MO2 ' & !038
     ,'AKAP' & !039
     ,'AKEP' & !040
     ,'BIOP' & !041
     ,'PHO ' & !042
     ,'ADD ' & !043
     ,'AROP' & !044
     ,'CBOP' & !045
     ,'OLN ' & !046
     ,'XO2 ' & !047
     ,'DMS ' & !048
     ,'NH3 ' & !049
     ,'CFC ' & !050
     ,'N2O ' & !050
   /)


  !Number of each specie
  INTEGER,PARAMETER :: O3  =001
  INTEGER,PARAMETER :: H2O2=002
  INTEGER,PARAMETER :: NO  =003
  INTEGER,PARAMETER :: NO2 =004
  INTEGER,PARAMETER :: NO3 =005
  INTEGER,PARAMETER :: N2O5=006
  INTEGER,PARAMETER :: HONO=007
  INTEGER,PARAMETER :: HNO3=008
  INTEGER,PARAMETER :: HNO4=009
  INTEGER,PARAMETER :: SO2 =010
  INTEGER,PARAMETER :: SULF=011
  INTEGER,PARAMETER :: CO  =012
  INTEGER,PARAMETER :: CO2 =013
  INTEGER,PARAMETER :: N2  =014
  INTEGER,PARAMETER :: O2  =015
  INTEGER,PARAMETER :: H2O =016
  INTEGER,PARAMETER :: H2  =017
  INTEGER,PARAMETER :: O3P =018
  INTEGER,PARAMETER :: O1D =019
  INTEGER,PARAMETER :: HO  =020
  INTEGER,PARAMETER :: HO2 =021
  INTEGER,PARAMETER :: CH4 =022
  INTEGER,PARAMETER :: ETH =023
  INTEGER,PARAMETER :: ALKA=024
  INTEGER,PARAMETER :: ALKE=025
  INTEGER,PARAMETER :: BIO =026
  INTEGER,PARAMETER :: ARO =027
  INTEGER,PARAMETER :: HCHO=028
  INTEGER,PARAMETER :: ALD =029
  INTEGER,PARAMETER :: KET =030
  INTEGER,PARAMETER :: CRBO=031
  INTEGER,PARAMETER :: ONIT=032
  INTEGER,PARAMETER :: PAN =033
  INTEGER,PARAMETER :: OP1 =034
  INTEGER,PARAMETER :: OP2 =035
  INTEGER,PARAMETER :: ORA1=036
  INTEGER,PARAMETER :: ORA2=037
  INTEGER,PARAMETER :: MO2 =038
  INTEGER,PARAMETER :: AKAP=039
  INTEGER,PARAMETER :: AKEP=040
  INTEGER,PARAMETER :: BIOP=041
  INTEGER,PARAMETER :: PHO =042
  INTEGER,PARAMETER :: ADD =043
  INTEGER,PARAMETER :: AROP=044
  INTEGER,PARAMETER :: CBOP=045
  INTEGER,PARAMETER :: OLN =046
  INTEGER,PARAMETER :: XO2 =047
  INTEGER,PARAMETER :: DMS =048
  INTEGER,PARAMETER :: NH3 =049
  INTEGER,PARAMETER :: CFC =050
  INTEGER,PARAMETER :: N2O =051

!     HENRYS LAW COEFFICIENTS
!     Henrys law coefficient
!     [KH298]=mole/(l atm)
!     Referencias em R. Sander (1999)
!     Compilation of Henry Law Constants
!     for Inorganic and Organic Species
!     of Potential Importance in
!     Environmental Chemistry (Version 3)
!     http://www.henrys-law.org
!     * indica artigos nao encontrados nesse endereço eletronico
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: hstar=(/&
    1.10E-2              ,   & ! O3 - 001
    8.30E+4              ,   & ! H2O2 - 002
    1.90E-3              ,   & ! NO - 003
    1.20E-2              ,   & ! NO2 - 004
    6.1E-01              ,   & ! NO3 - 005
    2.1E+00              ,   & ! N2O5 - 006
    5.00E+1              ,   & ! HONO - 007
    2.10E+5              ,   & ! HNO3 - 008
    1.20E+4              ,   & ! HNO4 - 009
    1.40E+0              ,   & ! SO2 - 010
    2.10E+5              ,   & ! SULF - 011
    9.90E-4              ,   & ! CO - 012
    3.6E-02              ,   & ! CO2 - 013
    6.1E-04              ,   & ! N2 - 014
    1.3E-03              ,   & ! O2 - 015
    0.0E+00              ,   & ! H2O - 016
    7.8E-04              ,   & ! H2 - 017
    0.00E+0              ,   & ! O3P - 018
    0.00E+0              ,   & ! O1D - 019
    3.00E+1              ,   & ! HO - 020
    5.70E+3              ,   & ! HO2 - 021
    1.40E-3              ,   & ! CH4 - 022
    1.90E-3              ,   & ! ETH - 023
    1.00E-3              ,   & ! ALKA - 024
    5.00E-3              ,   & ! ALKE - 025
    2.80E-2              ,   & ! BIO - 026
    1.73E-1              ,   & ! ARO - 027
    3.20E+3              ,   & ! HCHO - 028
    1.40E+1              ,   & ! ALD - 029
    3.00E+1              ,   & ! KET - 030
    2.1E+05              ,   & ! CRBO - 031
    1.00E+0              ,   & ! ONIT - 032
    3.60E+0              ,   & ! PAN - 033
    3.10E+2              ,   & ! OP1 - 034
    3.40E+2              ,   & ! OP2 - 035
    8.90E+3              ,   & ! ORA1 - 036
    4.10E+3              ,   & ! ORA2 - 037
    2.00E+3              ,   & ! MO2 - 038
    0.0E+00              ,   & ! AKAP - 039
    0.0E+00              ,   & ! AKEP - 040
    0.0E+00              ,   & ! BIOP - 041
    0.0E+00              ,   & ! PHO - 042
    0.0E+00              ,   & ! ADD - 043
    0.0E+00              ,   & ! AROP - 044
    1.14E+1              ,   & ! CBOP - 045
    0.0E+00              ,   & ! OLN - 046
    0.0E+00              ,   & ! XO2 - 047
    5.6E-01              ,   & ! DMS - 048
    1.05E+06              ,   & ! NH3 - 048
    -1.                  ,   & ! CFC - 048
    2.4E-02                  & ! N2O - 051
    /)


!     -DH/R (for temperature correction)
!     [-DH/R]=K
!     Referencias em R. Sander (1999)
!     Compilation of Henry Law Constants
!     for Inorganic and Organic Species
!     of Potential Importance in
!     Environmental Chemistry (Version 3)
!     http://www.henrys-law.org
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: dhr=(/&
    2400.         ,   & ! O3 - 001
    7400.         ,   & ! H2O2 - 002
    1400.         ,   & ! NO - 003
    2500.         ,   & ! NO2 - 004
    2000.         ,   & ! NO3 - 005
    3400.         ,   & ! N2O5 - 006
    4900.         ,   & ! HONO - 007
    8700.         ,   & ! HNO3 - 008
    6900.         ,   & ! HNO4 - 009
    2900.         ,   & ! SO2 - 010
    0.            ,   & ! SULF - 011
    1300.         ,   & ! CO - 012
    2200.         ,   & ! CO2 - 013
    1300.         ,   & ! N2 - 014
    1500.         ,   & ! O2 - 015
    0.            ,   & ! H2O - 016
    500.          ,   & ! H2 - 017
    0.            ,   & ! O3P - 018
    0.            ,   & ! O1D - 019
    4500.         ,   & ! HO - 020
    5900.         ,   & ! HO2 - 021
    1600.         ,   & ! CH4 - 022
    2300.         ,   & ! ETH - 023
    2700.         ,   & ! ALKA - 024
    3000.         ,   & ! ALKE - 025
    0.            ,   & ! BIO - 026
    4045.         ,   & ! ARO - 027
    6800.         ,   & ! HCHO - 028
    5600.         ,   & ! ALD - 029
    4600.         ,   & ! KET - 030
    5300.         ,   & ! CRBO - 031
    5800.         ,   & ! ONIT - 032
    6500.         ,   & ! PAN - 033
    5200.         ,   & ! OP1 - 034
    6000.         ,   & ! OP2 - 035
    5700.         ,   & ! ORA1 - 036
    6300.         ,   & ! ORA2 - 037
    6600.         ,   & ! MO2 - 038
    0.            ,   & ! AKAP - 039
    0.            ,   & ! AKEP - 040
    0.            ,   & ! BIOP - 041
    0.            ,   & ! PHO - 042
    0.            ,   & ! ADD - 043
    0.            ,   & ! AROP - 044
    0.            ,   & ! CBOP - 045
    0.            ,   & ! OLN - 046
    0.            ,   & ! XO2 - 047
    3500.         ,   & ! DMS - 048
    4200.         ,   & ! NH3 - 048
    -1.           ,   & ! CFC - 048
    2700.             & ! N2O - 048
    /)


  REAL,PARAMETER,DIMENSION(nspecies_HL) :: weight=(/&
    48.  ,   & ! O3 - 001
    34.  ,   & ! H2O2 - 002
    30.  ,   & ! NO - 003
    46.  ,   & ! NO2 - 004
    62.  ,   & ! NO3 - 005
    108. ,   & ! N2O5 - 006
    47.  ,   & ! HONO - 007
    63.  ,   & ! HNO3 - 008
    79.  ,   & ! HNO4 - 009
    64.  ,   & ! SO2 - 010
    98.  ,   & ! SULF - 011
    28.  ,   & ! CO - 012
    44.  ,   & ! CO2 - 013
    28.  ,   & ! N2 - 014
    32.  ,   & ! O2 - 015
    18.  ,   & ! H2O - 016
    2.   ,   & ! H2 - 017
    16.  ,   & ! O3P - 018
    16.  ,   & ! O1D - 019
    17.  ,   & ! HO - 020
    33.  ,   & ! HO2 - 021
    16.  ,   & ! CH4 - 022
    30.  ,   & ! ETH - 023
    61.6 ,   & ! ALKA - 024
    33.0 ,   & ! ALKE - 025
    68.  ,   & ! BIO - 026
    97.9 ,   & ! ARO - 027
    30.  ,   & ! HCHO - 028
    44.  ,   & ! ALD - 029
    72.  ,   & ! KET - 030
    68.6 ,   & ! CRBO - 031
    119. ,   & ! ONIT - 032
    122. ,   & ! PAN - 033
    48.  ,   & ! OP1 - 034
    62.  ,   & ! OP2 - 035
    46.  ,   & ! ORA1 - 036
    60.  ,   & ! ORA2 - 037
    47.  ,   & ! MO2 - 038
    102. ,   & ! AKAP - 039
    88.4 ,   & ! AKEP - 040
    117. ,   & ! BIOP - 041
    107. ,   & ! PHO - 042
    107. ,   & ! ADD - 043
    151. ,   & ! AROP - 044
    85.4 ,   & ! CBOP - 045
    136. ,   & ! OLN - 046
    44.  ,   & ! XO2 - 047
    62.13,   & ! DMS - 048
    17.03,   & ! NH3 - 048
    -1.  ,   & ! CFC - 048
    44.      & ! CFC - 048
   /)


!    ACID DISSOCIATION CONSTANT AT 298K
!     [mole/liter of liquid water]
!     Referencias: Barth et al. JGR 112, D13310 2007
!     Martell and Smith, 1976, Critical stability
!     vol1-4 Plenum Press New York
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: ak0=(/&
    0.00E+00     ,   & ! O3 - 001
    2.20E-12     ,   & ! H2O2 - 002
    0.00E+00     ,   & ! NO - 003
    0.00E+00     ,   & ! NO2 - 004
    0.00E+00     ,   & ! NO3 - 005
    0.00E+00     ,   & ! N2O5 - 006
    7.10E-04     ,   & ! HONO - 007
    1.54E+01     ,   & ! HNO3 - 008
    0.00E+00     ,   & ! HNO4 - 009
    1.30E-02     ,   & ! SO2 - 010
    1.00E-02     ,   & ! SULF - 011
    0.00E+00     ,   & ! CO - 012
    4.50E-07     ,   & ! CO2 - 013
    0.00E+00     ,   & ! N2 - 014
    0.00E+00     ,   & ! O2 - 015
    0.00E+00     ,   & ! H2O - 016
    0.00E+00     ,   & ! H2 - 017
    0.00E+00     ,   & ! O3P - 018
    0.00E+00     ,   & ! O1D - 019
    0.00E+00     ,   & ! HO - 020
    3.50E-05     ,   & ! HO2 - 021
    0.00E+00     ,   & ! CH4 - 022
    0.00E+00     ,   & ! ETH - 023
    0.00E+00     ,   & ! ALKA - 024
    0.00E+00     ,   & ! ALKE - 025
    0.00E+00     ,   & ! BIO - 026
    0.00E+00     ,   & ! ARO - 027
    0.00E+00     ,   & ! HCHO - 028
    0.00E+00     ,   & ! ALD - 029
    0.00E+00     ,   & ! KET - 030
    0.00E+00     ,   & ! CRBO - 031
    0.00E+00     ,   & ! ONIT - 032
    0.00E+00     ,   & ! PAN - 033
    0.00E+00     ,   & ! OP1 - 034
    0.00E+00     ,   & ! OP2 - 035
    1.80E-04     ,   & ! ORA1 - 036
    1.75E-05     ,   & ! ORA2 - 037
    0.00E+00     ,   & ! MO2 - 038
    0.00E+00     ,   & ! AKAP - 039
    0.00E+00     ,   & ! AKEP - 040
    0.00E+00     ,   & ! BIOP - 041
    0.00E+00     ,   & ! PHO - 042
    0.00E+00     ,   & ! ADD - 043
    0.00E+00     ,   & ! AROP - 044
    0.00E+00     ,   & ! CBOP - 045
    0.00E+00     ,   & ! OLN - 046
    0.00E+00     ,   & ! XO2 - 047
    0.00E+00     ,   & ! DMS - 048
    0.00E+00     ,   & ! NH3 - 049
    0.00E+00     ,   & ! NH3 - 049
    0.00E+00         & ! CFC - 050
   /)

!     Temperature correction factor for
!     acid dissociation constants
!     [K]
!     Referencias: Barth et al. JGR 112, D13310 2007
  REAL,PARAMETER,DIMENSION(nspecies_HL) :: dak=(/&
    0.         ,   & ! O3 - 001
    -3700.     ,   & ! H2O2 - 002
    0.         ,   & ! NO - 003
    0.         ,   & ! NO2 - 004
    0.         ,   & ! NO3 - 005
    0.         ,   & ! N2O5 - 006
    0.         ,   & ! HONO - 007
    0.         ,   & ! HNO3 - 008
    0.         ,   & ! HNO4 - 009
    2000.      ,   & ! SO2 - 010
    0.         ,   & ! SULF - 011
    0.         ,   & ! CO - 012
    -1000.     ,   & ! CO2 - 013
    0.         ,   & ! N2 - 014
    0.         ,   & ! O2 - 015
    0.         ,   & ! H2O - 016
    0.         ,   & ! H2 - 017
    0.         ,   & ! O3P - 018
    0.         ,   & ! O1D - 019
    0.         ,   & ! HO - 020
    0.         ,   & ! HO2 - 021
    0.         ,   & ! CH4 - 022
    0.         ,   & ! ETH - 023
    0.         ,   & ! ALKA - 024
    0.         ,   & ! ALKE - 025
    0.         ,   & ! BIO - 026
    0.         ,   & ! ARO - 027
    0.         ,   & ! HCHO - 028
    0.         ,   & ! ALD - 029
    0.         ,   & ! KET - 030
    0.         ,   & ! CRBO - 031
    0.         ,   & ! ONIT - 032
    0.         ,   & ! PAN - 033
    0.         ,   & ! OP1 - 034
    0.         ,   & ! OP2 - 035
    -1500.     ,   & ! ORA1 - 036
    0.         ,   & ! ORA2 - 037
    0.         ,   & ! MO2 - 038
    0.         ,   & ! AKAP - 039
    0.         ,   & ! AKEP - 040
    0.         ,   & ! BIOP - 041
    0.         ,   & ! PHO - 042
    0.         ,   & ! ADD - 043
    0.         ,   & ! AROP - 044
    0.         ,   & ! CBOP - 045
    0.         ,   & ! OLN - 046
    0.         ,   & ! XO2 - 047
    0.         ,   & ! DMS - 048
    0.         ,   & ! NH3 - 049
    0.         ,   & ! NH3 - 049
    0.             & ! CFC - 050
    /)


!EOP
!-------------------------------------------------------------------------
!  Begin
       found = 0
loop2: DO l = 1,nspecies_HL
        IF(TRIM(spc_name(l)) == TRIM(name)) then
          c1  = hstar(l)
          c2  =   dhr(l)
          c3  =   ak0(l)
          c4  =   dak(l)
          found = 1
          EXIT loop2
        ENDIF
       enddo loop2
       IF(found == 0) then
          c1  = notfound
          c2  = notfound
          c3  = notfound
          c4  = notfound
       ENDIF

       __RETURN__(__SUCCESS__)
   end subroutine get_HenrysLawCts

!==================================================================================
!BOP
! !IROUTINE: NIthermo

   subroutine NIthermo (km, klid, cdt, grav, delp, rhoa, tmpu, rh, fMassHNO3, fMassAir, &
                        SO4, NH3, NO3an1, NH4a, xhno3, &
                        NI_pno3aq, NI_pnh4aq, NI_pnh3aq, rc)


! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km    ! total model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt   ! model time step [s]
   real, intent(in)    :: grav  ! gravity [m s-2]
   real, dimension(:,:,:), intent(in)  :: delp  ! pressure thickness [Pa]
   real, dimension(:,:,:), intent(in)  :: rhoa   ! Layer air density [kg m-3]
   real, dimension(:,:,:), intent(in)  :: tmpu   ! Layer temperature [K]
   real, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [0-1]
   real, intent(in)  :: fMassHNO3   ! gram molecular weight of HNO3
   real, intent(in)  :: fMassAir    ! gram molecular weight of air

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout)  :: SO4    ! Sulphate aerosol [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: NH3    ! Ammonia (NH3, gas phase) [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: NO3an1 ! Nitrate size bin 001 [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: NH4a   ! Ammonium ion (NH4+, aerosol phase) [kg kg-1]
   real, dimension(:,:,:), intent(inout)  :: xhno3  ! buffer for NITRATE_HNO3 [kg m-2 sec-1]
   real, pointer, dimension(:,:), intent(inout) :: NI_pno3aq ! Nitrate Production from Aqueous Chemistry [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout) :: NI_pnh4aq ! Ammonium Production from Aqueous Chemistry [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout) :: NI_pnh3aq ! Ammonia Change from Aqueous Chemistry [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc                   ! Error return code:


! !DESCRIPTION: Prepares variables and calls the RPMARES (thermodynamics module)
!
! !REVISION HISTORY:
!
! Aug2020 E.Sherman - Refactored for process library
!

! !Local Variables
   real   :: fmmr_to_conc
   real(kind=DP) :: SO4_, GNO3, GNH3, RH_, TEMP, ASO4, AHSO4, AH2O, ANO3, ANH4
   integer :: k, j, i

   integer :: status

!EOP
!-------------------------------------------------------------------------
!  Begin...

   do k = klid, km
    do j = 1, ubound(tmpu,2)
     do i = 1, ubound(tmpu,1)

!     Conversion of mass mixing ratio to concentration (ug m-3)
      fmmr_to_conc = 1.e9 * rhoa(i,j,k)

!     Unit conversion for input to thermodynamic module
!     Per grid box call to RPMARES thermodynamic module
!     We do not presently treat chemistry of sulfate completely,
!     hence we ignore terms for ASO4, AHSO4, AH2O, and we do
!     not update SO4 on output from RPMARES.
!     At present we are importing HNO3 from offline file, so we
!     do not update on return.
      SO4_  = 1.d-32
      SO4_  = max(1.d-32,SO4(i,j,k) * fmmr_to_conc)
      GNO3  = max(1.d-32,xhno3(i,j,k) * fMassHNO3 / fMassAir * fmmr_to_conc)
      GNH3  = max(1.d-32,NH3(i,j,k)  * fmmr_to_conc)
      RH_   = rh(i,j,k)
      TEMP  = tmpu(i,j,k)
      ASO4  = 1.d-32
      AHSO4 = 1.d-32
      ANO3  = max(1.d-32,NO3an1(i,j,k) * fmmr_to_conc)
      AH2O  = 1.d-32
      ANH4  = max(1.d-32,NH4a(i,j,k) * fmmr_to_conc)

      call RPMARES (  SO4_, GNO3,  GNH3, RH_,  TEMP, &
                      ASO4, AHSO4, ANO3, AH2O, ANH4, __RC__ )

!     Diagnostic terms
      if(associated(NI_pno3aq)) &
       NI_pno3aq(i,j) = NI_pno3aq(i,j) &
        + (ANO3 / fmmr_to_conc - NO3an1(i,j,k)) &
          * delp(i,j,k)/grav/cdt
      if(associated(NI_pnh4aq)) &
       NI_pnh4aq(i,j) = NI_pnh4aq(i,j) &
        + (ANH4 / fmmr_to_conc - NH4a(i,j,k)) &
          * delp(i,j,k)/grav/cdt
      if(associated(NI_pnh3aq)) &
       NI_pnh3aq(i,j) = NI_pnh3aq(i,j) &
        + (GNH3 / fmmr_to_conc - NH3(i,j,k)) &
          * delp(i,j,k)/grav/cdt

!     Unit conversion back on return from thermodynamic module
      NH3(i,j,k)    = GNH3 / fmmr_to_conc
      NO3an1(i,j,k) = ANO3 / fmmr_to_conc
      NH4a(i,j,k)   = ANH4 / fmmr_to_conc
      xhno3(i,j,k) = max(1.d-32, GNO3 * fMassAir / fMassHNO3 / fmmr_to_conc)

     enddo
    enddo
   enddo


   __RETURN__(__SUCCESS__)
   end subroutine NIthermo

!==================================================================================
!BOP
! !IROUTINE: RPMARES

   subroutine RPMARES( SO4,  GNO3,  GNH3, RH,   TEMP, &
                       ASO4, AHSO4, ANO3, AH2O, ANH4, rc )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   real(kind=DP) :: SO4              ! Total sulfate in micrograms / m**3
   real(kind=DP) :: GNO3             ! Gas-phase nitric acid in micrograms / m**3
   real(kind=DP) :: GNH3             ! Gas-phase ammonia in micrograms / m**3
   real(kind=DP) :: RH               ! Fractional relative humidity
   real(kind=DP) :: TEMP             ! Temperature in Kelvin
   real(kind=DP) :: ASO4             ! Aerosol sulfate in micrograms / m**3
   real(kind=DP) :: AHSO4            ! Aerosol bisulfate in micrograms / m**3
   real(kind=DP) :: ANO3             ! Aerosol nitrate in micrograms / m**3
   real(kind=DP) :: AH2O             ! Aerosol liquid water content water in
                                     !   micrograms / m**3
   real(kind=DP) :: ANH4             ! Aerosol ammonium in micrograms / m**3

! !OUTPUT PARAMETERS:

   integer, intent(out) :: rc


! !DESCRIPTION:
!   ARES calculates the chemical composition of a sulfate/nitrate/
!   ammonium/water aerosol based on equilibrium thermodynamics.
!
!   This code considers two regimes depending upon the molar ratio
!   of ammonium to sulfate.
!
!   For values of this ratio less than 2,the code solves a cubic for
!   hydrogen ion molality, H+,  and if enough ammonium and liquid
!   water are present calculates the dissolved nitric acid. For molal
!   ionic strengths greater than 50, nitrate is assumed not to be present.
!
!   For values of the molar ratio of 2 or greater, all sulfate is assumed
!   to be ammonium sulfate and a calculation is made for the presence of
!   ammonium nitrate.
!
!   The Pitzer multicomponent approach is used in subroutine ACTCOF to
!   obtain the activity coefficients. Abandoned -7/30/97 FSB
!
!   The Bromley method of calculating the multicomponent activity coefficients
!    is used in this version 7/30/97 SJR/FSB
!
!   The calculation of liquid water
!   is done in subroutine water. Details for both calculations are given
!   in the respective subroutines.
!
!   Based upon MARS due to
!   P. Saxena, A.B. Hudischewskyj, C. Seigneur, and J.H. Seinfeld,
!   Atmos. Environ., vol. 20, Number 7, Pages 1471-1483, 1986.
!
!   and SCAPE due to
!   Kim, Seinfeld, and Saxeena, Aerosol Sience and Technology,
!   Vol 19, number 2, pages 157-181 and pages 182-198, 1993.
!
! NOTE: All concentrations supplied to this subroutine are TOTAL
!       over gas and aerosol phases

!
! !REVISION HISTORY:
!
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   11/10/87  Received the first version of the MARS code
!   S.Roselle   12/30/87  Restructured code
!   S.Roselle   2/12/88   Made correction to compute liquid-phase
!                         concentration of H2O2.
!   S.Roselle   5/26/88   Made correction as advised by SAI, for
!                         computing H+ concentration.
!   S.Roselle   3/1/89    Modified to operate with EM2
!   S.Roselle   5/19/89   Changed the maximum ionic strength from
!                         100 to 20, for numerical stability.
!   F.Binkowski 3/3/91    Incorporate new method for ammonia rich case
!                         using equations for nitrate budget.
!   F.Binkowski 6/18/91   New ammonia poor case which
!                         omits letovicite.
!   F.Binkowski 7/25/91   Rearranged entire code, restructured
!                         ammonia poor case.
!   F.Binkowski 9/9/91    Reconciled all cases of ASO4 to be output
!                         as SO4--
!   F.Binkowski 12/6/91   Changed the ammonia defficient case so that
!                         there is only neutralized sulfate (ammonium
!                         sulfate) and sulfuric acid.
!   F.Binkowski 3/5/92    Set RH bound on AWAS to 37 % to be in agreement
!                          with the Cohen et al. (1987)  maximum molality
!                          of 36.2 in Table III.( J. Phys Chem (91) page
!                          4569, and Table IV p 4587.)
!   F.Binkowski 3/9/92    Redid logic for ammonia defficient case to remove
!                         possibility for denomenator becoming zero;
!                         this involved solving for H+ first.
!                         Note that for a relative humidity
!                          less than 50%, the model assumes that there is no
!                          aerosol nitrate.
!   F.Binkowski 4/17/95   Code renamed  ARES (AeRosol Equilibrium System)
!                          Redid logic as follows
!                         1. Water algorithm now follows Spann & Richardson
!                         2. Pitzer Multicomponent method used
!                         3. Multicomponent practical osmotic coefficient
!                            use to close iterations.
!                         4. The model now assumes that for a water
!                            mass fraction WFRAC less than 50% there is
!                            no aerosol nitrate.
!   F.Binkowski 7/20/95   Changed how nitrate is calculated in ammonia poor
!                         case, and changed the WFRAC criterion to 40%.
!                         For ammonium to sulfate ratio less than 1.0
!                         all ammonium is aerosol and no nitrate aerosol
!                         exists.
!   F.Binkowski 7/21/95   Changed ammonia-ammonium in ammonia poor case to
!                         allow gas-phase ammonia to exist.
!   F.Binkowski 7/26/95   Changed equilibrium constants to values from
!                         Kim et al. (1993)
!   F.Binkowski 6/27/96   Changed to new water format
!   F.Binkowski 7/30/97   Changed to Bromley method for multicomponent
!                         activity coefficients. The binary activity
!                         coefficients
!                         are the same as the previous version
!   F.Binkowski 8/1/97    Changed minimum sulfate from 0.0 to 1.0e-6 i.e.
!                         1 picogram per cubic meter
!   F.Binkowski 2/23/98   Changes to code made by Ingmar Ackermann to
!                         deal with precision problems on workstations
!                         incorporated in to this version.  Also included
!                         are his improved descriptions of variables.
!  F. Binkowski 8/28/98   changed logic as follows:
!                         If iterations fail, initial values of nitrate
!                          are retained.
!                         Total mass budgets are changed to account for gas
!                         phase returned.
!  F.Binkowski 10/01/98   Removed setting RATIO to 5 for low to
!                         to zero sulfate sulfate case.
!  F.Binkowski 01/10/2000 reconcile versions
!
!  F.Binkowski 05/17/2000 change to logic for calculating RATIO
!  F.Binkowski 04/09/2001 change for very low values of RATIO,
!                         RATIO < 0.5, no iterative calculations are done
!                         in low ammonia case a MAX(1.0e-10, MSO4) IS
!                         applied, and the iteration count is
!                         reduced to fifty for each iteration loop.
!  R. Yantosca 09/25/2002 Bundled into "rpmares_mod.f".  Declared all REALs
!                          as REAL*8's.  Cleaned up comments.  Also now force
!                          double precision explicitly with "D" exponents.
!  P. Le Sager and        Bug fix for low ammonia case -- prevent floating
!  R. Yantosca 04/10/2008  point underflow and NaN's.
!  S. Steenrod 04/15/2010 Modified to include into GMI model
!  E. Sherman  08/06/2020 Moved to GOCART2G process library

! !Local Variables
  !=================================================================
  ! PARAMETERS and their descriptions:
  !=================================================================
  ! Molecular weights
   real(kind=DP), PARAMETER :: MWNO3  = 62.0049d0                ! NO3
   real(kind=DP), PARAMETER :: MWHNO3 = 63.01287d0               ! HNO3
   real(kind=DP), PARAMETER :: MWSO4  = 96.0576d0                ! SO4
   real(kind=DP), PARAMETER :: MWNH3  = 17.03061d0               ! NH3
   real(kind=DP), PARAMETER :: MWNH4  = 18.03858d0               ! NH4

   ! Minimum value of sulfate aerosol concentration
   real(kind=DP), PARAMETER :: MINSO4 = 1.0d-6 / MWSO4

   ! Minimum total nitrate cncentration
   real(kind=DP), PARAMETER :: MINNO3 = 1.0d-6 / MWNO3

   ! Force a minimum concentration
   real(kind=DP), PARAMETER :: FLOOR  = 1.0d-30

   ! Tolerances for convergence test.  NOTE: We now have made these
   ! parameters so they don't lose their values (phs, bmy, 4/10/08)
   real(kind=DP), PARAMETER :: TOLER1 = 0.00001d0
   real(kind=DP), PARAMETER :: TOLER2 = 0.001d0

   ! Limit to test for zero ionic activity (phs, bmy, 4/10/08)
   real(kind=DP), PARAMETER :: EPS    = 1.0d-30

   !=================================================================
   ! SCRATCH LOCAL VARIABLES and their descriptions:
   !=================================================================

   INTEGER :: IRH              ! Index set to percent relative humidity
   INTEGER :: NITR             ! Number of iterations for activity
                               !   coefficients
   INTEGER :: NNN              ! Loop index for iterations
   INTEGER :: NR               ! Number of roots to cubic equation for
                               ! H+ ciaprecision
   real(kind=DP)  :: A0        ! Coefficients and roots of
   real(kind=DP)  :: A1        ! Coefficients and roots of
   real(kind=DP)  :: A2        ! Coefficients and roots of
   REAL    :: AA               ! Coefficients and discriminant for
                               ! quadratic equation for ammonium nitrate
   real(kind=DP)  :: BAL       ! internal variables ( high ammonia case)
   real(kind=DP)  :: BB        ! Coefficients and discriminant for
                               !   quadratic equation for ammonium nitrate
   real(kind=DP)  :: BHAT      ! Variables used for ammonia solubility
   real(kind=DP)  :: CC        ! Coefficients and discriminant for
                               !   quadratic equation for ammonium nitrate
   real(kind=DP)  :: CONVT     ! Factor for conversion of units
   real(kind=DP)  :: DD        ! Coefficients and discriminant for
                               !   quadratic equation for ammonium nitrate
   real(kind=DP)  :: DISC      ! Coefficients and discriminant for
                               !   quadratic equation for ammonium nitrate
   real(kind=DP)  :: EROR      ! Relative error used for convergence test
   real(kind=DP)  :: FNH3      ! "Free ammonia concentration", that
                               !   which exceeds TWOSO4
   real(kind=DP)  :: GAMAAB    ! Activity Coefficient for (NH4+,
                               !   HSO4-)GAMS( 2,3 )
   real(kind=DP)  :: GAMAAN    ! Activity coefficient for (NH4+, NO3-)
                               !   GAMS( 2,2 )
   real(kind=DP)  :: GAMAHAT   ! Variables used for ammonia solubility
   real(kind=DP)  :: GAMANA    ! Activity coefficient for (H+ ,NO3-)
                               !   GAMS( 1,2 )
   real(kind=DP)  :: GAMAS1    ! Activity coefficient for (2H+, SO4--)
                               !   GAMS( 1,1 )
   real(kind=DP)  :: GAMAS2    ! Activity coefficient for (H+, HSO4-)
                               !   GAMS( 1,3 )
   real(kind=DP)  :: GAMOLD    ! used for convergence of iteration
   real(kind=DP)  :: GASQD     ! internal variables ( high ammonia case)
   real(kind=DP)  :: HPLUS     ! Hydrogen ion (low ammonia case) (moles
                               !   / kg water)
   real(kind=DP)  :: K1A       ! Equilibrium constant for ammonia to
                               !   ammonium
   real(kind=DP)  :: K2SA      ! Equilibrium constant for
                               !   sulfate-bisulfate (aqueous)
   real(kind=DP)  :: K3        ! Dissociation constant for ammonium
                               !   nitrate
   real(kind=DP)  :: KAN       ! Equilibrium constant for ammonium
                               !   nitrate (aqueous)
   real(kind=DP)  :: KHAT      ! Variables used for ammonia solubility
   real(kind=DP)  :: KNA       ! Equilibrium constant for nitric acid
                               !   (aqueous)
   real(kind=DP)  :: KPH       ! Henry's Law Constant for ammonia
   real(kind=DP)  :: KW        ! Equilibrium constant for water
                               !  dissociation
   real(kind=DP)  :: KW2       ! Internal variable using KAN
   real(kind=DP)  :: MAN       ! Nitrate (high ammonia case) (moles /
                               !   kg water)
   real(kind=DP)  :: MAS       ! Sulfate (high ammonia case) (moles /
                               !   kg water)
   real(kind=DP)  :: MHSO4     ! Bisulfate (low ammonia case) (moles /
                               !   kg water)
   real(kind=DP)  :: MNA       ! Nitrate (low ammonia case) (moles / kg
                               !   water)
   real(kind=DP)  :: MNH4      ! Ammonium (moles / kg water)
   real(kind=DP)  :: MOLNU     ! Total number of moles of all ions
   real(kind=DP)  :: MSO4      ! Sulfate (low ammonia case) (moles / kg
                               !   water)
   real(kind=DP)  :: PHIBAR    ! Practical osmotic coefficient
   real(kind=DP)  :: PHIOLD    ! Previous value of practical osmotic
                               !   coefficient used for convergence of
                               !   iteration
   real(kind=DP)  :: RATIO     ! Molar ratio of ammonium to sulfate
   real(kind=DP)  :: RK2SA     ! Internal variable using K2SA
   real(kind=DP)  :: RKNA      ! Internal variables using KNA
   real(kind=DP)  :: RKNWET    ! Internal variables using KNA
   real(kind=DP)  :: RR1
   real(kind=DP)  :: RR2
   real(kind=DP)  :: STION     ! Ionic strength
   real(kind=DP)  :: T1        ! Internal variables for temperature
                               !   corrections
   real(kind=DP)  :: T2        ! Internal variables for temperature
                               !   corrections
   real(kind=DP)  :: T21       ! Internal variables of convenience (low
                               !   ammonia case)
   real(kind=DP)  :: T221      ! Internal variables of convenience (low
                               !   ammonia case)
   real(kind=DP)  :: T3        ! Internal variables for temperature
                               !   corrections
   real(kind=DP)  :: T4        ! Internal variables for temperature
                               !   corrections
   real(kind=DP)  :: T6        ! Internal variables for temperature
                               !   corrections
   real(kind=DP)  :: TNH4      ! Total ammonia and ammonium in
                               !   micromoles / meter ** 3
   real(kind=DP)  :: TNO3      ! Total nitrate in micromoles / meter ** 3
   !-----------------------------------------------------------------------
   ! Prior to 4/10/08:
   ! Now make these PARAMETERS instead of variables (bmy, 4/10/08)
   !real(kind=DP)  :: TOLER1           ! Tolerances for convergence test
   !real(kind=DP)  :: TOLER2           ! Tolerances for convergence test
   !-----------------------------------------------------------------------
   real(kind=DP)  :: TSO4      ! Total sulfate in micromoles / meter ** 3
   real(kind=DP)  :: TWOSO4    ! 2.0 * TSO4  (high ammonia case) (moles
                               !   / kg water)
   real(kind=DP)  :: WFRAC     ! Water mass fraction
   real(kind=DP)  :: WH2O      ! Aerosol liquid water content (internally)
                               !   micrograms / meter **3 on output
                               !   internally it is 10 ** (-6) kg (water)
                               !   / meter ** 3
                               !   the conversion factor (1000 g = 1 kg)
                               !   is applied for AH2O output
   real(kind=DP)  :: WSQD      ! internal variables ( high ammonia case)
   real(kind=DP)  :: XNO3      ! Nitrate aerosol concentration in
                               ! micromoles / meter ** 3
   real(kind=DP)  :: XXQ       ! Variable used in quadratic solution
   real(kind=DP)  :: YNH4      ! Ammonium aerosol concentration in
                               !  micromoles / meter** 3
   real(kind=DP)  :: ZH2O      ! Water variable saved in case ionic
                               !  strength too high.
   real(kind=DP)  :: ZSO4      ! Total sulfate molality - mso4 + mhso4
                               !  (low ammonia case) (moles / kg water)
   real(kind=DP)  :: CAT( 2 )  ! Array for cations (1, H+); (2, NH4+)
                               !  (moles / kg water)
   real(kind=DP)  :: AN ( 3 )  ! Array for anions (1, SO4--); (2,
                               !   NO3-); (3, HSO4-)  (moles / kg water)
   real(kind=DP)  :: CRUTES( 3 )      ! Coefficients and roots of
   real(kind=DP)  :: GAMS( 2, 3 )     ! Array of activity coefficients
   real(kind=DP)  :: TMASSHNO3        ! Total nitrate (vapor and particle)
   real(kind=DP)  :: GNO3_IN, ANO3_IN
   character (len=75) :: err_msg

   integer :: status

!EOP
!-------------------------------------------------------------------------
!  Begin...

      rc = __SUCCESS__

      ! For extremely low relative humidity ( less than 1% ) set the
      ! water content to a minimum and skip the calculation.
      IF ( RH .LT. 0.01 ) THEN
         AH2O = FLOOR
         RETURN
      ENDIF

      ! total sulfate concentration
      TSO4 = MAX( FLOOR, SO4 / MWSO4  )
      ASO4 = SO4

      !Cia models3 merge NH3/NH4 , HNO3,NO3 here
      !c *** recommended by Dr. Ingmar Ackermann

      ! total nitrate
      TNO3      = MAX( 0.0d0, ( ANO3 / MWNO3 + GNO3 / MWHNO3 ) )

      ! total ammonia
      TNH4      = MAX( 0.0d0, ( GNH3 / MWNH3 + ANH4 / MWNH4 )  )

      GNO3_IN   = GNO3
      ANO3_IN   = ANO3
      TMASSHNO3 = MAX( 0.0d0, GNO3 + ANO3 )

      ! set the  molar ratio of ammonium to sulfate
      RATIO = TNH4 / TSO4

      ! validity check for negative concentration
      IF ( TSO4 < 0.0d0 .OR. TNO3 < 0.0d0 .OR. TNH4 < 0.0d0 ) THEN
          !$omp critical (G2G_proc_7)
          PRINT*, 'TSO4 : ', TSO4
          PRINT*, 'TNO3 : ', TNO3
          PRINT*, 'TNH4 : ', TNH4
          !$omp end critical (G2G_proc_7)


!.sds          CALL GEOS_CHEM_STOP
          err_msg = 'negative concen problem in RPMARES - TSO4, TNO3, TNH4:'
          call PrintError  &
     &      (err_msg, .true., 0, 0, 0, 2, TSO4, TNO3, __RC_NO_OPT__)
!     &      (err_msg, .true., 0, 0, 0, 2, TSO4, TNO3, rc=status)
      ENDIF

      ! now set humidity index IRH as a percent
      IRH = NINT( 100.0 * RH )

      ! now set humidity index IRH as a percent
      IRH = MAX(  1, IRH )
      IRH = MIN( 99, IRH )

      !=================================================================
      ! Specify the equilibrium constants at  correct temperature.
      ! Also change units from ATM to MICROMOLE/M**3 (for KAN, KPH, and K3 )
      ! Values from Kim et al. (1993) except as noted.
      ! Equilibrium constant in Kim et al. (1993)
      !   K = K0 exp[ a(T0/T -1) + b(1+log(T0/T)-T0/T) ], T0 = 298.15 K
      !   K = K0 EXP[ a T3 + b T4 ] in the code here.
      !=================================================================
      CONVT = 1.0d0 / ( 0.082d0 * TEMP )
      T6    = 0.082d-9 *  TEMP
      T1    = 298.0d0 / TEMP
      T2    = LOG( T1 )
      T3    = T1 - 1.0d0
      T4    = 1.0d0 + T2 - T1

      !=================================================================
      ! Equilibrium Relation
      !
      ! HSO4-(aq)         = H+(aq)   + SO4--(aq)  ; K2SA
      ! NH3(g)            = NH3(aq)               ; KPH
      ! NH3(aq) + H2O(aq) = NH4+(aq) + OH-(aq)    ; K1A
      ! HNO3(g)           = H+(aq)   + NO3-(aq)   ; KNA
      ! NH3(g) + HNO3(g)  = NH4NO3(s)             ; K3
      ! H2O(aq)           = H+(aq)   + OH-(aq)    ; KW
      !=================================================================
      KNA  = 2.511d+06 *  EXP(  29.17d0 * T3 + 16.83d0 * T4 ) * T6
      K1A  = 1.805d-05 *  EXP(  -1.50d0 * T3 + 26.92d0 * T4 )
      K2SA = 1.015d-02 *  EXP(   8.85d0 * T3 + 25.14d0 * T4 )
      KW   = 1.010d-14 *  EXP( -22.52d0 * T3 + 26.92d0 * T4 )
      KPH  = 57.639d0  *  EXP(  13.79d0 * T3 -  5.39d0 * T4 ) * T6
      !K3   =  5.746E-17 * EXP( -74.38 * T3 + 6.12  * T4 ) * T6 * T6
      KHAT =  KPH * K1A / KW
      KAN  =  KNA * KHAT

      ! Compute temperature dependent equilibrium constant for NH4NO3
      ! (from Mozurkewich, 1993)
      K3 = EXP( 118.87d0  - 24084.0d0 / TEMP -  6.025d0  * LOG( TEMP ) )

      ! Convert to (micromoles/m**3) **2
      K3     = K3 * CONVT * CONVT

      WH2O   = 0.0d0
      STION  = 0.0d0
!.sds      AH2O   = 0.0d0
      AH2O   = FLOOR

      MAS    = 0.0d0
      MAN    = 0.0d0
      HPLUS  = 0.0d0
      !--------------------------------------------------------------
      ! Prior to 4/10/08:
      ! Now make these parameters so that they won't lose their
      ! values. (phs, bmy, 4/10/08)
      !TOLER1 = 0.00001d0
      !TOLER2 = 0.001d0
      !--------------------------------------------------------------
      NITR   = 0
      NR     = 0
      GAMAAN = 1.0d0
      GAMOLD = 1.0d0

      ! If there is very little sulfate and  nitrate
      ! set concentrations to a very small value and return.
      IF ( ( TSO4 .LT. MINSO4 ) .AND. ( TNO3 .LT. MINNO3 ) ) THEN
         ASO4  = MAX( FLOOR, ASO4  )
         AHSO4 = MAX( FLOOR, AHSO4 ) ! [rjp, 12/12/01]
         ANO3  = MAX( FLOOR, ANO3  )
         ANH4  = MAX( FLOOR, ANH4  )
         WH2O  = FLOOR
         AH2O  = FLOOR
         GNH3  = MAX( FLOOR, GNH3  )
         GNO3  = MAX( FLOOR, GNO3  )

         RETURN
      ENDIF
      !=================================================================
      ! High Ammonia Case
      !=================================================================
      IF ( RATIO .GT. 2.0d0 ) THEN

         GAMAAN = 0.1d0

         ! Set up twice the sulfate for future use.
         TWOSO4 = 2.0d0 * TSO4
         XNO3   = 0.0d0
         YNH4   = TWOSO4

         ! Treat different regimes of relative humidity
         !
         ! ZSR relationship is used to set water levels. Units are
         !  10**(-6) kg water/ (cubic meter of air)
         !  start with ammomium sulfate solution without nitrate

         CALL AWATER( IRH, TSO4, YNH4, TNO3, AH2O ) !**** note TNO3
         WH2O = 1.0d-3 * AH2O

         ASO4 = TSO4   * MWSO4
         ! In sulfate poor case, Sulfate ion is preferred
         ! Set bisulfate equal to zero [rjp, 12/12/01]
         AHSO4 = 0.0d0
         ANO3  = 0.0d0
         ANH4  = YNH4 * MWNH4
         WFRAC = AH2O / ( ASO4 + ANH4 +  AH2O )

        !IF ( WFRAC .EQ. 0.0 )  RETURN   ! No water
        IF ( WFRAC .LT. 0.2d0 ) THEN

           ! "dry" ammonium sulfate and ammonium nitrate
           ! compute free ammonia
           FNH3 = TNH4 - TWOSO4
           CC   = TNO3 * FNH3 - K3

           ! check for not enough to support aerosol
           IF ( CC .LE. 0.0d0 ) THEN
              XNO3 = 0.0d0
           ELSE
              AA   = 1.0d0
              BB   = -( TNO3 + FNH3 )
              DISC = BB * BB - 4.0d0 * CC

              ! Check for complex roots of the quadratic
              ! set retain initial values of nitrate and RETURN
              ! if complex roots are found
              IF ( DISC .LT. 0.0d0 ) THEN
                 XNO3  = 0.0d0
                 AH2O  = 1000.0d0 * WH2O
                 YNH4  = TWOSO4
                 ASO4  = TSO4 * MWSO4
                 AHSO4 = 0.0d0
                 ANH4  = YNH4 * MWNH4
                 GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
                 GNO3  = GNO3_IN
                 ANO3  = ANO3_IN
                 RETURN
              ENDIF

              ! to get here, BB .lt. 0.0, CC .gt. 0.0 always
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN ( 1.0d0, BB ) * DD )


              ! Since both roots are positive, select smaller root.
              XNO3 = MIN( XXQ / AA, CC / XXQ )

           ENDIF                ! CC .LE. 0.0

           AH2O  = 1000.0d0 * WH2O
           YNH4  = TWOSO4 + XNO3
           ASO4  = TSO4 * MWSO4
           AHSO4 = FLOOR
           ANO3  = XNO3 * MWNO3
           ANH4  = YNH4 * MWNH4
           GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 )  )
           GNO3  = MAX( FLOOR, ( TMASSHNO3 - ANO3 ) )
           RETURN
        ENDIF                  ! WFRAC .LT. 0.2

        ! liquid phase containing completely neutralized sulfate and
        ! some nitrate.  Solve for composition and quantity.
        MAS    = TSO4 / WH2O
        MAN    = 0.0d0
        XNO3   = 0.0d0
        YNH4   = TWOSO4
        PHIOLD = 1.0d0

        !===============================================================
        ! Start loop for iteration
        !
        ! The assumption here is that all sulfate is ammonium sulfate,
        ! and is supersaturated at lower relative humidities.
        !===============================================================
        DO NNN = 1, 50 ! loop count reduced 0409/2001 by FSB

           NITR  = NNN
           GASQD = GAMAAN * GAMAAN
           WSQD  = WH2O * WH2O
           KW2   = KAN * WSQD / GASQD
           AA    = 1.0 - KW2
           BB    = TWOSO4 + KW2 * ( TNO3 + TNH4 - TWOSO4 )
           CC    = -KW2 * TNO3 * ( TNH4 - TWOSO4 )

           ! This is a quadratic for XNO3 [MICROMOLES / M**3]
           ! of nitrate in solution
           DISC = BB * BB - 4.0d0 * AA * CC

           ! Check for complex roots, retain inital values and RETURN
           IF ( DISC .LT. 0.0 ) THEN
              XNO3  = 0.0d0
              AH2O  = 1000.0d0 * WH2O
              YNH4  = TWOSO4
              ASO4  = TSO4 * MWSO4
              AHSO4 = FLOOR     ! [rjp, 12/12/01]
              ANH4  = YNH4 * MWNH4
              GNH3  = MWNH3 * MAX( FLOOR, (TNH4 - YNH4 ) )
              GNO3  = GNO3_IN
              ANO3  = ANO3_IN
              RETURN
           ENDIF

           ! Deal with degenerate case (yoj)
           IF ( AA .NE. 0.0d0 ) THEN
              DD  = SQRT( DISC )
              XXQ = -0.5d0 * ( BB + SIGN( 1.0d0, BB ) * DD )
              RR1 = XXQ / AA
              RR2 = CC / XXQ

              ! choose minimum positve root
              IF ( ( RR1 * RR2 ) .LT. 0.0d0 ) THEN
                 XNO3 = MAX( RR1, RR2 )
              ELSE
                 XNO3 = MIN( RR1, RR2 )
              ENDIF
           ELSE
              XNO3 = - CC / BB  ! AA equals zero here.
           ENDIF

           XNO3 = MIN( XNO3, TNO3 )

           ! This version assumes no solid sulfate forms (supersaturated )
           ! Now update water
           CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

           ! ZSR relationship is used to set water levels. Units are
           ! 10**(-6) kg water/ (cubic meter of air).  The conversion
           ! from micromoles to moles is done by the units of WH2O.
           WH2O = 1.0d-3 * AH2O

           ! Ionic balance determines the ammonium in solution.
           MAN  = XNO3 / WH2O
           MAS  = TSO4 / WH2O
           MNH4 = 2.0d0 * MAS + MAN
           YNH4 = MNH4 * WH2O

           ! MAS, MAN and MNH4 are the aqueous concentrations of sulfate,
           ! nitrate, and ammonium in molal units (moles/(kg water) ).
           STION    = 3.0d0 * MAS + MAN
           CAT( 1 ) = 0.0d0
           CAT( 2 ) = MNH4
           AN ( 1 ) = MAS
           AN ( 2 ) = MAN
           AN ( 3 ) = 0.0d0
!           CALL ACTCOF ( CAT, AN, GAMS, MOLNU, PHIBAR )
           CALL ACTCOF ( CAT, AN, GAMS )
           GAMAAN = GAMS( 2, 2 )

           ! Use GAMAAN for convergence control
           EROR   = ABS( GAMOLD - GAMAAN ) / GAMOLD
           GAMOLD = GAMAAN

           ! Check to see if we have a solution
           IF ( EROR .LE. TOLER1 ) THEN
              ASO4  = TSO4 * MWSO4
              AHSO4 = 0.0d0       ! [rjp, 12/12/01]
              ANO3  = XNO3 * MWNO3
              ANH4  = YNH4 * MWNH4
              GNO3  = MAX( FLOOR, ( TMASSHNO3  - ANO3 ) )
              GNH3  = MWNH3 * MAX( FLOOR, ( TNH4 - YNH4 ) )
              AH2O  = 1000.0d0 * WH2O
              RETURN
           ENDIF

        ENDDO

        ! If after NITR iterations no solution is found, then:
        ! FSB retain the initial values of nitrate particle and vapor
        ! note whether or not convert all bisulfate to sulfate
        ASO4  = TSO4 * MWSO4
        AHSO4 = FLOOR
        XNO3  = TNO3 / MWNO3
        YNH4  = TWOSO4
        ANH4  = YNH4 * MWNH4

        CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

        GNO3  = GNO3_IN
        ANO3  = ANO3_IN
        GNH3  = MAX( FLOOR, MWNH3 * (TNH4 - YNH4 ) )
        RETURN

      !================================================================
      ! Low Ammonia Case
      !
      ! Coded by Dr. Francis S. Binkowski 12/8/91.(4/26/95)
      ! modified 8/28/98
      ! modified 04/09/2001
      !
      ! All cases covered by this logic
      !=================================================================
      ELSE

         WH2O = 0.0d0
         CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )
         WH2O = 1.0d-3 * AH2O
         ZH2O = AH2O

         ! convert 10**(-6) kg water/(cubic meter of air) to micrograms
         ! of water per cubic meter of air (1000 g = 1 kg)
         ! in sulfate rich case, preferred form is HSO4-
         !ASO4 = TSO4 * MWSO4
         ASO4  = FLOOR          ![rjp, 12/12/01]
         AHSO4 = TSO4 * MWSO4   ![rjp, 12/12/01]
         ANH4  = TNH4 * MWNH4
         ANO3  = ANO3_IN
         GNO3  = TMASSHNO3 - ANO3
         GNH3  = FLOOR

         !==============================================================
         ! *** Examine special cases and return if necessary.
         !
         ! FSB For values of RATIO less than 0.5 do no further
         ! calculations.  The code will cycle and still predict the
         ! same amount of ASO4, ANH4, ANO3, AH2O so terminate early
         ! to swame computation
         !==============================================================
         IF ( RATIO .LT. 0.5d0 ) RETURN ! FSB 04/09/2001

         ! Check for zero water.
         IF ( WH2O .EQ. 0.0d0 ) RETURN
         ZSO4 = TSO4 / WH2O

         ! ZSO4 is the molality of total sulfate i.e. MSO4 + MHSO4
         ! do not solve for aerosol nitrate for total sulfate molality
         ! greater than 11.0 because the model parameters break down
         !### IF ( ZSO4 .GT. 11.0 ) THEN
         !IF ( ZSO4 .GT. 9.0 ) THEN ! 18 June 97
         !IF ( ZSO4 .GT. 9.d0 ) THEN ! H. Bian 24 June 2015
         IF ( ZSO4 .GT. 9.00 ) THEN ! H. Bian 24 June 2015
            RETURN
         ENDIF
         IF ( ZSO4 .GT. 0.1d0 .and. TEMP .le. 220.d0) THEN ! H. Bian 24 June 2015
            RETURN
         ENDIF

         ! *** Calculation may now proceed.
         !
         ! First solve with activity coeffs of 1.0, then iterate.
         PHIOLD = 1.0d0
         GAMANA = 1.0d0
         GAMAS1 = 1.0d0
         GAMAS2 = 1.0d0
         GAMAAB = 1.0d0
         GAMOLD = 1.0d0

         ! All ammonia is considered to be aerosol ammonium.
         MNH4 = TNH4 / WH2O

         ! MNH4 is the molality of ammonium ion.
         YNH4 = TNH4

         ! loop for iteration
         DO NNN = 1, 50    ! loop count reduced 04/09/2001 by FSB
            NITR = NNN

            ! set up equilibrium constants including activities
            ! solve the system for hplus first then sulfate & nitrate
            RK2SA  = K2SA * GAMAS2 * GAMAS2 / (GAMAS1 * GAMAS1 * GAMAS1)
            RKNA   = KNA / ( GAMANA * GAMANA )
            RKNWET = RKNA * WH2O
            T21    = ZSO4 - MNH4
            T221   = ZSO4 + T21

            ! set up coefficients for cubic
            A2 = RK2SA + RKNWET - T21
            A1 = RK2SA * RKNWET - T21 * ( RK2SA + RKNWET ) &
     &           - RK2SA * ZSO4 - RKNA * TNO3
            A0 = - (T21 * RK2SA * RKNWET &
     &           + RK2SA * RKNWET * ZSO4 + RK2SA * RKNA * TNO3 )

            CALL CUBIC ( A2, A1, A0, NR, CRUTES, __RC_NO_OPT__ )
!            CALL CUBIC ( A2, A1, A0, NR, CRUTES )

            ! Code assumes the smallest positive root is in CRUTES(1)
            HPLUS = CRUTES( 1 )
            BAL   = HPLUS **3 + A2 * HPLUS**2 + A1 * HPLUS + A0

            ! molality of sulfate ion
            MSO4  = RK2SA * ZSO4 / ( HPLUS + RK2SA )

            ! molality of bisulfate ion
            ! MAX added 04/09/2001 by FSB
            MHSO4 = MAX( 1.0d-10, ZSO4 - MSO4 )

            ! molality of nitrate ion
            MNA   = RKNA * TNO3 / ( HPLUS + RKNWET )
            MNA   = MAX( 0.0d0, MNA )
            MNA   = MIN( MNA, TNO3 / WH2O )
            XNO3  = MNA * WH2O
            ANO3  = MNA * WH2O * MWNO3
            GNO3  = MAX( FLOOR, TMASSHNO3 - ANO3 )
            ASO4  = MSO4 * WH2O * MWSO4 ![rjp, 12/12/01]
            AHSO4 = MHSO4 * WH2O * MWSO4 ![rjp, 12/12/01]

            ! Calculate ionic strength
            STION = 0.5d0 * ( HPLUS + MNA + MNH4 + MHSO4 + 4.0d0 * MSO4)
            ! Update water
            CALL AWATER ( IRH, TSO4, YNH4, XNO3, AH2O )

            ! Convert 10**(-6) kg water/(cubic meter of air) to micrograms
            ! of water per cubic meter of air (1000 g = 1 kg)
            WH2O     = 1.0d-3 * AH2O
            CAT( 1 ) = HPLUS
            CAT( 2 ) = MNH4
            AN ( 1 ) = MSO4
            AN ( 2 ) = MNA
            AN ( 3 ) = MHSO4

            CALL ACTCOF ( CAT, AN, GAMS )

            GAMANA = GAMS( 1, 2 )
            GAMAS1 = GAMS( 1, 1 )
            GAMAS2 = GAMS( 1, 3 )
            GAMAAN = GAMS( 2, 2 )

            !------------------------------------------------------------
            ! Add robustness: now check if GAMANA or GAMAS1 is too small
            ! for the division in RKNA and RK2SA. If they are, return w/
            ! original values: basically replicate the procedure used
            ! after the current DO-loop in case of no-convergence
            ! (phs, bmy, rjp, 4/10/08)
            !--------------------------------------------------------------
            IF ( ( ABS( GAMANA ) < EPS ) .OR. ( ABS( GAMAS1 ) < EPS ) ) THEN

               ! Reset to original values
               ANH4  = TNH4 * MWNH4
               GNH3  = FLOOR
               GNO3  = GNO3_IN
               ANO3  = ANO3_IN
               ASO4  = TSO4 * MWSO4
               AHSO4 = FLOOR

               ! Update water
               CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )

               ! Exit this subroutine
               RETURN
            ENDIF

            GAMAHAT = ( GAMAS2 * GAMAS2 / ( GAMAAB * GAMAAB ) )
            BHAT = KHAT * GAMAHAT
            !### EROR = ABS ( ( PHIOLD - PHIBAR ) / PHIOLD )
            !### PHIOLD = PHIBAR
            EROR = ABS ( GAMOLD - GAMAHAT ) / GAMOLD
            GAMOLD = GAMAHAT
            ! return with good solution
            IF ( EROR .LE. TOLER2 ) THEN
               RETURN
            ENDIF

         ENDDO

         ! after NITR iterations, failure to solve the system
         ! convert all ammonia to aerosol ammonium and return input
         ! values of NO3 and HNO3
         ANH4 = TNH4 * MWNH4
         GNH3 = FLOOR
         GNO3 = GNO3_IN
         ANO3 = ANO3_IN
         ASO4 = TSO4 * MWSO4    ! [rjp, 12/17/01]
         AHSO4= FLOOR           ! [rjp, 12/17/01]

         CALL AWATER ( IRH, TSO4, TNH4, TNO3, AH2O )

         RETURN

      ENDIF                     ! ratio .gt. 2.0

      ! Return to calling program

   end subroutine RPMARES

!------------------------------------------------------------------------------

      SUBROUTINE AWATER( IRHX, MSO4, MNH4, MNO3, WH2O )
!
!******************************************************************************
! NOTE!!! wh2o is returned in micrograms / cubic meter
!         mso4,mnh4,mno3 are in microMOLES / cubic meter
!
!  This  version uses polynomials rather than tables, and uses empirical
! polynomials for the mass fraction of solute (mfs) as a function of water
! activity
!   where:
!
!            mfs = ms / ( ms + mw)
!             ms is the mass of solute
!             mw is the mass of water.
!
!  Define y = mw/ ms
!
!  then  mfs = 1 / (1 + y)
!
!    y can then be obtained from the values of mfs as
!
!             y = (1 - mfs) / mfs
!
!
!     the aerosol is assumed to be in a metastable state if the rh is
!     is below the rh of deliquescence, but above the rh of crystallization.
!
!     ZSR interpolation is used for sulfates with x ( the molar ratio of
!     ammonium to sulfate in eh range 0 <= x <= 2, by sections.
!     section 1: 0 <= x < 1
!     section 2: 1 <= x < 1.5
!     section 3: 1.5 <= x < 2.0
!     section 4: 2 <= x
!     In sections 1 through 3, only the sulfates can affect the amount of water
!     on the particles.
!     In section 4, we have fully neutralized sulfate, and extra ammonium which
!     allows more nitrate to be present. Thus, the ammount of water is
!     calculated
!     using ZSR for ammonium sulfate and ammonium nitrate. Crystallization is
!     assumed to occur in sections 2,3,and 4. See detailed discussion below.
!
! definitions:
!     mso4, mnh4, and mno3 are the number of micromoles/(cubic meter of air)
!      for sulfate, ammonium, and nitrate respectively
!     irhx is the relative humidity (%)
!     wh2o is the returned water amount in micrograms / cubic meter of air
!     x is the molar ratio of ammonium to sulfate
!     y0,y1,y1.5, y2 are the water contents in mass of water/mass of solute
!     for pure aqueous solutions with x equal 1, 1.5, and 2 respectively.
!     y3 is the value of the mass ratio of water to solute for
!     a pure ammonium nitrate  solution.
!
!
!     coded by Dr. Francis S. Binkowski, 4/8/96.
!
! *** modified 05/30/2000
!     The use of two values of mfs at an ammonium to sulfate ratio
!     representative of ammonium sulfate led to an minor inconsistancy
!     in nitrate behavior as the ratio went from a value less than two
!     to a value greater than two and vice versa with either ammonium
!     held constant and sulfate changing, or sulfate held constant and
!     ammonium changing. the value of Chan et al. (1992) is the only value
!     now used.
!
! *** Modified 09/25/2002
!     Ported into "rpmares_mod.f".  Now declare all variables with REAL*8.
!     Also cleaned up comments and made cosmetic changes.  Force double
!     precision explicitly with "D" exponents.
!******************************************************************************
!
      ! Arguments
      INTEGER           :: IRHX
      REAL*8            :: MSO4, MNH4, MNO3, WH2O

      ! Local variables
      INTEGER           :: IRH
      REAL*8            :: TSO4,  TNH4,  TNO3,  X,      AW,     AWC
      REAL*8            :: MFS0,  MFS1,  MFS15, Y
      REAL*8            :: Y0,    Y1,    Y15,   Y2,     Y3,     Y40
      REAL*8            :: Y140,  Y1540, YC,    MFSSO4, MFSNO3

      ! Molecular weight parameters
      REAL*8, PARAMETER :: MWSO4  = 96.0636d0
      REAL*8, PARAMETER :: MWNH4  = 18.0985d0
      REAL*8, PARAMETER :: MWNO3  = 62.0649d0
      REAL*8, PARAMETER :: MW2    = MWSO4 + 2.0d0 * MWNH4
      REAL*8, PARAMETER :: MWANO3 = MWNO3 + MWNH4

      !=================================================================
      ! The polynomials use data for aw as a function of mfs from Tang
      ! and Munkelwitz, JGR 99: 18801-18808, 1994.  The polynomials were
      ! fit to Tang's values of water activity as a function of mfs.
      !
      ! *** coefficients of polynomials fit to Tang and Munkelwitz data
      !     now give mfs as a function of water activity.
      !=================================================================
      REAL*8 :: C1(4)  = (/ 0.9995178d0,  -0.7952896d0, &
     &                      0.99683673d0, -1.143874d0 /)

      REAL*8 :: C15(4) = (/ 1.697092d0, -4.045936d0, &
     &                      5.833688d0, -3.463783d0 /)

      !=================================================================
      ! The following coefficients are a fit to the data in Table 1 of
      !    Nair & Vohra, J. Aerosol Sci., 6: 265-271, 1975
      !      data c0/0.8258941, -1.899205, 3.296905, -2.214749 /
      !
      ! New data fit to data from
      !       Nair and Vohra J. Aerosol Sci., 6: 265-271, 1975
      !       Giaque et al. J.Am. Chem. Soc., 82: 62-70, 1960
      !       Zeleznik J. Phys. Chem. Ref. Data, 20: 157-1200
      !=================================================================
      REAL*8 :: C0(4)  =  (/ 0.798079d0, -1.574367d0, &
     &                       2.536686d0, -1.735297d0 /)

      !=================================================================
      ! Polynomials for ammonium nitrate and ammonium sulfate are from:
      ! Chan et al.1992, Atmospheric Environment (26A): 1661-1673.
      !=================================================================
      REAL*8 :: KNO3(6) = (/  0.2906d0,   6.83665d0, -26.9093d0, &
     &                       46.6983d0, -38.803d0,    11.8837d0 /)

      REAL*8 :: KSO4(6) = (/   2.27515d0, -11.147d0,   36.3369d0, &
     &                       -64.2134d0,   56.8341d0, -20.0953d0 /)

      !=================================================================
      ! AWATER begins here!
      !=================================================================

      ! Check range of per cent relative humidity
      IRH  = IRHX
      IRH  = MAX( 1, IRH )
      IRH  = MIN( IRH, 100 )

      ! Water activity = fractional relative humidity
      AW   = DBLE( IRH ) / 100.0d0
      TSO4 = MAX( MSO4 , 0.0d0 )
      TNH4 = MAX( MNH4 , 0.0d0 )
      TNO3 = MAX( MNO3 , 0.0d0 )
      X    = 0.0d0

      ! If there is non-zero sulfate calculate the molar ratio
      ! otherwise check for non-zero nitrate and ammonium
      IF ( TSO4 .GT. 0.0d0 ) THEN
         X = TNH4 / TSO4
      ELSE
         IF ( TNO3 .GT. 0.0d0 .AND. TNH4 .GT. 0.0d0 ) X = 10.0d0
      ENDIF

      ! *** begin screen on x for calculating wh2o
      IF ( X .LT. 1.0d0 ) THEN
         MFS0 = nh3_POLY4( C0, AW )
         MFS1 = nh3_POLY4( C1, AW )
         Y0   = ( 1.0d0 - MFS0 ) / MFS0
         Y1   = ( 1.0d0 - MFS1 ) / MFS1
         Y    = ( 1.0d0 - X    ) * Y0 + X * Y1

      ELSE IF ( X .LT. 1.5d0 ) THEN

         IF ( IRH .GE. 40 ) THEN
            MFS1  = nh3_POLY4( C1,  AW )
            MFS15 = nh3_POLY4( C15, AW )
            Y1    = ( 1.0d0 - MFS1  ) / MFS1
            Y15   = ( 1.0d0 - MFS15 ) / MFS15
            Y     = 2.0d0 * ( Y1 * ( 1.5d0 - X ) + Y15 *( X - 1.0d0 ) )

         !==============================================================
         ! Set up for crystalization
         !
         ! Crystallization is done as follows:
         !
         ! For 1.5 <= x, crystallization is assumed to occur
         ! at rh = 0.4
         !
         ! For x <= 1.0, crystallization is assumed to occur at an
         ! rh < 0.01, and since the code does not allow ar rh < 0.01,
         ! crystallization is assumed not to occur in this range.
         !
         ! For 1.0 <= x <= 1.5 the crystallization curve is a straignt
         ! line from a value of y15 at rh = 0.4 to a value of zero at
         ! y1. From point B to point A in the diagram.  The algorithm
         ! does a double interpolation to calculate the amount of
         ! water.
         !
         !        y1(0.40)               y15(0.40)
         !         +                     + Point B
         !
         !
         !
         !
         !         +--------------------+
         !       x=1                   x=1.5
         !      Point A
         !==============================================================
         ELSE

            ! rh along the crystallization curve.
            AWC = 0.80d0 * ( X - 1.0d0 )
            Y   = 0.0d0

            ! interpolate using crystalization curve
            IF ( AW .GE. AWC ) THEN
               MFS1  = nh3_POLY4( C1,  0.40d0 )
               MFS15 = nh3_POLY4( C15, 0.40d0 )
               Y140  = ( 1.0d0 - MFS1  ) / MFS1
               Y1540 = ( 1.0d0 - MFS15 ) / MFS15
               Y40   = 2.0d0 * ( Y140  * ( 1.5d0 - X ) + &
     &                           Y1540 * ( X - 1.0d0 ) )

               ! Y along crystallization curve
               YC   = 2.0d0 * Y1540 * ( X - 1.0d0 )
               Y    = Y40 - (Y40 - YC) * (0.40d0 - AW) / (0.40d0 - AWC)
            ENDIF
         ENDIF

      ELSE IF ( X .LT. 2.0d0 ) then               ! changed 12/11/2000 by FSB
         Y = 0.0D0

         IF ( IRH .GE. 40 ) THEN
            MFS15  = nh3_POLY4( C15, AW )
            !MFS2  = nh3_POLY4( C2,  AW )
            Y15    = ( 1.0d0 - MFS15 ) / MFS15
            !y2    = ( 1.0d0 - MFS2  ) / MFS2
            MFSSO4 = nh3_POLY6( KSO4, AW )             ! Changed 05/30/2000 by FSB
            Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
            Y      = 2.0d0 * (Y15 * (2.0d0 - X) + Y2 * (X - 1.5d0) )
         ENDIF

      ELSE                                 ! 2.0 <= x changed 12/11/2000 by FSB

         !==============================================================
         ! Regime where ammonium sulfate and ammonium nitrate are
         ! in solution.
         !
         ! following cf&s for both ammonium sulfate and ammonium nitrate
         ! check for crystallization here. their data indicate a 40%
         ! value is appropriate.
         !==============================================================
         Y2 = 0.0d0
         Y3 = 0.0d0

         IF ( IRH .GE. 40 ) THEN
            MFSSO4 = nh3_POLY6( KSO4, AW )
            MFSNO3 = nh3_POLY6( KNO3, AW )
            Y2     = ( 1.0d0 - MFSSO4 ) / MFSSO4
            Y3     = ( 1.0d0 - MFSNO3 ) / MFSNO3

         ENDIF

      ENDIF                     ! end of checking on x

      !=================================================================
      ! Now set up output of WH2O
      ! WH2O units are micrograms (liquid water) / cubic meter of air
      !=================================================================
      IF ( X .LT. 2.0D0 ) THEN  ! changed 12/11/2000 by FSB

         WH2O =  Y * ( TSO4 * MWSO4 + MWNH4 * TNH4 )

      ELSE

         ! this is the case that all the sulfate is ammonium sulfate
         ! and the excess ammonium forms ammonum nitrate
         WH2O =   Y2 * TSO4 * MW2 + Y3 * TNO3 * MWANO3

      ENDIF

      ! Return to calling program
      END SUBROUTINE AWATER

!------------------------------------------------------------------------------

      FUNCTION nh3_POLY4( A, X ) RESULT( Y )

      ! Arguments
      REAL*8, INTENT(IN) :: A(4), X

      ! Return value
      REAL*8             :: Y

      !=================================================================
      ! nh3_POLY4 begins here!
      !=================================================================
      Y = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) )))

      ! Return to calling program
      END FUNCTION nh3_POLY4

!------------------------------------------------------------------------------

      FUNCTION nh3_POLY6( A, X ) RESULT( Y )

      ! Arguments
      REAL*8, INTENT(IN) :: A(6), X

      ! Return value
      REAL*8             :: Y

      !=================================================================
      ! nh3_POLY6 begins here!
      !=================================================================
      Y = A(1) + X * ( A(2) + X * ( A(3) + X * ( A(4) +  &
     &           X * ( A(5) + X * ( A(6)  )))))

      ! Return to calling program
      END FUNCTION nh3_POLY6

!------------------------------------------------------------------------------

      SUBROUTINE CUBIC( A2, A1, A0, NR, CRUTES, rc )
!      SUBROUTINE CUBIC( A2, A1, A0, NR, CRUTES )

!
!******************************************************************************
! Subroutine to find the roots of a cubic equation / 3rd order polynomial
! Formulae can be found in numer. recip.  on page 145
!   kiran  developed  this version on 25/4/1990
!   Dr. Francis S. Binkowski modified the routine on 6/24/91, 8/7/97
! ***
! *** modified 2/23/98 by fsb to incorporate Dr. Ingmar Ackermann's
!     recommendations for setting a0, a1,a2 as real*8 variables.
!
! Modified by Bob Yantosca (10/15/02)
! - Now use upper case / white space
! - force double precision with "D" exponents
! - updated comments / cosmetic changes
! - now call ERROR_STOP from "error_mod.f" to stop the run safely
!******************************************************************************
!
      ! Arguments
      INTEGER           :: NR
      REAL*8            :: A2, A1, A0
      REAL*8            :: CRUTES(3)

      integer, intent(out) :: rc

      ! Local variables
      REAL*8            :: QQ,    RR,    A2SQ,  THETA, DUM1, DUM2
      REAL*8            :: PART1, PART2, PART3, RRSQ,  PHI,  YY1
      REAL*8            :: YY2,   YY3,   COSTH, SINTH
      REAL*8, PARAMETER :: ONE    = 1.0d0
      REAL*8, PARAMETER :: SQRT3  = 1.732050808d0
      REAL*8, PARAMETER :: ONE3RD = 0.333333333d0
      ! !LOCAL VARIABLES:
      character (len=75) :: err_msg

      integer :: status

      !=================================================================
      ! CUBIC begins here!
      !=================================================================
      rc = __SUCCESS__

      A2SQ = A2 * A2
      QQ   = ( A2SQ - 3.d0*A1 ) / 9.d0
      RR   = ( A2*( 2.d0*A2SQ - 9.d0*A1 ) + 27.d0*A0 ) / 54.d0

      ! CASE 1 THREE REAL ROOTS or  CASE 2 ONLY ONE REAL ROOT
      DUM1 = QQ * QQ * QQ
      RRSQ = RR * RR
      DUM2 = DUM1 - RRSQ

      IF ( DUM2 .GE. 0.d0 ) THEN

         ! Now we have three real roots
         PHI = SQRT( DUM1 )

         IF ( ABS( PHI ) .LT. 1.d-20 ) THEN
            CRUTES(1) = 0.0d0
            CRUTES(2) = 0.0d0
            CRUTES(3) = 0.0d0
            NR        = 0
!.sds no such module - what is ours?
!.sds            CALL ERROR_STOP( 'PHI < 1d-20', 'CUBIC (rpmares_mod.f)' )
            !$omp critical (G2G_proc_8)
            print *,'PHI < 1d-20 in  CUBIC (rpmares_mod.f)'
            !$omp end critical (G2G_proc_8)
            err_msg = 'PHI < 1d-20 in  CUBIC (rpmares_mod.f):'
            call PrintError  &
     &         (err_msg, .true., 0, 0, 0, 0, 0.0d0, 0.0d0, __RC_NO_OPT__)

         ENDIF

         THETA = ACOS( RR / PHI ) / 3.0d0
         COSTH = COS( THETA )
         SINTH = SIN( THETA )

         ! Use trig identities to simplify the expressions
         ! Binkowski's modification
         PART1     = SQRT( QQ )
         YY1       = PART1 * COSTH
         YY2       = YY1 - A2/3.0d0
         YY3       = SQRT3 * PART1 * SINTH
         CRUTES(3) = -2.0d0*YY1 - A2/3.0d0
         CRUTES(2) = YY2 + YY3
         CRUTES(1) = YY2 - YY3

         ! Set negative roots to a large positive value
         IF ( CRUTES(1) .LT. 0.0d0 ) CRUTES(1) = 1.0d9
         IF ( CRUTES(2) .LT. 0.0d0 ) CRUTES(2) = 1.0d9
         IF ( CRUTES(3) .LT. 0.0d0 ) CRUTES(3) = 1.0d9

         ! Put smallest positive root in crutes(1)
         CRUTES(1) = MIN( CRUTES(1), CRUTES(2), CRUTES(3) )
         NR        = 3

      ELSE

         ! Now here we have only one real root
         PART1     = SQRT( RRSQ - DUM1 )
         PART2     = ABS( RR )
         PART3     = ( PART1 + PART2 )**ONE3RD
         CRUTES(1) = -SIGN(ONE,RR) * ( PART3 + (QQ/PART3) ) - A2/3.D0
         CRUTES(2) = 0.D0
         CRUTES(3) = 0.D0
         NR        = 1

      ENDIF

      ! Return to calling program
      END SUBROUTINE CUBIC

!------------------------------------------------------------------------------

       SUBROUTINE ACTCOF( CAT, AN, GAMA, MOLNU, PHIMULT )
!
!******************************************************************************
!
! DESCRIPTION:
!
!  This subroutine computes the activity coefficients of (2NH4+,SO4--),
!  (NH4+,NO3-),(2H+,SO4--),(H+,NO3-),AND (H+,HSO4-) in aqueous
!  multicomponent solution, using Bromley's model and Pitzer's method.
!
! REFERENCES:
!
!   Bromley, L.A. (1973) Thermodynamic properties of strong electrolytes
!     in aqueous solutions.  AIChE J. 19, 313-320.
!
!   Chan, C.K. R.C. Flagen, & J.H.  Seinfeld (1992) Water Activities of
!     NH4NO3 / (NH4)2SO4 solutions, Atmos. Environ. (26A): 1661-1673.
!
!   Clegg, S.L. & P. Brimblecombe (1988) Equilibrium partial pressures
!     of strong acids over saline solutions - I HNO3,
!     Atmos. Environ. (22): 91-100
!
!   Clegg, S.L. & P. Brimblecombe (1990) Equilibrium partial pressures
!     and mean activity and osmotic coefficients of 0-100% nitric acid
!     as a function of temperature,   J. Phys. Chem (94): 5369 - 5380
!
!   Pilinis, C. and J.H. Seinfeld (1987) Continued development of a
!     general equilibrium model for inorganic multicomponent atmospheric
!     aerosols.  Atmos. Environ. 21(11), 2453-2466.
!
!
!
!
! ARGUMENT DESCRIPTION:
!
!     CAT(1) : conc. of H+    (moles/kg)
!     CAT(2) : conc. of NH4+  (moles/kg)
!     AN(1)  : conc. of SO4-- (moles/kg)
!     AN(2)  : conc. of NO3-  (moles/kg)
!     AN(3)  : conc. of HSO4- (moles/kg)
!     GAMA(2,1)    : mean molal ionic activity coeff for (2NH4+,SO4--)
!     GAMA(2,2)    :  "    "     "       "       "    "  (NH4+,NO3-)
!     GAMA(2,3)    :  "    "     "       "       "    "  (NH4+. HSO4-)
!     GAMA(1,1)    :  "    "     "       "       "    "  (2H+,SO4--)
!     GAMA(1,2)    :  "    "     "       "       "    "  (H+,NO3-)
!     GAMA(1,3)    :  "    "     "       "       "    "  (H+,HSO4-)
!     MOLNU   : the total number of moles of all ions.
!     PHIMULT : the multicomponent paractical osmotic coefficient.
!
! REVISION HISTORY:
!      Who       When        Detailed description of changes
!   ---------   --------  -------------------------------------------
!   S.Roselle   7/26/89   Copied parts of routine BROMLY, and began this
!                         new routine using a method described by Pilinis
!                         and Seinfeld 1987, Atmos. Envirn. 21 pp2453-2466.
!   S.Roselle   7/30/97   Modified for use in Models-3
!   F.Binkowski 8/7/97    Modified coefficients BETA0, BETA1, CGAMA
!   R.Yantosca  9/25/02   Ported into "rpmares_mod.f" for GEOS-CHEM.  Cleaned
!                         up comments, etc.  Also force double precision by
!                         declaring REALs as REAL*8 and by using "D" exponents.
!******************************************************************************
      ! Error codes



      !=================================================================
      ! PARAMETERS and their descriptions:
      !=================================================================
      INTEGER, PARAMETER :: NCAT = 2         ! number of cation
      INTEGER, PARAMETER :: NAN  = 3         ! number of anions
      REAL*8,  PARAMETER :: XSTAT0 = 0       ! Normal, successful completion
      REAL*8,  PARAMETER :: XSTAT1 = 1       ! File I/O error
      REAL*8,  PARAMETER :: XSTAT2 = 2       ! Execution error
      REAL*8,  PARAMETER :: XSTAT3 = 3       ! Special  error

      !=================================================================
      ! ARGUMENTS and their descriptions
      !=================================================================
      REAL*8, optional   :: MOLNU            ! tot # moles of all ions
      REAL*8, optional   :: PHIMULT          ! multicomponent paractical
                                             !   osmotic coef
      REAL*8             :: CAT(NCAT)        ! cation conc in moles/kg (input)
      REAL*8             :: AN(NAN)          ! anion conc in moles/kg (input)
      REAL*8             :: GAMA(NCAT,NAN)   ! mean molal ionic activity coefs

      !=================================================================
      ! SCRATCH LOCAL VARIABLES and their descriptions:
      !=================================================================
      INTEGER            :: IAN              ! anion indX
      INTEGER            :: ICAT             ! cation indX
      REAL*8             :: FGAMA            !
      REAL*8             :: I                ! ionic strength
      REAL*8             :: R                !
      REAL*8             :: S                !
      REAL*8             :: TA               !
      REAL*8             :: TB               !
      REAL*8             :: TC               !
      REAL*8             :: TEXPV            !
      REAL*8             :: TRM              !
      REAL*8             :: TWOI             ! 2*ionic strength
      REAL*8             :: TWOSRI           ! 2*sqrt of ionic strength
      REAL*8             :: ZBAR             !
      REAL*8             :: ZBAR2            !
      REAL*8             :: ZOT1             !
      REAL*8             :: SRI              ! square root of ionic strength
      REAL*8             :: F2(NCAT)         !
      REAL*8             :: F1(NAN)          !
      REAL*8             :: BGAMA (NCAT,NAN) !
      REAL*8             :: X     (NCAT,NAN) !
      REAL*8             :: M     (NCAT,NAN) ! molality of each electrolyte
      REAL*8             :: LGAMA0(NCAT,NAN) ! binary activity coefficients
      REAL*8             :: Y     (NAN,NCAT) !
      REAL*8             :: BETA0 (NCAT,NAN) ! binary activity coef parameter
      REAL*8             :: BETA1 (NCAT,NAN) ! binary activity coef parameter
      REAL*8             :: CGAMA (NCAT,NAN) ! binary activity coef parameter
      REAL*8             :: V1    (NCAT,NAN) ! # of cations in electrolyte
                                             !   formula
      REAL*8             :: V2    (NCAT,NAN) ! # of anions in electrolyte
                                             !   formula
      ! absolute value of charges of cation
      REAL*8             :: ZP(NCAT) = (/ 1.0d0, 1.0d0 /)

      ! absolute value of charges of anion
      REAL*8             :: ZM(NAN)  = (/ 2.0d0, 1.0d0, 1.0d0 /)

      ! Character values.
      CHARACTER(LEN=120)      :: XMSG  = ' '
!      CHARACTER(LEN=16), SAVE :: PNAME = ' driver program name'

      !================================================================
      ! *** Sources for the coefficients BETA0, BETA1, CGAMA
      ! (1,1);(1,3)  - Clegg & Brimblecombe (1988)
      ! (2,3)        - Pilinis & Seinfeld (1987), cgama different
      ! (1,2)        - Clegg & Brimblecombe (1990)
      ! (2,1);(2,2)  - Chan, Flagen & Seinfeld (1992)
      !================================================================

      ! now set the basic constants, BETA0, BETA1, CGAMA
      DATA BETA0(1,1) /2.98d-2/,      BETA1(1,1) / 0.0d0/,  &
     &     CGAMA(1,1) /4.38d-2/                                 ! 2H+SO4-

      DATA BETA0(1,2) /  1.2556d-1/,  BETA1(1,2) / 2.8778d-1/,  &
     &     CGAMA(1,2) / -5.59d-3/                               ! HNO3

      DATA BETA0(1,3) / 2.0651d-1/,   BETA1(1,3) / 5.556d-1/,  &
     &     CGAMA(1,3) /0.0d0/                                   ! H+HSO4-

      DATA BETA0(2,1) / 4.6465d-2/,   BETA1(2,1) /-0.54196d0/,  &
     &     CGAMA(2,1) /-1.2683d-3/                              ! (NH4)2SO4

      DATA BETA0(2,2) /-7.26224d-3/,  BETA1(2,2) /-1.168858d0/,  &
     &     CGAMA(2,2) / 3.51217d-5/                             ! NH4NO3

      DATA BETA0(2,3) / 4.494d-2/,    BETA1(2,3) / 2.3594d-1/,  &
     &     CGAMA(2,3) /-2.962d-3/                               ! NH4HSO4

      DATA V1(1,1), V2(1,1) / 2.0d0, 1.0d0 /     ! 2H+SO4-
      DATA V1(2,1), V2(2,1) / 2.0d0, 1.0d0 /     ! (NH4)2SO4
      DATA V1(1,2), V2(1,2) / 1.0d0, 1.0d0 /     ! HNO3
      DATA V1(2,2), V2(2,2) / 1.0d0, 1.0d0 /     ! NH4NO3
      DATA V1(1,3), V2(1,3) / 1.0d0, 1.0d0 /     ! H+HSO4-
      DATA V1(2,3), V2(2,3) / 1.0d0, 1.0d0 /     ! NH4HSO4

      !=================================================================
      ! ACTCOF begins here!
      !=================================================================

      ! Compute ionic strength
      I = 0.0d0
      DO ICAT = 1, NCAT
         I = I + CAT( ICAT ) * ZP( ICAT ) * ZP( ICAT )
      ENDDO

      DO IAN = 1, NAN
         I = I + AN( IAN ) * ZM( IAN ) * ZM( IAN )
      ENDDO

      I = 0.5d0 * I

      ! check for problems in the ionic strength
      IF ( I .EQ. 0.0d0 ) THEN

         DO IAN  = 1, NAN
         DO ICAT = 1, NCAT
            GAMA( ICAT, IAN ) = 0.0d0
         ENDDO
         ENDDO

         XMSG = 'Ionic strength is zero...returning zero activities'
         !CALL M3WARN ( PNAME, 0, 0, XMSG )
         RETURN

      ELSE IF ( I .LT. 0.0d0 ) THEN
         XMSG = 'Ionic strength below zero...negative concentrations'
         !CALL M3EXIT ( PNAME, 0, 0, XMSG, XSTAT1 )
      ENDIF

      ! Compute some essential expressions
      SRI    = SQRT( I )
      TWOSRI = 2.0d0 * SRI
      TWOI   = 2.0d0 * I
      TEXPV  = 1.0d0 - EXP( -TWOSRI ) * ( 1.0d0 + TWOSRI - TWOI )
      R      = 1.0d0 + 0.75d0 * I
      S      = 1.0d0 + 1.5d0  * I
      ZOT1   = 0.511d0 * SRI / ( 1.0d0 + SRI )

      ! Compute binary activity coeffs
      FGAMA = -0.392d0 * ( ( SRI / ( 1.0d0 + 1.2d0 * SRI )  &
     &      + ( 2.0d0 / 1.2d0 ) * LOG( 1.0d0 + 1.2d0 * SRI ) ) )

      DO ICAT = 1, NCAT
      DO IAN  = 1, NAN

         BGAMA( ICAT, IAN ) = 2.0d0 * BETA0( ICAT, IAN )  &
     &        + ( 2.0d0 * BETA1( ICAT, IAN ) / ( 4.0d0 * I ) )  &
     &        * TEXPV

         ! Compute the molality of each electrolyte for given ionic strength
         M( ICAT, IAN ) = ( CAT( ICAT )**V1( ICAT, IAN )  &
     &                   *   AN( IAN )**V2( ICAT, IAN ) )**( 1.0d0  &
     &                   / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) ) )

         ! Calculate the binary activity coefficients
         LGAMA0( ICAT, IAN ) = ( ZP( ICAT ) * ZM( IAN ) * FGAMA  &
     &        + M( ICAT, IAN )  &
     &        * ( 2.0d0 * V1( ICAT, IAN ) * V2( ICAT, IAN )  &
     &        / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )  &
     &        * BGAMA( ICAT, IAN ) )  &
     &        + M( ICAT, IAN ) * M( ICAT, IAN )  &
     &        * ( 2.0d0 * ( V1( ICAT, IAN )  &
     &        * V2( ICAT, IAN ) )**1.5d0  &
     &        / ( V1( ICAT, IAN ) + V2( ICAT, IAN ) )  &
     &        * CGAMA( ICAT, IAN ) ) ) / 2.302585093d0

      ENDDO
      ENDDO

      ! prepare variables for computing the multicomponent activity coeffs
      DO IAN = 1, NAN
      DO ICAT = 1, NCAT
         ZBAR           = ( ZP( ICAT ) + ZM( IAN ) ) * 0.5d0
         ZBAR2          = ZBAR * ZBAR
         Y( IAN, ICAT ) = ZBAR2 * AN( IAN ) / I
         X( ICAT, IAN ) = ZBAR2 * CAT( ICAT ) / I
      ENDDO
      ENDDO

      DO IAN = 1, NAN
         F1( IAN ) = 0.0d0
         DO ICAT = 1, NCAT
            F1( IAN ) = F1( IAN ) + X( ICAT, IAN ) * LGAMA0( ICAT, IAN )  &
     &                + ZOT1 * ZP( ICAT ) * ZM( IAN ) * X( ICAT, IAN )
         ENDDO
      ENDDO

      DO ICAT = 1, NCAT
         F2( ICAT ) = 0.0d0
         DO IAN = 1, NAN
            F2( ICAT ) = F2( ICAT ) + Y( IAN, ICAT ) * LGAMA0(ICAT, IAN)  &
     &                 + ZOT1 * ZP( ICAT ) * ZM( IAN ) * Y( IAN, ICAT )
         ENDDO
      ENDDO

      ! now calculate the multicomponent activity coefficients
      DO IAN  = 1, NAN
      DO ICAT = 1, NCAT

         TA  = -ZOT1 * ZP( ICAT ) * ZM( IAN )
         TB  = ZP( ICAT ) * ZM( IAN ) / ( ZP( ICAT ) + ZM( IAN ) )
         TC  = ( F2( ICAT ) / ZP( ICAT ) + F1( IAN ) / ZM( IAN ) )
         TRM = TA + TB * TC

         IF ( TRM .GT. 30.0d0 ) THEN
            GAMA( ICAT, IAN ) = 1.0d+30
         ELSE
            GAMA( ICAT, IAN ) = 10.0d0**TRM
         ENDIF

      ENDDO
      ENDDO

      ! Return to calling program
      END SUBROUTINE ACTCOF

!------------------------------------------------------------------------------

!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: PrintError
!
! !INTERFACE:
!
      subroutine PrintError  &
        (err_msg, err_do_stop, err_num_ints, err_int1, err_int2,  &
         err_num_reals, err_real1, err_real2, rc)
!
      implicit none
!
! !INPUT PARAMETERS:
!!     err_msg       : error message to be printed out
!!     err_do_stop   : do stop on error?
!!     err_num_ints  : number of integers to be printed out (0, 1, or 2)
!!     err_int1      : integer 1 to print out
!!     err_int2      : integer 2 to print out
!!     err_num_reals : number of reals to be printed out (0, 1, or 2)
!!     err_real1     : real 1 to print out
!!     err_real2     : real 2 to print out
      character (len=*), intent(in) :: err_msg
      logical          , intent(in) :: err_do_stop
      integer          , intent(in) :: err_num_ints
      integer          , intent(in) :: err_int1
      integer          , intent(in) :: err_int2
      integer          , intent(in) :: err_num_reals
      real*8           , intent(in) :: err_real1
      real*8           , intent(in) :: err_real2

      integer, intent(out) :: rc
!
! !DESCRIPTION:
!  Output error messages, and exits if requested.
!
! !AUTHOR:
!  Jules Kouatchou
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
      rc = __SUCCESS__
!BOC
      !$omp critical (G2G_proc_9)
      Write (6,*)
      Write (6,*) &
        '--------------------------------------------------------------'

      Write (6,*) '!! ' // Trim (err_msg)

      if (err_num_ints == 1) then
         Write (6,*) '   ', err_int1
      else if (err_num_ints == 2) then
         Write (6,*) '   ', err_int1, err_int2
      end if

      if (err_num_reals == 1) then
         Write (6,*) '   ', err_real1
      else if (err_num_reals == 2) then
         Write (6,*) '   ', err_real1, err_real2
      end if

      Write (6,*) &
        '--------------------------------------------------------------'
      Write (6,*)
      !$omp end critical (G2G_proc_9)

      if (err_do_stop) then
        rc = __FAIL__
        return
      end if

      return

      end subroutine PrintError

!==================================================================================
!BOP
!
! !IROUTINE:  Chem_UtilResVal --- returns resolution dependent value
!
! !INTERFACE:
!
   function Chem_UtilResVal( im_World, jm_World, res_value, rc ) result (val)

! !USES:

   implicit NONE

   real :: val                                ! resolution dependent value

! !INPUT/OUTPUT PARAMETERS:
   integer, intent(in) :: im_World, jm_World  ! number of global grid cells
   real,    intent(in) :: res_value(:)        ! array with the resolution dependent values:
                                              ! the 'a', 'b', ..., 'e' resolution values have
                                              ! indexes 1, 2, ..., 5.

! !OUTPUT PARAMETERS:
   integer, intent(inout) :: rc               ! return code


! !DESCRIPTION:
!
! !REVISION HISTORY:
!
! 13 Feb2012   Anton Darmenov  First crack.
! 25 Oct2012   Anton Darmenov  Added support for FV3 resolutions.
! 19 Aug2020   E. Sherman - moved from Chem_UtilMod.F90 to process library
!
!EOP
!-------------------------------------------------------------------------
       character(len=*), parameter :: Iam = 'Chem_UtilResVal'

       integer            :: i_res
       integer, parameter :: res_a = 1  ! 'a' to 'e' resolution indexes
       integer, parameter :: res_b = 2  !
       integer, parameter :: res_c = 3  !
       integer, parameter :: res_d = 4  !
       integer, parameter :: res_e = 5  !
       integer, parameter :: res_f = 6  !

       i_res = 0

       if ((im_World < 1) .or. (jm_World < 1)) then
!           call die(Iam, 'incorrect model resolution')
           print*,'GOCART2G_Process::Chem_UtilResVal - incorrect model resolution'
           return
       end if

       if (jm_World == 6*im_World) then
           if (im_World <= 24) then
               i_res = res_a
           else if (im_World <=  48) then
               i_res = res_b
           else if (im_World <=  90) then
               i_res = res_c
           else if (im_World <= 180) then
               i_res = res_d
           else if (im_World <= 360) then
               i_res = res_e
           else if (im_World <= 720) then
               i_res = res_f
           else
               i_res = res_f
           end if
       else
           if ((im_World <= 72) .and. (jm_World <= 46)) then
               i_res = res_a
           else if ((im_World <=  144) .and. (jm_World <=  91)) then
               i_res = res_b
           else if ((im_World <=  288) .and. (jm_World <= 181)) then
               i_res = res_c
           else if ((im_World <=  576) .and. (jm_World <= 361)) then
               i_res = res_d
           else if ((im_World <= 1152) .and. (jm_World <= 721)) then
               i_res = res_e
           else if ((im_World <= 2304) .and. (jm_World <=1441)) then
               i_res = res_f
           else
               i_res = res_f
           end if


       end if

       if ((i_res < 1) .or. (i_res > size(res_value))) then
           val = 0.0
           rc  = __FAIL__
       else
           val = res_value(i_res)
           rc  = __SUCCESS__
       end if

   end function Chem_UtilResVal

!==================================================================================

   function Chem_UtilIdow(nymd) result (idow)
     implicit NONE
     integer, intent(in) :: nymd
     integer :: idow ! day of the week: Sun=1, Mon=2, etc.
     integer :: y, m, d
     integer, parameter :: t(0:11) = (/ 0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4 /)
     y = nymd / 10000
     m = (nymd - y*10000)/100
     d = nymd - (y*10000 + m*100)
     if ( m<3 ) then
        y = y - 1
     end if
     idow = 1+mod(y + y/4 - y/100 + y/400 + t(m-1) + d,7)
     return
   end function Chem_UtilIdow

   function Chem_UtilCdow(nymd) result (cdow)
     implicit NONE
     integer, intent(in) :: nymd
     character(len=3) :: cdow ! day of the week: Sun, Mon, etc.
     character(len=3) :: cday(7) = (/ 'Sun','Mon', 'Tue', 'Wed', 'Thu', 'Fri', 'Sat' /)
     cdow = cday(Chem_UtilIdow(nymd))
     return
   end function Chem_UtilCdow


!==================================================================================

!BOP
!
! !ROUTINE:  Chem_BiomassDiurnal - Applies diurnal cycle to biomass emissions.
!
! !INTERFACE:
     subroutine Chem_BiomassDiurnal ( Eout, Ein, lons, lats, nhms, cdt)

! !USES:

  IMPLICIT NONE

! !ARGUMENTS:

       real, intent(out)   :: Eout(:,:) ! Emissions valid at NHMS
       real, intent(in)    :: Ein(:,:)  ! Daily-mean emissions
       real, intent(in)    :: lons(:,:) ! Latitudes in degrees
       real, intent(in)    :: lats(:,:) ! Latitudes in degrees
       integer, intent(in) :: nhms
       real, intent(in)    :: cdt       ! time step in seconds

! !DESCRIPTION:
!
!      Applies diurnal cycle to biomass emissions.
!
! !DESCRIPTION:
!
!  This module implements assorted odds & ends for fvChem.
!
! !REVISION HISTORY:
!
!  13nov2009  da Silva  First crack.
!  19Aug2020  E. Sherman - moved from Chem_UtilMod.F90 to process library
!
!EOP
!-------------------------------------------------------------------------

!      Hardwired diurnal cycle (multiplied by 100)
!      These numbers were derived from GOES-12
!      fire counts for 2003-2007.
!      -------------------------------------------
       integer, parameter :: N = 240
       real,    parameter :: DT = 86400. / N

!      Apply flat diurnal cycle for boreal forests as a
!      temporary solution to prevent very high aerosol
!      optical depth during the day
       real,    parameter :: Boreal(N) = 1.0
!      real,    parameter :: Boreal(N) = &
!      (/ 0.0277, 0.0292, 0.0306, 0.0318, 0.0327, 0.0335, &
!         0.0340, 0.0342, 0.0341, 0.0338, 0.0333, 0.0326, &
!         0.0316, 0.0305, 0.0292, 0.0278, 0.0263, 0.0248, &
!         0.0233, 0.0217, 0.0202, 0.0187, 0.0172, 0.0158, &
!         0.0145, 0.0133, 0.0121, 0.0110, 0.0100, 0.0091, &
!         0.0083, 0.0075, 0.0068, 0.0062, 0.0056, 0.0051, &
!         0.0046, 0.0042, 0.0038, 0.0035, 0.0032, 0.0030, &
!         0.0028, 0.0026, 0.0025, 0.0024, 0.0024, 0.0024, &
!         0.0024, 0.0026, 0.0027, 0.0030, 0.0033, 0.0036, &
!         0.0041, 0.0046, 0.0052, 0.0060, 0.0069, 0.0079, &
!         0.0090, 0.0104, 0.0119, 0.0137, 0.0157, 0.0180, &
!         0.0205, 0.0235, 0.0268, 0.0305, 0.0346, 0.0393, &
!         0.0444, 0.0502, 0.0565, 0.0634, 0.0711, 0.0794, &
!         0.0884, 0.0982, 0.1087, 0.1201, 0.1323, 0.1453, &
!         0.1593, 0.1742, 0.1900, 0.2069, 0.2249, 0.2439, &
!         0.2642, 0.2858, 0.3086, 0.3329, 0.3587, 0.3860, &
!         0.4149, 0.4455, 0.4776, 0.5115, 0.5470, 0.5840, &
!         0.6227, 0.6628, 0.7043, 0.7470, 0.7908, 0.8355, &
!         0.8810, 0.9271, 0.9735, 1.0200, 1.0665, 1.1126, &
!         1.1580, 1.2026, 1.2460, 1.2880, 1.3282, 1.3664, &
!         1.4023, 1.4356, 1.4660, 1.4933, 1.5174, 1.5379, &
!         1.5548, 1.5679, 1.5772, 1.5826, 1.5841, 1.5818, &
!         1.5758, 1.5661, 1.5529, 1.5365, 1.5169, 1.4944, &
!         1.4693, 1.4417, 1.4119, 1.3801, 1.3467, 1.3117, &
!         1.2755, 1.2383, 1.2003, 1.1616, 1.1225, 1.0832, &
!         1.0437, 1.0044, 0.9653, 0.9265, 0.8882, 0.8504, &
!         0.8134, 0.7771, 0.7416, 0.7070, 0.6734, 0.6407, &
!         0.6092, 0.5787, 0.5493, 0.5210, 0.4939, 0.4680, &
!         0.4433, 0.4197, 0.3974, 0.3763, 0.3565, 0.3380, &
!         0.3209, 0.3051, 0.2907, 0.2777, 0.2662, 0.2561, &
!         0.2476, 0.2407, 0.2352, 0.2313, 0.2289, 0.2279, &
!         0.2283, 0.2300, 0.2329, 0.2369, 0.2417, 0.2474, &
!         0.2536, 0.2602, 0.2670, 0.2738, 0.2805, 0.2869, &
!         0.2927, 0.2979, 0.3024, 0.3059, 0.3085, 0.3101, &
!         0.3107, 0.3102, 0.3087, 0.3061, 0.3026, 0.2983, &
!         0.2931, 0.2871, 0.2806, 0.2735, 0.2659, 0.2579, &
!         0.2497, 0.2412, 0.2326, 0.2240, 0.2153, 0.2066, &
!         0.1979, 0.1894, 0.1809, 0.1726, 0.1643, 0.1562, &
!         0.1482, 0.1404, 0.1326, 0.1250, 0.1175, 0.1101, &
!         0.1028, 0.0956, 0.0886, 0.0818, 0.0751, 0.0687 /)
       real,    parameter :: NonBoreal(N) = &
       (/ 0.0121, 0.0150, 0.0172, 0.0185, 0.0189, 0.0184, &
          0.0174, 0.0162, 0.0151, 0.0141, 0.0133, 0.0126, &
          0.0121, 0.0117, 0.0115, 0.0114, 0.0114, 0.0116, &
          0.0120, 0.0126, 0.0133, 0.0142, 0.0151, 0.0159, &
          0.0167, 0.0174, 0.0180, 0.0184, 0.0187, 0.0189, &
          0.0190, 0.0190, 0.0191, 0.0192, 0.0192, 0.0193, &
          0.0194, 0.0194, 0.0193, 0.0192, 0.0190, 0.0187, &
          0.0185, 0.0182, 0.0180, 0.0178, 0.0177, 0.0176, &
          0.0174, 0.0172, 0.0169, 0.0166, 0.0162, 0.0158, &
          0.0153, 0.0149, 0.0144, 0.0138, 0.0132, 0.0126, &
          0.0118, 0.0109, 0.0101, 0.0092, 0.0085, 0.0081, &
          0.0080, 0.0083, 0.0091, 0.0102, 0.0117, 0.0135, &
          0.0157, 0.0182, 0.0210, 0.0240, 0.0273, 0.0308, &
          0.0345, 0.0387, 0.0432, 0.0483, 0.0540, 0.0606, &
          0.0683, 0.0775, 0.0886, 0.1022, 0.1188, 0.1388, &
          0.1625, 0.1905, 0.2229, 0.2602, 0.3025, 0.3500, &
          0.4031, 0.4623, 0.5283, 0.6016, 0.6824, 0.7705, &
          0.8650, 0.9646, 1.0676, 1.1713, 1.2722, 1.3662, &
          1.4491, 1.5174, 1.5685, 1.6014, 1.6173, 1.6200, &
          1.6150, 1.6082, 1.6040, 1.6058, 1.6157, 1.6353, &
          1.6651, 1.7045, 1.7513, 1.8024, 1.8541, 1.9022, &
          1.9429, 1.9738, 1.9947, 2.0072, 2.0132, 2.0141, &
          2.0096, 1.9994, 1.9829, 1.9604, 1.9321, 1.8977, &
          1.8562, 1.8052, 1.7419, 1.6646, 1.5738, 1.4734, &
          1.3693, 1.2676, 1.1724, 1.0851, 1.0052, 0.9317, &
          0.8637, 0.8004, 0.7414, 0.6862, 0.6348, 0.5871, &
          0.5434, 0.5037, 0.4682, 0.4368, 0.4097, 0.3864, &
          0.3667, 0.3499, 0.3355, 0.3231, 0.3123, 0.3029, &
          0.2944, 0.2862, 0.2773, 0.2670, 0.2547, 0.2402, &
          0.2238, 0.2061, 0.1882, 0.1712, 0.1562, 0.1434, &
          0.1332, 0.1251, 0.1189, 0.1141, 0.1103, 0.1071, &
          0.1043, 0.1018, 0.0996, 0.0979, 0.0968, 0.0964, &
          0.0966, 0.0970, 0.0973, 0.0970, 0.0959, 0.0938, &
          0.0909, 0.0873, 0.0831, 0.0784, 0.0732, 0.0676, &
          0.0618, 0.0565, 0.0521, 0.0491, 0.0475, 0.0473, &
          0.0480, 0.0492, 0.0504, 0.0514, 0.0519, 0.0521, &
          0.0520, 0.0517, 0.0513, 0.0510, 0.0507, 0.0507, &
          0.0508, 0.0512, 0.0515, 0.0518, 0.0519, 0.0518, &
          0.0513, 0.0506, 0.0496, 0.0482, 0.0465, 0.0443, &
          0.0418, 0.0387, 0.0351, 0.0310, 0.0263, 0.0214 /)

!      Fixed normalization factors; a more accurate normalization would take
!      in consideration longitude and time step
!      ---------------------------------------------------------------------
       real*8 :: fBoreal, fNonBoreal
       real :: fDT

       integer :: hh, mm, ss, ndt, i, j, k
       integer :: NN
       real :: secs, secs_local, aBoreal, aNonBoreal, alpha

!                              -----

!      Normalization factor depends on timestep
!      ----------------------------------------
       fBoreal = 0.0
       fNonBoreal = 0.0
       NN = 0
       ndt = max(1,nint(cdt/DT))

       do k = 1, N, ndt
          NN = NN + 1
          fBoreal    = fBoreal    + Boreal(k)
          fNonBoreal = fNonBoreal + NonBoreal(k)
       end do

       fBoreal    = fBoreal / NN
       fnonBoreal = fnonBoreal / NN


!      Find number of secs since begining of the day (GMT)
!      ---------------------------------------------------
       hh = nhms/10000
       mm = (nhms - 10000*hh) / 100
       ss = nhms - 10000*hh - 100*mm
       secs = 3600.*hh + 60.*mm + ss

!      Apply factors depending on latitude
!      -----------------------------------
       do j = lbound(Ein,2), ubound(Ein,2)
         do i = lbound(Ein,1), ubound(Ein,1)

!            Find corresponding index in hardwired diurnal cycle
!            240 = 24 * 60 * 60 secs / 360 deg
!            ---------------------------------------------------
             secs_local = secs + 240. * lons(i,j)
             k = 1 + mod(nint(secs_local/DT),N)
             if ( k < 1 ) k = N + k

!            Apply diurnal cycle
!            -------------------
             aBoreal = Boreal(k) / fBoreal
             aNonBoreal = NonBoreal(k) / fNonBoreal

                if ( lats(i,j) >= 50. ) then
                   Eout(i,j) = aBoreal    * Ein(i,j)
                else if ( lats(i,j) >= 30. ) then
                   alpha = (lats(i,j) - 30. ) / 20.
                   Eout(i,j) = (1-alpha) * aNonBoreal * Ein(i,j) + &
                                  alpha  * aBoreal    * Ein(i,j)
                else
                   Eout(i,j) = aNonBoreal * Ein(i,j)
                end if
          end do
       end do

     end subroutine Chem_BiomassDiurnal
!==================================================================================

   subroutine ReadPointEmissions( nymd, filename, nPts, vLat, vLon, vBase, vTop, vEmis, vStart, vEnd, unusable, label, rc)
      integer, intent(in)            :: nymd
      character(*), intent(in) :: filename
      integer, intent(out)           :: nPts
      real, allocatable, dimension(:), intent(out)    :: vLat, vLon, vTop, vBase, vEmis
      integer, allocatable, dimension(:), intent(out) :: vStart, vEnd

      type(KeywordEnforcer), optional, intent(in) :: unusable
      character(*), optional, intent(in) :: label
      integer, optional, intent(out) :: rc

      ! Local arguments
      type(EmissionReader) :: reader
      character(:), allocatable :: label_
      real, allocatable :: table(:,:)
      integer :: nCols
      integer :: status, status1, status2, status3

      if (present(label)) then
         label_ = trim(label)
      else
         label_ = 'source'
      end if

      reader = EmissionReader()
      !$omp critical (process1)
      call reader%open(filename, rc=status1)
      table = reader%read_table(label=label_, rc=status2)
      call reader%close(rc=status3)
      !$omp end critical (process1)
      __VERIFY__(status1)
      __VERIFY__(status2)
      __VERIFY__(status3)

      nCols = size(table,1)
      nPts = size(table,2)
      vStart = spread(-1, 1, nPts)
      vEnd = spread(-1, 1, nPts)

      vLat  = table(1,:)
      vLon  = table(2,:)
      vEmis = table(3,:)
      vBase = table(4,:)
      vTop  = table(5,:)
      if (nCols >= 6) vStart = table(6,:)
      if (nCols >= 7) vEnd = table(7,:)

      where(vStart < 0) vStart = 000000
      where(vEnd < 0)   vEnd   = 240000
      !call reader%close()

      __RETURN__(__SUCCESS__)
   end subroutine ReadPointEmissions

!==================================================================================
   subroutine open(this, filename, rc)
      class(EmissionReader), intent(inout) :: this
      character(*), intent(in) :: filename
      integer, optional, intent(out) :: rc

      integer :: status

      __ASSERT__(.not. allocated(this%unit))
      allocate(this%unit)

      open(newunit=this%unit, file=filename,  &
           form='formatted', access = 'sequential', status='old', &
           action='read', __IOSTAT__)

      __RETURN__(__SUCCESS__)
   end subroutine open


   subroutine close(this, rc)
      class(EmissionReader), intent(inout) :: this
      integer, optional, intent(out) :: rc

      integer :: status

      __ASSERT__(allocated(this%unit))
      close(this%unit, __IOSTAT__)
      deallocate(this%unit)

      __RETURN__(__SUCCESS__)
   end subroutine close


   subroutine rewind_reader(this, rc)
      class(EmissionReader), intent(in) :: this
      integer, optional, intent(out) :: rc

      integer :: status

      __ASSERT__(allocated(this%unit))
      rewind(this%unit, __IOSTAT__)

      __RETURN__(__SUCCESS__)
   end subroutine rewind_reader

   function get_dims(this, label, rc) result(dims)
      integer :: dims(2)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: label
      integer, optional, intent(out) :: rc

      integer :: status
      logical :: eof
      character(:), allocatable :: line
      integer :: n_words

      call this%rewind(__RC__)
      call this%scan_to_label(label, __RC__)

      dims = 0
      do
         line = this%next_line(eof=eof, __RC__)
         __ASSERT__(.not. eof)
         if (this%is_end_marker(line)) exit

         dims(2) = dims(2) + 1

         n_words = this%count_words(line)
         dims(1) = max(dims(1), n_words)
      end do

      __RETURN__(__SUCCESS__)
   end function get_dims

   integer function count_words(this, line) result(n_words)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: line

      integer :: idx, i0

      n_words = 0
      i0 = 0
      do
         ! scan to start of next word
         idx = verify(line(i0+1:), ' ')

         n_words = n_words + 1
         i0 = i0 + idx

         ! scan to end of current word
         idx = index(line(i0+1:), ' ')
         i0 = i0 + idx
         if (idx == 0) exit

      end do

      return
   end function count_words

   logical function is_end_marker(this, line)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: line

      is_end_marker = (line == '::')

   end function is_end_marker

   function read_table(this, label, rc) result(table)
      class(EmissionReader), intent(in) :: this
      real, allocatable :: table(:,:)
      character(*), intent(in) :: label
      integer, optional, intent(out) :: rc

      integer :: i, j
      integer :: dims(2)
      integer :: status
      logical :: eof
      character(:), allocatable :: line

      dims = this%get_dims(label, __RC__)
      call this%scan_to_label(label, __RC__)

      associate (n_words => dims(1), n_lines => dims(2))
        allocate(table(n_words, n_lines), __STAT__)

        do j = 1, n_lines
           line = this%next_line(eof=eof)
           __ASSERT__(.not. eof)

           read(line,*, iostat=status) (table(i,j),i=1,n_words)
           __VERIFY__(status)
        end do

      end associate

      __RETURN__(__SUCCESS__)
   end function read_table

   function next_line(this, eof, rc) result(line)
      character(:), allocatable :: line
      class(EmissionReader), intent(in) :: this
      logical, intent(out) :: eof
      integer, optional, intent(out) :: rc

      integer, parameter :: MAX_LINE_LEN=1024
      character(len=MAX_LINE_LEN) :: buffer
      integer :: idx
      integer :: status

      eof = .false.
      do

         read(this%unit,'(a)', iostat=status) buffer
         if (status == IOSTAT_END) then
            eof = .true.
            __RETURN__(__SUCCESS__)
         end if
         __VERIFY__(status)

         idx = index(buffer, '#')
         if (idx == 0) idx = len(buffer)

         line = trim(buffer(:idx-1))
         if (line /= '')  exit

      end do

      __RETURN__(__SUCCESS__)
   end function next_line

   subroutine scan_to_label(this, label, rc)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: label
      integer, optional, intent(out) :: rc

      integer :: status
      logical :: eof
      character(:), allocatable :: line

      call this%rewind(__RC__)
      do
         line = this%next_line(eof=eof, __RC__)
         if (line == label // '::') exit
      end do

      __RETURN__(__SUCCESS__)
   end subroutine scan_to_label
 end module GOCART2G_Process
