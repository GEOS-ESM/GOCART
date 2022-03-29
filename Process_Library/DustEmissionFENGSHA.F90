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
