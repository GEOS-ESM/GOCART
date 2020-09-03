#include "unused_dummy.H"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  DustEmissionMod.F90 --- Calculate the dust emissions
!
! !INTERFACE:
!

   module  DustEmissionMod

! !USES:

   use Chem_Mod
   use Chem_ConstMod, only: grav        ! Constants !
   use Chem_UtilMod
   use MAPL_ConstantsMod, only : MAPL_PI

   use m_mpout

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  DustEmissionDEAD
   PUBLIC  DustEmissionGOCART
   PUBLIC  DustEmissionSimple
   PUBLIC  MAM_DustEmissionGOCART
   PUBLIC  MAM_DustEmission
   PUBLIC  KokSizeDistribution    


  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

!
! !DESCRIPTION:
!
!  This module implements the Dust Emission calculations
!
! !REVISION HISTORY:
!
!  29Dec2009 Colarco    First crack!
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
!
! DEAD-based dust emission scheme (Zender et al., JGR, 2003)
! Pass in a grid of wind speeds and surface characteristics and
! returns the total dust emission mass flux [kg m-2 s-1].  
! Emissions need to be scaled by source function (ie., fraction of
! grid cell emitting) and distributed over a particle size distribution.
!
! !IROUTINE:  DustEmissionDEAD - Compute the dust emissions
!
! !INTERFACE:
!

   subroutine DustEmissionDEAD( i1, i2, j1, j2, km, &
                                gwettop, oro, ustar, u10m, v10m, &
                                emissions, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km
   real, pointer, dimension(:,:) :: gwettop, oro, ustar, u10m, v10m

! !OUTPUT PARAMETERS:

   real     ::  emissions(i1:i2,j1:j2)             ! Local emission
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'DustEmissionsDEAD'

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
!  29Dec2009, Colarco - Modifications to change calling
!  10Oct2007, Nowottnick/Colarco - Implement simplified Zender source
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j
   real, parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
   real            ::  u_thresh0        ! dry bed, non-salting saltation threshold [m s-1]
   real            ::  w10m

!  Variables and parameters specific to Zender source implementation
   real            ::  fd               ! drag partitioning eff. factor
                     ! dust effective diameter for monomodal soil
   real, parameter ::  soil_diameter = 75.e-6
                     ! soil grain particle density [kg m-3]
   real, parameter ::  soil_density  = 2650.
   real, parameter ::  zoms = 33.e-6     ! smooth roughness length [m]
   real, parameter ::  zom = 100.e-6     ! roughness length [m]
   real, parameter ::  mclay = 0.2       ! mass fraction of clay
   real            ::  fw               ! water eff. factor
   real            ::  wt               ! Threshold water content
   real            ::  u_thresh_drag
   real            ::  u_thresh_drag_water
   real            ::  ustars           ! Modified fric. vel for Owen Effect
   real            ::  k_z              ! Log profile for nonsaltating
   real            ::  wt10m            ! Threshold 10m wind speed
   real, parameter ::  cs=2.61          ! Constant for horiz flux
   real            ::  rat              ! Ratio of threshold to friction vel
   real            ::  horiz_flux       ! Horizontal Mass Flux
   real            ::  vert_flux        ! Vertical Mass Flux
   real            ::  alpha            ! Vertical Flux Conversion Factor [m-1]

   _UNUSED_DUMMY(km)

!  Initialize local variables
!  --------------------------
   emissions(:,:) = 0.

!  Calculate drag partitioning efficiency to represent sink of atmos momentum 
!  to nonerodible roughness elements
!  Assumes constant smooth roughness length zoms and constant roughness
!  length for momentum transfer zom
!  Zender 2003 eq. 3
   fd = ( 1.0 - &
          ( log( zom/zoms ) / &
            log( 0.35 * ( (0.1/zoms)**0.8 ) ) &
          ) &
        ) ** (-1.)

!  Calculate the threshold water content following Fecan et al. [1999]
!  Assumes a globally constant mass fraction of clay in soil, mclay
!  Zender 2003 eq. 5 (implicit a = 1)
   wt = 0.17 * mclay + 0.14 * mclay**2

!  Compute parameter for scaling horizontal to vertical mass flux,
!  alpha = [1/m], Zender 2003 eq. 11
   alpha = 100. * exp( (13.4*mclay - 6.)*log(10.) )

!  Calculate constant k for nonsaltating profile
!  Gillette et al. [1998] eq. 3
   k_z = 0.4 / log(10./zom)

!  Calculate the threshold velocity of wind erosion [m/s] for each radius
!  for a dry soil, as in Marticorena et al. [1997], eq. 1
!  Assumptions are a constant surface level air density (air_dens)
!  and we compute threshold for a single particle size/density (soil_diameter
!  and soil_density), equivalent to choosing a monomodal bed of potentially
!  saltating particles.
   u_thresh0 = 0.13 * sqrt(soil_density*grav*soil_diameter/air_dens) &
                    * sqrt(1.+6.e-7/(soil_density*grav*soil_diameter**2.5)) &
        / sqrt(1.928*(1331.*(100.*soil_diameter)**1.56+0.38)**0.092 - 1.)

!  Apply the drag partitioning correction
!  Equivalent to Marticorena [1997] eq. 4, where our fd = 1/(their)feff
   u_thresh_drag = u_thresh0*fd


!  Spatially dependent part of calculation
!  ---------------------------------------
   do j = j1, j2
    do i = i1, i2

     if ( oro(i,j) /= LAND ) cycle ! only over LAND gridpoints

     w10m = sqrt(u10m(i,j)**2.+v10m(i,j)**2.)

!    Modify the threshold depending on soil moisture as in Fecan 1999
!    Zender 2003 eq. 6
     if (gwettop(i,j) <= wt) then
       fw = 1.0
     else       
       fw = sqrt( 1. + 1.21 * (100. * (gwettop(i,j)-wt) )**0.68)
     endif
     u_thresh_drag_water = u_thresh_drag*fw 
       
!    Modify friction velocity for Owen Effect
!    Assumption of stable atmospheric profile to go from saltation
!    wind speed to equivalent threshold at z = 10m
!    Gillette et al. [1998] eq. 3
     wt10m = u_thresh_drag_water/k_z 
     if (w10m >= wt10m) then
       ustars = ustar(i,j) + 0.003*((w10m-wt10m)**2)  
     else
       ustars = ustar(i,j)
     endif 
      
!    Calculate the horizontal mass flux of dust [kg m-1 s-1]
!    Marticorena et al. 1997 eq. 5 
!    Note: differs from Zender et al. 2003 eq. 10
     rat = u_thresh_drag_water / ustars
     if (rat < 1.0) then
       horiz_flux = cs * air_dens * ustars**3 /grav * &
                     (1 - rat**2) * (1+rat)

     else
       horiz_flux = 0.0
     endif

!    Calculate the vertical mass flux of dust and scale to source [kg m-2 s-1]
     vert_flux = alpha * horiz_flux

     emissions(i,j) = vert_flux

    end do
   end do

   rc = 0

   end subroutine DustEmissionDEAD



!
! GOCART-based dust emission scheme (modified from Ginoux et al., JGR, 2001)
! Pass in a grid of wind speeds, surface characteristics, and an aerosol
! particle radius and returns the total dust emission mass flux [kg m-2 s-1].  
! Emissions need to be scaled by source function (ie., fraction of
! grid cell emitting), tuning factor, and fractional soil content of particle
! size.
!
! !IROUTINE:  DustEmissionGOCART - Compute the dust emissions
!
! !INTERFACE:
!

   subroutine DustEmissionGOCART( i1, i2, j1, j2, km, radius, &
                                  fraclake, gwettop, oro, u10m, v10m, &
                                  emissions, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km
   real                :: radius                   ! particle radius [m]
   real, pointer, dimension(:,:) :: fraclake, gwettop, oro, u10m, v10m

! !OUTPUT PARAMETERS:

   real     ::  emissions(i1:i2,j1:j2)             ! Local emission
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
   character(len=*), parameter :: myname = 'DustEmissionsGOCART'

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
!  29Dec2009, Colarco - Modifications to change calling
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j
   real, parameter ::  air_dens = 1.25  ! Air density = 1.25 kg m-3
   real, parameter ::  soil_density  = 2650.  ! km m-3
   real            ::  diameter         ! dust effective diameter [m]
   real            ::  u_thresh0
   real            ::  u_thresh
   real            ::  w10m

   _UNUSED_DUMMY(km)

!  Initialize local variables
!  --------------------------
   emissions(:,:) = 0.

!  Calculate the threshold velocity of wind erosion [m/s] for each radius
!  for a dry soil, as in Marticorena et al. [1997].
!  The parameterization includes the air density which is assumed 
!  = 1.25 kg m-3 to speed the calculation.  The error in air density is
!  small compared to errors in other parameters.
   diameter = 2. * radius
   u_thresh0 = 0.13 * sqrt(soil_density*grav*diameter/air_dens) &
                    * sqrt(1.+6.e-7/(soil_density*grav*diameter**2.5)) &
           / sqrt(1.928*(1331.*(100.*diameter)**1.56+0.38)**0.092 - 1.)

!if(mapl_am_i_root())print*,'DustEmissionGOCART u_thresh0 = ', u_thresh0

!  Spatially dependent part of calculation
!  ---------------------------------------
   do j = j1, j2
    do i = i1, i2

     if ( oro(i,j) /= LAND ) cycle ! only over LAND gridpoints

     w10m = sqrt(u10m(i,j)**2.+v10m(i,j)**2.)

!     This should give emissions equal to about 200 Tg month-1
!      if(gcDU%src(i,j) .lt. 1.) then
!       DU_emis(n)%data2d(i,j) = &
!           4.3064e-8*gcDU%sfrac(n)*gcDU%src(i,j)
!       w_c%qa(n1+n-1)%data3d(i,j,km) =  w_c%qa(n1+n-1)%data3d(i,j,km) &
!                         + DU_emis(n)%data2d(i,j)*cdt*grav/w_c%delp(i,j,km)
!      endif

!    Modify the threshold depending on soil moisture as in Ginoux et al. [2001]
     if(gwettop(i,j) .lt. 0.5) then
      u_thresh = amax1(0.,u_thresh0* &
       (1.2+0.2*alog10(max(1.e-3,gwettop(i,j)))))

       if(w10m .gt. u_thresh) then     
!       Emission of dust [kg m-2 s-1]
        emissions(i,j) = &
            (1.-fraclake(i,j)) * w10m**2. * (w10m-u_thresh)

!if(mapl_am_i_root())print*,'DustEmissionGOCART u_thresh = ', u_thresh
!if(mapl_am_i_root())print*,'DustEmissionGOCART w10m     = ', w10m

       endif
      endif

     end do   ! i
    end do    ! j

   rc = 0

   end subroutine DustEmissionGOCART


! Simple version of GOCART like scheme, where a threshold value is provided
! along with a wind speed and the emissions are returned.  No size information
! is implied.
!
! !IROUTINE:  DustEmissionSimple - Compute the dust emissions
!
! !INTERFACE:
!

   subroutine DustEmissionSimple( i1, i2, j1, j2, ut0, &
                                  fraclake, gwettop, oro, w, &
                                  tsoil, tsoil_freezing, emissions, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in)           :: i1, i2, j1, j2
   real                          :: ut0            ! dry threshold speed [m s-1]
   real, pointer, dimension(:,:) :: fraclake, &    ! fraction of grid cell in lake
                                    gwettop,  &    ! soil moisture
                                    oro,      &    ! orography flag
                                    w,        &    ! wind speed [m s-1]
                                    tsoil          ! soil surface temperature [K]
   real                          :: tsoil_freezing ! soil freezing point [K]

! !OUTPUT PARAMETERS:

   real                          :: emissions(i1:i2,j1:j2)  ! Local emission [kg m-2 s-1]
   integer, intent(out)          :: rc          ! Error return code:
                                                !  0 - all is well
                                                !  1 - 
   character(len=*), parameter :: myname = 'DustEmissionsSimple'

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
!  19Nov2012, Colarco - Introduced
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j
   integer  ::  n1, n2
   real     ::  emis(i1:i2,j1:j2)       ! Local bin emission
   real     ::  ut
   real     ::  w10m

!  Initialize local variables
!  --------------------------
   emissions(:,:) = 0.

!  Spatially dependent part of calculation
!  ---------------------------------------
   do j = j1, j2
    do i = i1, i2

     if ( oro(i,j) /= LAND ) cycle ! only over LAND gridpoints

!    Modify the threshold depending on soil moisture as in Ginoux et al. [2001]
     if ((gwettop(i,j) .lt. 0.5) .and. (tsoil(i,j) .gt. tsoil_freezing)) then
      ut = max(0.0, ut0 * (1.2 + 0.2*alog10(max(1.0e-3, gwettop(i,j)))) / &
                          (1.2 + 0.2*alog10(1.0e-3)))

       if(w(i,j) .gt. ut) then     
!       Emission of dust [kg m-2 s-1]
        emissions(i,j) = &
            (1.-fraclake(i,j)) * w(i,j)**2. * (w(i,j)-ut)

       endif
      endif

     end do   ! i
    end do    ! j

   rc = 0

   end subroutine DustEmissionSimple


   subroutine MAM_DustEmissionGOCART( i1, i2, j1, j2, km, &
                                      fraclake, gwettop, oro, u10m, v10m, &
                                      emission_total, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km
   real, pointer, dimension(:,:) :: fraclake, gwettop, oro, u10m, v10m

! !OUTPUT PARAMETERS:

   real, intent(inout)  :: emission_total(i1:i2,j1:j2)  ! Local emission
   integer, intent(out) :: rc    ! Error return code:
                                 !  0 - all is well
                                 !  1 - 

   character(len=*), parameter :: myname = 'MAM_DustEmissionGOCART'

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
!  11Oct2011, Darmenov - For now use the GOCART emission scheme to 
!                        calculate the total emission
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
  integer :: n
  real    :: emission(i1:i2, j1:j2)
  double precision :: emission_total_(i1:i2, j1:j2)

  integer, parameter :: nbins = 5
  real,    parameter :: radius(nbins) = (/0.73, 1.4, 2.4, 4.5, 8.0/) * 1e-6 ! [m]
  real,    parameter :: src_fraction(nbins) = (/0.1, 0.25, 0.25, 0.25, 0.25/)

  emission_total  = 0.0    ! total emission
  emission_total_ = 0.0d0  ! total emission: double precision

  do n = 1, nbins
      emission = 0.0
      call DustEmissionGOCART( i1, i2, j1, j2, km, radius(n), &
                               fraclake, gwettop, oro, u10m, v10m, &
                               emission, rc )

      emission_total_ = emission_total_ +  src_fraction(n) * emission
  end do

  emission_total = emission_total_

  end subroutine MAM_DustEmissionGOCART


  subroutine MAM_DustEmission( i1, i2, j1, j2, km, &
                               rLow, rUp, &
                               emission_bulk, &
                               emission_mass, emission_num, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km
   real, intent(in)    :: rLow, rUp      ! Dry particle bin edge radii [m]
   real, intent(in)    :: emission_bulk(i1:i2, j1:j2)

! !OUTPUT PARAMETERS:

   real, intent(inout)  :: emission_mass(i1:i2, j1:j2)    ! Local emission
   real, intent(inout)  :: emission_num (i1:i2, j1:j2)
   integer, intent(out) :: rc    ! Error return code:
                                 !  0 - all is well
                                 !  1 - 

   character(len=*), parameter :: myname = 'MAM_DustEmission'

! !DESCRIPTION: Computes the dust emissions for one time step
!
! !REVISION HISTORY:
!
!  11Oct2011, Darmenov - MAM
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer :: n
   real    :: pi

!! Initial size distribution of dust, D'Almeida[1987]
!  integer, parameter :: nmodes = 3
!  real, parameter    :: MMD(nmodes) = (/0.832, 4.82, 19.38/) * 1e-6   ! mass median diameter, [m]
!  real, parameter    :: sigma_g(nmodes) = (/2.1, 1.9, 1.6/)           ! geometric standard deviation
!  real, parameter    :: weight_mass(nmodes) = (/0.036, 0.957, 0.007/) ! mass weights

!! Initial size distribution of dust, Claquin[1999]  
!  integer, parameter :: nmodes = 3
!  real, parameter    :: MMD(nmodes) = (/0.011, 2.54, 42.10/) * 1e-6   ! mass median diameter, [m]
!  real, parameter    :: sigma_g(nmodes) = (/1.89, 2.0, 2.13/)         ! geometric standard deviation
!  real, parameter    :: weight_mass(nmodes) = (/2.6e-6, 0.78, 0.22/)  ! mass weights

!! Initial size distribution of dust, Alfaro and Gomes[2001], and Foret et al[2006] 
   integer, parameter :: nmodes = 3
   real, parameter    :: MMD(nmodes) = (/1.5, 6.7, 14.2/) * 1e-6       ! mass median diameter, [m]
   real, parameter    :: sigma_g(nmodes) = (/1.7, 1.6, 1.5/)           ! geometric standard deviation
   real, parameter    :: weight_mass(nmodes) = (/0.02, 0.27, 0.71/)    ! mass weights
   
   real, parameter    :: soil_density = 2650.0  ! km m-3

   real               :: NMD(nmodes)
   real               :: weight_number(nmodes), w_n(nmodes)
   real               :: integral_mass(nmodes), integral_number(nmodes)

   _UNUSED_DUMMY(km)

!  Initialize local variables
!  --------------------------
   rc = 0

   pi = 4 * atan(1.0)

   NMD = MMD * exp(-3 * log(sigma_g)**2)

   w_n = weight_mass / (pi/6 * soil_density * NMD**3 * exp(4.5*log(sigma_g)**2))
   weight_number = w_n / sum(w_n)

   ! compute the size distribution integrals 
   do n = 1, nmodes
       integral_mass(n)   = lognormal_integral(2*rLow, 2*rUp, MMD(n), sigma_g(n))
       integral_number(n) = lognormal_integral(2*rLow, 2*rUp, NMD(n), sigma_g(n))
   end do
    
   emission_mass = emission_bulk * sum(weight_mass * integral_mass)
   emission_num  = emission_bulk * sum(w_n * integral_number)

!  emission_mass(:,:) = 0.0
!  emission_num (:,:) = 0.0 
!
!  do n = 1, nmodes
!      emission_mass = emission_mass + emission_bulk * weight_mass(n) * integral_mass(n)
!      emission_num  = emission_num  + emission_bulk * w_n(n) * integral_number(n)
!  end do

   contains
       real function lognormal_cdf(x, median, sigma) result(cdf)
          implicit none
          real, intent(in) :: x, median, sigma
          real             :: erf

          cdf = 0.5 * (1 + erf(log(x/median) / (sqrt(2.0) * log(sigma))))
       end function lognormal_cdf

       real function lognormal_integral(x1, x2, median, sigma) result(integral)
          implicit none
          real, intent(in) :: x1, x2, median, sigma

          integral = lognormal_cdf(x2, median, sigma) - lognormal_cdf(x1, median, sigma)
       end  function lognormal_integral

   end subroutine MAM_DustEmission



! ------------------------------------------------------------------------------
! Algorithm to compute mass fraction of dust emissions in particle bin of radius
! range spanned by rlow and rup [m] based on brittle fragmentation theory, after
! Kok PNAS (2011), his equations 5 & 6
!
! !IROUTINE:  KokSizeDistribution -- pass mass fraction size bins at dust emission
!
! !INTERFACE:
!

   subroutine KokSizeDistribution( r, rlow, rup, dm, dn, rhod, rhog)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   real, intent(in), dimension(:)  :: r, rlow, rup  ! bin center, lower, and
                                                    ! upper radius [m]
   real, intent(in), optional      :: rhod, rhog    ! dust and group densities [kg m-3]

! !OUTPUT PARAMETERS:
   real, intent(out), dimension(:) :: dm            ! mass fraction of each bin
   real, intent(out), dimension(:), optional :: dn  ! number fraction of each bin

! !DESCRIPTION: Compute the dust emission mass fraction per size bin based on 
!               Kok PNAS (2011) brittle fragmentation theory
!
! !REVISION HISTORY:
!
!  13Apr2017, Colarco - Introduced
!  01Oct2018, Colarco - For CARMA it is possible for mixed group dust element
!               to have different density than mixed group "pc" element, which
!               requires a correction here in terms of the radius of a "pure"
!               dust particle with same mass as the mixed group particle.
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer         :: nbin, ibin
   real            :: erf
   real, parameter :: ds = 3.4e-6, &                ! soil median diameter [m]
                      sigma = 3., &                 ! width of soil distribution
                      lambda = 12.e-6, &            ! side crack propagation length
                      cv = 12.62e-6, &              ! some constants
                      cn = 0.9539e-6
   real            :: rho_dust, rho_grp             ! dust and group densities [kg m-3]

!  parameters for a fine sub-bin distribution
   integer, parameter :: nbin_ = 1000
   real               :: r_(nbin_), dr_(nbin_), dm_(nbin_), dn_(nbin_)
   real               :: rmrat_, rmin_

!  Set the particle densities that go with the radii passed in
   rho_dust = 2650.   ! nominal dust particle density
   if(present(rhod)) rho_dust = rhod
   rho_grp  = rho_dust
   if(present(rhog)) rho_grp  = rhog

!  Get number of major bins
   nbin = size(r)

!  Loop over the major size bins
   do ibin = 1, nbin

    rmrat_ = (rup(ibin)/rlow(ibin))**(3./nbin_)
    rmin_  = rlow(ibin)*((1.+rmrat_)/2.)**(1./3.)
    call carmabins(rmin_, rmrat_, r_, dr_)

!   If the radii passed in are for a "mixed" group with a different density
!   than "pure" dust we need to determine the radius of the "pure" dust
!   particle with the same mass as the "mixed" group particle of specified
!   radius, which is simply:
    r_  = r_*(rho_grp/rho_dust)**(1./3.)
    dr_ = dr_*(rho_grp/rho_dust)**(1./3.)

!   Calculate the mass fraction from Kok PNAS 2011 Eqn. 6
    dm_ = ( dr_/r_ ) * (2.*r_/cv) * &
           (1. + erf(alog(2.*r_/ds)/sqrt(2.)/alog(sigma)))*exp(-(2.*r_/lambda)**3)
    dm(ibin) = sum(dm_)
if(MAPL_AM_I_ROOT()) print *, 'DUST USING THESE BIN MASS FRACTIONS:', ibin, r(ibin), dm(ibin)

!   Calculate the number fraction from Kok PNAS 2011 Eqn. 5
    if(present(dn)) then
     dn_ = ( dr_/r_ ) * (1./cn/(2.*r_)**2) * &
            (1. + erf(alog(2.*r_/ds)/sqrt(2.)/alog(sigma)))*exp(-(2.*r_/lambda)**3)
     dn(ibin) = sum(dn_)
    endif
   enddo

   contains
    subroutine carmabins(rmin, rmrat, r, dr)
!   mimics the CARMA size bin routine, provides onlly r and dr and assumes
!   particle density = 1.
    real, intent(in) :: rmin, rmrat
    real, intent(out), dimension(:) :: r, dr

    integer         :: nbin, ibin
    real, parameter :: pi = MAPL_PI
    real            :: rmassmin, vrfact, rmass

    nbin = size(r)

    rmassmin = 4./3.*pi*rmin**3
    vrfact   = ( (3./2./pi / (rmrat+1.))**(1./3.))*(rmrat**(1./3.) - 1.)
    do ibin = 1, nbin
     rmass    = rmassmin*rmrat**(ibin-1)
     r(ibin)  = (rmass/(4./3.*pi))**(1./3.)
     dr(ibin) = vrfact*rmass**(1./3.)
    enddo

    end subroutine carmabins

   end subroutine KokSizeDistribution

end module 
