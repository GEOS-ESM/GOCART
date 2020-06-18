#include "unused_dummy.H"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  SmokeEmissionMod 
!
! This routine emits smoke particles as a partile size distribution.
! The user has control over the vertical distrubtion of the smoke
!  emission as well.
!
! Jamison A. Smith
!
! !INTERFACE:

Module SmokeEmissionMod

! !USES:

use Chem_Mod                     ! Chemistry Base Class
use Chem_ConstMod, only: grav    ! Constants
use Chem_UtilMod                 ! I/O
use m_mpout                      ! new to GEOS-5
use MAPL                     ! To test for rootproc

Implicit None

! !PUBLIC TYPES:

! This private statement keeps definitions from used mods from clashing
!  with definition in other mods, like CARMA_GridComp, that use this mod.

PRIVATE

! !PUBLIC MEMBER FUNCTIONS:

Public SmokeEmission

! !DESCRIPTION:
!
! !REVISION HISTORY:
!
!  2011Jan21 Adding standard deviation of MINX altitude observations, JAS
!  2010May06 Adjusted spatial mask for Siberian fire, JAS
!  2010Apr01 Final form, JAS
!  2009Feb03 J. A. Smith converting SM_GridCompMod (GEOS-4)
!             to SmokeEmissionMod (GEOS-5)
!  2007Aug28 J. A. Smith converting BC_GridCompMod (GOCART)
!             to SM_GridCompMod (CARMA)
!  
!EOP
!-------------------------------------------------------------------------

CONTAINS

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  SmokeEmission - Adds Smoke from passed argument emission
!                             Smoke is distributed across several size bins
!                             Smoke is also distributed vertically
!
! !INTERFACE:

Subroutine SmokeEmission ( i1, i2, j1, j2, km, &
                           nymd, nhms, dtime, &
                           rhoa, ple, & 
                           emission, &
                           profile, &
                           binfirst, binlast, &
                           r, & 
                           dr, & 
                           rhop, &
                           z_bl, &
                           z_lo, &
                           lon, lat, &
                           qa, rc )

! !USES:

Implicit None

! !INPUT PARAMETERS:

integer, intent(in) :: i1, i2, j1, j2, km         ! spatial indices for bounds
                                                  !  of subdomain
integer, intent(in) :: nymd, nhms                 ! date and time
real, intent(in)    :: dtime                      ! timestep
real, pointer, dimension(:,:,:)  :: rhoa, ple     ! air density in grid box
                                                  !  center and pressure at
                                                  !  grid box edges
real, pointer, dimension(:,:)  :: emission        ! emission rate
real, pointer, dimension(:,:,:)  :: profile       ! observed altitudes
integer, intent(in) :: binfirst, binlast          ! bounds for mass bins

! Are the indicies of r and dr going to be binfirst and binlast?
!  maybe check this with lbound and ubound?

real, pointer, dimension(:)    :: r    ! radius for mass grid read from rc
                                       !  file and computed via setupbins
real, pointer, dimension(:)    :: dr   ! bin width consistent with mass grid

real, pointer, dimension(:)    :: rhop ! density of smoke particles         

real, pointer :: z_bl(:,:)    ! altitude at the top of the boundary layer
real, pointer :: z_lo(:,:,:)  ! altitude of the lower edge of grid boxes

real, pointer, dimension(:,:)    :: lon  ! longitude
real, pointer, dimension(:,:)    :: lat  ! latitude

! !OUTPUT PARAMETERS:

integer, intent(out) :: rc    ! Error return code:
                              !  0 - all is well
                              !  1 - badness

! INOUT PARAMETERS

type( Chem_Array ), pointer  :: qa(:)  ! mixing ratio of smoke bins

! !DESCRIPTION: Computes the smoke emission for one timestep
!               Reads emission rate from a file
!               Distributes mass across several size bins
!               Distributes aerosol vertically as well, can emit into the
!                boundary layer or in elevated layers as seems more appropriate
!                for the Siberian fires of May 2003

! !REVISION HISTORY:
!
!  2010Apr01, J. A. Smith, final form
!  2010Feb02, J. A. Smith, adapting GEOS-4 SM_GridComp to GEOS-5
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables

character(len=255), parameter :: myname = 'SmokeEmission'

integer :: i, j, k       ! counter indices for spatial coords

integer            :: ibin        ! counter index for mass bins
integer, parameter :: nelem = 1   ! # of CARMA elements
integer, parameter :: ngroup = 1  ! # of CARMA groups

! vars for particle size distribution

real, parameter :: pi = 3.14159265358979323846

real :: sigma, rm        ! characteristics of log-normal size dist.
real :: factor           ! coefficients for size distribution

real :: dist_mass, dndlnr  ! total mass of size distribution with 1 particle
                           !  per cubic meter, the size distribution function
      
real, dimension(binfirst:binlast) :: dm, mass_frac  ! mass in each bin for size
                                                    !  dist., frac of total
                                                    !  mass that lies in each
                                                    !  bin

! vars for vertical distribtion of smoke plumes

logical :: elevated     ! Will want to run most smoke in the BL and
                        !  Siberian smoke elevated 

real :: xlo, xhi, ylo, yhi    ! bounds of lat/lon for Siberian spatial mask

!real :: ptop, pbot      ! pressure at top and bottom of a grid box

!integer, parameter      :: nlayer = 1   ! # of elevated smoke plume layers
!integer                 :: ilayer       ! counter index
!real, dimension(nlayer) :: p_layer      ! pressure levels (mbar)
                                        !  of elevated smoke plumes 

real, dimension(i1:i2,j1:j2)      :: sum_frac     ! sum of air dens for
                                                  !  smoke plume layers
                                                  !  for a column, used to
                                                  !  put equal mass mixing
                                                  !  ratios in each layer
real, dimension(i1:i2,j1:j2,1:km) :: alt_frac     ! fraction of smoke mass
                                                  !  emitted into a press.
                                                  !  level for a column

! MINX altitudes for Siberia coresponding to the dates of May 2003

real, dimension(31) :: MINXalt = &
 (/ -9999.99, 1164.85, 1946.29, 1983.82, 2102.06, 1329.31, 1987.82 & 
  , 1860.88, 1791.34, 4373.53, -9999.99, 3239.38, 1446.67, 2514.30 & 
  , 2422.74, 1923.87, 2946.53, 2021.50, 2341.06, 3191.82, 1629.18 & 
  , 3025.59, 1895.57, 1568.98, 1247.66, 1904.46, 1243.16, 1948.45 & 
  , 1159.08, 1744.56, 1534.30 /)

real, dimension(31) :: MINXdev = &
 (/ -9999.99, 982.614, 548.590, 839.832, 738.142, 506.484, 532.365 & 
  , 526.686, 393.708, 227.842, -9999.99, 282.753, 947.079, 551.947 & 
  , 597.896, 673.824, 1228.36, 833.620, 538.828, 1172.63, 770.214 & 
  , 1347.20, 655.367, 411.847, 669.359, 521.719, 293.923, 469.504 & 
  , 430.608, 494.765, 655.285 /)

integer :: date
real :: altitude           ! altitude to be used for a particular day
logical :: use_profile     ! If MINX observations exist, use them
real :: z

real, parameter :: smallval = 1.0e-37          ! FPE testing

_UNUSED_DUMMY(nymd)
_UNUSED_DUMMY(nhms)

!--
!  Initialize local variables

rc = 0

!if ( MAPL_AM_I_ROOT () ) then
!  write(*,*) 'Hey! I am SmokeEmission. Here me roar!'
!endif

!--
! Generate a log-normal size distribution for emitting the smoke bins
!
! rc file and setupbins are now being used to generate mass grid
! rc file uses cgs units

! Particle size distribution a la Westphal and Toon, 1991
!
! Log-normal characteristics

sigma = 1.45
rm = 1.e-5           ! cm

! The mass of the entire distribution

dist_mass = 0.

!if ( MAPL_AM_I_ROOT () ) then
!  write(*,*) 'JAS: Check 1'
!  print *, binfirst, binlast
!  print *, pi, rm, sigma, log(sigma)
!  print *, rhop, ubound(rhop,1), lbound(rhop,1)
!  print *, r, ubound(r,1), lbound(r,1)
!endif

do ibin = binfirst, binlast

  ! Size distribution with 1 particle / m**3

  dndlnr = 1. / ( 2. * pi ) ** 0.5 / log( sigma ) &
           * exp( - ( &
                      log( r(ibin) / rm ) ** 2 &
                      / &
                      ( 2 * log( sigma ) ** 2 ) &
                    ) &
                )

  ! Mass in each bin

  factor =  4. / 3. * pi * rhop(ibin)
  dm(ibin) = dndlnr * dr(ibin) * factor * r(ibin) ** 2

  ! Mass of the entire distribution

  dist_mass = dist_mass + dm(ibin)

enddo  ! loop over bins

!if ( MAPL_AM_I_ROOT () ) then
!  write(*,*) 'JAS: Check 2'
!endif

if ( dist_mass < smallval ) then
  write(*,*) 'JAS: dist_mass too small', dist_mass
  dist_mass = 42.0
  rc = 42
else

  ! Calc frac of mass of distribution within each bin

  do ibin = binfirst, binlast
    mass_frac(ibin) = dm(ibin) / dist_mass
  enddo

endif

!if ( MAPL_AM_I_ROOT () ) then
!  write(*,*) 'JAS: Check 3'
!endif

!--
! lon and lat vals for identifying Siberian region

xlo = 80. 
xhi = 140.
ylo = 44.
yhi = 64.

!--
! This code segment distributes the smoke emission into mulitple elevated
!  layers or into the boundary layer.
! Emitted mass is proportional to the air density in that layer, so that
!  plume parcels have the same mass mixing ratio of smoke aerosol.

! Elevated layers

!p_layer(1) = 756.  ! mbar
!p_layer(2) = 669.
!p_layer(3) = 614.
!p_layer(4) = 511.

! May want to shift these emission layers up a few mbars for sensitivity
! Trajectories indicate that Siberian smoke comes from pressure levels of
!  400 mbar or less

!p_layer = p_layer - 200.
  
!p_layer = p_layer * 1.e2    ! mbar => Pa

!--
! Compute the fraction of emission to emit at each level

sum_frac(:,:) = 0.
alt_frac(:,:,:) = 0.

date = date - 20030500

!--
! Is this code obsolete? This code was used when I didn't have a MINX
!  profile, only the mean and standard deviation

altitude = -9999.99
if ( date == 1 .or. date == 11 ) then
 altitude = -9999.99                   ! No valid plume height observations
else
 if ( date > 0 .and. date < 32 ) then
   altitude = MINXalt(date)            ! Use mean MINX obs for alt
 endif
endif

!--

!if ( MAPL_AM_I_ROOT () ) then
!  write(*,*) 'JAS: Check 4'
!endif

do j = j1, j2
  do i = i1, i2

    if ( lon(i,j) .gt. xlo .and. lon(i,j) .lt. xhi &
         .and. lat(i,j) .gt. ylo .and. lat(i,j) .lt. yhi ) then

      ! Siberian

!      elevated = .false.
      elevated = .true.

    else

      ! Other sources

      elevated = .false.
!      emission(i,j) = 0.

    endif   ! spatial mask for Siberian fire

    ! Lowest valid altitude on Earth is ~ -500.0 m

    if ( elevated .eqv. .true. .AND. altitude > -500.0 ) then

      use_profile = .false.       ! will set to true if valid obs are found

      do k = km, 1, -1

        ! ple has bounds of 0 and km, ple(:,:,0) is the top of the model
        !                             ple(:,:,km) is the surface

        !ptop = ple(i,j,k-1)
        !pbot = ple(i,j,k)

        ! JAS -- careful b/c two smoke layers may lie in the same model layer
        
        ! do ilayer = 1, nlayer
        !   if( p_layer(ilayer) .ge. ptop .and. p_layer(ilayer) .lt. pbot ) then
        !     alt_frac(i,j,k) = alt_frac(i,j,k) + rhoa(i,j,k)
        !     sum_frac(i,j) = sum_frac(i,j) + rhoa(i,j,k)
        !   endif
        ! enddo  ! ilayer

        ! JAS, shifting to geometric vertical coordinate for use with MINX
        !
        ! z_lo has bounds of 0 and km, z_lo(;,:,0) is the top of the model
        !                              z_lo(:,:,km) is the surface

        !if ( altitude >= z_lo(i,j,k) .AND. altitude < z_lo(i,j,k-1) ) then
        !  alt_frac(i,j,k) = rhoa(i,j,k)
        !  sum_frac(i,j) = rhoa(i,j,k)
        !else
        !  alt_frac(i,j,k) = 0.0
        !endif

        !--
        ! By ending up in this branch, we have a date with valid MINX
        !  observations
        !
        ! Not every Siberian column has valid MINX observations, however
        !
        ! Need to check for valid observations within each column.
        ! If they exist, then use the profile
        ! If they don't exist, spread the smoke evenly around the
        !  mean altitude for the region such that the smoke is
        !  confined to +/- 1 standard deviation from the regional mean

        if ( profile(i,j,k) > smallval ) use_profile = .true.

      enddo  ! k

      do k = km, 1, -1
        if ( use_profile .eqv. .true. ) then
          alt_frac(i,j,k) = rhoa(i,j,k) * profile(i,j,k)
          sum_frac(i,j) = sum_frac(i,j) + alt_frac(i,j,k)
        else
          z = 0.5 * ( z_lo(i,j,k) + z_lo(i,j,k-1) )
          if ( z >= MINXalt(date) - MINXdev(date) .and. &
               z <  MINXalt(date) + MINXdev(date) ) then
            alt_frac(i,j,k) = rhoa(i,j,k)
            sum_frac(i,j) = sum_frac(i,j) + alt_frac(i,j,k)
          endif
        endif
      enddo  ! k

    else 

!     or emit into the boundary layer

!     z_lo has bounds of 0 and km, z_lo(;,:,0) is the top of the model
!                                  z_lo(:,:,km) is the surface

      do k = km, 1, -1

        if ( z_lo(i,j,k) .lt. z_bl(i,j) ) then
          alt_frac(i,j,k) = rhoa(i,j,k)
          sum_frac(i,j) = sum_frac(i,j) + alt_frac(i,j,k)
        endif
 
      enddo  ! k

    endif  ! elevated or boundary layer

    !--
    !
    ! Caculate alt_frac

    ! protect against divide-by-zero error


    if ( sum_frac(i,j) > smallval ) then

      do k = km, 1, -1
        alt_frac(i,j,k) = alt_frac(i,j,k) / sum_frac(i,j)
      enddo

    else

      ! if sum_frac = 0., force all emission into surface layer
        
      alt_frac(i,j,:) = 0.0
      alt_frac(i,j,km) = 1.0

    endif

  enddo  ! i
enddo  ! j 

!--------------------------------------
!
! Update mixing ratios of CARMA species

!if ( MAPL_AM_I_ROOT () ) then
!  write(*,*) 'JAS: Check 5'
!  write(*,*) 
!  write(*,*) dtime, grav
!  write(*,*) 
!  write(*,*) 'Bounds of mass_frac and qa'
!  write(*,*) binfirst, binlast
!  write(*,*) lbound ( mass_frac, 1 ) , ubound ( mass_frac, 1 )
!  write(*,*) lbound ( qa, 1 ), ubound ( qa, 1 )
!  write(*,*) 
!  write(*,*) 'Bounds of emission'
!  write(*,*) i1, i2, j1, j2
!  write(*,*) lbound ( emission, 1 ) , ubound ( emission, 1 )
!  write(*,*) lbound ( emission, 2 ) , ubound ( emission, 2 )
!  write(*,*) 
!  write(*,*) 'Bounds of alt_frac'
!  write(*,*) i1, i2, j1, j2, 1, km
!  write(*,*) lbound ( alt_frac, 1 ) , ubound ( alt_frac, 1 )
!  write(*,*) lbound ( alt_frac, 2 ) , ubound ( alt_frac, 2 )
!  write(*,*) lbound ( alt_frac, 3 ) , ubound ( alt_frac, 3 )
!  write(*,*) 
!  write(*,*) 'Bounds of ple'
!  write(*,*) i1, i2, j1, j2, 1, km
!  write(*,*) lbound( ple, 1 ), ubound( ple, 1 )
!  write(*,*) lbound( ple, 2 ), ubound( ple, 2 )
!  write(*,*) lbound( ple, 3 ), ubound( ple, 3 )
!
!endif

do i = i1, i2
do j = j1, j2
do k = km, 1, -1

!  if ( MAPL_AM_I_ROOT () ) then
!    write(*,*) 'JAS: Check 6'
!    write(*,*) i, j, k, ple(i,j,k), ple(i,j,k-1)
!  endif

  if ( ( ple(i,j,k) - ple(i,j,k-1 ) ) < smallval ) then
    write(*,*) 'JAS: delta ple is too small' &
               , k, k-1, ple(i,j,k) - ple(i,j,k-1)
    rc = 43
  endif
  do ibin = binfirst, binlast

!    if ( MAPL_AM_I_ROOT () ) then
!      write(*,*) 'JAS: Check 7'
!      write(*,*) ibin, lbound ( qa, 1 ), ubound ( qa, 1 )
!      write(*,*) 
!      write(*,*) 'Checking qa(ibin)%data3d'
!      write(*,*) ibin, i, j, k
!      write(*,*) i, lbound ( qa(ibin)%data3d, 1 ), ubound ( qa(ibin)%data3d, 1 )
!      write(*,*) j, lbound ( qa(ibin)%data3d, 2 ), ubound ( qa(ibin)%data3d, 2 )
!      write(*,*) k, lbound ( qa(ibin)%data3d, 3 ), ubound ( qa(ibin)%data3d, 3 )
!      write(*,*)
!      write(*,*) 'Checking values'
!      write(*,*) qa(ibin)%data3d(i,j,k), mass_frac(ibin)
!      write(*,*) emission(i,j), alt_frac(i,j,k)
!      write(*,*) ple(i,j,k), ple(i,j,k-1)
!      write(*,*) 'Ending Check 7'
!    endif

      qa(ibin)%data3d(i,j,k) = qa(ibin)%data3d(i,j,k) &
                             + mass_frac(ibin) * emission(i,j) &
                             * alt_frac(i,j,k) &
                             * dtime * grav &
                             / ( ple(i,j,k) - ple(i,j,k-1) )
  enddo  ! ibin

enddo  ! k
enddo  ! j
enddo  ! i

!if ( MAPL_AM_I_ROOT () ) then
!  write(*,*) 'I am SmokeEmission, and I''m gonna take a nap.'
!endif

rc = 0

End Subroutine SmokeEmission

End Module SmokeEmissionMod

!---------------------------------------------------------------------------
!
! JAS: I'm storing old comments and code down here at the end of the routine

#define JAS
#ifndef JAS

!   Note: I (J. A. Smith) am modifying this routine extensively because
!         I want to move it from GEOS-4 to GEOS-5. The organization is
!         changing a lot as I try to simplify this thing and get away
!         from the old structure of having three routines to
!         1) Initialize
!         2) Run
!         3) Finalize
!         My impression is that the initialize and finalize steps are
!         occuring somewhere else in the code like the coupler routines
!         written by Colarco and Bardeen. - JAS 2010Feb04

!   use Chem_StateMod                         ! Chemistry State
!   use m_inpak90                             ! Resource file management
!   use m_die, only: die

!   real, parameter          :: SM_rhop = 1350.   ! MKS units for now
!   real, dimension(binlast-binfirst+1)    :: SM_radius
!   real, dimension(binlast-binfirst+1)    :: SM_mass
!   real                     :: rmrat
!   real                     :: dr
!   real                     :: vrfact

!real :: rhop = 1.35
!real, dimension(binfirst:binlast)    :: r
!real, dimension(binfirst:binlast)    :: dr

!   logical, parameter :: elevated = .false.
!   logical, parameter :: elevated = .true.

! + JAS -- I do not understand the next section
!
! Update emissions/production if necessary (daily)
!  -----------------------------------------------
!   if( gcSM%nymd .ne. nymd )then
!
!   Biomass Burning -- select on known inventories
!   ----------------------------------------------
!
!   Daily files (e.g., MODIS) or GFED v.2 (1997 - 2005 valid)
!
!    if(  index( gcSM%BiomassBurningFile,'%' ) .gt. 0 .or. &
!         index( gcSM%BiomassBurningFile,'gfed' ) .gt. 0 )then
!       nymd1 = nymd
!       nhms1 = 120000
!
!   Assume GFED climatology or Martin (Duncan) climatology
!
!    else
!       nymd1 = 1971*10000 + mod ( nymd, 10000 )  ! assumes 1971
!!       nymd1 = nymd
!       nhms1 = 120000
!    end if
!
!   Get the information from the biomass burning file
!
! fileread command went here orginally
!
!    gcSM%nymd = nymd
!
!   endif  ! gcSM%nymd .ne. nymd
!
! - JAS -- I do not understand the above section

! JAS, I ought to be using carma_bins
! JAS, carma_bins got lost in converstion to GEOS-5

! JAS, looks like I'm ignoring the .rc file for now
!      originally, we were using MKS -- now, rc file uses cgs
!      seeing how I'm ignoring the rc file for now, I can stay w/ MKS as
!        long as I'm consistent

!  if ( NBIN .eq. 16 ) then
!    SM_radius(1) = 25.e-9                 ! meter(m)
!    rmrat = 2.
!  elseif ( NBIN .eq. 8 ) then
!    SM_radius(1) = 25.e-9                 ! m
!    rmrat = exp( 15. * alog(2.) / 7. )    ! = 4.416358
!  else
!    print *, '*** Error ***'
!    print *, 'Chem_Shared/SmokeEmissionMod.F90: code not ready for this NBIN'
!    print *, 'NBIN = ', NBIN
!    rc = 2
!    return
!  endif

!  SM_mass(1) = factor * SM_radius(1) ** 3
!  do ibin = 2, nbin
!    SM_mass(ibin) = rmrat * SM_mass(ibin-1)
!    SM_radius(ibin) = ( SM_mass(ibin) / factor ) ** ( 1. / 3. )
!  enddo

! For bin width calc

!  vrfact = ( 3. / 2. / pi / ( rmrat + 1. ) ) ** ( 1. / 3. ) &
!             * ( rmrat ** ( 1. / 3. ) - 1. )

! Bin width

!    dr = vrfact &
!         * ( 4. / 3. * pi * SM_radius(ibin) ** 3 ) ** ( 1. / 3. )

!  k1 = lbound( z_lo, 3 )
!  k2 = ubound( z_lo, 3 )
!  print *, 'z_lo of k1', z_lo(i1,j1,k1), 'k1', k1 
!  print *, 'z_lo of k2', z_lo(i1,j1,k2), 'k2', k2 
!  k1 = lbound( ple, 3 )
!  k2 = ubound( ple, 3 )
!  print *, 'ple of k1', ple(i1,j1,k1), 'k1', k1 
!  print *, 'ple of k2', ple(i1,j1,k2), 'k2', k2 

!print *, 'i j', i, j, 'lon lat', lon(i,j), lat(i,j)

! JAS, dbg
!        if( i .eq. 1 .and. j .eq. 1 ) print *,'k = ', k, ' Ps @ k = ', &
!         ple(i,j,k), ple(i,j,k-1)

! JAS -- I am coopting Pete's code for finding the alitudes of various
!        layers to do my elevated layers code

!  Now update the tracer mixing ratios with the aerosol sources
!  ------------------------------------------------------------
!   p0 = ps

K_LOOP: do k = km, 1, -1

!    print *, 'SM_Emissions: getting emissions for layer ', k

!   First determine emissions for this layer
!   ----------------------------------------
    do j = j1, j2
     do i = i1, i2

!!! Gotta use ple now; find out the bounds of this !!!

      p1 = p0(i,j) - w_c%delp(i,j,k)

! !!! JAS: code for 4 elevated layers emitted
!

! Based on old vtracer code which is still here but commented

!           srcBiomass(i,j) = 0.
!
! JAS This next part looks like code for tracing the time of emission
!
! JAS if p is in between good src altitudes then
!     n = mod(date,8) + 1
!     MUST comment out previous vtracer stuff so that extra emissions
!     don't occur
!
!          n = 0
!          if ( ( p0(i,j) .ge. 600.e2 .and. p1 .lt. 600.e2 ) &
!               .or. ( p0(i,j) .ge. 400.e2 .and. p1 .lt. 400.e2 ) &
!               .or. ( p0(i,j) .ge. 300.e2 .and. p1 .lt. 300.e2 ) )then
!            n = mod( date + 5, 16 ) / 2 + 1
!          endif
!
! Emission of the tracers will be prescribed as 1 ppbm per time step.
!  This will permit me to see how many time steps air resides in an
!  emitting column.
!
!          if ( n .gt. 0 .and. n .le. nbins ) then
!            w_c%qa(nbeg+n-1)%data3d(i,j,k) = &
!              w_c%qa(nbeg+n-1)%data3d(i,j,k) + 1.e-9
!          endif

! JAS This next part looks like my poor attempt at emitting in 3 elevated
!     layers. This thing split the mass equally among the layers, making
!     the mixing ratio increase with altitude of the layer
!
!          if ( ( p0(i,j) .ge. 600.e2 .and. p1 .lt. 600.e2 ) &
!               .or. ( p0(i,j) .ge. 400.e2 .and. p1 .lt. 400.e2 ) &
!               .or. ( p0(i,j) .ge. 300.e2 .and. p1 .lt. 300.e2 ) )then
!
!            srcBiomass(i,j) = 1. / 3. *eBiomass*gcSM%biomass_src(i,j)
!
!          endif

! !!! JAS: end of elevated layers code

!     Update pressure of lower level
!     ------------------------------
      p0(i,j) = p1

      end do ! i
     end do  ! j

   end do K_LOOP

#endif

! endif of C-preprocessor statement to comment out old code
