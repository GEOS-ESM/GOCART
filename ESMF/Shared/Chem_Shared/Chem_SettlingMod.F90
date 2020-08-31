#include "unused_dummy.H"
!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_SettlingMod --- Gravitional Sedimentation & Settling Speed
!
! !INTERFACE:
!

   module  Chem_SettlingMod

! !USES:

   use Chem_Mod
   use Chem_ConstMod, only: grav        ! Constants !
   use Chem_UtilMod

   use m_mpout

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  Chem_Settling
   PUBLIC  Chem_SettlingSimple
   PUBLIC  Chem_CalcVsettle

!
! !DESCRIPTION:
!
!  This module implements the (pre-ESMF) DU Grid Component. 
!
! !REVISION HISTORY:
!
!  16Sep2003 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_Settling - 
!
! !INTERFACE:
!

   subroutine Chem_Settling ( i1, i2, j1, j2, km, nbeg, nend, nbins, flag, &
                              radiusInp, rhopInp, cdt, w_c, tmpu, rhoa, &
                              hsurf, hghte, fluxout, rc, &
                              vsettleOut, correctionMaring )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, nbeg, nend, nbins
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt
   real, pointer, dimension(:)     :: radiusInp, rhopInp
   real, pointer, dimension(:,:)   :: hsurf
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, hghte

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields
   type(Chem_Array), pointer        :: fluxout(:) ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 
!  Optionally output the settling velocity calculated
   type(Chem_Array), pointer, optional, dimension(:)  :: vsettleOut

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
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, iit, n
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]
   real :: vsettle(i1:i2,j1:j2,km)               ! fall speed [m s-1]
   real*8 :: cmass_before(i1:i2,j1:j2)
   real*8 :: cmass_after(i1:i2,j1:j2)
   real :: diff_coef                 ! Brownian diffusion coefficient [m2 s-1]
   real*8 :: qdel(i1:i2,j1:j2), qsrc(i1:i2,j1:j2), dp(i1:i2,j1:j2), dpm1(i1:i2,j1:j2)
   real :: dz(i1:i2,j1:j2,km)                    ! layer thickness [m]
   real*8 :: dzd(i1:i2,j1:j2,km), vsd(i1:i2,j1:j2,km), qa(i1:i2,j1:j2,km)

!  The following parameters relate to the swelling of seasalt like particles
!  following Fitzgerald, Journal of Applied Meteorology, 1975.
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2

!  parameter from Gerber 1985 (units require radius in cm, see rcm)
   real :: rcm
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424
!  parameters for ammonium sulfate
   real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428


!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

!
   real :: sat, rrat
   real :: radius, rhop   ! particle radius and density passed to
                          ! fall velocity calculation
   real    :: minTime, qmin, qmax
   integer :: nSubSteps, dk, ijl
   real*8  :: dt_settle, g

   _UNUSED_DUMMY(nend)
   rc = 0

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   g = grav

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

!  Loop over the number of dust bins
   do n = 1, nbins

    qa = w_c%qa(nbeg+n-1)%data3d

    radius = radiusInp(n)
    rhop = rhopInp(n)

!   Reset a (large) minimum time to cross a grid cell in settling
    minTime = cdt

    if( associated(fluxout(n)%data2d) ) fluxout(n)%data2d(i1:i2,j1:j2) = 0.0
    cmass_before(:,:) = 0.d0
    cmass_after(:,:) = 0.d0

!   If radius le 0 then get out of loop
    if(radius .le. 0.) cycle

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

!      Find the column dry mass before sedimentation
       cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k)/g * w_c%delp(i,j,k)

!      Adjust the particle size for relative humidity effects
       sat = max(w_c%rh(i,j,k),tiny(1.0)) ! to avoid zero FPE

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
        radius = alpha * radiusInp(n)**beta
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 2) then   ! Gerber
        sat = min(0.995,sat)
        rcm = radiusInp(n)*100.
        radius = 0.01 * (   c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                          + rcm**3.)**(1./3.)
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 3) then   
!       Gerber parameterization for Ammonium Sulfate
        sat = min(0.995,sat)
        rcm = radiusInp(n)*100.
        radius = 0.01 * (   SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                      + rcm**3.)**(1./3.)
        rrat = (radiusInp(n)/radius)**3.
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       elseif(flag .eq. 4) then
!       Petters and Kreidenweis (ACP2007) parameterization
        sat = min(0.99,sat)
        radius = (radiusInp(n)**3 * (1+1.19*sat/(1-sat)))**(1./3.)
        rrat = (radiusInp(n)/radius)**3
        rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
       endif

!      Calculate the settling velocity
       call Chem_CalcVsettle(radius, rhop, rhoa(i,j,k), &
                        tmpu(i,j,k), diff_coef, vsettle(i,j,k))
      end do
     end do
    end do

    if(present(correctionMaring)) then
     if ((correctionMaring) .and. (radiusInp(n) .le. (0.5*diameterMaring))) then
       vsettle = max(1.0e-9, vsettle - v_upwardMaring)
     endif
    endif

    vsd = vsettle

    if(present(vsettleOut)) then
     vsettleOut(n)%data3d = vsettle
    endif

!   Determine global max/min time to cross grid cell
    call pmaxmin ( 'Chem_Settling: dt', dz(i1:i2,j1:j2,1:km)/vsettle(i1:i2,j1:j2,1:km), &
                                        qmin, qmax, ijl, km, 0. )
    minTime = min(minTime,qmin)


!   Now, how many iterations do we need to do?
    if ( minTime < 0 ) then
         nSubSteps = 0
         call mpout_log(myname,'no Settling because minTime = ', minTime )
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

     qdel = qa(i1:i2,j1:j2,1)*dt_settle*vsd(i1:i2,j1:j2,1)/dzd(i1:i2,j1:j2,1)
     qa(i1:i2,j1:j2,1) = qa(i1:i2,j1:j2,1) - qdel

     do k = 2, km
      dp   = w_c%delp(i1:i2,j1:j2,k)
      dpm1 = w_c%delp(i1:i2,j1:j2,k-1)
      qsrc = qdel * dpm1 / dp
      qdel = qa(i1:i2,j1:j2,k)*dt_settle*vsd(i1:i2,j1:j2,k)/dzd(i1:i2,j1:j2,k)
      qa(i1:i2,j1:j2,k) = qa(i1:i2,j1:j2,k) - qdel + qsrc
     enddo

!     An alternative accumulator approach to computing the outgoing flux
!     if( associated(fluxout(n)%data2d) ) then
!        fluxout(n)%data2d = fluxout(n)%data2d + qdel * pdog/grav / dt_settle
!     endif

    end do  ! iit

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k)/ g * w_c%delp(i,j,k)
      enddo
     enddo
    enddo

    if( associated(fluxout(n)%data2d) ) then
       fluxout(n)%data2d(i1:i2,j1:j2) &
        = (cmass_before(i1:i2,j1:j2) - cmass_after(i1:i2,j1:j2))/cdt
    endif

    w_c%qa(nbeg+n-1)%data3d = qa

   end do   ! n

 end  subroutine Chem_Settling


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_SettlingSimple - support for single bin settling call
!
! !INTERFACE:
!

   subroutine Chem_SettlingSimple ( i1, i2, j1, j2, km, ibin, flag, &
                              radiusInp, rhopInp, cdt, w_c, tmpu, rhoa, &
                              hsurf, hghte, fluxout, rc, &
                              vsettleOut, correctionMaring )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   integer, intent(in) :: i1, i2, j1, j2, km, ibin
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt, radiusInp, rhopInp
   real, pointer, dimension(:,:)   :: hsurf
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, hghte

! !OUTPUT PARAMETERS:

   type(Chem_Bundle), intent(inout) :: w_c        ! Chemical tracer fields
   type(Chem_Array), pointer        :: fluxout    ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 
!  Optionally output the settling velocity calculated
   type(Chem_Array), pointer, optional  :: vsettleOut

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
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   integer  ::  i, j, k, iit
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]
   real :: vsettle(i1:i2,j1:j2,km)               ! fall speed [m s-1]
   real*8 :: cmass_before(i1:i2,j1:j2)
   real*8 :: cmass_after(i1:i2,j1:j2)
   real :: diff_coef                 ! Brownian diffusion coefficient [m2 s-1]
   real*8 :: qdel(i1:i2,j1:j2), qsrc(i1:i2,j1:j2), dp(i1:i2,j1:j2), dpm1(i1:i2,j1:j2)
   real :: dz(i1:i2,j1:j2,km)                    ! layer thickness [m]
   real*8 :: dzd(i1:i2,j1:j2,km), vsd(i1:i2,j1:j2,km), qa(i1:i2,j1:j2,km)

!  The following parameters relate to the swelling of seasalt like particles
!  following Fitzgerald, Journal of Applied Meteorology, 1975.
   real, parameter :: epsilon = 1.   ! soluble fraction of deliqeuscing particle
   real, parameter :: alphaNaCl = 1.35
   real :: alpha, alpha1, alpharat, beta, theta, f1, f2

!  parameter from Gerber 1985 (units require radius in cm, see rcm)
   real :: rcm
   real, parameter :: c1=0.7674, c2=3.079, c3=2.573e-11, c4=-1.424
!  parameters for ammonium sulfate
   real, parameter :: SU_c1=0.4809, SU_c2=3.082, SU_c3=3.110e-11, SU_c4=-1.428


!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

!
   real :: sat, rrat
   real :: radius, rhop   ! particle radius and density passed to
                          ! fall velocity calculation
   real    :: minTime, qmin, qmax
   integer :: nSubSteps, dk, ijl
   real*8  :: dt_settle, g

   rc = 0

   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )
   g = grav

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

   qa = w_c%qa(ibin)%data3d

    radius = radiusInp
    rhop = rhopInp

!   Reset a (large) minimum time to cross a grid cell in settling
    minTime = cdt

    if( associated(fluxout%data2d) ) fluxout%data2d(i1:i2,j1:j2) = 0.0
    cmass_before(:,:) = 0.d0
    cmass_after(:,:) = 0.d0

!   If radius le 0 then get out of loop
    if(radius .le. 0.) return

    do k = 1, km
     do j = j1, j2
      do i = i1, i2

!      Find the column dry mass before sedimentation
       cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k)/g * w_c%delp(i,j,k)

!      Adjust the particle size for relative humidity effects
       sat = max(w_c%rh(i,j,k),tiny(1.0)) ! to avoid zero FPE

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
       call Chem_CalcVsettle(radius, rhop, rhoa(i,j,k), &
                        tmpu(i,j,k), diff_coef, vsettle(i,j,k))
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
     if(associated(vsettleOut%data3d)) vsettleOut%data3d = vsettle
    endif

!   Determine global max/min time to cross grid cell
    call pmaxmin ( 'Chem_Settling: dt', dz(i1:i2,j1:j2,1:km)/vsettle(i1:i2,j1:j2,1:km), &
                                        qmin, qmax, ijl, km, 0. )
    minTime = min(minTime,qmin)


!   Now, how many iterations do we need to do?
    if ( minTime < 0 ) then
         nSubSteps = 0
         call mpout_log(myname,'no Settling because minTime = ', minTime )
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

     qdel = qa(i1:i2,j1:j2,1)*dt_settle*vsd(i1:i2,j1:j2,1)/dzd(i1:i2,j1:j2,1)
     qa(i1:i2,j1:j2,1) = qa(i1:i2,j1:j2,1) - qdel

     do k = 2, km
      dp   = w_c%delp(i1:i2,j1:j2,k)
      dpm1 = w_c%delp(i1:i2,j1:j2,k-1)
      qsrc = qdel * dpm1 / dp
      qdel = qa(i1:i2,j1:j2,k)*dt_settle*vsd(i1:i2,j1:j2,k)/dzd(i1:i2,j1:j2,k)
      qa(i1:i2,j1:j2,k) = qa(i1:i2,j1:j2,k) - qdel + qsrc
     enddo

    end do  ! iit

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k)/ g * w_c%delp(i,j,k)
      enddo
     enddo
    enddo

    if( associated(fluxout%data2d) ) then
       fluxout%data2d(i1:i2,j1:j2) &
        = (cmass_before(i1:i2,j1:j2) - cmass_after(i1:i2,j1:j2))/cdt
    endif

    w_c%qa(ibin)%data3d = qa

 end  subroutine Chem_SettlingSimple


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Chem_CalcVsettle - Calculate the aerosol settling velocity
!
! !INTERFACE:
!

   subroutine Chem_CalcVsettle ( radius, rhop, rhoa, tmpu, &
                                 diff_coef, vsettle )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real, intent(in)    :: radius              ! Particle radius [m]
   real, intent(in)    :: rhop                ! Particle density [kg m-3]
   real, intent(in)    :: rhoa                ! Layer air density [kg m-3]
   real, intent(in)    :: tmpu                ! Layer temperature [K]

! !OUTPUT PARAMETERS:

   real, intent(out)   :: diff_coef               ! Brownian diffusion 
                                                  ! coefficient [m2 s-1]
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
   real*8 rmu                       ! Dynamic viscosity [kg m-1 s-1]
   real*8 vt                        ! Thermal velocity of air molecule [m s-1]
   real*8 rmfp                      ! Air molecule mean free path [m]
   real*8 bpm                       ! Cunningham slip correction factor
   real*8 rkn                       ! Knudsen number
   real*8 re, x, y                  ! reynold's number and parameters
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]
   real, parameter :: pi = 3.141529265

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


   end subroutine Chem_CalcVsettle

   end module Chem_SettlingMod
