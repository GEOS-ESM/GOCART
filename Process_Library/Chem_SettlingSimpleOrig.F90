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
        dp, dpm1, qsrc

!EOP
!-------------------------------------------------------------------------
!  Begin

   gravDP = grav

   i2 = ubound(hghte,1)
   j2 = ubound(hghte,2)

   hsurf => hghte(i1:i2,j1:j2,km)

!  Allocate arrays
!  ---------------
   allocate(dz, mold=rhoa);
   allocate(dzd(i2,j2,km), vsd(i2,j2,km), qa(i2,j2,km), vsettle(i2,j2,km), qa_temp(i2,j2,km))
   allocate(cmass_before(i2,j2), cmass_after(i2,j2), qdel(i2,j2), dp(i2,j2), &
            dpm1(i2,j2), qsrc(i2,j2))

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
      dp   = delp(i1:i2,j1:j2,k)
      dpm1 = delp(i1:i2,j1:j2,k-1)
      qsrc = qdel * dpm1 / dp
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
