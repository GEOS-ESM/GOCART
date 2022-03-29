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
   real(kind=DP) :: qdel, qsrc, d_p, dpm1

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
   allocate(dz, mold=rhoa);
   allocate(dzd(i2,j2,km), vsd(i2,j2,km), qa(i2,j2,km), vsettle(i2,j2,km), qa_temp(i2,j2,km))
   allocate(cmass_before(i2,j2), cmass_after(i2,j2))

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
               dpm1 = delp(i,j,k-1)
               qsrc = qdel * dpm1 / d_p
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
