
!-------------------------------------------------------------------------
!
! !MODULE:  GriddedEmission.F90 --- Calculate the gridded emissions
!
! !INTERFACE:
!

   module  GOCART2G_Process

! !USES:
!  Only instrinsic fortran types and functions are allowed.
   use Chem_MieTableMod2G

   implicit none

! !PUBLIC TYPES:
!
   private

!
! !PUBLIC MEMBER FUNCTIONS:
!
   public DustEmissionGOCART2G
   public DistributePointEmission
   public updatePointwiseEmissions
   public Chem_Settling2G
   public Chem_Settling2Gorig
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

   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   integer, parameter     :: DP=kind(1.0d0)


!
! !DESCRIPTION:
!
!  This module implements the gridded emission calculations
!
! !REVISION HISTORY:
!
!  11Feb2020  E.Sherman, A.da Silva, T.Clune, A.Darmenov -  First attempt at refactor
!
!EOP
!-------------------------------------------------------------------------
CONTAINS

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
   real, pointer, dimension(:,:), intent(in) :: fraclake, gwettop, oro, u10m, v10m, du_src
   real, intent(in) :: Ch_DU, grav

! !OUTPUT PARAMETERS:
   real  ::  emissions(:,:)    ! Local emission
   integer, intent(out) :: rc                   ! Error return code:


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
   real, allocatable ::  emissions_(:,:)

!EOP
!-------------------------------------------------------------------------
!  Begin

!  Initialize local variables
!  --------------------------
   emissions(:,:) = 0.
   rc = 824

!  Get dimensions
!  ---------------
   nbins = size(radius)
   dims = shape(u10m)
   i1 = 1; j1 = 1
   i2 = dims(1); j2 = dims(2)

   allocate(emissions_(i2,j2))

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

      emissions_(:,:) = 0.

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
                  emissions_(i,j) = (1.-fraclake(i,j)) * w10m**2. * (w10m-u_thresh)
               endif
            endif !(gwettop(i,j) .lt. 0.5)
         end do ! i
      end do ! j
      emissions = emissions + (Ch_DU * du_src * emissions_)
    end do ! n
 
   rc=0

   end subroutine DustEmissionGOCART2G

!==================================================================================
!BOP

! !IROUTINE: updatePointwiseEmissions

  subroutine updatePointwiseEmissions (km, pBase, pTop, pEmis, nPts, pStart, &
                                       pEnd, airdens, delp, GRAV, area, &
                                       iPoint, jPoint, nhms, emissions_point, rc)

    implicit none

!   !ARGUMENTS:
    integer,                              intent(in   )  :: km
    real, dimension(:),                   intent(in   )  :: pBase
    real, dimension(:),                   intent(in   )  :: pTop
    real, dimension(:),                   intent(in   )  :: pEmis
    integer,                              intent(in   )  :: nPts
    integer, dimension(:),                intent(in   )  :: pStart
    integer, dimension(:),                intent(in   )  :: pEnd
    real, dimension(:,:,:),               intent(in   )  :: airdens
    real, dimension(:,:,:),               intent(in   )  :: delp
    real,                                 intent(in   )  :: GRAV
    real, dimension(:,:),                 intent(in   )  :: area
    integer, pointer, dimension(:),       intent(in   )  :: iPoint
    integer, pointer, dimension(:),       intent(in   )  :: jPoint
    integer,                              intent(in   )  :: nhms
    real, dimension(:,:,:),               intent(inout)  :: emissions_point
    integer, optional,                    intent(  out)  :: rc

!   !Local
    real, dimension(km)          :: point_column_emissions
    integer                           :: n, i, j

!   Description: Returns 3D array of pointwise emissions.
!
!   Revision History:
!EOP
!-----------------------------------------------------------------------------
!    Begin...

        do  n = 1, nPts
            i = iPoint(n)
            j = jPoint(n)
            if( i<1 .OR. j<1 ) cycle    ! Point emission not in this sub-domain
            ! Emissions not occurring in current time step
            if(nhms < pStart(n) .or. nhms >= pEnd(n)) cycle
            call DistributePointEmission(km, delp(i,j,:), airdens(i,j,:), pBase(n), &
                                         pTop(n), GRAV, pEmis(n), area(i,j), &
                                         point_column_emissions, rc)

            emissions_point(i,j,:) = point_column_emissions
        end do



!    RETURN_(ESMF_SUCCESS)

  end subroutine updatePointwiseEmissions

!==================================================================================
!BOP
! !IROUTINE: DistributePointEmissions

! !INTERFACE:

   subroutine DistributePointEmission(km, delp, rhoa, z_bot, z_top, &
                                      GRAV, emissions_point, area, &
                                      point_column_emissions, rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)             :: km
   real, dimension(:), intent(in)  :: delp, rhoa
   real, intent(in)                :: z_bot, z_top
   real,               intent(in)  :: GRAV, emissions_point, area


! !OUTPUT PARAMETERS:
   real, dimension(:), intent(out) ::  point_column_emissions

   integer, intent(out) :: rc                                 ! Error return code:


! !DESCRIPTION: Distributes piont emissions
!
! !REVISION HISTORY:
! ??? P. Colarco
! 16March2020 E.Sherman - Moved from original DU_GridCompMod.F90 and generalized
!
! !Locals
   integer :: k
   integer :: k_bot, k_top
   real    :: z_
   real, dimension(km) :: z, dz, w_

!EOP
!-------------------------------------------------------------------------
! Begin

!   find level height
    z = 0.0
    z_= 0.0

    do k = km, 1, -1
       dz(k) = delp(k)/rhoa(k)/grav
       z_    = z_ + dz(k)
       z(k)  = z_
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

!   distribute emissions in the vertical
!    point_column_emissions(:) = (w_ / sum(w_)) * emissions_point
    point_column_emissions(:) = ((w_ / sum(w_)) * emissions_point) / area

    rc = 0

    end subroutine DistributePointEmission
!==================================================================================

!BOP
!
! !IROUTINE:  Chem_Settling - 
!
! !INTERFACE:
!
   subroutine Chem_Settling2G (km, flag, int_qa, grav, delp, &
                               radiusInp, rhopInp, cdt, tmpu, rhoa, &
                               rh, hghte, fluxout, rc, &
                               vsettleOut, correctionMaring)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: grav     ! mapl_grav
   real, intent(in)    :: cdt
   real, dimension(:), intent(in)  :: radiusInp, rhopInp
   real, pointer, dimension(:,:,:) :: tmpu, rhoa, rh, hghte, delp
   real, dimension(:,:,:,:), intent(inout) :: int_qa

! !OUTPUT PARAMETERS:
   real, pointer, dimension(:,:,:)  :: fluxout ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 
!  Optionally output the settling velocity calculated
   real, pointer, optional, dimension(:,:,:,:)  :: vsettleOut

!  Optionally correct the settling velocity following Maring et al, 2003
   logical, optional, intent(in)    :: correctionMaring

   character(len=*), parameter :: myname = 'Chem_Settling2G'

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop) 
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the 
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
!  05May2020  Sherman   Refactor - only uses intrinsic fortran types/functions
!  15May2019  Darmenov  Refactor and speed up code
!  17Sep2004  Colarco   Strengthen sedimentation flux out at surface
!                       by setting removal to be valid from middle of
!                       surface layer
!  06Nov2003  Colarco   Based on Ginoux
!  23Jan2003  da Silva  Standardization
!

! !Local Variables
   integer  ::  i, j, k, n, dk

   integer         :: i1=1, i2, j1=1, j2, nbins
   integer         :: dims(3)

!   integer, parameter     :: DP=kind(1.0d0)
   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

   real, dimension(:,:,:), allocatable  :: vsettle   ! fall speed [m s-1]
   real, dimension(:,:,:), allocatable  :: radius, rhop ! output for ParticleSwelling
   real, dimension(:,:,:), allocatable  :: dz, qa

   real(kind=DP), dimension(:,:), allocatable  :: cmass_before, cmass_after
   
   real, allocatable, dimension(:,:)    :: hsurf

!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

!   real, parameter :: one_over_g = 1.0/grav

   real, dimension(:,:,:), allocatable  :: qa_temp
   real(kind=DP) ::  gravDP
   rc = 0

!EOP
!-------------------------------------------------------------------------
!  Begin

!  Get dimensions
!  ---------------
   nbins = size(radiusInp)
   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

   gravDP = grav

!  Allocate arrays
!  ---------------
   allocate(dz(i1:i2,j1:j2,km), vsettle(i1:i2,j1:j2,km), radius(i1:i2,j1:j2,km), qa(i1:i2,j1:j2,km), &
            rhop(i1:i2,j1:j2,km), qa_temp(i1:i2,j1:j2,km))
   allocate(cmass_before(i1:i2,j1:j2), cmass_after(i1:i2,j1:j2), hsurf(i1:i2,j1:j2))

   hsurf = hghte(i1:i2,j1:j2,km)

!print*,'sum(tmpu) = ',sum(tmpu)
!print*,'sum(rhoa) = ',sum(rhoa)
!print*,'sum(rh) = ',sum(rh)
!print*,'sum(hghte) = ',sum(hghte)
!print*,'sum(delp) = ',sum(delp)


!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1
  
!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo

!  Loop over the number of dust bins
   do n = 1, nbins

    ! alias
    qa(:,:,:) = int_qa(:,:,:,n)

!    if( associated(fluxout(n)%data2d) ) fluxout(n)%data2d(i1:i2,j1:j2) = 0.0
    if(associated(fluxout)) fluxout(:,:,n) = 0.0

    cmass_before(:,:) = 0.d0
    cmass_after(:,:)  = 0.d0

!   If radius le 0 then get out of loop
    if(radiusInp(n) .le. 0.) cycle

!   Find the column dry mass before sedimentation
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k)/grav * delp(i,j,k)
      enddo
     enddo
    enddo
!    cmass_before = sum(qa * delp * 1.0/grav, 3)

!   Particle swelling
    call ParticleSwelling(i1, i2, j1, j2, km, rh, radiusInp(n), rhopInp(n), radius, rhop, flag)

!print*,'sum(radius) = ',sum(radius)
!print*,'sum(rhop) = ',sum(rhop)


!   Settling velocity of the wet particle
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       call Chem_CalcVsettle2G(radius(i,j,k), rhop(i,j,k), rhoa(i,j,k), &
                               tmpu(i,j,k), grav, vsettle(i,j,k))
      end do
     end do
    end do


    if(present(correctionMaring)) then
     if (correctionMaring) then
       vsettle = max(1.0e-9, vsettle - v_upwardMaring)
     endif
    endif

!print*,'sum(vsettle) = ',sum(vsettle)

    if(present(vsettleOut)) then
!     vsettleOut(n)%data3d = vsettle
     vsettleOut(:,:,:,n) = vsettle
    endif

qa_temp = qa
!print*,'before sum(qa) = ',sum(qa)

!   Time integration
    call SettlingSolver(i1, i2, j1, j2, km, cdt, delp, dz, vsettle, qa)

!print*,'after sum(qa)  = ',sum(qa)


!if (sum(qa_temp) < sum(qa)) then
!   print*,'< qa_temp = ',sum(qa_temp), ' : qa = ', sum(qa)
!end if

!if (sum(qa_temp) == sum(qa)) then
!   print*,'== qa_temp = ',sum(qa_temp), ' : qa = ', sum(qa)
!end if

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = 1, km
     do j = j1, j2
      do i = i1, i2
       cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k)/grav * delp(i,j,k)
      enddo
     enddo
    enddo
!    cmass_after = sum(qa * delp * 1.0/grav, 3)

    if( associated(fluxout) ) then
       fluxout(i1:i2,j1:j2,n) = (cmass_before(i1:i2,j1:j2) - cmass_after(i1:i2,j1:j2))/cdt

!if(sum(cmass_before) < sum(cmass_after)) then
!  print*,'< cmass_before = ',sum(cmass_before), ' : cmass_after = ', sum(cmass_after)
!endif

!if(sum(cmass_before) ==  sum(cmass_after)) then
! print*,'== cmass_before = ',sum(cmass_before), ' : cmass_after = ', sum(cmass_after)
!endif

    endif

!print*,'sum(cmass_before) = ',sum(cmass_before)
!print*,'sum(cmass_after) = ',sum(cmass_after)

!print*,'sum(fluxout(:,:,n)) = ',sum(fluxout(:,:,n))

    int_qa(:,:,:,n) = qa
   end do   ! n

 end subroutine Chem_Settling2G





!==================================================================================

!BOP
!
! !IROUTINE:  Chem_CalcVsettle - Calculate the aerosol settling velocity
!
! !INTERFACE:
   subroutine Chem_CalcVsettle2G ( radius, rhop, rhoa, tmpu, grav, &
                                 vsettle )
! !USES:
  implicit NONE

! !INPUT PARAMETERS:
   real, intent(in)    :: radius              ! Particle radius [m]
   real, intent(in)    :: rhop                ! Particle density [kg m-3]
   real, intent(in)    :: rhoa                ! Layer air density [kg m-3]
   real, intent(in)    :: tmpu                ! Layer temperature [K]
   real, intent(in)    :: grav

! !OUTPUT PARAMETERS:
   real, intent(out)   :: vsettle                 ! Layer fall speed [m s-1]

   character(len=*), parameter :: myname = 'Chem_CalcVsettle2G'

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

!-------------------------------------------------------------------------
!  Begin

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
   end subroutine Chem_CalcVsettle2G

!==================================================================================
   
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

   subroutine ParticleSwelling(i1, i2, j1, j2, km, rh, radius_dry, rhop_dry, radius, rhop, flag)
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
!----------------------------------------------------------------------

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
!----------------------------------------------------------------------

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
!----------------------------------------------------------------------
   
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
!----------------------------------------------------------------------
   
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
! !IROUTINE: Chem_Settling2G

   subroutine Chem_Settling2Gorig (km, flag, int_qa, grav, delp, &
                               radiusInp, rhopInp, cdt, tmpu, rhoa, &
                               rh, hghte, fluxout, vsettleOut, correctionMaring, rc)

! !USES:
   implicit none

! !INPUT PARAMETERS:
   integer,    intent(in)    :: km
   integer,    intent(in)    :: flag           ! flag to control particle swelling (see note)
   real,       intent(in)    :: grav           !mapl_grav
   real,       intent(in)    :: cdt
   real, dimension(:), intent(in)  :: radiusInp, rhopInp
   real, dimension(:,:,:,:), intent(inout) :: int_qa(:,:,:,:)
!   real, dimension(:,:,:), intent(in) :: tmpu, rhoa, rh, hghte, delp
   real, pointer, dimension(:,:,:), intent(in) :: tmpu, rhoa, rh, hghte, delp
   logical, optional, intent(in)    :: correctionMaring


! !OUTPUT PARAMETERS:

   real, pointer, dimension(:,:,:), intent(inout)  :: fluxout
   real, dimension(:,:,:,:), optional, intent(out)  :: vsettleOut

   integer,   intent(out)   :: rc

! !DESCRIPTION: Gravitational settling of aerosol between vertical
!               layers.  Assumes input radius in [m] and density (rhop) 
!               in [kg m-3]. If flag is set, use the Fitzgerald 1975 (flag = 1)
!               or Gerber 1985 (flag = 2) parameterization to update the 
!               particle radius for the calculation (local variables radius
!               and rhop).
!
! !REVISION HISTORY:
!
! 11Feb2020 E.Sherman - First attempt at refactor
!

!  !Local
   integer, parameter     :: DP=kind(1.0d0)
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

   integer         :: i1=1, i2, j1=1, j2, nbins
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
!   real, allocatable, dimension(:,:)   :: hsurf
   real, pointer, dimension(:,:)   :: hsurf

   real, dimension(:,:,:), allocatable  :: vsettle   ! fall speed [m s-1]
   real(kind=DP), dimension(:,:,:), allocatable   :: dzd, vsd, qa, qa_temp
   real(kind=DP), dimension(:,:), allocatable  :: cmass_before, cmass_after, qdel, &
        d_p, dpm1, qsrc


!EOP
!-------------------------------------------------------------------------

!  Get dimensions
!  ---------------
   nbins = size(radiusInp)
   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)
   ijl  = ( i2 - i1 + 1 ) * ( j2 - j1 + 1 )

!print*,'nbins = ',nbins
!print*,'i2 = ',i2
!print*,'j2 = ',j2

   gravDP = grav

   hsurf => hghte(i1:i2,j1:j2,km)

!  Allocate arrays
!  ---------------
   allocate(dz, mold=rhoa); 
   allocate(dzd(i2,j2,km), vsd(i2,j2,km), qa(i2,j2,km), vsettle(i2,j2,km), qa_temp(i2,j2,km))
   allocate(cmass_before(i2,j2), cmass_after(i2,j2), qdel(i2,j2), d_p(i2,j2), &
            dpm1(i2,j2), qsrc(i2,j2))

!#if 0
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

!print*,'CHEM_SETTLING TEST A'

!  Loop over the number of dust bins
   do n = 1, nbins

!    qa = w_c%qa(nbeg+n-1)%data3d
    qa(:,:,:) = int_qa(:,:,:,n)

    radius = radiusInp(n)
    rhop = rhopInp(n)

!   Reset a (large) minimum time to cross a grid cell in settling
    minTime = cdt

    if(associated(fluxout)) fluxout(:,:,n) = 0.0
    cmass_before(:,:) = 0.d0
    cmass_after(:,:) = 0.d0

!   If radius le 0 then get out of loop
    if(radius .le. 0.) cycle
    do k = 1, km
       do j = j1, j2
          do i = i1, i2
!            Find the column dry mass before sedimentation
             cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k)/gravDP * delp(i,j,k)

!            Adjust the particle size for relative humidity effects
             sat = max(rh(i,j,k),tiny(1.0)) ! to avoid zero FPE

!            Fitzgerald
             select case(flag)
             case(1)
               if (sat >= 0.80) then 
!             if(flag .eq. 1 .and. sat .ge. 0.80) then
!            parameterization blows up for RH > 0.995, so set that as max
!            rh needs to be scaled 0 - 1
                sat = min(0.995,sat)
!               Calculate the alpha and beta parameters for the wet particle
!               relative to amonium sulfate
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
!               radius is the radius of the wet particle
                radius = alpha * radiusInp(n)**beta
                rrat = (radiusInp(n)/radius)**3.
                rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
               end if
              case(2)
!             elseif(flag .eq. 2) then   ! Gerber
                sat = min(0.995,sat)
                rcm = radiusInp(n)*100.
                radius = 0.01 * (c1*rcm**c2 / (c3*rcm**c4-alog10(sat)) &
                                 + rcm**3.)**(1./3.)
                rrat = (radiusInp(n)/radius)**3.
                rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
              case(3)
!             elseif(flag .eq. 3) then
!               Gerber parameterization for Ammonium Sulfate
                sat = min(0.995,sat)
                rcm = radiusInp(n)*100.
                radius = 0.01 * (SU_c1*rcm**SU_c2 / (SU_c3*rcm**SU_c4-alog10(sat)) &
                                 + rcm**3.)**(1./3.)
                rrat = (radiusInp(n)/radius)**3.
                rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
              case(4)
!             elseif(flag .eq. 4) then
!               Petters and Kreidenweis (ACP2007) parameterization
                sat = min(0.99,sat)
                radius = (radiusInp(n)**3 * (1+1.19*sat/(1-sat)))**(1./3.)
                rrat = (radiusInp(n)/radius)**3
                rhop = rrat*rhopInp(n) + (1.-rrat)*rhow
!             endif
              end select

!            Calculate the settling velocity
             call Chem_CalcVsettle2Gorig(radius, rhop, rhoa(i,j,k), tmpu(i,j,k), &
                                     grav, diff_coef, vsettle(i,j,k))
          end do !do i
       end do !do j
    end do !do k

!print*,'CHEM_SETTLING TEST2'

    if(present(correctionMaring)) then
       if ((correctionMaring) .and. (radiusInp(n) .le. (0.5*diameterMaring))) then
          vsettle = max(1.0e-9, vsettle - v_upwardMaring)
       endif
    endif

    vsd = vsettle

    if(present(vsettleOut)) then
!       vsettleOut(n)%data3d = vsettle
       vsettleOut(:,:,:,n) = vsettle
    endif

!   Determine global min time to cross grid cell
    qmin = minval(dz/vsettle)

!    call pmaxmin ( 'Chem_Settling: dt', dz(i1:i2,j1:j2,1:km)/vsettle(i1:i2,j1:j2,1:km), &
!                                        qmin, qmax, ijl, km, 0. )
 
    minTime = min(minTime,qmin)

!   Now, how many iterations do we need to do?
    if ( minTime < 0 ) then
         nSubSteps = 0
!REVISIT THIS!!!
!----------------------
!         call mpout_log(myname,'no Settling because minTime = ', minTime ) 
!---------------------
    else if(minTime .ge. cdt) then
       nSubSteps = 1
       dt_settle = cdt
    else
       nSubSteps = cdt/minTime+1
       dt_settle = cdt/nSubSteps
    endif

!print*,'CHEM_SETTLING TEST3'

qa_temp = qa

!   Loop over sub-timestep
    do iit = 1, nSubSteps

!      Try a simple forward Euler scheme
       qdel = qa(i1:i2,j1:j2,1)*dt_settle*vsd(i1:i2,j1:j2,1)/dzd(i1:i2,j1:j2,1)
       qa(i1:i2,j1:j2,1) = qa(i1:i2,j1:j2,1) - qdel

       do k = 2, km
          d_p   = delp(i1:i2,j1:j2,k)
          dpm1 = delp(i1:i2,j1:j2,k-1)
          qsrc = qdel * dpm1 / d_p
          qdel = qa(i1:i2,j1:j2,k)*dt_settle*vsd(i1:i2,j1:j2,k)/dzd(i1:i2,j1:j2,k)
          qa(i1:i2,j1:j2,k) = qa(i1:i2,j1:j2,k) - qdel + qsrc
       enddo

!commented out in original Chem_Settling
!!     An alternative accumulator approach to computing the outgoing flux
!!     if( associated(fluxout(n)%data2d) ) then
!!        fluxout(n)%data2d = fluxout(n)%data2d + qdel * pdog/grav / dt_settle
!!     endif

    end do  ! iit

!if (sum(qa_temp) < sum(qa)) then
!   print*,'qa_temp = ',sum(qa_temp), ' : qa = ', sum(qa)
!end if

!   Find the column dry mass after sedimentation and thus the loss flux
!    do k = 1, km
!       do j = j1, j2
!          do i = i1, i2
!             cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k)/ gravDP * delp(i,j,k)
!          enddo
!       enddo
!    enddo

    cmass_after = sum(qa / gravDP * delp,3)

    if( associated(fluxout) ) then
       fluxout(:,:,n) = (cmass_before - cmass_after)/cdt
    endif

!    w_c%qa(nbeg+n-1)%data3d = qa
    int_qa(:,:,:,n) = qa
   end do   ! n

   rc=0

!#endif
   end subroutine Chem_Settling2Gorig

!=========================================================================================

!BOP
!
! !IROUTINE:  Chem_CalcVsettle2G - Calculate the aerosol settling velocity
!
! !INTERFACE:
!

   subroutine Chem_CalcVsettle2Gorig ( radius, rhop, rhoa, tmpu, grav, &
                                 diff_coef, vsettle )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   real, intent(in)    :: radius              ! Particle radius [m]
   real, intent(in)    :: rhop                ! Particle density [kg m-3]
   real, intent(in)    :: rhoa                ! Layer air density [kg m-3]
   real, intent(in)    :: tmpu                ! Layer temperature [K]
   real, intent(in)    :: grav                ! MAPL_GRAV
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

! !Local Variables
   integer, parameter :: DP=kind(1.0d0)
   real(kind=DP) :: rmu                   ! Dynamic viscosity [kg m-1 s-1]
   real(kind=DP) :: vt                    ! Thermal velocity of air molecule [m s-1]
   real(kind=DP) :: rmfp                  ! Air molecule mean free path [m]
   real(kind=DP) :: bpm                   ! Cunningham slip correction factor
   real(kind=DP) :: rkn                   ! Knudsen number
   real(kind=DP) :: re, x, y              ! reynold's number and parameters
   real, parameter :: kb = 1.3807e-23 ! Boltzmann constant [kg m2 s-1 K-1 mol-1]
   real, parameter :: m_air = 4.8096e-26 ! Mass of <avg> air molecule [kg]
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

!============================================================================
!BOP
!
! !IROUTINE: DryDepositionGOCART - Calculate aerosol dry deposition for lowest layer
!
! !INTERFACE:
!

   subroutine DryDeposition ( km, tmpu, rhoa, hghte, oro, ustar, pblh, shflux, & 
                              von_karman, cpd, grav, z0h, drydepf, rc, &
                              radius, rhop, u10m, v10m, fraclake, gwettop )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km
   real, pointer, dimension(:,:,:) :: tmpu      ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa      ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: hghte     ! top of layer geopotential height [m]
   real, pointer, dimension(:,:)   :: oro       ! orography flag
   real, pointer, dimension(:,:)   :: ustar     ! friction speed [m s-1]
   real, pointer, dimension(:,:)   :: pblh      ! PBL height [m]
   real, pointer, dimension(:,:)   :: shflux    ! sfc. sens. heat flux [W m-2]
   real, intent(in)                :: von_karman 
   real, intent(in)                :: cpd       
   real, intent(in)                :: grav      ! gravity
   real, pointer, dimension(:,:)   :: z0h       ! rough height, sens. heat [m]

! !OUTPUT PARAMETERS:
   real, intent(inout)        :: drydepf(:,:)     ! Deposition frequency [s-1]
   integer, intent(out)          :: rc          ! Error return code:

! !OPTIONAL PARAMETERS:
!  If these parameters are provided we compute a "resuspension" term as
!  if the particles are lifted like dust
   real, optional                            :: radius    ! particle radius [m]
   real, optional                            :: rhop      ! particle density [kg m-3]
   real, pointer, dimension(:,:), optional   :: u10m      ! 10-m u-wind component [m s-1]
   real, pointer, dimension(:,:), optional   :: v10m      ! 10-m v-wind component [m s-1]
   real, pointer, dimension(:,:), optional   :: fraclake  ! fraction covered by water
   real, pointer, dimension(:,:), optional   :: gwettop   ! fraction soil moisture



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
            obk(i2,j2))

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
!print*,'DU sum(obk) = ',sum(obk)
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

!print*,'DU sum(vdep) = ',sum(vdep)
!print*,'DU sum(dz) = ',sum(dz)

!     Save the dry deposition frequency for the chemical removal terms
!     in units of s-1
      drydepf(i,j) = max(0.,vdep(i,j) / dz(i,j))

    end do  ! i
   end do   ! j

!print*,'DU sum(drydepf) = ', sum(drydepf)


   rc = 0

   end subroutine DryDeposition

!====================================================================================
! !IROUTINE: ObukhovLength - Calculate the Obukhov length scale stability parameter
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
   real, intent(in) :: von_karman
   real, intent(in) :: cpd
   real, intent(in) :: grav
   real, dimension(i1:i2,j1:j2)  :: t         ! temperature [K]
   real, dimension(i1:i2,j1:j2)  :: rhoa      ! air density [kg m-3]
   real, pointer, dimension(:,:) :: ustar     ! friction speed [m s-1]
   real, pointer, dimension(:,:) :: shflux    ! sfc. sens. heat flux [W m-2]

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

! !IROUTINE: updatePointwiseEmissions
!#if 0
   subroutine WetRemovalGOCART2G ( km, n1, n2, bin_ind, cdt, aero_type, kin, grav, fwet, &
                                   aerosol, ple, tmpu, rhoa, pfllsan, pfilsan, &
                                   precc, precl, fluxout, rc )

! !USES:
  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km, n1, n2, bin_ind
   real, intent(in)    :: cdt
   character(len=*)    :: aero_type
   logical, intent(inout)  :: KIN ! true for aerosol
   real, intent(in)    :: grav
   real, intent(in)    :: fwet
   real, dimension(:,:,:), intent(inout) :: aerosol  ! internal state aerosol
   real, pointer, dimension(:,:,:)     :: ple, tmpu, rhoa
   real, pointer, dimension(:,:,:)     :: pfllsan, pfilsan
   real, pointer, dimension(:,:)       :: precc, precl
   real, pointer, dimension(:,:,:)       :: fluxout  ! tracer loss flux [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, intent(out)             :: rc          ! Error return code:
!   real, pointer, dimension(:,:)       :: fluxout  ! tracer loss flux [kg m-2 s-1]

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
!   real :: pdog(i1:i2,j1:j2,km)      ! air mass factor dp/g [kg m-2]
!   real :: delz(i1:i2,j1:j2,km)      ! box height  dp/g/rhoa [m]
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
!   real :: c_h2o(i1:i2,j1:j2,km), cldliq(i1:i2,j1:j2,km), cldice(i1:i2,j1:j2,km)
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

   rc=0

!EOP
!-----------------------------------------------------------------------------
!  Begin...

   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

!  Allocate arrays
   allocate(c_h2o(i2,j2,km), cldliq(i2,j2,km), cldice(i2,j2,km), pdog(i2,j2,km), &
            delz(i2,j2,km), dpfli(i2,j2,km))

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

   allocate(fd(km,nbins),stat=ios)
   allocate(dc(nbins),stat=ios)

!   if( associated(fluxout%data2d) ) fluxout%data2d(i1:i2,j1:j2) = 0.0
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
!        if (MAPL_AM_I_ROOT()) print *, 'stop in WetRemoval, need Kstar298 and H298_R'
        stop
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
     if(pac .le. 0.) goto 100
     pls = precl(i,j)
     pcv = precc(i,j)

!    Initialize the precipitation fields
     qls(:)  = 0.
     qcv(:)  = 0.
     fd(:,:) = 0.

!    Find the highest model layer experiencing rainout.  Assumes no
!    scavenging if T < 258 K
     LH = 0
     do k = 1, km
      if(dpfli(i,j,k) .gt. 0. ) then
       LH = k
       goto 15
      endif
     end do
 15  continue
     if(LH .lt. 1) goto 100

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
      if (qls(k) .gt. 0.) then
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
!        DC(n) = qa(n1+n-1)%data3d(i,j,k) * F * effRemoval *(1.-exp(-BT))
        DC(n) = aerosol(i,j,k) * F * effRemoval *(1.-exp(-BT))
        if (DC(n).lt.0.) DC(n) = 0.
!        qa(n1+n-1)%data3d(i,j,k) = qa(n1+n-1)%data3d(i,j,k)-DC(n)
        aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
!        if (qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) qa(n1+n-1)%data3d(i,j,k) = 1.0E-32
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
          goto 333
         end if
        end do

 333    continue
        F = F0_ls / (1. + F0_ls*B0_ls*XL_ls/(Qmx*cdt/Td_ls))
     ! if (MAPL_AM_I_ROOT()) then
     !    print *, 'hbianwdep WASHFmax =',F
     ! endif
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
!            DC(n) = qa(n1+n-1)%data3d(i,j,k) * F * (1.-exp(-BT))
            DC(n) = aerosol(i,j,k) * F * (1.-exp(-BT))
         else
!            DC(n) = qa(n1+n-1)%data3d(i,j,k) * F * WASHFRAC
            DC(n) = aerosol(i,j,k) * F * WASHFRAC
         endif
         if (DC(n).lt.0.) DC(n) = 0.
!         qa(n1+n-1)%data3d(i,j,k) = qa(n1+n-1)%data3d(i,j,k)-DC(n)
         aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
!         if (qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) &
         if (aerosol(i,j,k) .lt. 1.0E-32) &
          aerosol(i,j,k) = 1.0E-32
!         if( associated(fluxout%data2d) ) then
         if( associated(fluxout)) then
!          fluxout%data2d(i,j) = fluxout%data2d(i,j)+DC(n)*pdog(i,j,k)/cdt
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
!        effRemoval = qa(n1+n-1)%fwet
        effRemoval = fwet
!        DC(n) = qa(n1+n-1)%data3d(i,j,k) * F * effRemoval * (1.-exp(-BT))
        DC(n) = aerosol(i,j,k) * F * effRemoval * (1.-exp(-BT))
        if (DC(n).lt.0.) DC(n) = 0.
!        qa(n1+n-1)%data3d(i,j,k) = qa(n1+n-1)%data3d(i,j,k)-DC(n)
        aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
!        if (qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) qa(n1+n-1)%data3d(i,j,k) = 1.0E-32
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
          goto 444
         end if
        end do

 444    continue
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
!         DC(n) = qa(n1+n-1)%data3d(i,j,k) * F * (1.-exp(-BT))
         DC(n) = aerosol(i,j,k) * F * (1.-exp(-BT))
         if (DC(n).lt.0.) DC(n) = 0.
!         qa(n1+n-1)%data3d(i,j,k) = qa(n1+n-1)%data3d(i,j,k)-DC(n)
         aerosol(i,j,k) = aerosol(i,j,k)-DC(n)
!         if (qa(n1+n-1)%data3d(i,j,k) .lt. 1.0E-32) &
         if (aerosol(i,j,k) .lt. 1.0E-32) &
!          qa(n1+n-1)%data3d(i,j,k) = 1.0E-32
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
!           qa(n1+n-1)%data3d(i,j,k) = qa(n1+n-1)%data3d(i,j,k) + DC(n)
           aerosol(i,j,k) = aerosol(i,j,k) + DC(n)
!           qa(n1+n-1)%data3d(i,j,k) = max(qa(n1+n-1)%data3d(i,j,k),1.e-32)
           aerosol(i,j,k) = max(aerosol(i,j,k),1.e-32)
!          Adjust the flux out of the bottom of the layer
           Fd(k,n)  = Fd(k,n) - DC(n)*pdog(i,j,k)
          end do

        endif
       endif                                   ! if -moistq < 0
      endif
     end do  ! k


     do n = 1, nbins
!      if( associated(fluxout%data2d) ) then
      if( associated(fluxout)) then
!       fluxout%data2d(i,j) = fluxout%data2d(i,j)+Fd(km,n)/cdt
       fluxout(i,j,bin_ind) = fluxout(i,j,bin_ind)+Fd(km,n)/cdt
      endif
     end do

 100 continue
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
   real, pointer, dimension(:,:)     :: emissions_surface
   real, dimension(:,:,:,:), intent(inout) :: emissions
   real, dimension(:,:,:), intent(in) :: emissions_point

   real, dimension(:), intent(in)   :: sfrac
   integer, intent(in)  :: nPts
   integer, intent(in) :: km
   real, intent(in)  :: cdt
   real, intent(in) :: grav
   integer, intent(in) :: nbins
   real, pointer, dimension(:,:,:) :: delp
   real, pointer, dimension(:,:,:,:)  :: aero


! !OUTPUT PARAMETERS:
   integer, intent(out)             :: rc          ! Error return code:

! !DESCRIPTION: Updates internal state variables
!
! !REVISION HISTORY:
!
!  15May2020 - Sherman
!
! !Local Variables
   integer :: n, kmin


!EOP
!--------------------------------------------------------------------------------
!   Begin...

    rc = 0

    do n = 1, nbins
       emissions(:,:,km,n) = emissions_surface * sfrac(n)
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
!

   subroutine Aero_Compute_Diags ( mie_table, km, nbins, rlow, rup, channels, aerosol, grav, tmpu, rhoa, &
                                 rh, u, v, delp, sfcmass, colmass, mass, exttau, scatau,     &
                                 sfcmass25, colmass25, mass25, exttau25, scatau25, &
                                 fluxu, fluxv, conc, extcoef, scacoef,    &
                                 exttaufm, scataufm, angstrom, aerindx, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   type(Chem_Mie),  intent(in) :: mie_table        ! mie table
   integer, intent(in) :: km, nbins
   real, dimension(:), intent(in)   :: rlow   ! bin radii - low bounds
   real, dimension(:), intent(in)   :: rup    ! bin radii - upper bounds
   real, dimension(:), intent(in) :: channels
   real, pointer, dimension(:,:,:,:) :: aerosol     ! 
   real :: grav
   real, pointer, dimension(:,:,:) :: tmpu     ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa     ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: delp     ! 
   real, pointer, dimension(:,:,:) :: rh     ! 
   real, pointer, dimension(:,:,:) :: u        ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:) :: v        ! north-south wind [m s-1]


! !OUTPUT PARAMETERS:
!  Total mass
   real, pointer, dimension(:,:), intent(inout)  :: sfcmass   ! sfc mass concentration kg/m3
   real, pointer, dimension(:,:), intent(inout)  :: colmass   ! col mass density kg/m2
   real, pointer, dimension(:,:,:), intent(inout)  :: mass      ! 3d mass mixing ratio kg/kg
!  Total optical properties
   real, pointer, dimension(:,:), intent(inout)  :: exttau    ! ext. AOT at 550 nm
   real, pointer, dimension(:,:), intent(inout)  :: scatau    ! sct. AOT at 550 nm
   real, pointer, dimension(:,:), intent(inout)  :: sfcmass25 ! sfc mass concentration kg/m3 (pm2.5)
   real, pointer, dimension(:,:), intent(inout)  :: colmass25 ! col mass density kg/m2 (pm2.5)
   real, pointer, dimension(:,:,:), intent(inout)  :: mass25    ! 3d mass mixing ratio kg/kg (pm2.5)
   real, pointer, dimension(:,:), intent(inout)  :: exttau25  ! ext. AOT at 550 nm (pm2.5)
   real, pointer, dimension(:,:), intent(inout)  :: scatau25  ! sct. AOT at 550 nm (pm2.5)
   real, pointer, dimension(:,:),  intent(inout)  :: aerindx   ! TOMS UV AI
   real, pointer, dimension(:,:), intent(inout)  :: fluxu     ! Column mass flux in x direction
   real, pointer, dimension(:,:), intent(inout)  :: fluxv     ! Column mass flux in y direction
   real, pointer, dimension(:,:,:), intent(inout)  :: conc      ! 3d mass concentration, kg/m3
   real, pointer, dimension(:,:,:), intent(inout)  :: extcoef   ! 3d ext. coefficient, 1/m
   real, pointer, dimension(:,:,:), intent(inout)  :: scacoef   ! 3d scat.coefficient, 1/m
   real, pointer, dimension(:,:), intent(inout)  :: exttaufm  ! fine mode (sub-micron) ext. AOT at 550 nm
   real, pointer, dimension(:,:), intent(inout)  :: scataufm  ! fine mode (sub-micron) sct. AOT at 550 nm
   real, pointer, dimension(:,:), intent(inout)  :: angstrom  ! 470-870 nm Angstrom parameter
   integer, intent(out)             :: rc        ! Error return code:
                                                 !  0 - all is well
                                                 !  1 - 

! !DESCRIPTION: Calculates some simple 2d diagnostics from the dust fields
!
! !REVISION HISTORY:
!
!  16APR2004, Colarco
!  11MAR2010, Nowottnick  

! !Local Variables
   character(len=*), parameter :: myname = 'Aero_Compute_Diags'
   integer :: i, j, k, n, ios, nch, idx
   integer :: i1 =1, i2, j1=1, j2
   real :: ilam550, ilam470, ilam870
   real :: tau, ssa
   real :: fPMfm(nbins)  ! fraction of bin with particles diameter < 1.0 um
   real :: fPM25(nbins)  ! fraction of bin with particles diameter < 2.5 um
   character(len=255) :: qname
   logical :: do_angstrom
   real, dimension(:,:), allocatable :: tau470, tau870

!EOP
!-------------------------------------------------------------------------
!  Begin...
 
   rc = 0

!  Initialize local variables
!  --------------------------
   nch = size(channels)
   i2 = size(rhoa,1)
   j2 = size(rhoa,2)

!  Get the wavelength indices
!  --------------------------
!  Must provide ilam550 for AOT calculation
   ilam550 = 1.
   ilam470 = 0.
   ilam870 = 0.
   if(nch .gt. 1) then
    do i = 1, nch
     if ( channels(i) .ge. 5.49e-7 .and. &
          channels(i) .le. 5.51e-7) ilam550 = i
     if ( channels(i) .ge. 4.69e-7 .and. &
          channels(i) .le. 4.71e-7) ilam470 = i
     if ( channels(i) .ge. 8.69e-7 .and. &
          channels(i) .le. 8.71e-7) ilam870 = i
    enddo
   endif

   do_angstrom = .false.
!  If both 470 and 870 channels provided (and not the same) then
!  possibly will do Angstrom parameter calculation
   if(ilam470 .ne. 0. .and. &
      ilam870 .ne. 0. .and. &
      ilam470 .ne. ilam870) do_angstrom = .true.

!   if( associated(angstrom) .and. do_angstrom ) then
      allocate(tau470(i1:i2,j1:j2), tau870(i1:i2,j1:j2))
!   end if

!print*,'SS2G ilam550 = ', ilam550
!print*,'SS2G ilam470 = ', ilam470
!print*,'SS2G ilam870 = ', ilam870
!print*,'SS2G do_angstrom = ',do_angstrom

!  Compute the fine mode (sub-micron) and PM2.5 bin-wise fractions
!  ------------------------------------
   call Aero_Binwise_PM_Fractions(fPMfm, 0.50, rlow, rup, nbins)   ! 2*r < 1.0 um
   call Aero_Binwise_PM_Fractions(fPM25, 1.25, rlow, rup, nbins)   ! 2*r < 2.5 um

   if (associated(aerindx))  aerindx = 0.0  ! for now

!  Calculate the diagnostic variables if requested
!  -----------------------------------------------

!  Calculate the surface mass concentration
   if( associated(sfcmass) ) then
      sfcmass(i1:i2,j1:j2) = 0.
      do n = 1, nbins
         sfcmass(i1:i2,j1:j2) &
              =   sfcmass(i1:i2,j1:j2) &
              + aerosol(i1:i2,j1:j2,km,n)*rhoa(i1:i2,j1:j2,km)
      end do
   endif
   if( associated(sfcmass25) ) then
      sfcmass25(i1:i2,j1:j2) = 0.
      do n = 1, nbins
         sfcmass25(i1:i2,j1:j2) &
              =   sfcmass25(i1:i2,j1:j2) &
              + aerosol(i1:i2,j1:j2,km,n)*rhoa(i1:i2,j1:j2,km)*fPM25(n)
      end do
   endif

!  Calculate the aerosol column loading
   if( associated(colmass) ) then
      colmass(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        colmass(i1:i2,j1:j2) &
         =   colmass(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav
       end do
      end do
   endif
   if( associated(colmass25)) then
      colmass25(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        colmass25(i1:i2,j1:j2) &
         =   colmass25(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*fPM25(n)
       end do
      end do
   endif

!  Calculate the total mass concentration
   if( associated(conc) ) then
      conc(i1:i2,j1:j2,1:km) = 0.
      do n = 1, nbins
       conc(i1:i2,j1:j2,1:km) &
         =   conc(i1:i2,j1:j2,1:km) &
           + aerosol(i1:i2,j1:j2,1:km,n)*rhoa(i1:i2,j1:j2,1:km)
      end do
   endif

!  Calculate the total mass mixing ratio
   if( associated(mass) ) then
      mass(i1:i2,j1:j2,1:km) = 0.
      do n = 1, nbins
       mass(i1:i2,j1:j2,1:km) &
         =   mass(i1:i2,j1:j2,1:km) &
           + aerosol(i1:i2,j1:j2,1:km,n)
      end do
   endif
   if( associated(mass25) ) then
      mass25(i1:i2,j1:j2,1:km) = 0.
      do n = 1, nbins
       mass25(i1:i2,j1:j2,1:km) &
         =   mass25(i1:i2,j1:j2,1:km) &
           + aerosol(i1:i2,j1:j2,1:km,n)*fPM25(n)
      end do
   endif

!  Calculate the column mass flux in x direction
   if( associated(fluxu) ) then
      fluxu(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        fluxu(i1:i2,j1:j2) &
         =   fluxu(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
       end do
      end do
   endif

!  Calculate the column mass flux in y direction
   if( associated(fluxv) ) then
      fluxv(i1:i2,j1:j2) = 0.
      do n = 1, nbins
       do k = 1, km
        fluxv(i1:i2,j1:j2) &
         =   fluxv(i1:i2,j1:j2) &
           + aerosol(i1:i2,j1:j2,k,n)*delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
       end do
      end do
   endif

!  Calculate the extinction and/or scattering AOD
   if( associated(exttau) .or. associated(scatau) ) then

      if( associated(exttau)) exttau(i1:i2,j1:j2) = 0.
      if( associated(scatau)) scatau(i1:i2,j1:j2) = 0.

      if( associated(exttau25)) exttau25(i1:i2,j1:j2) = 0.
      if( associated(scatau)) scatau(i1:i2,j1:j2) = 0.

      if( associated(exttau25)) exttau25(i1:i2,j1:j2) = 0.
      if( associated(scatau25)) scatau25(i1:i2,j1:j2) = 0.

      if( associated(exttaufm)) exttaufm(i1:i2,j1:j2) = 0.
      if( associated(scataufm)) scataufm(i1:i2,j1:j2) = 0.

      if( associated(extcoef)) extcoef(i1:i2,j1:j2,1:km) = 0.
      if( associated(scacoef)) scacoef(i1:i2,j1:j2,1:km) = 0.

      do n = 1, nbins

!      Select the name for species
!       qname = trim(w_c%reg%vname(w_c%reg%i_DU+n-1))
!       idx = Chem_MieQueryIdx(gcDU%mie_tables,qname,rc)
!       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2
          call Chem_MieQuery(mie_table, n, ilam550, &
              aerosol(i,j,k,n)*delp(i,j,k)/grav, &
              rh(i,j,k), tau=tau, ssa=ssa)

!         Calculate the total ext. and scat. coefficients
          if( associated(extcoef) ) then
              extcoef(i,j,k) = extcoef(i,j,k) + &
                                      tau * (grav * rhoa(i,j,k) / delp(i,j,k))
          endif
          if( associated(scacoef) ) then
              scacoef(i,j,k) = scacoef(i,j,k) + &
                                      ssa * tau * (grav * rhoa(i,j,k) / delp(i,j,k))
          endif

!         Integrate in the vertical
          if( associated(exttau) ) exttau(i,j) = exttau(i,j) + tau
          if( associated(exttaufm)) &
                         exttaufm(i,j) = exttaufm(i,j) + tau*fPMfm(n)
          if( associated(exttau25)) &
                         exttau25(i,j) = exttau25(i,j) + tau*fPM25(n)

          if( associated(scatau) ) scatau(i,j) = scatau(i,j) + tau*ssa
          if( associated(scataufm) ) &
                         scataufm(i,j) = scataufm(i,j) + tau*ssa*fPMfm(n)
          if( associated(scatau25) ) &
                         scatau25(i,j) = scatau25(i,j) + tau*ssa*fPM25(n)

         enddo
        enddo
       enddo

      enddo  ! nbins

   endif

!  Calculate the 470-870 Angstrom parameter
   if( associated(angstrom) .and. do_angstrom ) then

      angstrom(i1:i2,j1:j2) = 0.
!     Set tau to small number by default
      tau470(i1:i2,j1:j2) = tiny(1.0)
      tau870(i1:i2,j1:j2) = tiny(1.0)

      do n = 1, nbins

!      Select the name for species
!       qname = trim(w_c%reg%vname(w_c%reg%i_DU+n-1))
!       idx = Chem_MieQueryIdx(gcDU%mie_tables,qname,rc)
!       if(rc .ne. 0) call die(myname, 'cannot find proper Mie table index')

       do k = 1, km
        do j = j1, j2
         do i = i1, i2

          call Chem_MieQuery(mie_table, n, ilam470, &
              aerosol(i,j,k,n)*delp(i,j,k)/grav, &
              rh(i,j,k), tau=tau)
          tau470(i,j) = tau470(i,j) + tau

          call Chem_MieQuery(mie_table, idx, ilam870, &
              aerosol(i,j,k,n)*delp(i,j,k)/grav, &
              rh(i,j,k), tau=tau)
          tau870(i,j) = tau870(i,j) + tau


         enddo
        enddo
       enddo

      enddo  ! nbins

!print*,'SS2G sum(tau470) = ',sum(tau470)
!print*,'SS2G sum(tau870) = ',sum(tau870)

      angstrom(i1:i2,j1:j2) = &
        -log(tau470(i1:i2,j1:j2)/tau870(i1:i2,j1:j2)) / &
         log(470./870.)
!print*,'SS2G sum(angstrom) = ', sum(angstrom)

   endif

   rc = 0

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
   real, pointer, dimension(:,:)        :: lats
   real, pointer, dimension(:,:)        :: lons          
   real                                :: radToDeg
          
! !OUTPUT PARAMETERS:
   integer, optional, intent(out) :: rc
!EOP

! !Local Variables
    real                :: dummylon
    integer             :: i, j
!EOP
!-------------------------------------------------------------------------
!  Begin...

   rc = 0

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
  
   rc = 0

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
            gweibull(i,j)  = (c / wm(i,j))**3.41d0 * igamma(a,x)
         end if
      end do ! i
   end do ! j
   endif

   deallocate(wm)

   end subroutine weibullDistribution

!=====================================================================================

 DOUBLE PRECISION function igamma(A, X)
!----------------------------------------------------------------------- 
! incomplete Gamma function
!----------------------------------------------------------------------- 
 IMPLICIT NONE
 double precision, intent(in) ::        A
 DOUBLE PRECISION, INTENT(IN) ::      X
! LOCAL VARIABLE
 DOUBLE PRECISION :: XAM, GIN,  S, R, T0
 INTEGER K
        XAM=-X+A*LOG(X)
        IF (XAM.GT.700.0.OR.A.GT.170.0) THEN
           WRITE(*,*)'IGAMMA: a and/or x too large, X = ', X
           WRITE(*,*) 'A = ', A
           STOP

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
   subroutine SeasaltEmission ( rLow, rUp, method, u10m, v10m, ustar, &
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

   real, intent(in)                  :: rLow, rUp   ! Dry particle bin edge radii [um]
   real, intent(in)                  :: u10m(:,:)   ! 10-meter eastward wind [m s-1]
   real, intent(in)                  :: v10m(:,:)   ! 10-m northward wind [m s-1]
   real, target, intent(in)        :: ustar(:,:)  ! friction velocity [m s-1]
   integer, intent(in)               :: method      ! Algorithm to use

! !INOUTPUT PARAMETERS:

!   real, pointer, dimension(:,:) :: memissions      ! Mass Emissions Flux [kg m-2 s-1]
!   real, pointer, dimension(:,:) :: nemissions      ! Number Emissions Flux [# m-2 s-1]
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
   real, parameter    :: pi = 3.1415       ! ratio of circumference to diameter of circle
   integer, parameter :: nr = 10                    ! Number of (linear) sub-size bins

   character(len=*), parameter :: myname = 'SeasaltEmission'

!EOP
!-------------------------------------------------------------------------
!  Begin...

   rc = 0

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
      print *, 'SeasaltEmission missing algorithm method'
      rc = 1
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

   rc = 0

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
! !IROUTINE: weibullDistribution - Compute the wet radius of sea salt particle
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

   rc = 0

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
   integer, intent(in)   :: rhFlag      ! 1 (Fitzgerald, 1975)
                                        ! 2 (Gerber, 1985)
   real, dimension(:,:), intent(in)  :: rh    ! relative humidity [0-1]
   real, dimension(:,:), intent(in)  :: dz    ! surface layer height [m]
   real, dimension(:,:), intent(in)  :: ustar ! surface velocity scale [m s-1]
   real, dimension(:,:), intent(in)  :: airdens
   real, dimension(:,:), intent(in)  :: t  ! temperature [k]
   real, intent(in)  :: grav    ! gravity
   real, intent(in)  :: karman  ! Von Karman constant


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


!EOP
!------------------------------------------------------------------------------------
!  Begin..

   rc = 0
   fhoppel = 1.0
   allocate(vsettle, mold=rh)

   do j = 1, ubound(rh,2)
      do i = 1, ubound(rh,1)
         call wetRadius (radius, rhop, rh(i,j), rhFlag, &
                         radius_wet, rhop_wet, rc)
         if (rc /= 0) return
         call Chem_CalcVsettle2Gorig (radius_wet, rhop_wet, airdens(i,j), t(i,j), &
                                      GRAV, diff_coef, vsettle(i,j))
         fhoppel(i,j) = (10./dz(i,j)) ** (vsettle(i,j)/KARMAN/ustar(i,j))
      end do
   end do


   deallocate(vsettle)

   end subroutine hoppelCorrection





 end module GOCART2G_Process
