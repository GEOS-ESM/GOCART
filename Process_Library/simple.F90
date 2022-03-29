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
#include "DustAerosolDistributionKok.F90"
!===============================================================================
!BOP
!
! !IROUTINE: soilMoistureConvertVol2Grav - volumetric to gravimetric soil moisture
!
! !INTERFACE:
#include "soilMoistureConvertVol2Grav.F90"
!===============================================================================
!BOP
!
! !IROUTINE: moistureCorrectionFecan - Correction factor for Fecan soil moisture
!
! !INTERFACE:
#include "moistureCorrectionFecan.F90"
!===============================================================================
!BOP
!
! !IROUTINE: DustFluxV2HRatioMB95 - vertical-to-horizontal dust flux ratio (MB95)
!
! !INTERFACE:
#include "DustFluxV2HRatioMB95.F90"
!==================================================================================
!BOP
!
! !IROUTINE: DustEmissionFENGSHA - Compute dust emissions using NOAA/ARL FENGSHA model
!
! !INTERFACE:
#include "DustEmissionFENGSHA.F90"
!==================================================================================
!BOP
! !IROUTINE: DustEmissionGOCART2G

#include "DustEmissionGOCART2G.F90"
!==================================================================================
!BOP
! !IROUTINE: DustEmissionK14

#include "DustEmissionK14.F90"
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

#include "updatePointwiseEmissions.F90"
!==================================================================================
!BOP
! !IROUTINE: DistributePointEmissions

! !INTERFACE:

#include "DistributePointEmission.F90"

!==================================================================================
!BOP
! !IROUTINE: Chem_SettlingSimple

#include "Chem_SettlingSimple.F90"
!==================================================================================
!BOP
! !IROUTINE: Chem_Settling

#include "Chem_Settling.F90"

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

#include "Chem_Settling2Gorig.F90"
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

#include "Chem_SettlingSimpleOrig.F90"

!============================================================================
!BOP
!
! !IROUTINE: DryDeposition - Calculate aerosol dry deposition for lowest layer
!
! !INTERFACE:
!
#include "DryDeposition.F90"
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
#include "WetRemovalGOCART2G.F90"
!=============================================================================
!BOP

! !IROUTINE: UpdateAerosolState
#include "UpdateAerosolState.F90"
!==============================================================================

!BOP
!
! !IROUTINE:  Aero_Compute_Diags - Calculate aerosol diagnostics
!
! !INTERFACE:
#include "Aero_Compute_Diags.F90"!====================================================================

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
#include "deepLakesMask.F90"
!========================================================================================
!BOP
!
! !IROUTINE:  jeagleSSTcorrection - Apply SST correction following Jaegle et al. 2011
!
! !INTERFACE:
#include "jeagleSSTcorrection.F90"
!=====================================================================================
!BOP
!
! !IROUTINE: weibullDistribution - Apply a Weibull distribution to emissions wind speeds
!
! !INTERFACE:
#include "weibullDistribution.F90"
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
#include "SeasaltEmission.F90"

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
#include "wetRadius.F90"
!===============================================================================

!BOP
!
! !IROUTINE: hoppelCorrection
!
! !INTERFACE:
#include "hoppelCorrection.F90"
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

#include "CAEmission.F90"
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
#include "phobicTophilic.F90"

!============================================================================
!BOP
!
! !IROUTINE: NIheterogenousChemOpt
!
! !INTERFACE:
#include "NIheterogenousChem.F90"
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

#include "SulfateDistributeEmissions.F90"
!==================================================================================
!BOP
! !IROUTINE: DMSemission

#include "DMSemission.F90"

!==================================================================================
!BOP
! !IROUTINE: SUvolcanicEmissions

#include "SUvolcanicEmissions.F90"
!==================================================================================
!BOP
! !IROUTINE: SulfateUpdateOxidants

#include "SulfateUpdateOxidants.F90"
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

#include "SU_Wet_Removal.F90"
!==================================================================================
!BOP
! !IROUTINE: SU_Compute_Diags

#include "SU_Compute_Diags.F90"
!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver

#include "SulfateChemDriver.F90"
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

#include "get_HenrysLawCts.F90"
!==================================================================================
!BOP
! !IROUTINE: NIthermo

#include "NIthermo.F90"
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
          PRINT*, 'TSO4 : ', TSO4
          PRINT*, 'TNO3 : ', TNO3
          PRINT*, 'TNH4 : ', TNH4


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
            print *,'PHI < 1d-20 in  CUBIC (rpmares_mod.f)'
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
#include "Chem_UtilResVal.F90"
!==================================================================================

#include "Chem_UtilIdow.F90"
#include "Chem_UtilCdow.F90"

!==================================================================================

!BOP
!
! !ROUTINE:  Chem_BiomassDiurnal - Applies diurnal cycle to biomass emissions.
!
! !INTERFACE:
#include "Chem_BiomassDiurnal.F90"!==================================================================================

#include "ReadPointEmissions.F90"
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
!      print*,__FILE__,__LINE__, ' found label'

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
