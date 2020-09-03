!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  DryDepositionMod --- Aerosol Turbulent Deposition Module
!
! !INTERFACE:
!

   module  DryDepositionMod

! !USES:

   use Chem_Mod
   use Chem_ConstMod, only: grav, von_karman, cpd        ! Constants !

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  DryDepositionGOCART

!
! !DESCRIPTION:
!
!  This module implements various dry deposition schemes
!
! !REVISION HISTORY:
!
!  08Jan2010 - Colarco, first crack
!
!EOP
  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

!-------------------------------------------------------------------------
CONTAINS

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

   subroutine ObukhovLength ( i1, i2, j1, j2, &
                              t, rhoa, shflux, ustar, &
                              obk )

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2
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
   end subroutine ObukhovLength



!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: DryDepositionGOCART - Calculate aerosol dry deposition for lowest layer
!
! !INTERFACE:
!

   subroutine DryDepositionGOCART ( i1, i2, j1, j2, km, &
                                    tmpu, rhoa, hghte, oro, ustar, &
                                    pblh, shflux, z0h, drydepf, rc, &
                                    radius, rhop, u10m, v10m, fraclake, gwettop )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km
   real, pointer, dimension(:,:,:) :: tmpu      ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa      ! air density [kg m-3]
   real, pointer, dimension(:,:,:) :: hghte     ! top of layer geopotential height [m]
   real, pointer, dimension(:,:)   :: oro       ! orography flag
   real, pointer, dimension(:,:)   :: ustar     ! friction speed
   real, pointer, dimension(:,:)   :: pblh      ! PBL height [m]
   real, pointer, dimension(:,:)   :: shflux    ! sfc. sens. heat flux [W m-2]
   real, pointer, dimension(:,:)   :: z0h       ! rough height, sens. heat [m]

! !OUTPUT PARAMETERS:
   real             :: drydepf(i1:i2,j1:j2)     ! Deposition frequency [s-1]
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 
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
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   character(len=*), parameter :: myname = 'DryDepositionGOCART'
   integer :: i, j
   real, parameter :: rhow = 1000.      ! density of water [kg m-3]
   real, parameter :: coll_size = 0.002 ! collector size [m]
   real :: dz(i1:i2,j1:j2)              ! lowest layer thickness
   real :: rmu(i1:i2,j1:j2)             ! dynamic viscosity [kg m-1 s-1]
   real :: Ra(i1:i2,j1:j2)              ! aerodynamic resistance
   real :: Rs(i1:i2,j1:j2)              ! surface resistance
   real :: vdep(i1:i2,j1:j2)            ! Deposition speed [m s-1]
   real :: obk(i1:i2,j1:j2)             ! Obukhov Length [m]

   real*8 Rttl        ! total surface resistance

   real*8 R2, w10m, u_thresh0
   real*8 vds, vdsmax, czh
   real*8 :: frac, cz, psi_h, eps, logmfrac, z0h_min, z0h_
   real*8 :: one = 1.0, zero = 0.0

!  Initialize local variables
!  --------------------------

!  Calculate the viscosity and thickness of the surface level
   dz = hghte(:,:,km-1) - hghte(:,:,km)
   rmu = 1.8325e-5*(416.16/(tmpu(i1:i2,j1:j2,km)+120.)) &
                  *(tmpu(i1:i2,j1:j2,km)/296.16)**1.5

   z0h_min = 100. * tiny(1.0)  ! because sometimes we may get z0h=0.

!  =========================================================================
!  Calculate the Obukhov length scale
   call ObukhovLength ( i1, i2, j1, j2, &
                        tmpu(:,:,km), rhoa(:,:,km), shflux, ustar, &
                        obk )


!  =========================================================================
!  Aerodynamic Resistance
!  psi_h and Ra are equations 2, 4-5 of Walcek et al. 1986 Atmospheric Environment
!
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

!  =======================================================================
!  Surface Resistance term for aerosols
!  Rs formulation from eqn. 15 - 18 of Walcek et al. 1986 Atmospheric Environment

!  Loop over space
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

   rc = 0

   end subroutine DryDepositionGOCART

   end module DryDepositionMod
