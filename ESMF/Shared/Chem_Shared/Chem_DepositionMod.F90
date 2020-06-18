!-------------------------------------------------------------------------
!         NASA/GSFC, Data Assimilation Office, Code 910.3, GEOS/DAS      !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE:  Chem_DepositionMod --- Aerosol Turbulent Deposition Module
!
! !INTERFACE:
!

   module  Chem_DepositionMod

! !USES:

   use Chem_Mod
   use Chem_ConstMod, only: grav, von_karman, cpd        ! Constants !
!!!   use Chem_UtilMod

   implicit none

! !PUBLIC TYPES:
!
   PRIVATE

!
! !PUBLIC MEMBER FUNCTIONS:
!

   PUBLIC  Chem_Deposition

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
  real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

!-------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 900.3     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: Chem_Deposition - Calculate aerosol dry deposition for lowest layer
!
! !INTERFACE:
!

   subroutine Chem_Deposition ( i1, i2, j1, j2, km, nbeg, nend, nbins, cdt, w_c, &
                                radius, rhop, &
                                tmpu, rhoa, hsurf, hghte, oro, ustar, &
                                u10m, v10m, fraclake, gwettop, pblh, shflux, &
                                z0h, fluxout, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: i1, i2, j1, j2, km, nbins, nbeg, nend
   real, intent(in)    :: cdt
   real, pointer, dimension(:,:,:) :: tmpu      ! temperature [K]
   real, pointer, dimension(:,:,:) :: rhoa      ! air density [kg m-3]
   real, pointer, dimension(:,:)   :: hsurf     ! surface geopotential height [m]
   real, pointer, dimension(:,:,:) :: hghte     ! top of layer geopotential height [m]
   real, pointer, dimension(:,:)   :: oro       ! orography flag
   real, pointer, dimension(:,:)   :: fraclake  ! fraction covered by water
   real, pointer, dimension(:,:)   :: gwettop   ! fraction soil moisture
   real, pointer, dimension(:,:)   :: ustar     ! friction speed
   real, pointer, dimension(:,:)   :: u10m      ! 10-m u-wind component [m s-1]
   real, pointer, dimension(:,:)   :: v10m      ! 10-m v-wind component [m s-1]
   real, pointer, dimension(:,:)   :: pblh      ! PBL height [m]
   real, pointer, dimension(:,:)   :: shflux    ! sfc. sens. heat flux [W m-2]
   real, pointer, dimension(:,:)   :: z0h       ! rough height, sens. heat [m]
   real, pointer, dimension(:)     :: radius    ! particle radius [m]
   real, pointer, dimension(:)     :: rhop      ! particle density [kg m-3]

! !OUTPUT PARAMETERS:
   type(Chem_Bundle), intent(inout) :: w_c
   type(Chem_Array), pointer :: fluxout(:)      ! Mass lost by deposition
                                                      ! to surface, kg/m2/s
   integer, intent(out)             :: rc          ! Error return code:
                                                   !  0 - all is well
                                                   !  1 - 

! !DESCRIPTION: Calculates the deposition velocity and removal for the
!               lowest model layer
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!
!EOP
!-------------------------------------------------------------------------

! !Local Variables
   character(len=*), parameter :: myname = 'Chem_Deposition'
   character(len=255) :: NAME
   integer :: i, j, k, n, im, dk
   real, parameter :: rhow = 1000.   ! density of water [kg m-3]
   real, parameter :: coll_size = 0.002 ! collector size [m]
   real :: pm(i1:i2,j1:j2)           ! pressure [Pa]
   real :: dz(i1:i2,j1:j2)           ! lowest layer thickness
   real :: rmu(i1:i2,j1:j2)          ! dynamic viscosity [kg m-1 s-1]
   real :: Ra(i1:i2,j1:j2)         ! aerodynamic resistance
   real :: Rs(i1:i2,j1:j2)          ! surface resistance
   real :: vdep(i1:i2,j1:j2)        ! Deposition speed [m s-1]
   real :: drydepf(i1:i2,j1:j2)     ! Deposition frequency [s-1]
   real :: outflux(i1:i2,j1:j2)  

   real qmin, qmax

   real*8 Rttl        ! total surface resistance
   real*8 dc
   real*8 diff_coef, vsettle
   real*8 Sc          ! Schmidt number
   real*8 Eb          ! Brownian diffusion collection efficiency
   real*8 St          ! Stokes number
   real*8 Ein         ! Interception collection efficiency
   real*8 Eim         ! Impaction collection efficiency
   real*8 alpha, gamma

   real*8 R1, R2, w10m, u_thresh0
   real*8 obk, vds, vdsmax, czh, factor
   real*8 frac, cz, psi_h, eps, logmfrac, z0h_min, z0h_, r8_cdt
   real*8 :: one = 1.0, zero = 0.0

   r8_cdt = cdt

!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1

!  Initialize local variables
!  --------------------------
   im = w_c%grid%im

!  Calculate the pressure, air density, viscosity, and thickness of the
!  surface level
   pm = 0.5*w_c%delp(i1:i2,j1:j2,1)
   do k = 2, km
    pm = pm + 0.5*(w_c%delp(i1:i2,j1:j2,k)+w_c%delp(i1:i2,j1:j2,k-1))
   end do
   dz = hghte(:,:,km+dk) - hsurf(:,:)

!!!   call pmaxmin ( 'dep dz', dz, qmin, qmax, i2-i1+1, j2-j1+1, 1. )
    
   rmu = 1.8325e-5*(416.16/(tmpu(i1:i2,j1:j2,km)+120.)) &
                  *(tmpu(i1:i2,j1:j2,km)/296.16)**1.5

   z0h_min = 100. * tiny(1.0)  ! because sometimes we may get z0h=0.

!  =========================================================================
!  Aerodynamic Resistance
!  psi_h and Ra are equations 2, 4-5 of Walcek et al. 1986 Atmospheric Environment
!  obk supposedly from Wesely and Hicks 1977 J. Air Poll. Control Assoc.
!
   do j = j1, j2
    do i = i1, i2

!     Calculate the Monin-Obhukov length:
!            -Air denisity * Cp * T(surface) * Ustar^3
!     OBK = -------------------------------------------
!                 vK * g * Sensible heat flux
!     vK = 0.4               von Karman constant
!     Cp = 1000 J kg-1 K-1   specific heat of air at constant pressure
!     If OBK < 0 the air is unstable; if OBK > 0 the air is stable
!     For sensible heat flux of zero OBK goes to infinity (set to 1.e5)
      if(shflux(i,j) .eq. 0.) then
       obk = 1.e5
      else
       obk =   -rhoa(i,j,km)*cpd*tmpu(i,j,km)*ustar(i,j)**3. &
             / (von_karman*grav*shflux(i,j))
      endif

      cz = dz(i,j) / 2.
      frac = cz / obk
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
!  obk supposedly from Wesely and Hicks 1977 J. Air Poll. Control Assoc.

!  Loop over the bins
   do n = 1, nbins

    NAME = trim(w_c%reg%vname(nbeg+n-1))

    outflux = 0.0

!   Loop over space
    do j = j1, j2
     do i = i1, i2

!     Calculate the Monin-Obhukov length:
!            -Air denisity * Cp * T(surface) * Ustar^3
!     OBK = -------------------------------------------
!                 vK * g * Sensible heat flux
!     vK = 0.4               von Karman constant
!     Cp = 1000 J kg-1 K-1   specific heat of air at constant pressure
!     If OBK < 0 the air is unstable; if OBK > 0 the air is stable
!     For sensible heat flux of zero OBK goes to infinity (set to 1.e5)
      if(shflux(i,j) .eq. 0.) then
       obk = 1.e5
      else
       obk =   -rhoa(i,j,km)*cpd*tmpu(i,j,km)*ustar(i,j)**3. &
             / (von_karman*grav*shflux(i,j))
      endif

!     Calculate the surface resistance term
      vds = 0.002*ustar(i,j)
!     Set to small value of vds if ustar too small
      vds = max(vds, 0.002 * 0.00001)
      if(obk .lt. 0.) vds = vds*(1.+(-300./obk)**0.6667)
      czh = pblh(i,j)/obk
      if(czh .lt. -30.) vds = 0.0009*ustar(i,j)*(-czh)**0.6667
!     vdsMax is from Table 2 of Walcek et al. 1986
!     There are actually seasonal and regionally varying values,
!     but for most of the world a value of 1.0 cm s-1 is used.
      vdsMax = 0.01

      Rs(i,j) = 1./min(vds,vdsmax)

      if(Rs(i,j) .gt. 9999.) Rs(i,j) = 9999.
      if(Rs(i,j) .lt. 1.)    Rs(i,j) = 1.

!     If doing dust over land, possibly re-emit
      R2 = 1.
      if(trim(NAME(1:2)) .eq. 'du') then

!      Calculate the threshold velocity for dust emissions
       u_thresh0 = 0.13 * sqrt(rhop(n)*grav*2.*radius(n)/rhoa(i,j,km)) &
                        * sqrt(1.+6.e-7/(rhop(n)*grav*(2.*radius(n))**2.5)) &
              / sqrt(1.928*(1331.*(100.*2.*radius(n))**1.56+0.38)**0.092 - 1.)
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
      drydepf(i,j) = vdep(i,j) / dz(i,j)

!     Update the mass mixing ratio and compute the flux out
      factor = (1.0d0-exp(-drydepf(i,j)*r8_cdt))
      dc = max(zero,w_c%qa(n+nbeg-1)%data3d(i,j,km)*factor)
      w_c%qa(n+nbeg-1)%data3d(i,j,km) = w_c%qa(n+nbeg-1)%data3d(i,j,km) - dc

      outflux(i,j) = dc * w_c%delp(i,j,km)/grav/cdt

     end do  ! i
    end do   ! j

!   Diagnostic output if requested
    if( associated(fluxout(n)%data2d) ) fluxout(n)%data2d = outflux

   end do    ! n

   rc = 0

   end subroutine Chem_Deposition

   end module Chem_DepositionMod
