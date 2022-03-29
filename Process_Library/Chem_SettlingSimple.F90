   subroutine Chem_SettlingSimple ( km, klid, flag, cdt, grav, &
                                    radiusInp, rhopInp, int_qa, tmpu, &
                                    rhoa, rh, hghte, delp, fluxout,  &
                                    vsettleOut, correctionMaring, rc)

! !USES:

  implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in)    :: km     ! total model levels
   integer, intent(in)    :: klid   ! index for pressure lid
   integer, intent(in) :: flag     ! flag to control particle swelling (see note)
   real, intent(in)    :: cdt
   real, intent(in)    :: grav   ! gravity [m/sec^2]
   real, intent(in)  :: radiusInp  ! particle radius [microns]
   real, intent(in)  :: rhopInp    ! soil class density [kg/m^3]
   real, dimension(:,:,:), intent(inout) :: int_qa  ! aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(in)  :: tmpu   ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in)  :: rhoa   ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in)  :: rh     ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in)  :: hghte  ! geopotential height [m]
   real, pointer, dimension(:,:,:), intent(in)  :: delp   ! pressure level thickness [Pa]

! !OUTPUT PARAMETERS:

   real, pointer, dimension(:,:), intent(inout)  :: fluxout ! Mass lost by settling
                                                  ! to surface, kg/m2/s
   integer, optional, intent(out)             :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -
!  Optionally output the settling velocity calculated
   real, pointer, optional, dimension(:,:,:)  :: vsettleOut

!  Optionally correct the settling velocity following Maring et al, 2003
   logical, optional, intent(in)    :: correctionMaring

   character(len=*), parameter :: myname = 'SettlingSimple'

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
   integer  ::  i, j, k, n, dk

   real, parameter ::  rhow = 1000.  ! Density of water [kg m-3]

!  parameters from Maring et al, 2003
   real, parameter :: v_upwardMaring = 0.33e-2   ! upward velocity, [m s-1]
   real, parameter :: diameterMaring = 7.30e-6   ! particle diameter, [m]

   integer         :: i1=1, i2, j1=1, j2
   integer         :: dims(3)

   real, pointer, dimension(:,:)   :: hsurf
   real(kind=DP), dimension(:,:), allocatable  :: cmass_before, cmass_after
   real, allocatable    :: dz(:,:,:)
   real, dimension(:,:,:), allocatable  :: radius, rhop, qa
   real, dimension(:,:,:), allocatable  :: vsettle   ! fall speed [m s-1]
   real ::  ONE_OVER_G
   integer :: status

!EOP
!-------------------------------------------------------------------------

!  Get dimensions
!  ---------------
   dims = shape(rhoa)
   i2 = dims(1); j2 = dims(2)

   ONE_OVER_G = 1.0/grav

   hsurf => hghte(i1:i2,j1:j2,km)

   allocate(dz(i2,j2,km), radius(i2,j2,km), rhop(i2,j2,km), vsettle(i2,j2,km), qa(i2,j2,km))
   allocate(cmass_before(i2,j2), cmass_after(i2,j2))
   cmass_before = 0.d0
   cmass_after = 0.d0

   qa = int_qa

   if(associated(fluxout)) fluxout(:,:) = 0.0

!  Handle the fact that hghte may be in the range [1,km+1] or [0,km]
!  -----------------------------------------------------------------
   dk = lbound(hghte,3) - 1  ! This is either 0 or 1

!  Layer thickness from hydrostatic equation
   k = km
   dz(:,:,k) = hghte(:,:,k+dk)-hsurf(:,:)
   do k = km-1, 1, -1
    dz(:,:,k) = hghte(:,:,k+dk) - hghte(:,:,k+dk+1)
   enddo

!  If radius le 0 then get out
   if(radiusInp .le. 0.) then
      status = 100
      __RETURN__(STATUS)
   end if

!   Find the column dry mass before sedimentation
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             cmass_before(i,j) = cmass_before(i,j) + qa(i,j,k) * delp(i,j,k) * ONE_OVER_G
          enddo
       enddo
    enddo

!   Particle swelling
    call ParticleSwelling(i1, i2, j1, j2, km, rh, radiusInp, rhopInp, radius, rhop, flag)

!   Settling velocity of the wet particle
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             call Chem_CalcVsettle(radius(i,j,k), rhop(i,j,k), rhoa(i,j,k), &
                                   tmpu(i,j,k), vsettle(i,j,k), grav)
          end do
       end do
    end do

    if(present(correctionMaring)) then
       if (correctionMaring) then
          vsettle = max(1.0e-9, vsettle - v_upwardMaring)
       endif
    endif

    if(present(vsettleOut)) then
       vsettleOut = vsettle
    endif

!   Time integration
    call SettlingSolver(i1, i2, j1, j2, km, cdt, delp, dz, vsettle, qa)

!   Find the column dry mass after sedimentation and thus the loss flux
    do k = klid, km
       do j = j1, j2
          do i = i1, i2
             cmass_after(i,j) = cmass_after(i,j) + qa(i,j,k) * delp(i,j,k) * ONE_OVER_G
          enddo
       enddo
    enddo

    if( associated(fluxout) ) then
       fluxout(:,:) = (cmass_before - cmass_after)/cdt
    endif

    int_qa = qa

   __RETURN__(__SUCCESS__)

   end subroutine Chem_SettlingSimple
