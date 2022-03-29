   subroutine SU_Compute_Diags ( km, klid, rmed, sigma, rhop, grav, pi, nSO4, mie, &
                                 wavelengths_profile, wavelengths_vertint, &
                                 tmpu, rhoa, delp, ple, tropp,rh, u, v, &
                                 DMS, SO2, SO4, MSA, &
                                 dmssfcmass, dmscolmass, &
                                 msasfcmass, msacolmass, &
                                 so2sfcmass, so2colmass, &
                                 so4sfcmass, so4colmass, &
                                 exttau, stexttau,scatau, stscatau,so4mass, so4conc, extcoef, &
                                 scacoef, angstrom, fluxu, fluxv, sarea, snum, rc )

! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km    ! number of model levels
   integer,    intent(in)    :: klid   ! index for pressure lid
   real, intent(in)    :: rmed  ! mean radius [um]
   real, intent(in)    :: sigma ! Sigma of lognormal number distribution
   real, intent(in)    :: rhop  ! dry particle density [kg m-3]
   real, intent(in)    :: grav  ! gravity [m/sec]
   real, intent(in)    :: pi    ! pi constant
   integer, intent(in) :: nSO4  ! index of SO4 relative to other internal variables
   type(GOCART2G_Mie), intent(in) :: mie   ! mie table
   real, dimension(:), intent(in)  :: wavelengths_profile
   real, dimension(:), intent(in)  :: wavelengths_vertint
   real, pointer, dimension(:,:,:), intent(in) :: tmpu    ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa    ! air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: delp    ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: ple   ! level edge air pressure [Pa]
   real, pointer, dimension(:,:), intent(in)   :: tropp ! tropopause pressure [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: rh      ! relative humidity [1]
   real, pointer, dimension(:,:,:), intent(in) :: u       ! east-west wind [m s-1]
   real, pointer, dimension(:,:,:), intent(in) :: v       ! north-south wind [m s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: DMS  ! dimethyl sulfide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO2  ! sulfer dioxide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: SO4  ! sulfate aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout)  :: MSA  ! methanesulphonic acid [kg/kg]
   real, pointer, dimension(:,:),   intent(inout)  :: dmssfcmass ! sfc mass concentration [kg/m3]
   real, pointer, dimension(:,:),   intent(inout)  :: dmscolmass ! col mass density [kg/m2]
   real, pointer, dimension(:,:),   intent(inout)  :: msasfcmass ! sfc mass concentration [kg/m3]
   real, pointer, dimension(:,:),   intent(inout)  :: msacolmass ! col mass density [kg/m2]
   real, pointer, dimension(:,:),   intent(inout)  :: so2sfcmass ! sfc mass concentration [kg/m3]
   real, pointer, dimension(:,:),   intent(inout)  :: so2colmass ! col mass density [kg/m2]
   real, pointer, dimension(:,:),   intent(inout)  :: so4sfcmass ! sfc mass concentration [kg/m3]
   real, pointer, dimension(:,:),   intent(inout)  :: so4colmass ! col mass density [kg/m2]
   real, pointer, dimension(:,:,:), intent(inout)  :: exttau     ! ext. AOT at 550 nm
   real, pointer, dimension(:,:,:), intent(inout)  :: stexttau   ! Stratosphere ext. AOT at 550 nm
   real, pointer, dimension(:,:,:), intent(inout)  :: scatau     ! sct. AOT at 550 nm
   real, pointer, dimension(:,:,:), intent(inout)  :: stscatau   ! Stratosphere sct. AOT at 550 nm
   real, pointer, dimension(:,:,:), intent(inout)  :: so4mass    ! 3D sulfate mass mr
   real, pointer, dimension(:,:,:), intent(inout)  :: so4conc    ! 3D mass concentration, [kg/m3]
   real, pointer, dimension(:,:,:,:), intent(inout)  :: extcoef    ! 3D ext. coefficient, [1/m]
   real, pointer, dimension(:,:,:,:), intent(inout)  :: scacoef    ! 3D scat.coefficient, [1/m]
   real, pointer, dimension(:,:),   intent(inout)  :: angstrom   ! 470-870 nm Angstrom parameter
   real, pointer, dimension(:,:),   intent(inout)  :: fluxu      ! Column mass flux in x direction
   real, pointer, dimension(:,:),   intent(inout)  :: fluxv      ! Column mass flux in y direction
   real, pointer, dimension(:,:,:), intent(inout)  :: sarea      ! Sulfate surface area density [m2 m-3]
   real, pointer, dimension(:,:,:), intent(inout)  :: snum       ! Sulfate number density [# m-2]
   integer, optional, intent(out)   :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -


! !DESCRIPTION: Calculates some simple 2d diagnostics from the SU fields
!
! !REVISION HISTORY:
!
!  16APR2004, Colarco
!  29july2020, E.Sherman - refactored for process library

! !Local Variables
   integer :: i, j, k, w, i1=1, j1=1, i2, j2, status
   real, dimension(:,:,:), allocatable :: tau, ssa
   real, dimension(:,:), allocatable :: tau470, tau870
   integer    :: ilam470, ilam870
   logical :: do_angstrom
   real :: rh_, gf, rwet, svol


!EOP
!-------------------------------------------------------------------------
!  Begin
   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

   allocate(tau470(i1:i2,j1:j2), tau870(i1:i2,j1:j2))

!  Get the wavelength indices
!  --------------------------

   ilam470 = mie%getChannel(4.70e-7)
   if(ilam470 <= 0) ilam470 = 0

   ilam870 = mie%getChannel(8.70e-7)
   if(ilam870 <= 0) ilam870 = 0

!  Determine if going to do Angstrom parameter calculation
!  -------------------------------------------------------
   do_angstrom = .false.
!  If both 470 and 870 channels provided (and not the same) then
!  possibly will do Angstrom parameter calculation
   if(ilam470 .ne. 0 .and. &
      ilam870 .ne. 0 .and. &
      ilam470 .ne. ilam870) do_angstrom = .true.


!  Calculate the diagnostic variables if requested
!  -----------------------------------------------

!  Calculate the surface mass concentration
   if( associated(so4sfcmass) ) then
      so4sfcmass(i1:i2,j1:j2) = 0.
      so4sfcmass(i1:i2,j1:j2) &
       =  SO4(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( associated(so2sfcmass) ) then
      so2sfcmass(i1:i2,j1:j2) = 0.
      so2sfcmass(i1:i2,j1:j2) &
       =   SO2(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( associated(dmssfcmass) ) then
      dmssfcmass(i1:i2,j1:j2) = 0.
      dmssfcmass(i1:i2,j1:j2) &
       =   DMS(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif
   if( associated(msasfcmass) .and. associated(MSA)) then
      msasfcmass(i1:i2,j1:j2) = 0.
      msasfcmass(i1:i2,j1:j2) &
       =   MSA(i1:i2,j1:j2,km)*rhoa(i1:i2,j1:j2,km)
   endif


!  Initialize the diagnostic variables
!  -----------------------------------

!  Calculate the column loading
   if( associated(so4colmass) ) then
      so4colmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       so4colmass(i1:i2,j1:j2) &
        =   so4colmass(i1:i2,j1:j2) &
          + SO4(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( associated(so2colmass) ) then
      so2colmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       so2colmass(i1:i2,j1:j2) &
        =   so2colmass(i1:i2,j1:j2) &
          + SO2(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( associated(dmscolmass) ) then
      dmscolmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       dmscolmass(i1:i2,j1:j2) &
        =   dmscolmass(i1:i2,j1:j2) &
          + DMS(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif
   if( associated(msacolmass) .and. associated(MSA)) then
      msacolmass(i1:i2,j1:j2) = 0.
      do k = klid, km
       msacolmass(i1:i2,j1:j2) &
        =   msacolmass(i1:i2,j1:j2) &
          + MSA(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav
      enddo
   endif


!  Calculate the mass concentration of sulfate
   if( associated(so4conc) ) then
      so4conc(i1:i2,j1:j2,1:km) = 0.
      so4conc(i1:i2,j1:j2,1:km) = SO4(i1:i2,j1:j2,1:km)*rhoa(i1:i2,j1:j2,1:km)
   endif

!  Mass mixing ratio of sulfate
   if( associated(so4mass) ) then
      so4mass(i1:i2,j1:j2,1:km) = 0.
      so4mass(i1:i2,j1:j2,1:km) = SO4(i1:i2,j1:j2,1:km)
   endif

!  Calculate the column mass flux in x direction
   if( associated(fluxu) ) then
      fluxu(i1:i2,j1:j2) = 0.
       do k = klid, km
        fluxu(i1:i2,j1:j2) &
         =   fluxu(i1:i2,j1:j2) &
           + SO4(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav*u(i1:i2,j1:j2,k)
       end do
   endif

!  Calculate the column mass flux in y direction
   if( associated(fluxv) ) then
      fluxv(i1:i2,j1:j2) = 0.
       do k = klid, km
        fluxv(i1:i2,j1:j2) &
         =   fluxv(i1:i2,j1:j2) &
           + SO4(i1:i2,j1:j2,k)*delp(i1:i2,j1:j2,k)/grav*v(i1:i2,j1:j2,k)
       end do
   endif

!  Calculate the extinction and/or scattering AOD
   allocate(tau(i1:i2,j1:j2,km), source = 0.)
   allocate(ssa(i1:i2,j1:j2,km), source = 0.)
   if( associated(extcoef) .or. associated(scacoef) ) then

      if (associated(extcoef)) extcoef = 0.
      if (associated(scacoef)) scacoef = 0.

      do w = 1, size(wavelengths_profile)
         call mie%Query(wavelengths_profile(w), 1, & ! Only SO4 exists in the MieTable, so its index is 1
                        SO4*delp/grav, rh,         &
                        tau=tau, ssa=ssa, __RC__)

!         Calculate the total ext. and scat. coefficients
         if( associated(extcoef) ) then
              extcoef(:,:,:,w) = extcoef(:,:,:,w) + &
                              tau * (grav * rhoa / delp)
         endif
         if( associated(scacoef) ) then
              scacoef(:,:,:,w) = scacoef(:,:,:,w) + &
                              ssa * tau * (grav * rhoa / delp)
         endif
      enddo
   endif

   if( associated(exttau) .or. associated(stexttau) .or. &
       associated(scatau) .or. associated(stscatau)) then

      if (associated(exttau)) exttau = 0.
      if (associated(stexttau)) stexttau = 0.
      if (associated(scatau)) scatau = 0.
      if (associated(stscatau)) stscatau = 0.

      do w = 1, size(wavelengths_vertint)
         call mie%Query(wavelengths_vertint(w), 1,  & ! Only SO4 exists in the MieTable, so its index is 1
                        SO4*delp/grav, rh,          &
                        tau=tau, ssa=ssa, __RC__)

         do k = klid, km
!           Integrate in the vertical
            if ( associated(exttau) ) then
               exttau(:,:,w) = exttau(:,:,w) + tau(:,:,k)
            endif

            if (associated(stexttau) ) then
               where (ple(:,:,k) .le. tropp) 
                  stexttau(:,:,w) = stexttau(:,:,w) + tau(:,:,k)
               elsewhere(ple(:,:,k-1) .lt. tropp) 
                 stexttau(:,:,w)  = stexttau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)
               endwhere
            endif

            if ( associated(scatau) ) then
               scatau(:,:,w) = scatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
            endif

            if ( associated(stscatau) ) then
               where (ple(:,:,k) .le. tropp) 
                  stscatau(:,:,w) = stscatau(:,:,w) + tau(:,:,k)*ssa(:,:,k)
               elsewhere(ple(:,:,k-1) .lt. tropp) 
                  stscatau(:,:,w) = stscatau(:,:,w) + log(tropp/ple(:,:,k-1))/log(ple(:,:,k)/ple(:,:,k-1))*tau(:,:,k)*ssa(:,:,k)
               endwhere
           endif
         enddo
      enddo
   endif

!  Calculate the 470-870 Angstrom parameter
   if( associated(angstrom) .and. do_angstrom ) then

      angstrom(i1:i2,j1:j2) = 0.
!     Set tau to small number by default
      tau470(i1:i2,j1:j2) = tiny(1.0)
      tau870(i1:i2,j1:j2) = tiny(1.0)

      call mie%Query(4.70E-7,  1,       & ! Only SO4 exists in the MieTable, so its index is 1
                     SO4*delp/grav, rh, &
                     tau=tau, __RC__)
      do k = klid, km
         tau470 = tau470 + tau(:,:,k)
      enddo

      call mie%Query(8.70E-7, 1,       &
                     SO4*delp/grav,rh, &
                     tau=tau, __RC__)
      do k = klid, km
         tau870 = tau870 + tau(:,:,k)
      enddo

!      enddo  ! nbins
      angstrom(i1:i2,j1:j2) = &
        -log(tau470(i1:i2,j1:j2)/tau870(i1:i2,j1:j2)) / &
         log(470./870.)
   endif

!  Calculate the sulfate surface area density [m2 m-3], possibly for use in
!  StratChem or other component.  Assumption here is a specified effective
!  radius (gcSU%radius for sulfate) and standard deviation of lognormal
!  distribution.  Hydration is by grid box provided RH and is follows Petters
!  and Kreeidenweis (ACP2007)
   if(associated(sarea) .or. associated(snum)) then
!        rmed   = w_c%reg%rmed(n1+nSO4-1)                    ! median radius, m
        if(rmed > 0.) then
!         sigma  = w_c%reg%sigma(n1+nSO4-1)                  ! width of lognormal distribution
         do k = klid, km
         do j = j1, j2
          do i = i1, i2
           rh_ = min(0.95,rh(i,j,k))
           gf = (1. + 1.19*rh_/(1.-rh_) )                   ! ratio of wet/dry volume, eq. 5
           rwet = rmed * gf**(1./3.)                      ! wet effective radius, m
!          Wet particle volume m3 m-3
           svol = SO4(i,j,k) * rhoa(i,j,k) / rhop * gf
!          Integral of lognormal surface area m2 m-3
           if(associated(sarea)) sarea(i,j,k) = 3./rwet*svol*exp(-5./2.*alog(sigma)**2.)
!          Integral of lognormal number density # m-3
           if(associated(snum)) snum(i,j,k) = svol / (rwet**3) * exp(-9/2.*alog(sigma)**2.) * 3./4./pi
          enddo
         enddo
        enddo
       endif
   endif
   deallocate(tau,ssa)
   __RETURN__(__SUCCESS__)
   end subroutine SU_Compute_Diags
