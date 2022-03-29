   subroutine SulfateChemDriver (km, klid, cdt, PI, radToDeg, von_karman, &
                                 airMolWght, nAvogadro, cpd, grav, &
                                 fMassMSA, fMassDMS, fMassSO2, fMassSO4, &
                                 nymd, nhms, lonRad, latRad, &
                                 dms, so2, so4, msa, &
                                 nDMS, nSO2, nSO4, nMSA, &
                                 xoh, xno3, xh2o2, h2o2_init, &
                                 delp, tmpu, cloud, rhoa, hghte, &
                                 ustar, shflux, oro, pblh, z0h, &
                                 SU_dep, SU_PSO2, SU_PMSA, &
                                 SU_PSO4, SU_PSO4g, SU_PSO4aq, &     ! 2d diagnostics
                                 pso2, pmsa, pso4, pso4g, pso4aq, drydepositionfrequency, & ! 3d diagnostics
                                 rc)


! !USES:
   implicit NONE

! !INPUT PARAMETERS:
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: PI     ! pi constnat
   real, intent(in)    :: radToDeg ! radians to degree conversion
   real, intent(in)    :: von_karman ! Von Karman constant [unitless]
   real, intent(in)    :: nAvogadro  ! Avogadro's number [1/kmol]
   real, intent(in)    :: airMolWght ! molecular weight of air [kg/kmol]
   real, intent(in)    :: cpd
   real, intent(in)    :: grav   ! gravity [m/sec]
   real, intent(in)    :: fMassMSA, fMassDMS, fMassSO2, fMassSO4 ! gram molecular weights of species
   integer, intent(in) :: nymd   ! model year month day
   integer, intent(in) :: nhms   ! model hour mintue second
   real, dimension(:,:), intent(in) :: lonRad   ! model grid lon [radians]
   real, dimension(:,:), intent(in) :: latRad   ! model grid lat [radians]
   real, dimension(:,:,:), intent(inout) :: dms  ! dimethyl sulfide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: so2  ! sulfer dioxide [kg/kg]
   real, dimension(:,:,:), intent(inout) :: so4  ! sulfate aerosol [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: msa  ! methanesulphonic acid [kg/kg]
   integer, intent(in) :: nDMS, nSO2, nSO4, nMSA ! index position of sulfates
   real, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]
   real, pointer, dimension(:,:,:), intent(in) :: tmpu   ! temperature [K]
   real, dimension(:,:,:), intent(in) :: cloud  ! cloud fraction for radiation [1]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa   ! layer air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: hghte  ! top of layer geopotential height [m]
   real, pointer, dimension(:,:), intent(in)   :: ustar  ! surface velocity scale [m/sec]
   real, pointer, dimension(:,:), intent(in)   :: shflux ! sensible heat flux from turbulence [w/m^2]
   real, pointer, dimension(:,:), intent(in)   :: oro    ! land-ocean-ice mask
   real, pointer, dimension(:,:), intent(in)   :: pblh   ! planetary boundary layer height [m]
   real, pointer, dimension(:,:), intent(in)   :: z0h    ! surface roughness for heat [m]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: xoh, xno3, xh2o2 ! OH, NO3, H2O2 respectievly [kg/kg]
   real, dimension(:,:,:) :: h2o2_init ! private H2O2 that is saved and used to initialize [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO2 ! SO2 Prod from DMS Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PMSA ! MSA Prod from DMS Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO4 ! SO4 Prod from All SO2 Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO4g ! SO4 Prod from Gaseous SO2 Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:), intent(inout)   :: SU_PSO4aq ! SO4 Prod from Aqueous SO2 Oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso2 ! SO2 Prod from DMS oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pmsa ! MSA Prod from DMS oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso4 ! SO4 Prod from all SO2 oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso4g ! SO4 Prod from gaseous SO2 oxidation [kg m-2 s-1]
   real, pointer, dimension(:,:,:), intent(inout) :: pso4aq ! SO4 Prod from aqueous SO2 oxidation [kg m-2 s-1]
   real, dimension(:,:), allocatable, intent(out) :: drydepositionfrequency

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc         ! Error return code:
                                                  !  0 - all is well
                                                  !  1 -

! !DESCRIPTION: Updates the SU concentration due to chemistry
!  The SU grid component is currently established with 4 different
!  species (bins) following this convection:
!   1) DMS
!   2) SO2
!   3) SO4
!   4) MSA
!  Accordingly we have 4 chemical cycles to follow through, which are
!  sub-subroutines under this one.
!  The chemistry is a function of OH, NO3, and H2O2 concentrations
!  as well as DMS, SO2, SO4, MSA concentrations.  It is also a function
!  of solar zenith angle and temperature.  We pass in temperature.  SZA
!  will be a function of time of day and lat/lon.  For now we simply add
!  this to the grid component before calculating it.  I bet this is
!  somewhere else in the model.

!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!  30july2020 E.Sherman - ported to process library
!

! !Local Variables
   real, dimension(:,:), allocatable :: cossza, sza
   integer :: k, jday, i2, j2
   real, dimension(:,:,:), allocatable :: pSO2_DMS, pMSA_DMS, pSO4g_SO2, pSO4aq_SO2
   real    :: xhour
   integer :: status


!EOP
!-------------------------------------------------------------------------
!  Begin

   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

   allocate(drydepositionfrequency, mold=oro)
   allocate(cossza, mold=oro)
   allocate(sza, mold=oro)

!  Reset the production terms
   allocate(pSO2_DMS, mold=tmpu)
   allocate(pMSA_DMS, mold=tmpu)
   allocate(pSO4g_SO2, mold=tmpu)
   allocate(pSO4aq_SO2, mold=tmpu)
   pSO2_DMS = 0.
   pMSA_DMS = 0.
   pSO4g_SO2 = 0.
   pSO4aq_SO2 = 0.

   if( associated(su_pSO2) )  su_pSO2   = 0.
   if( associated(su_pMSA) )  su_pMSA   = 0.
   if( associated(su_pSO4) )  su_pSO4   = 0.
   if( associated(su_pSO4g) )  su_pSO4g  = 0.
   if( associated(su_pSO4aq) )  su_pSO4aq = 0.
   if( associated(pSO2) )     pSO2   = 0.
   if( associated(pMSA) )     pMSA   = 0.
   if( associated(pSO4) )     pSO4   = 0.
   if( associated(pSO4g) )    pSO4g  = 0.
   if( associated(pSO4aq) )   pSO4aq = 0.


!  Find the cossza
!  ----------------------------------
   jday = idaynum(nymd)
   xhour = (  real(nhms/10000)*3600. &
            + real(mod(nhms,10000)/100)*60. &
            + real(mod(nhms,100)) &
           ) / 3600.

   call szangle (jday, xhour, lonRad, latRad, PI, radToDeg, sza, cossza, i2, j2)
!  Reset the dry deposition fluxes & frequencies
   if( associated(su_dep) ) su_dep = 0.0

   call DryDeposition ( km, tmpu, rhoa, hghte, oro, ustar, pblh, shflux, &
                        von_karman, cpd, grav, z0h, drydepositionfrequency, __RC__)


!  Now call the chemistry packages...
!  ----------------------------------
!  DMS source and oxidation to SO2 and MSA
   call SulfateChemDriver_DMS (km, klid, cdt, airMolWght, nAvogadro, cpd,&
                               fMassMSA, fMassDMS, fMassSO2, &
                               dms, nDMS, xoh, xno3, &
                               cossza, tmpu, rhoa, &
                               pSO2_DMS, pMSA_DMS, SU_dep, &
                               __RC__)

   if( associated(pSO2) )  pSO2 = pSO2_DMS
   if( associated(su_pSO2)) then
     do k = klid, km
      su_pSO2(:,:) = su_pSO2(:,:) + pSO2_DMS(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

   if( associated(pMSA) )  pMSA = pMSA_DMS
   if( associated(su_pMSA)) then
     do k = klid, km
      su_pMSA(:,:) = su_pMSA(:,:) + pMSA_DMS(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

!  SO2 source and oxidation to SO4
   call SulfateChemDriver_SO2 (km, klid, cdt, airMolWght, nAvogadro, cpd, grav, &
                               fMassSO4, fMassSO2, &
                               so2, nSO2, xoh, xh2o2, &
                               tmpu, rhoa, delp, oro, cloud, drydepositionfrequency, &
                               pSO2_DMS, pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                               __RC__)

   if( associated(pSO4g) )  pSO4g = pSO4g_SO2
   if( associated(su_pSO4g)) then
     do k = klid, km
      su_pSO4g(:,:) = su_pSO4g(:,:) + pSO4g_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

   if( associated(pSO4aq) )  pSO4aq = pSO4aq_SO2
   if( associated(su_pSO4aq)) then
     do k = klid, km
      su_pSO4aq(:,:) = su_pSO4aq(:,:) + pSO4aq_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

   if( associated(pSO4) ) pSO4 = pSO4g_SO2 + pSO4aq_SO2
   if( associated(su_pSO4)) then
     do k = klid, km
      su_pSO4(:,:) = su_pSO4(:,:) + pSO4g_SO2(:,:,k)*delp(:,:,k)/grav &
                     + pSO4aq_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

!  SO4 source and loss
   call SulfateChemDriver_SO4 (km, klid, cdt, grav, so4, nSO4, delp, &
                               drydepositionfrequency, pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                               __RC__)

!  MSA source and loss
   if( associated(msa)) then
      call SulfateChemDriver_MSA (km, klid, cdt, grav, msa, nMSA, delp, &
                                  drydepositionfrequency, pMSA_DMS, SU_dep, &
                                  __RC__)
   end if

!  Save the h2o2 value after chemistry
   h2o2_init = xh2o2

   __RETURN__(__SUCCESS__)
   end subroutine SulfateChemDriver 
