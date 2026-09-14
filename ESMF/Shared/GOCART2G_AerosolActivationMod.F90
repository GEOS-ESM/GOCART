#include "MAPL_Generic.h"

!=============================================================================
! Aerosol activation properties for GOCART2G.
!
! This module implements aerosol properties required for aerosol activation 
!used through the ESMF callback interface and method.
!=============================================================================
module GOCART2G_AerosolActivationMod

   use ESMF
   use MAPL
   use Chem_AeroGeneric

   implicit none
   private

   public :: setup_aerosol_activation
   public :: aerosol_activation_properties

contains

  subroutine setup_aerosol_activation(aero, grid, CF, rc)

    implicit none

    type(ESMF_State),  intent(inout) :: aero
    type(ESMF_Grid),   intent(inout) :: grid
    type(ESMF_Config), intent(inout) :: CF
    integer,           intent(out)   :: rc

    integer                                :: n_modes
    integer, parameter                     :: n_gocart_modes = 14
    integer, parameter                     :: n_mamnet_modes = 7
    character(len=ESMF_MAXSTR)             :: aero_aci_modes(n_gocart_modes)
    character(len=ESMF_MAXSTR)             :: mamnet_aci_modes(n_mamnet_modes)
    real                                   :: maxclean, ccntuning
    logical                                :: use_mamnet

    __Iam__('setup_aerosol_activation')

!   Begin adding necessary aerosol cloud interaction information
!   ------------------------------------------------------------
    aero_aci_modes =  (/'du001    ', 'du002    ', 'du003    ', &
                        'du004    ', 'du005    ',              &
                        'ss001    ', 'ss002    ', 'ss003    ', &
                        'sulforg01', 'sulforg02', 'sulforg03', &
                        'bcphilic ', 'ocphilic ', 'brcphilic'/)

    mamnet_aci_modes = (/'NUM_A_ACC', 'NUM_A_AIT', 'NUM_A_CDU', 'NUM_A_CSS', &
                         'NUM_A_FDU', 'NUM_A_FSS', 'NUM_A_PCM'/)

    call ESMF_ConfigGetAttribute(CF, use_mamnet, default=.FALSE., label='USE_MAMNET:', __RC__)

    if (use_mamnet) then
       n_modes = size(mamnet_aci_modes)
       call ESMF_AttributeSet(aero, name='number_of_aerosol_modes', value=n_modes, __RC__)
       call ESMF_AttributeSet(aero, name='aerosol_modes', itemcount=n_modes, valuelist=mamnet_aci_modes, __RC__)
       call WRITE_PARALLEL ('Using MAMnet')

    else
       n_modes = size(aero_aci_modes)
       call ESMF_AttributeSet(aero, name='number_of_aerosol_modes', value=n_modes, __RC__)
       call ESMF_AttributeSet(aero, name='aerosol_modes', itemcount=n_modes, valuelist=aero_aci_modes, __RC__)
    end if

    ! max mixing ratio before switching to "polluted" size distributions
    call ESMF_ConfigGetAttribute(CF, maxclean, default=1.0e-9, label='MAXCLEAN:', __RC__)
    call ESMF_AttributeSet(aero, name='max_q_clean', value=maxclean, __RC__)

    call ESMF_ConfigGetAttribute(CF, CCNtuning, default=1.8, label='CCNTUNING:', __RC__)
    call ESMF_AttributeSet(aero, name='ccn_tuning', value=CCNtuning, __RC__)

!   Add variables to AERO state
    call add_aero (aero, label='air_temperature', label2='T', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_land_type', label2='FRLAND', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='width_of_aerosol_mode', label2='SIGMA', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_number_concentration', label2='NUM', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_dry_size', label2='DGN', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_density', label2='density', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='aerosol_hygroscopicity', label2='KAPPA', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_dust_aerosol', label2='FDUST', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_soot_aerosol', label2='FSOOT', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='fraction_of_organic_aerosol', label2='FORGANIC', grid=grid, typekind=MAPL_R4, __RC__)

!   Attach aerosol activation method. Keep the callback label and signature unchanged.
    call ESMF_MethodAdd(aero, label='aerosol_activation_properties', userRoutine=aerosol_activation_properties, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine setup_aerosol_activation

  subroutine aerosol_activation_properties(state, rc)

    implicit none

!   Arguments
!   ---------
    type(ESMF_State)     :: state
    integer, intent(out) :: rc

!   Local
!   ---------
    character(len=ESMF_MAXSTR)      :: mode              ! mode name
    character(len=ESMF_MAXSTR)      :: mode_             ! lowercase mode name

    type(ESMF_State)                :: child_state

    real, dimension(:,:,:), pointer :: ple               ! pressure at the edges of model layers
    real, dimension(:,:,:), pointer :: temperature       ! air temperature
    real, dimension(:,:,:), allocatable :: air_density   ! air density calculated on the fly
    real, dimension(:,:),   pointer :: f_land            ! fraction of land type in a grid cell

    real, dimension(:,:,:), pointer :: f                 ! correction factor for sea salt

    real, dimension(:,:,:), allocatable :: q             ! aerosol mass mixing ratio
    real, dimension(:,:,:,:), pointer   :: ptr_4d        ! aerosol mass mixing ratio (temporary)
    real, dimension(:,:,:), pointer     :: ptr_3d        ! aerosol mass mixing ratio (temporary)

    real, dimension(:,:,:), pointer :: num               ! number concentration of aerosol particles
    real, dimension(:,:,:), pointer :: diameter          ! dry size of aerosol
    real, dimension(:,:,:), pointer :: sigma             ! width of aerosol mode
    real, dimension(:,:,:), pointer :: density           ! density of aerosol
    real, dimension(:,:,:), pointer :: hygroscopicity    ! hygroscopicity of aerosol
    real, dimension(:,:,:), pointer :: f_dust            ! fraction of dust aerosol
    real, dimension(:,:,:), pointer :: f_soot            ! fraction of soot aerosol
    real, dimension(:,:,:), pointer :: f_organic         ! fraction of organic aerosol

    real                            :: max_clean          ! max mixing ratio before considered polluted
    real                            :: ccn_tuning         ! tunes conversion factors for sulfate

    character(len=ESMF_MAXSTR)      :: fld_name

    integer                         :: i2, j2, km
    integer                         :: b, i, j, n, aerosol_bin
    integer                         :: ple_k0
    integer                         :: varNameLen

    character (len=ESMF_MAXSTR), allocatable  :: itemList(:), aeroList(:)
    type (ESMF_StateItem_Flag), allocatable   :: itemTypes(:)

!   auxilliary parameters
!   ---------------------
    real, parameter :: densSO4 = 1700.0
    real, parameter :: densORG = 1600.0
    real, parameter :: densSS  = 2200.0
    real, parameter :: densDU  = 1700.0
    real, parameter :: densBC  = 1600.0
    real, parameter :: densOC  =  900.0
    real, parameter :: densBR  =  900.0

    real, parameter :: k_SO4   = 0.65
    real, parameter :: k_ORG   = 0.20
    real, parameter :: k_SS    = 1.28
    real, parameter :: k_DU    = 0.0001
    real, parameter :: k_BC    = 0.0001
    real, parameter :: k_OC    = 0.0001
    real, parameter :: k_BR    = 0.0001

    real, parameter :: rd_air  = 287.04   ! dry-air gas constant [J kg-1 K-1]

!   PySR/MAMnet emulator unit-conversion and numerical guard constants
    real, parameter :: mam_mass_kgkg_to_ugkg     = 1.0e9
    real, parameter :: mam_mass_ugkg_to_kgkg     = 1.0e-9
    real, parameter :: mam_air_kgm3_to_mgm3      = 1.0e6
    real, parameter :: mam_min_input_mass_kgkg   = 1.0e-30
    real, parameter :: mam_min_air_density_kgm3  = 1.0e-3
    real, parameter :: mam_max_output_mass_kgkg  = 1.0e-3
    real, parameter :: mam_max_number_m3         = 1.0e12
    real, parameter :: mam_min_diameter_m        = 1.0e-9
    real, parameter :: mam_max_diameter_m        = 20.0e-6
    real, parameter :: mam_min_kappa             = 0.0
    real, parameter :: mam_max_kappa             = 2.0

    integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015

    __Iam__('GOCART2G::aerosol_activation_properties')

!   Begin...

!   Get list of child states within state and add to aeroList
!   ---------------------------------------------------------
    call ESMF_StateGet (state, itemCount=n, __RC__)
    allocate (itemList(n), __STAT__)
    allocate (itemTypes(n), __STAT__)
    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes, __RC__)

    b=0
    do i = 1, n
       if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= 'NI')) then
          b = b + 1
       end if
    end do

    allocate (aeroList(b), __STAT__)

    j = 1
    do i = 1, n
       if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= 'NI')) then
          aeroList(j) = trim(itemList(i))
          j = j + 1
       end if
    end do

!   Aerosol mode
!   ------------
    call ESMF_AttributeGet(state, name='aerosol_mode', value=mode, __RC__)

!   Land fraction
!   -------------
    call ESMF_AttributeGet(state, name='fraction_of_land_type', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_land, trim(fld_name), __RC__)

!   Pressure at layer edges
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

!   Temperature
!   -----------
    call ESMF_AttributeGet(state, name='air_temperature', value=fld_name, __RC__)
    call MAPL_GetPointer(state, temperature, trim(fld_name), __RC__)

    i2 = ubound(temperature, 1)
    j2 = ubound(temperature, 2)
    km = ubound(temperature, 3)

!   Activation activation properties
!   --------------------------------
    call ESMF_AttributeGet(state, name='aerosol_number_concentration', value=fld_name, __RC__)
    call MAPL_GetPointer(state, num, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_dry_size', value=fld_name, __RC__)
    call MAPL_GetPointer(state, diameter, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='width_of_aerosol_mode', value=fld_name, __RC__)
    call MAPL_GetPointer(state, sigma, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_density', value=fld_name, __RC__)
    call MAPL_GetPointer(state, density, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_hygroscopicity', value=fld_name, __RC__)
    call MAPL_GetPointer(state, hygroscopicity, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_dust_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_dust, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_soot_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_soot, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_organic_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_organic, trim(fld_name), __RC__)

!   Sea salt scaling fctor
!   ----------------------
    call ESMF_AttributeGet(state, name='max_q_clean', value=max_clean, __RC__)
    call ESMF_AttributeGet(state, name='ccn_tuning', value=ccn_tuning, __RC__)

!   Aerosol mass mixing ratios
!   --------------------------
    mode_ = trim(mode)
    mode_ = ESMF_UtilStringLowerCase(mode_, __RC__)

    if (index(mode_, 'num_a_') == 1) then

!      Compute dry-air density.
       ple_k0 = lbound(ple, 3)
       allocate(air_density(1:i2,1:j2,km), __STAT__)
       do aerosol_bin = 1, km
          air_density(:,:,aerosol_bin) = 0.5 * (ple(:,:,ple_k0+aerosol_bin-1) + &
                                               ple(:,:,ple_k0+aerosol_bin)) / &
                                         (rd_air * max(temperature(:,:,aerosol_bin), 180.0))
       end do

       call mamnet_emulator_(mode_, state, aeroList, temperature, air_density, &
                             num, diameter, sigma, hygroscopicity, density, &
                             f_dust, f_soot, f_organic, 1, i2, 1, j2, km, __RC__)
       deallocate(air_density, itemList, itemTypes, aeroList, __STAT__)
       RETURN_(ESMF_SUCCESS)
    end if

    ! these modes are only called is mamnet is not active
    allocate(q(i2,j2,km),  __STAT__)
    q = 0.0

    if (index(mode_, 'du00') > 0) then ! Dust
       ! dust is mapped one-to-one
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'DU') > 0) then
             read (mode_(3:len(mode_)),*) aerosol_bin
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             call MAPL_GetPointer(child_state, ptr_4d, 'DU', __RC__)
             q = q + ptr_4d(:,:,:,aerosol_bin)
             ptr_3d => ptr_4d(:,:,:,aerosol_bin)

             hygroscopicity = k_DU
             density = densDU
          end if
       end do

    else if (index(mode_, 'ss00') > 0) then ! Sea Salt
       ! compute the total mass mixing ratio and impose a tri-modal size distribution
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'SS') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             call MAPL_GetPointer(child_state, ptr_4d, 'SS', __RC__)
             do j = 1, ubound(ptr_4d, 4)
               q = q + ptr_4d(:,:,:,j)
               ptr_3d => ptr_4d(:,:,:,j)
             end do

             hygroscopicity = k_SS
             density = densSS
          end if
       end do

    else if (index(mode_, 'sulforg') > 0) then ! Sulfate
       hygroscopicity = 0.0
       density = 0.0

       do i = 1, size(aeroList)
          if (index(aeroList(i), 'SU') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             call MAPL_GetPointer(child_state, ptr_3d, 'SO4', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_SO4 * ptr_3d + hygroscopicity
             density = densSO4 * ptr_3d + density
          end if

          if (index(aeroList(i), 'CA.oc') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.oc, or CA.oc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_ORG * ptr_3d + hygroscopicity
             density = densORG * ptr_3d + density
          end if

       end do

          where (q > 2.0e-12 .and. hygroscopicity > tiny(0.0))
             hygroscopicity = hygroscopicity / q
             hygroscopicity = max(0.001, hygroscopicity)

             density = density / q
             density = min(max(density, densORG), densSO4)
          elsewhere
             hygroscopicity = k_SO4
             density = densSO4
          end where

    else if (index(mode_, 'bcphilic') > 0) then ! Black Carbon
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'CA.bc') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.bc, or CA.bc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_BC
             density = densBC
          end if
       end do

    else if (index(mode_, 'ocphilic') > 0) then ! Organic Carbon
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'CA.oc') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.oc, or CA.oc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_OC
             density = densOC
          end if
       end do

    else if (index(mode_, 'brcphilic') > 0) then ! Organic Carbon
       do i = 1, size(aeroList)
          if (index(aeroList(i), 'CA.br') > 0) then
             call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)
             varNameLen = len_trim(aeroList(i))
!            the '5' refers to '_AERO', which we want to remove to get the CA component name (e.g. CA.bc, or CA.bc.data)
             varNameLen = varNameLen - 5
             call MAPL_GetPointer(child_state, ptr_3d, aeroList(i)(1:varNameLen)//'philic', __RC__)
             q = q + ptr_3d
             hygroscopicity = k_BR
             density = densBR
          end if
       end do

    end if !(index(mode_, 'du00') > 0) then

!   Obtain aerosol activation properties of this aerosol mode
!   ---------------------------------------------------------
    call aap_(mode,               &
              q,                  &
              num,                &
              diameter,           &
              sigma,              &
              f_dust,             &
              f_soot,             &
              f_organic,          &
              density,            &
              ptr_3d,             &
              1, i2, 1, j2, km, &
              __RC__)

    deallocate(q, __STAT__)

    RETURN_(ESMF_SUCCESS)

   contains

    subroutine aap_(mode, q, num, diameter, sigma, f_dust, f_soot, f_organic, dens_, q_, &
                    i1, i2, j1, j2, km, rc)

     implicit none

     integer, intent(in) :: i1, i2                                  ! dimension bounds
     integer, intent(in) :: j1, j2                                  ! ... // ..
     integer, intent(in) :: km                                      ! ... // ..

     character(len=*),  intent(in )               :: mode           ! name of aerosol mode
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: q              ! aerosol mass mixing ratio, kg kg-1
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: q_             ! auxiliary mass
     real, intent(in),  dimension(i1:i2,j1:j2,km) :: dens_          ! density

     real, intent(out), dimension(i1:i2,j1:j2,km) :: num            ! number concentration of aerosol particles
     real, intent(out), dimension(i1:i2,j1:j2,km) :: diameter       ! dry size of aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: sigma          ! width of aerosol mode
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_dust         ! fraction of dust aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_soot         ! fraction of soot aerosol
     real, intent(out), dimension(i1:i2,j1:j2,km) :: f_organic      ! fraction of organic aerosol

     integer, intent(out) :: rc                                     ! return code

     ! local
     integer :: STATUS
     character(len=ESMF_MAXSTR) :: mode_
     character(len=ESMF_MAXSTR) :: Iam = 'GOCART::aerosol_activation_properties::aap_()'

     integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015

     integer            :: kinx
     real               :: fmassaux, fmassclean
     real, dimension(3) :: TPI, DPGI, SIGI
     real, dimension(3) :: TPIclean, DPGIclean, SIGIclean
     real, dimension(i1:i2,j1:j2,km) :: qaux
      !real, parameter    :: max_clean = 5.0e-7  !max mixing ratio before considered polluted

     mode_ = trim(mode)
     mode_ = ESMF_UtilStringLowerCase(mode_, __RC__)

     num       = 0.0
     diameter  = 1.0e-9
     sigma     = log(2.0)
     f_dust    = 0.0
     f_soot    = 0.0
     f_organic = 0.0

     qaux=q !this corrects a bug

     if (index(mode_, 'ss00') > 0) then
         TPI  (1) = 230e6          ! num fraction (reduced 091015)
         DPGI (1) = 0.02e-6        ! modal diameter (m)
         SIGI (1) = log(1.6)       ! geometric dispersion (sigma_g)
         ! accumulation
         TPI  (2) = 60.0e6         ! total concentration (# m-3)
         DPGI (2) = 0.071e-6       ! modal diameter (m)
         SIGI (2) = log(2.0)       ! geometric dispersion (sigma_g)
         ! coarse
         TPI  (3) = 3.1e6          ! total concentration (# m-3)
         DPGI (3) = 0.62e-6        ! modal diameter (m)
         SIGI (3) = log(2.7)       ! geometric dispersion (sigma_g)

         fmassaux = 0.0
         do kinx = 1, 3
             fmassaux = (TPI(kinx)*densSS*MAPL_PI*exp(4.5*SIGI(kinx)*SIGI(kinx))*DPGI(kinx)*DPGI(kinx)*DPGI(kinx))/6.0 + fmassaux
         end do
     end if

     if (index(mode_, 'sulforg0') > 0) then
         TPI  (1) = 1.06e11        ! num fraction
         DPGI (1) = .014e-6        ! modal diameter (m)
         SIGI (1) = log(1.8)       ! geometric dispersion (sigma_g)
         ! accumulation
         TPI  (2) = 3.2e10         ! total concentration (# m-3)
         DPGI (2) = 0.054e-6       ! modal diameter (m)
         SIGI (2) = log(2.16)      ! geometric dispersion (sigma_g)
         !coarse
         TPI  (3) = 5.4e6          ! total concentration (# m-3)
         DPGI (3) = 0.86e-6        ! modal diameter (m)
         SIGI (3) = log(2.21)      ! geometric dispersion (sigma_g)

         fmassaux = 0.0
         do kinx = 1, 3
             ! density is multiplied below since this is a case of a 3-d field
             fmassaux = (TPI(kinx)*MAPL_PI*exp(4.5*SIGI(kinx)*SIGI(kinx))*DPGI(kinx)*DPGI(kinx)*DPGI(kinx))/6.0 + fmassaux
         end do

         ! clean continental polluted plus org
         ! fine
         TPIclean  (1) = 1.0e9      ! total concentration (# m-3)
         DPGIclean (1) = 0.016e-6   ! modal diameter (m)
         SIGIclean (1) = log(1.6)   ! geometric dispersion (sigma_g)
         ! accumulation
         TPIclean  (2) = 8.0e8      ! total concentration (# m-3)
         DPGIclean (2) = 0.067e-6   ! modal diameter (m)
         SIGIclean (2) = log(2.1)   ! geometric dispersion (sigma_g)
         !Coarse
         TPIclean  (3) = 2.0e6      ! total concentration (# m-3)
         DPGIclean (3) = 0.93e-6    ! modal diameter (m)
         SIGIclean (3) = log(2.2)   ! geometric dispersion (sigma_g)

         fmassclean= 0.0
         do kinx = 1, 3
             fmassclean = (TPIclean(kinx)*MAPL_PI*exp(4.5*SIGIclean(kinx)*SIGIclean(kinx))*DPGIclean(kinx)*DPGIclean(kinx)*DPGIclean(kinx))/6.0 + fmassclean  !
         end do
     end if

     select case(mode_)

     case ('du001')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 1.46e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du002')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 2.80e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du003')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 4.80e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du004')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 9.0e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('du005')
         sigma    = log(1.8)
         f_dust   = 1.0
         diameter = 16.0e-6
         num      = q / ((MAPL_PI/6.0) * densDU * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('ss001')
         sigma    = SIGI(1)
         diameter = DPGI(1)
         num      = TPI(1) * q / fmassaux

     case ('ss002')
         sigma    = SIGI(2)
         diameter = DPGI(2)
         num      = TPI(2) * q / fmassaux

     case ('ss003')
         sigma    = SIGI(3)
         diameter = DPGI(3)
         num      = TPI(3) * q / fmassaux

     case ('sulforg01')  !different distributions for clean and polluted environments
         where (q > max_clean)
             sigma    = SIGI(1)
             diameter = DPGI(1)
             num      = TPI(1) * qaux*ccn_tuning / (dens_*fmassaux)             ! only sulfate  mass
         elsewhere
             sigma    = SIGIclean(1)
             diameter = DPGIclean(1)
             num      = TPIclean(1) * qaux*ccn_tuning / (dens_*fmassclean)      ! only sulfate
         end where

     case ('sulforg02')
         where (q > max_clean)
             sigma    = SIGI(2)
             diameter = DPGI(2)
             num      = TPI(2) * qaux*ccn_tuning / (dens_*fmassaux)            ! only sulfate mass
         elsewhere
             sigma    = SIGIclean(2)
             diameter = DPGIclean(2)
             num      = TPIclean(2) * qaux*ccn_tuning / (dens_*fmassclean)     ! only sulfate
         end where

     case ('sulforg03')
         where (q > max_clean)
             sigma    = SIGI(3)
             diameter = DPGI(3)
             num      = TPI(3) * qaux*ccn_tuning / (dens_*fmassaux)           ! only sulfate mass
         elsewhere
             sigma    = SIGIclean(3)
             diameter = DPGIclean(3)
             num      = TPIclean(3) * qaux*ccn_tuning / (dens_*fmassclean)    ! only sulfate
         end where

     case ('bcphilic')
         sigma    = log(2.0)
         f_soot   = 1.0
         diameter = 0.0118*2e-6
         num = q / ((MAPL_PI/6.0) * densBC * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('ocphilic')
         sigma     = log(2.2)
         f_organic = 1.0
         diameter  = 0.0212*2.0e-6
         num = q / ((MAPL_PI/6.0) * densOrg * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case ('brcphilic')
         sigma     = log(2.2)
         f_organic = 1.0
         diameter  = 0.0212*2.0e-6
         num = q / ((MAPL_PI/6.0) * densOrg * diameter*diameter*diameter * exp(4.5*sigma*sigma))

     case default
         __raise__(UNKNOWN_AEROSOL_MODE,"Unknown aerosol mode used in the GOCART aerosol activation properties method: "//trim(mode))

     end select


     RETURN_(ESMF_SUCCESS)

    end subroutine aap_


    subroutine mamnet_emulator_(mode, state, aeroList, temperature, air_density, num, diameter, sigma, &
                                hygroscopicity, density, f_dust, f_soot, f_organic, &
                                i1, i2, j1, j2, km, rc)

      implicit none

      integer, intent(in) :: i1, i2
      integer, intent(in) :: j1, j2
      integer, intent(in) :: km

      character(len=*), intent(in)                  :: mode
      type(ESMF_State), intent(inout)               :: state
      character(len=ESMF_MAXSTR), intent(in)        :: aeroList(:)
      real, intent(in), dimension(i1:i2,j1:j2,km)   :: temperature
      real, intent(in), dimension(i1:i2,j1:j2,km)   :: air_density

      real, intent(out), dimension(i1:i2,j1:j2,km)  :: num
      real, intent(out), dimension(i1:i2,j1:j2,km)  :: diameter
      real, intent(out), dimension(i1:i2,j1:j2,km)  :: sigma
      real, intent(out), dimension(i1:i2,j1:j2,km)  :: hygroscopicity
      real, intent(out), dimension(i1:i2,j1:j2,km)  :: density
      real, intent(out), dimension(i1:i2,j1:j2,km)  :: f_dust
      real, intent(out), dimension(i1:i2,j1:j2,km)  :: f_soot
      real, intent(out), dimension(i1:i2,j1:j2,km)  :: f_organic

      integer, intent(out) :: rc

      type(ESMF_State)             :: child_state
      character(len=ESMF_MAXSTR)   :: mode_, aeroName, baseName, baseSuffix
      integer                      :: i, varNameLen
      logical                      :: added, added_any

      real, parameter :: sig_acc = log(1.8)
      real, parameter :: sig_ait = log(1.6)
      real, parameter :: sig_cdu = log(1.8)
      real, parameter :: sig_css = log(2.0)
      real, parameter :: sig_fdu = log(1.8)
      real, parameter :: sig_fss = log(2.0)
      real, parameter :: sig_pcm = log(1.6)

!     PySR emulator native units:
!       aerosol mass inputs  : log10(microgram aerosol / kg air)
!       aerosol mass outputs : log10(microgram aerosol / kg air)
!       aerosol number output: log10(number / mg air)
!     GEOS/MAPL native units in this callback are kg kg-1 for mass and
!     number m-3 for aerosol number concentration.
      real, allocatable, dimension(:,:,:) :: DUsum, SSsum, OGsum, BCsum, SUsum, T, AIRDENS
      real, allocatable, dimension(:,:,:) :: modal_mass, modal_number
      real, allocatable, dimension(:,:,:) :: mass1, mass2, mass3, mass4, mass5, mass6

      __Iam__('GOCART2G::aerosol_activation_properties::mamnet_emulator_')

      mode_ = ESMF_UtilStringLowerCase(trim(mode), __RC__)

      allocate(DUsum(i1:i2,j1:j2,km), SSsum(i1:i2,j1:j2,km), OGsum(i1:i2,j1:j2,km), &
               BCsum(i1:i2,j1:j2,km), SUsum(i1:i2,j1:j2,km), T(i1:i2,j1:j2,km), &
               AIRDENS(i1:i2,j1:j2,km), modal_mass(i1:i2,j1:j2,km), &
               modal_number(i1:i2,j1:j2,km), mass1(i1:i2,j1:j2,km), &
               mass2(i1:i2,j1:j2,km), mass3(i1:i2,j1:j2,km), mass4(i1:i2,j1:j2,km), &
               mass5(i1:i2,j1:j2,km), mass6(i1:i2,j1:j2,km), __STAT__)

      DUsum = 0.0
      SSsum = 0.0
      OGsum = 0.0
      BCsum = 0.0
      SUsum = 0.0

!     Sum native GOCART species into the MAMnet aggregate inputs.
      do i = 1, size(aeroList)
         aeroName = ESMF_UtilStringUpperCase(trim(aeroList(i)), __RC__)
         call ESMF_StateGet(state, trim(aeroList(i)), child_state, __RC__)

         if (aeroName(1:min(2,len_trim(aeroName))) == 'DU') then
            call add_field_sum_(child_state, 'DU', DUsum, added)

         else if (aeroName(1:min(2,len_trim(aeroName))) == 'SS') then
            call add_field_sum_(child_state, 'SS', SSsum, added)

         else if (aeroName(1:min(2,len_trim(aeroName))) == 'SU') then
            call add_field_sum_(child_state, 'SO4', SUsum, added)
            if (.not. added) call add_field_sum_(child_state, 'SU', SUsum, added)

         else if (aeroName(1:min(2,len_trim(aeroName))) == 'OG' .or. &
                  aeroName(1:min(2,len_trim(aeroName))) == 'OC' .or. &
                  aeroName(1:min(2,len_trim(aeroName))) == 'BR' .or. &
                  index(aeroName, 'CA.OC') > 0 .or. index(aeroName, 'CA.BR') > 0) then
            added_any = .false.
            baseName = trim(aeroList(i))
            varNameLen = len_trim(baseName)
            if (varNameLen > 5) then
               baseSuffix = ESMF_UtilStringUpperCase(baseName(varNameLen-4:varNameLen), __RC__)
               if (baseSuffix == '_AERO') then
                  varNameLen = varNameLen - 5
                  baseName = baseName(1:varNameLen)
               end if
            end if
            call add_field_sum_(child_state, trim(baseName)//'philic', OGsum, added)
            added_any = added_any .or. added
            call add_field_sum_(child_state, trim(baseName)//'phobic', OGsum, added)
            added_any = added_any .or. added
            if (.not. added_any) then
               call add_field_sum_(child_state, 'OG', OGsum, added)
               if (.not. added) call add_field_sum_(child_state, 'OC', OGsum, added)
            end if

         else if (aeroName(1:min(2,len_trim(aeroName))) == 'BC' .or. &
                  index(aeroName, 'CA.BC') > 0) then
            added_any = .false.
            baseName = trim(aeroList(i))
            varNameLen = len_trim(baseName)
            if (varNameLen > 5) then
               baseSuffix = ESMF_UtilStringUpperCase(baseName(varNameLen-4:varNameLen), __RC__)
               if (baseSuffix == '_AERO') then
                  varNameLen = varNameLen - 5
                  baseName = baseName(1:varNameLen)
               end if
            end if
            call add_field_sum_(child_state, trim(baseName)//'philic', BCsum, added)
            added_any = added_any .or. added
            call add_field_sum_(child_state, trim(baseName)//'phobic', BCsum, added)
            added_any = added_any .or. added
            if (.not. added_any) then
               call add_field_sum_(child_state, 'BC', BCsum, added)
            end if
         end if
      end do

!      PySR/Mamnet preprocessing:
!       1. Sum native GEOS mass fields [kg kg-1].
!       2. Convert each aggregate mass to log10(microgram kg-1), with
!          non-positive values set to MIN_MASS = -20.0.
!       3. Standardize all seven inputs as (x - INPUT_MU) / INPUT_SIGMA.
!     Temperature is standardized in K and AIRDENS in kg m-3; neither is
!     log-transformed.
      call mamnet_preprocess_mass_(SUsum, -1.95, 3.11)
      call mamnet_preprocess_mass_(SSsum, -2.32, 3.46)
      call mamnet_preprocess_mass_(OGsum, -3.01, 3.34)
      call mamnet_preprocess_mass_(BCsum, -4.12, 3.24)
      call mamnet_preprocess_mass_(DUsum, -4.72, 5.61)
      T       = (temperature - 243.8) / 27.35
      AIRDENS = (air_density - 0.39) / 0.44

      num           = 0.0
      diameter      = 1.0e-9
      sigma         = sig_acc
      density       = densSO4
      hygroscopicity= k_SO4
      f_dust        = 0.0
      f_soot        = 0.0
      f_organic     = 0.0
      modal_mass    = 0.0
      modal_number  = 0.0
      mass1         = 0.0
      mass2         = 0.0
      mass3         = 0.0
      mass4         = 0.0
      mass5         = 0.0
      mass6         = 0.0

      !Apply MAMnet/PYSR correlations
      select case (mode_)
        case ('num_a_acc')
        mass1 = mamnet_mass_from_logugkg_( &
             AIRDENS * 0.14940897 + OGsum * 0.13861534 + 2.96296735254423 * SUsum - (SUsum**2 * 0.010278153 - &
             0.9066853) * (-0.47732595218884 * SUsum**2 + (3.0 * SSsum - SUsum + SUsum * SUsum**2) * &
             (-0.12791413)) * (SUsum * SUsum**2 * (-0.0023168945) - 0.36351326) - 2.368718 )
        mass2 = mamnet_mass_from_logugkg_( &
             log(mamnet_safe_exp_(9.99418968269206e-8 * OGsum**9 * SUsum**10 + 3.4575536 * OGsum + (DUsum - T) * &
             mamnet_safe_exp_(-OGsum + SSsum - 1.46074765921952 * mamnet_safe_exp_(-BCsum + DUsum + OGsum + SUsum &
             - T))) + 7.1715164e-8) - 3.5188847 )
        mass3 = mamnet_mass_from_logugkg_( &
             (SSsum - (SSsum / (SSsum + (-SSsum - 6.9668455)**4 + 2.97739 + 7.2665043) - 0.44077685) - &
             mamnet_safe_exp_((SSsum + 0.34078407) / mamnet_safe_exp_(AIRDENS + 0.30479166) - mamnet_safe_exp_(-T &
             - 0.85498846) - mamnet_safe_exp_(-SSsum**2 + SUsum - 0.71152645) - 0.44651955) - 1.2695731) * &
             3.539773 )
        mass4 = mamnet_mass_from_logugkg_( &
             0.1340005 - (OGsum + (OGsum - 4.2014384) * (-3.1285706) - (mamnet_safe_exp_((OGsum + 1.5630813) * &
             (-0.37928277)) + 0.7314927) - mamnet_safe_exp_(2.437156 + 0.1856344 * (BCsum - &
             mamnet_safe_exp_((OGsum - 1.1092931) * (AIRDENS * 0.64947796 + mamnet_safe_exp_((OGsum + 2.5312295) &
             * (-0.330968))) * 0.093144685 * DUsum + 0.409221))) + 0.65003765) )
        mass5 = mamnet_mass_from_logugkg_( &
             2.0 * BCsum + 0.37393135 * DUsum - (-(AIRDENS**3 * (-0.042047553) + BCsum) + ((AIRDENS * &
             (-0.0174103851089025 * BCsum**3 + (DUsum - 0.3569122) / 0.86332375) / (-0.6214044) + AIRDENS + BCsum &
             * DUsum) * 0.022352552 * OGsum * (-0.24613461 * BCsum - 1.1892506) + 2.6788719)**2) + 2.8679202 )
        mass6 = mamnet_mass_from_logugkg_( &
             ((-SUsum + (0.73333204 - (OGsum + (BCsum**3 * (-0.034107774) + DUsum * SSsum + DUsum) * 0.42104927)) &
             * (AIRDENS + 0.6659424)) * (OGsum + 6.592014) * 0.014408116 + 1.0507722) * (OGsum**3 * DUsum**2 * &
             (-0.00072262995) + DUsum + (OGsum - 1.5162371) * 2.2371898) )
        modal_number = &
             log(0.159724471828337 * mamnet_safe_exp_(3.2600024 * SSsum) + mamnet_safe_exp_(-AIRDENS + 4.473384 * &
             SUsum - 1.1625581 * mamnet_safe_exp_(-AIRDENS * (AIRDENS + SSsum) + SSsum - SUsum - 0.25958222 * T)) &
             + mamnet_safe_exp_(-0.501604104847096 * BCsum * OGsum + 0.501604104847096 * DUsum + 3.6520624 * &
             OGsum) + 0.1504267) - 1.1134043

        modal_number = mamnet_number_from_logmg_(modal_number, air_density)
        modal_mass = mass1 + mass2 + mass3 + mass4 + mass5 + mass6
        sigma = sig_acc

        where ((modal_mass .gt. 0.0) .and. (modal_number .gt. 0.0))
           hygroscopicity = (k_SO4*mass1 + k_ORG*mass2 + k_SS*mass3 + k_ORG*mass4 + k_BC*mass5 + k_SO4*mass6) / &
           modal_mass
           density        = (densSO4*mass1 + densORG*mass2 + densSS*mass3 + densORG*mass4 + densBC*mass5 + &
           densSO4*mass6) / modal_mass
           diameter       = mamnet_diameter_m_(modal_mass, modal_number, density, air_density, sig_acc)
           f_organic      = (mass2 + mass4) / modal_mass
           f_soot         = mass5 / modal_mass
        end where
        diameter = min(max(diameter, 0.056e-6), 0.26e-6)

        case ('num_a_ait')
        mass1 = mamnet_mass_from_logugkg_( &
             SUsum * 1.435992 + (0.9068899 - mamnet_safe_exp_(-mamnet_safe_exp_(SUsum + 6.573142) - 0.39771122)) &
             * (-AIRDENS + DUsum * mamnet_safe_exp_((-OGsum - (SUsum + 1.1855656)**2 - &
             mamnet_safe_exp_(mamnet_safe_exp_(0.53339916 - AIRDENS)) + 0.5632794) * 0.20213664) + 2.0 * SUsum) - &
             2.6875386 )
        mass2 = mamnet_mass_from_logugkg_( &
             ((mamnet_safe_exp_((2.3166714 - DUsum) * 0.094840094 * (OGsum - mamnet_safe_exp_(-OGsum + T - &
             mamnet_safe_exp_(T * (OGsum + 1.1159782) * (-0.18184726) - (BCsum - OGsum + 1.7534208)**2 + &
             3.811902) + 0.50736135))) - 3.4294894)**3 + 9.272862) * 0.6449888 )
        mass3 = mamnet_mass_from_logugkg_( &
             -(SSsum * (-0.3775216 * SSsum + mamnet_safe_exp_((SSsum * (AIRDENS + (SSsum + 1.2909489) * &
             (-0.42036253)) - 1.3459847) * SSsum * (-0.31506303) + 0.29448342) - 4.450231) + 0.7113204 * (AIRDENS &
             * (-0.65032005) + T) * mamnet_safe_exp_(SSsum * (SSsum + T + 0.363342 + 1.1393063) * (-0.3971348))) &
             - 7.14481 )
        mass4 = mamnet_mass_from_logugkg_( &
             -AIRDENS + 0.227597092073677 * BCsum**2 + 3.6720307 * BCsum - (BCsum - 0.03374244) - 2.7662373 - 1 / &
             (mamnet_safe_exp_(DUsum * ((0.99282354 - OGsum) * mamnet_safe_exp_(DUsum) + (mamnet_safe_exp_(-OGsum &
             + SSsum) + mamnet_safe_exp_(-SUsum - 3.8900654)) * 0.31177384 + 0.36878952)) + 0.095603265) )
        modal_number = &
             mamnet_safe_exp_((1.4672223 - SUsum) * ((4.2974687 - DUsum) * (AIRDENS * (AIRDENS * 0.03419681 * &
             (BCsum + (+ 0.3641984 * AIRDENS - 0.8241336) / (0.0782392 * 3.6440687)) - 0.091314726) + SUsum / &
             4.2350483) - 0.08502343 * AIRDENS * AIRDENS * SSsum / AIRDENS) / 3.7598991 + 1.7792255) - 3.0110688

        modal_number = mamnet_number_from_logmg_(modal_number, air_density)
        modal_mass = mass1 + mass2 + mass3 + mass4
        sigma = sig_ait

        where ((modal_mass .gt. 0.0) .and. (modal_number .gt. 0.0))
           hygroscopicity = (k_SO4*mass1 + k_ORG*mass2 + k_SS*mass3 + k_SO4*mass4) / modal_mass
           density        = (densSO4*mass1 + densORG*mass2 + densSS*mass3 + densSO4*mass4) / modal_mass
           diameter       = mamnet_diameter_m_(modal_mass, modal_number, density, air_density, sig_ait)
           f_organic      = mass2 / modal_mass
        end where
        diameter = min(max(diameter, 0.015e-6), 0.052e-6)

        case ('num_a_cdu')
        mass1 = mamnet_mass_from_logugkg_( &
             DUsum * (-2.0153341) + (SUsum - 0.49038178 * (OGsum - 4.442369) + 1.3438548) * &
             mamnet_safe_exp_(DUsum * 0.78873074) - 20.168962 * mamnet_safe_exp_(mamnet_safe_exp_(DUsum) * &
             (-2.7150674)) + mamnet_safe_exp_(DUsum * DUsum / (-1.6112665) - mamnet_safe_exp_((+ 3.4170017 * &
             DUsum * DUsum + SSsum) * 0.6617511)) * (-1.9663713) - 8.657404 )
        mass2 = mamnet_mass_from_logugkg_( &
             (mamnet_safe_exp_(1.6663667 + mamnet_safe_exp_((DUsum - (+ 0.1757562 * 0.0236223521982394 * DUsum**6 &
             + 0.034614235 / mamnet_safe_exp_(3 * AIRDENS**3))) * (-0.54991704)) * (-0.38928843)) - 8.1776085) * &
             (DUsum * DUsum * (-0.21759257) - (DUsum + mamnet_safe_exp_(2.2243273 - DUsum) * (-0.0015917345)) + &
             1.2363818) + 0.94540584 )
        mass3 = mamnet_mass_from_logugkg_( &
             DUsum + (1.9855924 - DUsum) * (-4.0251822) + ((DUsum + mamnet_safe_exp_(-DUsum * (DUsum + 2.1608353) &
             - 0.88563025) + 0.30679715)**3 + mamnet_safe_exp_(-DUsum - 0.5742986)) * (-0.41017267) + &
             mamnet_safe_exp_((BCsum - ((AIRDENS * (SSsum - 0.8071881) - DUsum) * 0.23728389 - 0.39312544) - &
             1.6215085)**2 * (-1.0332288)) )
        modal_number = &
             DUsum**2 * mamnet_safe_exp_(0.81642973 - (AIRDENS + mamnet_safe_exp_(mamnet_safe_exp_(DUsum - (BCsum &
             + 0.54987484)**2 * mamnet_safe_exp_(2 * (-SSsum + SUsum) * (-DUsum - OGsum + 2.0 * SUsum - (T - &
             0.12808016) + 0.44125444))))) * 705.837329873139 * mamnet_safe_exp_(-12 * DUsum)) - 3.0000176

        modal_number = mamnet_number_from_logmg_(modal_number, air_density)
        modal_mass = mass1 + mass2 + mass3
        sigma = sig_cdu

        where ((modal_mass .gt. 0.0) .and. (modal_number .gt. 0.0))
           hygroscopicity = (k_SO4*mass1 + k_DU*mass2 + k_SO4*mass3) / modal_mass
           density        = (densSO4*mass1 + densDU*mass2 + densSO4*mass3) / modal_mass
           diameter       = mamnet_diameter_m_(modal_mass, modal_number, density, air_density, sig_cdu)
           f_dust         = mass2 / modal_mass
        end where
        diameter = min(max(diameter, 0.59e-6), 2.75e-6)

        case ('num_a_css')
        mass1 = mamnet_mass_from_logugkg_( &
             -0.0436979524995172 * SSsum**3 - (SSsum * (-3.3817623) - 0.30655035) + (mamnet_safe_exp_(-SSsum - &
             7.3036537) + 2.04504) * (-4.1291785) + mamnet_safe_exp_((0.236116749768772 * AIRDENS**2 * (AIRDENS + &
             0.00398151749312282 * SSsum**6 + SSsum + SUsum**2)**2) * (-0.06323494) + mamnet_safe_exp_(AIRDENS * &
             0.17233828 + SUsum * 0.15341063)) )
        mass2 = mamnet_mass_from_logugkg_( &
             -0.0369266690709064 * SSsum**3 + SSsum * 4.6673074 - 3.247442 * mamnet_safe_exp_((SSsum * (AIRDENS - &
             SSsum) * (-0.2504734) - T - (AIRDENS * SSsum * (-0.38934025) + 1.521244)**2 * (+ 1 * &
             0.477358704437942 * SSsum**3 + SSsum + T * (-0.77617097)) + 1.5156318) * (-0.10526849)) )
        mass3 = mamnet_mass_from_logugkg_( &
             DUsum + 2.0 * SSsum + 0.13189164992721 * (+ 0.5101209 * DUsum + SSsum - 0.5286089)**2 + &
             mamnet_safe_exp_(((OGsum + (-DUsum + SSsum) * (OGsum - SSsum)) * 0.10700623 - 0.37593406)**2 * &
             (BCsum * DUsum * OGsum * 0.1375481 - 2.2478037) + 2.2914555) - 13.653886 )
        modal_number = &
             (mamnet_safe_exp_(SSsum - (1.2254843 - SSsum)**6) - 0.06454122) * mamnet_safe_exp_((SUsum - (SSsum + &
             (AIRDENS**2 + SUsum) * 0.30078012 - 1.7685149 + (-T + mamnet_safe_exp_(SUsum + 0.21698205 * (AIRDENS &
             + SUsum))) * 0.058395915 * (-SSsum - T + 2.2507513) / SSsum)**2) * 0.16670467) - 2.9979303

        modal_number = mamnet_number_from_logmg_(modal_number, air_density)
        modal_mass = mass1 + mass2 + mass3
        sigma = sig_css

        where ((modal_mass .gt. 0.0) .and. (modal_number .gt. 0.0))
           hygroscopicity = (k_SO4*mass1 + k_SS*mass2 + k_SO4*mass3) / modal_mass
           density        = (densSO4*mass1 + densSS*mass2 + densSO4*mass3) / modal_mass
           diameter       = mamnet_diameter_m_(modal_mass, modal_number, density, air_density, sig_css)
        end where
        diameter = min(max(diameter, 0.63e-6), 3.7e-6)

        case ('num_a_fdu')
        mass1 = mamnet_mass_from_logugkg_( &
             2.0 * DUsum + SUsum - ((OGsum - SUsum + (DUsum - mamnet_safe_exp_(-AIRDENS + SSsum * &
             (-0.23903924)))**2) * mamnet_safe_exp_(SUsum) / 5.2433176 + ((0.00797231176144782 * DUsum**2 * &
             SUsum**2) * (-2.5036454) - (2.0 * DUsum + SUsum - 2.8992734) + mamnet_safe_exp_(BCsum))**2 * &
             0.04901954) - 4.4872775 )
        mass2 = mamnet_mass_from_logugkg_( &
             DUsum * (11.759894 - mamnet_safe_exp_((-(DUsum + 2.829828) - 0.058263373) * mamnet_safe_exp_((DUsum &
             + 3.7056472) / 1.605885))) + 0.46684366 * (-AIRDENS * mamnet_safe_exp_((BCsum + &
             mamnet_safe_exp_(1.2944329 - 4.168569 * (DUsum + 2.4785476)) + 0.15922968) * (-0.4279211)) + DUsum * &
             ((DUsum - 0.46338895) * (-0.14047985) - 12.628727)) - 5.4927025 )
        mass3 = mamnet_mass_from_logugkg_( &
             mamnet_safe_exp_(((DUsum + 0.2967915 * OGsum * mamnet_safe_exp_((T - 1.3718457) * (AIRDENS * &
             (-0.34895086) + SUsum * (-0.4814743))) + ((-0.212636957940678 * DUsum**3) * (-AIRDENS + OGsum - &
             mamnet_safe_exp_((DUsum + SUsum + 0.37202525) * 1.6236217)) + SUsum) * 0.13209933) / 5.135033 - &
             0.6809525)**3 + 2.8903885) - 20.02814 )
        modal_number = &
             mamnet_safe_exp_(1.3534774 + mamnet_safe_exp_(1.2618353 - mamnet_safe_exp_((DUsum + (OGsum + T + &
             (SUsum - 0.2181561) * (-3.0058098)) * 0.11055958 / (((BCsum + T)**2 - (AIRDENS + SUsum - &
             0.27325228)) * (-0.25686136) - 1.0907141)**3) * 2.0652635)) * (-5.42722)) * mamnet_safe_exp_(DUsum * &
             1.5242821 - 1.7857378) - 2.9989681

        modal_number = mamnet_number_from_logmg_(modal_number, air_density)
        modal_mass = mass1 + mass2 + mass3
        sigma = sig_fdu

        where ((modal_mass .gt. 0.0) .and. (modal_number .gt. 0.0))
           hygroscopicity = (k_SO4*mass1 + k_DU*mass2 + k_SO4*mass3) / modal_mass
           density        = (densSO4*mass1 + densDU*mass2 + densSO4*mass3) / modal_mass
           diameter       = mamnet_diameter_m_(modal_mass, modal_number, density, air_density, sig_fdu)
           f_dust         = mass2 / modal_mass
        end where
        diameter = min(max(diameter, 0.14e-6), 0.62e-6)

        case ('num_a_fss')
        mass1 = mamnet_mass_from_logugkg_( &
             (1.4358984 - 0.04234827 * SSsum * SSsum) * (-BCsum + SUsum + (BCsum * (DUsum + SSsum + SUsum * &
             (-0.5049907) - mamnet_safe_exp_(AIRDENS * 0.9217091) - 1.0473211) + (1.78869626866225 * SSsum**2) * &
             mamnet_safe_exp_(SUsum) + T) * (-0.1326069) + 1.788232) + (-SSsum - 0.8101632) * (-2.524106) - &
             7.904183 )
        mass2 = mamnet_mass_from_logugkg_( &
             3.0 * SSsum + 0.00637642878327194 * (SSsum**2 + (-(-DUsum + mamnet_safe_exp_(1.2761903 + &
             mamnet_safe_exp_(-DUsum + T + (T - (AIRDENS * SUsum - 0.1587111)) * 3.585926 / (-0.77560794)) * &
             (-0.7753292))) + (-0.024844296634896 * SSsum**3 + SSsum)**2) * (-0.747859) - 18.497986)**2 - &
             4.9092274 )
        mass3 = mamnet_mass_from_logugkg_( &
             -AIRDENS + 2.0 * DUsum + SSsum + (OGsum + 1 / (2.4145953956997e-8 * OGsum**9 - 0.28634813)) * &
             2.1531038 + (SSsum - (-DUsum * SSsum + (-AIRDENS + (0.5193177 - (-SSsum * SSsum + &
             mamnet_safe_exp_(DUsum) * (-4.0822825))) * 0.3172893) * mamnet_safe_exp_(OGsum))) * 0.39999217 + &
             3.330142 )
        modal_number = &
             -(BCsum * 0.013463435 + mamnet_safe_exp_((0.21528454 - mamnet_safe_exp_(SSsum - 0.25585648 * (-SUsum &
             + T - 0.7192855 * (AIRDENS - 1.0624889)) * mamnet_safe_exp_((-AIRDENS - 0.97842515) * ((BCsum + T) * &
             (SUsum + 0.7201301) - 0.22462068) * 0.20248717)))**3) * 1.9161335) + mamnet_safe_exp_((SSsum - &
             0.7816974) * 2.1192024) - 1.1324198

        modal_number = mamnet_number_from_logmg_(modal_number, air_density)
        modal_mass = mass1 + mass2 + mass3
        sigma = sig_fss

        where ((modal_mass .gt. 0.0) .and. (modal_number .gt. 0.0))
           hygroscopicity = (k_SO4*mass1 + k_SS*mass2 + k_SO4*mass3) / modal_mass
           density        = (densSO4*mass1 + densSS*mass2 + densSO4*mass3) / modal_mass
           diameter       = mamnet_diameter_m_(modal_mass, modal_number, density, air_density, sig_fss)
        end where
        diameter = min(max(diameter, 0.095e-6), 0.56e-6)

        case ('num_a_pcm')
        mass1 = mamnet_mass_from_logugkg_( &
             -BCsum * SUsum + BCsum + (1.3283285 - BCsum) * (-0.084837966 * (3.6884565 - BCsum) * (SUsum + &
             mamnet_safe_exp_((4.2200356 - SSsum) * (DUsum + SUsum + 1.1826622) * 0.11968093)) + &
             mamnet_safe_exp_( T + (-AIRDENS + OGsum)**3) * (-0.53216743) - 1.9602302) - 0.4937047 )
        mass2 = mamnet_mass_from_logugkg_( &
             BCsum * 3.067212 - 0.3232528 * DUsum + (BCsum + 3.9767632)**2 * (mamnet_safe_exp_(OGsum * &
             0.21598318) - 1.1018325)**2 + mamnet_safe_exp_((BCsum + SSsum + (-BCsum + T * (T + 1.3359089) * &
             (OGsum - SSsum + 0.9639113))**2 / (-3.801517)) * 0.41341227) - 6.0948753 - 0.10709444 )
        modal_number = &
             0.27706329088824 * ((T + 0.36539873) * (AIRDENS - T - 0.8513434) - 0.9586623) * &
             mamnet_safe_exp_(OGsum) - 0.27706329088824 * mamnet_safe_exp_(0.8977992 * DUsum + SUsum) + &
             0.27706329088824 * mamnet_safe_exp_(OGsum**3 + (SSsum - 0.07689504)**3) + 0.826208660649664 * &
             log(mamnet_safe_exp_(4.9132433 * BCsum) + 0.10393974) - 1.11439971011352

        modal_number = mamnet_number_from_logmg_(modal_number, air_density)
        modal_mass = mass1 + mass2
        sigma = sig_pcm

        where ((modal_mass .gt. 0.0) .and. (modal_number .gt. 0.0))
           hygroscopicity = (k_ORG*mass1 + k_BC*mass2) / modal_mass
           density        = (densORG*mass1 + densBC*mass2) / modal_mass
           diameter       = mamnet_diameter_m_(modal_mass, modal_number, density, air_density, sig_pcm)
           f_organic      = mass1 / modal_mass
           f_soot         = mass2 / modal_mass
        end where
        diameter = min(max(diameter, 0.039e-6), 0.13e-6)

        case default
           __raise__(UNKNOWN_AEROSOL_MODE, "Unknown MAMnet emulator aerosol mode in aerosol activation properties: ")
      end select

      num = modal_number

!     Final numerical guards for model stability and MAMnet physical ranges.
      density        = min(max(density,        900.0), 2200.0)
      hygroscopicity = min(max(hygroscopicity, mam_min_kappa), mam_max_kappa)
      diameter       = min(max(diameter, mam_min_diameter_m), mam_max_diameter_m)
      f_dust         = min(max(f_dust,          0.0),    1.0)
      f_soot         = min(max(f_soot,          0.0),    1.0)
      f_organic      = min(max(f_organic,       0.0),    1.0)
      num            = min(max(num, 0.0), mam_max_number_m3)

      deallocate(DUsum, SSsum, OGsum, BCsum, SUsum, T, AIRDENS, modal_mass, modal_number, &
                 mass1, mass2, mass3, mass4, mass5, mass6, __STAT__)

      RETURN_(ESMF_SUCCESS)

    end subroutine mamnet_emulator_


      subroutine mamnet_preprocess_mass_(x, input_mu, input_sigma)
        real, intent(inout), dimension(:,:,:) :: x
        real, intent(in)                      :: input_mu
        real, intent(in)                      :: input_sigma

        where (x > 0.0)
           x = log10(x * mam_mass_kgkg_to_ugkg)
        elsewhere
           x = -20.0
        end where
        x = (x - input_mu) / input_sigma
      end subroutine mamnet_preprocess_mass_

      elemental real function mamnet_safe_exp_(x)
        real, intent(in) :: x

        if (x /= x) then
           mamnet_safe_exp_ = 0.0
        else if (x > 20.0) then
           mamnet_safe_exp_ = exp(20.0)
        else if (x < -20.0) then
           mamnet_safe_exp_ = 0.0
        else
           mamnet_safe_exp_ = exp(x)
        end if
      end function mamnet_safe_exp_

      elemental real function mamnet_mass_from_logugkg_(log_mass_ugkg)
        real, intent(in) :: log_mass_ugkg
        real             :: log_max_mass_ugkg

        log_max_mass_ugkg = log10(mam_max_output_mass_kgkg / mam_mass_ugkg_to_kgkg)

        if (log_mass_ugkg /= log_mass_ugkg) then
           mamnet_mass_from_logugkg_ = 0.0
        else if (log_mass_ugkg <= -30.0) then
           mamnet_mass_from_logugkg_ = 0.0
        else if (log_mass_ugkg >= log_max_mass_ugkg) then
           mamnet_mass_from_logugkg_ = mam_max_output_mass_kgkg
        else
           mamnet_mass_from_logugkg_ = mam_mass_ugkg_to_kgkg * 10.0**log_mass_ugkg
        end if
      end function mamnet_mass_from_logugkg_

      elemental real function mamnet_number_from_logmg_(log_number_mg, rho_air_kgm3)
        real, intent(in) :: log_number_mg
        real, intent(in) :: rho_air_kgm3
        real             :: log_number_m3

        if (log_number_mg /= log_number_mg .or. rho_air_kgm3 <= 0.0) then
           mamnet_number_from_logmg_ = 0.0
        else
           log_number_m3 = log_number_mg + &
                            log10(max(rho_air_kgm3, mam_min_air_density_kgm3) * mam_air_kgm3_to_mgm3)
           if (log_number_m3 >= log10(mam_max_number_m3)) then
              mamnet_number_from_logmg_ = mam_max_number_m3
           else if (log_number_m3 <= -30.0) then
              mamnet_number_from_logmg_ = 0.0
           else
              mamnet_number_from_logmg_ = 10.0**log_number_m3
           end if
        end if
      end function mamnet_number_from_logmg_

      elemental real function mamnet_diameter_m_(mass_kgkg, number_m3, rho_particle_kgm3, &
                                                 rho_air_kgm3, sigma_g)
        real, intent(in) :: mass_kgkg
        real, intent(in) :: number_m3
        real, intent(in) :: rho_particle_kgm3
        real, intent(in) :: rho_air_kgm3
        real, intent(in) :: sigma_g
        real             :: mass_conc_kgm3, width2

        if (mass_kgkg > 0.0 .and. number_m3 > 0.0 .and. &
            rho_particle_kgm3 > 0.0 .and. rho_air_kgm3 > 0.0) then
           mass_conc_kgm3 = mass_kgkg * rho_air_kgm3
           width2 = log(max(sigma_g, 1.0001)) * log(max(sigma_g, 1.0001))
           mamnet_diameter_m_ = ((6.0 * mass_conc_kgm3) / &
                                 (MAPL_PI * rho_particle_kgm3 * number_m3 * &
                                  mamnet_safe_exp_(4.5 * width2)))**(1.0/3.0)
        else
           mamnet_diameter_m_ = mam_min_diameter_m
        end if

        mamnet_diameter_m_ = min(max(mamnet_diameter_m_, mam_min_diameter_m), mam_max_diameter_m)
      end function mamnet_diameter_m_

      subroutine add_field_sum_(child_state_, field_name, accumulator, added)
        type(ESMF_State), intent(inout)                  :: child_state_
        character(len=*), intent(in)                     :: field_name
        real, intent(inout), dimension(:,:,:)            :: accumulator
        logical, intent(out)                             :: added

        real, pointer, dimension(:,:,:)                  :: ptr3d
        real, pointer, dimension(:,:,:,:)                :: ptr4d
        type(ESMF_Field)                                 :: field
        character(len=ESMF_MAXSTR), allocatable          :: itemNames(:)
        type(ESMF_StateItem_Flag), allocatable           :: itemTypes(:)
        integer                                          :: status_, nitems, item, rank, bin
        logical                                          :: found

        added = .false.
        found = .false.
        nullify(ptr3d)
        nullify(ptr4d)

!       Do not probe with MAPL_GetPointer for optional names/ranks.  MAPL
!       emits GetPointer.H failures even when RC is caught locally.  Check
!       the child state's item list and field rank first, then call the
!       matching GetPointer only for fields that actually exist.
        status_ = ESMF_SUCCESS
        call ESMF_StateGet(child_state_, itemCount=nitems, RC=status_)
        if (status_ /= ESMF_SUCCESS .or. nitems <= 0) return

        allocate(itemNames(nitems), itemTypes(nitems), stat=status_)
        if (status_ /= 0) return

        status_ = ESMF_SUCCESS
        call ESMF_StateGet(child_state_, itemNameList=itemNames, itemTypeList=itemTypes, RC=status_)
        if (status_ /= ESMF_SUCCESS) then
           deallocate(itemNames, itemTypes)
           return
        end if

        do item = 1, nitems
           if (itemTypes(item) == ESMF_StateItem_Field .and. trim(itemNames(item)) == trim(field_name)) then
              found = .true.
              exit
           end if
        end do
        deallocate(itemNames, itemTypes)

        if (.not. found) return

        status_ = ESMF_SUCCESS
        call ESMF_StateGet(child_state_, trim(field_name), field, RC=status_)
        if (status_ /= ESMF_SUCCESS) return

        status_ = ESMF_SUCCESS
        call ESMF_FieldGet(field, rank=rank, RC=status_)
        if (status_ /= ESMF_SUCCESS) return

        if (rank == 3) then
           status_ = ESMF_SUCCESS
           call MAPL_GetPointer(child_state_, ptr3d, trim(field_name), RC=status_)
           if (status_ == ESMF_SUCCESS .and. associated(ptr3d)) then
              accumulator = accumulator + ptr3d
              added = .true.
           end if

        else if (rank == 4) then
           status_ = ESMF_SUCCESS
           call MAPL_GetPointer(child_state_, ptr4d, trim(field_name), RC=status_)
           if (status_ == ESMF_SUCCESS .and. associated(ptr4d)) then
              do bin = 1, ubound(ptr4d, 4)
                 accumulator = accumulator + ptr4d(:,:,:,bin)
              end do
              added = .true.
           end if
        end if
      end subroutine add_field_sum_

  end subroutine aerosol_activation_properties

end module GOCART2G_AerosolActivationMod
