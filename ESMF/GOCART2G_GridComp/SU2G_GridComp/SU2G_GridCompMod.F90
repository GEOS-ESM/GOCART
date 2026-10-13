#include "MAPL.h"

!MODULE: SU2G_GridCompMod - GOCART refactoring of the SU gridded component

!INTERFACE:
module SU2G_GridCompMod

   !USES:
   use ESMF
   use MAPL
   use GOCART2G_MieMod
   use Chem_AeroGeneric
   use ReplenishAlarm
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process ! GOCART2G process library
   use GA_EnvironmentMod
   !$ use omp_lib

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData = 2

   integer, parameter :: DP = kind(1.0d0)

   ! gram molecular weights of species
   real, parameter :: fMassSulfur = 32., fMassSO2 = 64., fMassSO4 = 96.
   real, parameter :: fMassDMS = 62., fMassMSA = 96.

   real, parameter :: cpd = 1004.16

   ! relative position of sulfate tracers
   integer, parameter :: nDMS = 1, nSO2 = 2, nSO4 = 3, nMSA = 4
   real, parameter :: OCEAN = 0.0, LAND = 1.0, SEA_ICE = 2.0

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices

   !DESCRIPTION: This module implements GOCART's Sulfer (SU) Gridded Component.

   !REVISION HISTORY:
   ! 04January2024  Collow - Update to ChemSettling Call
   ! 08July2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

   ! Sulfer state
   type :: ThreadWorkspace
      integer :: nymd_last = -1 ! Previous nymd. Updated daily

      ! Degassing volcanoes
      integer :: nVolc = 0
      real, allocatable, dimension(:) :: vLat, vLon, vSO2, vElev, vCloud
      integer, allocatable, dimension(:) :: vStart, vEnd
      ! Explosive volcanoes
      integer :: nVolcE = 0
      real, allocatable, dimension(:) :: vLatE, vLonE, vSO2E, vElevE, vCloudE
      integer, allocatable, dimension(:) :: vStartE, vEndE

      ! Other point emissions of Sulfate (SO4)
      integer :: nPts = -1
      integer, allocatable, dimension(:) :: pstart, pend
      real, allocatable, dimension(:) :: pLat, pLon, pBase, pTop, pEmis
      logical :: firstRun = .true.
      integer :: nymd_oxidants = -1 ! Update the oxidant files?
      logical :: recycle_H2O2 = .false.
   end type ThreadWorkspace

   type, extends(GA_Environment) :: SU2G_GridComp
      integer :: myDOW = -1 ! my Day of the week: Sun=1, Mon=2,...,Sat=7
      logical :: diurnal_bb ! diurnal biomass burning
      real :: eAircraftFuel ! Aircraft emission factor: go from kg fuel to kg SO2
      real, allocatable :: aviation_layers(:) ! heights of the LTO, CDS and CRS layers
      real :: fSO4anth ! Fraction of anthropogenic emissions that are SO4
      !logical :: firstRun = .true.

      ! PRC: logic for GMI coupling
      logical :: using_GMI
      logical :: disable_emissions

      ! Special handling for volcanic emissions
      character(len=:), allocatable :: volcano_srcfilen_degassing
      character(len=:), allocatable :: volcano_srcfilen_explosive

      ! Workspace for point emissions
      logical :: doing_point_emissions = .false.
      character(len=:), allocatable :: point_emissions_srcfilen ! filename for pointwise emissions
      type(ThreadWorkspace), allocatable :: workspaces(:)
      type(ESMF_Alarm) :: alarm
   end type SU2G_GridComp

   type wrap_
      type(SU2G_GridComp), pointer :: PTR !=> null()
   end type wrap_

contains

   !IROUTINE: SetServices

   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, intent(out) :: rc ! return code

      !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
      ! the Initialize and Finalize services to generic versions. It also
      ! allocates our instance of a generic state and puts it in the
      ! gridded component (GC). Here we only set the two-stage run method
      ! and declare the data services.

      !REVISION HISTORY:
      ! 08July2020   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

      ! Locals
      character(len=ESMF_MAXSTR) :: COMP_NAME
      type(wrap_) :: wrap
      type(SU2G_GridComp), pointer :: self

      character(len=ESMF_MAXSTR) :: field_name

      integer :: i
      real :: DEFVAL
      logical :: data_driven = .true.
      integer :: num_threads

      integer :: status
      character(len=255) :: Iam = 'SetServices'

      ! Begin...

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(gc, NAME=COMP_NAME, _RC)
      Iam = trim(COMP_NAME) // '::' // Iam

      ! Wrap internal state for storing in GC
      allocate(self, _STAT)
      wrap%PTR => self

      num_threads = MAPL_get_num_threads()
      allocate(self%workspaces(0:num_threads - 1), _STAT)

      ! process generic config items
      call self%GA_Environment%load_from_config(gc, _RC)

      ! process SU-specific items
      call MAPL_GridCompGetResource(gc, "using_GMI", self%using_GMI, _RC)
      call MAPL_GridCompGetResource(gc, "disable_emissions", self%disable_emissions, _RC)
      call MAPL_GridCompGetResource(gc, "volcano_srcfilen_degassing", self%volcano_srcfilen_degassing, _RC)
      call MAPL_GridCompGetResource(gc, "volcano_srcfilen_explosive", self%volcano_srcfilen_explosive, _RC)
      call MAPL_GridCompGetResource(gc, "aircraft_fuel_emission_factor", self%eAircraftFuel, _RC)
      call MAPL_GridCompGetResource(gc, "so4_anthropogenic_fraction", self%fSO4anth, _RC)
      call MAPL_GridCompGetResource(gc, "aviation_vertical_layers", self%aviation_layers, _RC)
      _ASSERT(size(self%aviation_layers) == 4, "aviation_vertical_layers must have 4 values")

      call MAPL_GridCompGetResource(gc, "point_emissions_srcfilen", self%point_emissions_srcfilen, default='/dev/null', _RC)
      if ((index(self%point_emissions_srcfilen, '/dev/null') > 0)) then
         self%doing_point_emissions = .false. ! disable it if no file specified
      else
         self%doing_point_emissions = .true. ! we are good to go
      end if

      ! Is SU data driven?
      call determine_data_driven(COMP_NAME, data_driven, _RC)

      ! Set entry points
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)
      if (data_driven .neqv. .true.) then
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run2, _RC)
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run0, _RC)
      end if

      DEFVAL = 0.0

      ! Import and Internal states if data instance
      if (data_driven) then

         ! Pressure at layer edges
         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='PLE', &
              LONG_NAME='air_pressure', &
              UNITS='Pa', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationEdge, &
              RESTART=MAPL_RestartSkip, _RC)

         ! RH: is between 0 and 1
         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='RH2', &
              LONG_NAME='Rel_Hum_after_moist', &
              UNITS='1', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartSkip, _RC)

         call MAPL_AddInternalSpec(gc,&
              SHORT_NAME='SO4', &
              LONG_NAME='Sulphate aerosol', &
              UNITS='kg kg-1', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartOptional, &
         ! friendlyto='DYNAMICS:TURBULENCE:MOIST', &
              add2export=.true., _RC)

         call MAPL_AddImportSpec(gc,&
              SHORT_NAME='climSO4', &
              LONG_NAME='Sulphate aerosol', &
              UNITS='kg kg-1', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartOptional, _RC)

         do i = 1, self%nbins
            write(field_name, '(A, I0.3)') '', i
            ! dry deposition
            call MAPL_AddImportSpec(gc, &
                 SHORT_NAME='climSUDP' // trim(field_name), &
                 LONG_NAME='Sulfate dry deposition (bin ' // trim(field_name) // ')', &
                 UNITS='kg kg-1 s-1', &
                 DIMS=MAPL_DimsHorzOnly, &
                 VLOCATION=MAPL_VLocationCenter, &
                 RESTART=MAPL_RestartSkip, _RC)

            ! wet deposition
            call MAPL_AddImportSpec(gc, &
                 SHORT_NAME='climSUWT' // trim(field_name), &
                 LONG_NAME='Sulfate wet deposition (bin ' // trim(field_name) // ')', &
                 UNITS='kg kg-1 s-1', &
                 DIMS=MAPL_DimsHorzOnly, &
                 VLOCATION=MAPL_VLocationCenter, &
                 RESTART=MAPL_RestartSkip, _RC)

            ! gravitational settling
            call MAPL_AddImportSpec(gc, &
                 SHORT_NAME='climSUSD' // trim(field_name), &
                 LONG_NAME='Sulfate settling (bin ' // trim(field_name) // ')', &
                 UNITS='kg kg-1 s-1', &
                 DIMS=MAPL_DimsHorzOnly, &
                 VLOCATION=MAPL_VLocationCenter, &
                 RESTART=MAPL_RestartSkip, _RC)

            ! convective scavenging
            call MAPL_AddImportSpec(gc, &
                 SHORT_NAME='climSUSV' // trim(field_name), &
                 LONG_NAME='Sulfate convective scavenging (bin ' // trim(field_name) // ')', &
                 UNITS='kg kg-1 s-1', &
                 DIMS=MAPL_DimsHorzOnly, &
                 VLOCATION=MAPL_VLocationCenter, &
                 RESTART=MAPL_RestartSkip, _RC)
         end do
      end if ! (data_driven)

      if (.not. data_driven) then
         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='SU_H2O2', &
              LONG_NAME='source species', &
              UNITS='1', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartSkip, _RC)

         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='SU_OH', &
              LONG_NAME='source species', &
              UNITS='1', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartSkip, _RC)

         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='SU_NO3', &
              LONG_NAME='source species', &
              UNITS='1', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartSkip, _RC)
      end if

      if (self%using_GMI) then

         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='GMI_OH', &
              LONG_NAME='Hydroxyl_radical', &
              UNITS='mol/mol', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartSkip, _RC)

         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='GMI_H2O2', &
              LONG_NAME='Hydrogen_peroxide', &
              UNITS='mol/mol', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartSkip, _RC)

         call MAPL_AddImportSpec(gc, &
              SHORT_NAME='GMI_NO3', &
              LONG_NAME='Nitrogen_trioxide', &
              UNITS='mol/mol', &
              DIMS=MAPL_DimsHorzVert, &
              VLOCATION=MAPL_VLocationCenter, &
              RESTART=MAPL_RestartSkip, _RC)
      end if

      ! Import, Export, Internal states for computational instance
      if (.not. data_driven) then
#include "SU2G_Export___.h"
#include "SU2G_Import___.h"
#include "SU2G_Internal___.h"
         if (mapl_am_i_root()) then
            write(*, *) trim(Iam) // ": Wet removal scheme is " // trim(self%wet_removal_scheme)
            write(*, *) trim(Iam) // ": Settling scheme is " // trim(self%settling_scheme)
         end if
      end if

      ! This state holds fields needed by radiation
      call MAPL_AddExportSpec(gc, &
           SHORT_NAME=trim(COMP_NAME) // '_AERO', &
           LONG_NAME='aerosols_from_' // trim(COMP_NAME), &
           UNITS='kg kg-1', &
           DIMS=MAPL_DimsHorzVert, &
           VLOCATION=MAPL_VLocationCenter, &
           DATATYPE=MAPL_StateItem, _RC)

      ! This bundle is needed by surface for snow albedo modification
      ! by aerosol settling and deposition
      ! DEVELOPMENT NOTE - Change to StateItem in future
      call MAPL_AddExportSpec(gc, &
           SHORT_NAME=trim(COMP_NAME) // '_AERO_DP', &
           LONG_NAME='aerosol_deposition_from_' // trim(COMP_NAME), &
           UNITS='kg m-2 s-1', &
           DIMS=MAPL_DimsHorzOnly, &
           DATATYPE=MAPL_BundleItem, _RC)

      ! Store internal state in GC
      call ESMF_UserCompSetInternalState(gc, 'SU2G_GridComp', wrap, _RC)

      ! Set generic services
      call MAPL_GenericSetServices(gc, _RC)

      _RETURN(ESMF_SUCCESS)

   end subroutine SetServices

   !IROUTINE: Initialize

   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code

      !DESCRIPTION: This initializes SU Grid Component.

      !REVISION HISTORY:
      ! 08July2019   E.Sherman  First attempt at refactoring

      ! Locals
      character(len=ESMF_MAXSTR) :: COMP_NAME
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Grid) :: grid
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(ESMF_State) :: providerState
      type(ESMF_Config) :: cfg, universal_cfg
      type(ESMF_FieldBundle) :: Bundle_DP
      type(wrap_) :: wrap
      type(SU2G_GridComp), pointer :: self

      integer, allocatable :: mieTable_pointer(:)
      integer :: i, DIMS(3), km
      integer :: instance
      type(ESMF_Field) :: field, fld
      character(len=ESMF_MAXSTR) :: prefix, diurnal_bb
      real, pointer, dimension(:, :) :: lats
      real, pointer, dimension(:, :) :: lons
      real :: CDT ! chemistry timestep (secs)
      integer :: HDT ! model     timestep (secs)
      logical :: data_driven
      integer :: NUM_BANDS
      logical :: bands_are_present

      type(ESMF_Calendar) :: calendar
      type(ESMF_Time) :: currentTime
      type(ESMF_Time) :: ringTime
      type(ESMF_TimeInterval) :: ringInterval
      integer :: year, month, day, hh, mm, ss

      real, dimension(4) :: Vect_Hcts
      integer, allocatable, dimension(:) :: channels_
      integer :: nmom_
      character(len=ESMF_MAXSTR) :: file_
      logical :: file_exists
      integer :: status
      character(len=255) :: Iam = 'Initialize'

      ! Begin...

      ! Get the target components name and set-up traceback handle.
      call ESMF_GridCompGet(gc, grid=grid, NAME=COMP_NAME, config=universal_cfg, _RC)
      Iam = trim(COMP_NAME) // '::' // trim(Iam)

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      ! Get my internal private state
      call ESMF_UserCompGetInternalState(gc, 'SU2G_GridComp', wrap, status)
      _VERIFY(status)
      self => wrap%PTR

      ! Get dimensions
      call MAPL_GridGet(grid, localCellCountPerDim=DIMS, _RC)
      km = DIMS(3)
      self%km = km

      ! Get DTs
      call MAPL_GetResource(MAPL, HDT, label='RUN_DT:', _RC)
      call MAPL_GetResource(MAPL, CDT, label='GOCART2G_DT:', default=real(HDT), _RC)
      self%CDT = CDT

      ! Check whether to de-activate diurnal biomass burning (default is *on*)
      call ESMF_ConfigGetAttribute(universal_cfg, diurnal_bb, label='DIURNAL_BIOMASS_BURNING:', &
           default='YES', _RC)
      diurnal_bb = ESMF_UtilStringUpperCase(diurnal_bb, _RC)
      if (trim(diurnal_bb) == 'YES') then
         self%diurnal_bb = .true.
      else
         self%diurnal_bb = .false.
      end if

      ! Load resource file and get number of bins
      cfg = ESMF_ConfigCreate(_RC)
      inquire(file='SU2G_instance_' // trim(COMP_NAME) // '.rc', exist=file_exists)
      if (file_exists) then
         call ESMF_ConfigLoadFile(cfg, 'SU2G_instance_' // trim(COMP_NAME) // '.rc', _RC)
      else
         if (mapl_am_i_root()) print*, 'SU2G_instance_' // trim(COMP_NAME) // '.rc does not exist! loading' // &
              ' SU2G_instance_SU.rc instead'
         call ESMF_ConfigLoadFile(cfg, 'SU2G_instance_SU.rc', _RC)
      end if

      ! Call Generic Initialize
      call MAPL_GenericInitialize(gc, import, export, clock, _RC)

      ! Get parameters from generic state.
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, &
           lons=lons, &
           lats=lats, _RC)

      ! Is SU data driven?
      call determine_data_driven(COMP_NAME, data_driven, _RC)

      ! If this is a data component, the data is provided in the import
      ! state via ExtData instead of the actual GOCART children
      if (data_driven) then
         providerState = import
         prefix = 'clim'
      else
         providerState = export
         prefix = ''
      end if

      ! Add attribute information
      if (.not. data_driven) then
         call ESMF_StateGet(internal, 'DMS', field, _RC)
         call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), _RC)
         ! Vect_Hcts=(/0.56, 3500.0, 0.0, 0.0/)
         call get_HenrysLawCts('DMS', Vect_Hcts(1), Vect_Hcts(2), Vect_Hcts(3), Vect_Hcts(4), _RC)
         call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, _RC)

         call ESMF_StateGet(internal, 'SO2', field, _RC)
         call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(2), _RC)
         ! Vect_Hcts=(/1.4, 2900.0, 1.3E-02, 2000.0/)
         call get_HenrysLawCts('SO2', Vect_Hcts(1), Vect_Hcts(2), Vect_Hcts(3), Vect_Hcts(4), _RC)
         call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, _RC)

         call ESMF_StateGet(internal, 'MSA', field, rc=status)
         if (status == 0) then
            call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(4), _RC)
         end if
      end if

      ! Fill AERO State with SO4
      call ESMF_StateGet(export, trim(COMP_NAME) // '_AERO', aero, _RC)
      call ESMF_StateGet(export, trim(COMP_NAME) // '_AERO_DP', Bundle_DP, _RC)

      call ESMF_StateGet(internal, 'SO4', field, _RC)
      call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(3), _RC)
      fld = MAPL_FieldCreate(field, 'SO4', _RC)
      call MAPL_StateAdd(aero, fld, _RC)

      if (data_driven) then
         instance = instanceData

         ! (SO4; only aerosol component; bin 003)
         ! Dry deposition
         call append_to_bundle('SUDP003', providerState, prefix, Bundle_DP, _RC)

         ! Wet deposition (Convective scavenging)
         call append_to_bundle('SUSV003', providerState, prefix, Bundle_DP, _RC)

         ! Wet deposition
         call append_to_bundle('SUWT003', providerState, prefix, Bundle_DP, _RC)

         ! Gravitational Settling
         call append_to_bundle('SUSD003', providerState, prefix, Bundle_DP, _RC)
      else
         instance = instanceComputational

         ! Dry deposition
         call append_to_bundle('SUDP', providerState, prefix, Bundle_DP, _RC)

         ! Wet deposition (Convective scavenging)
         call append_to_bundle('SUSV', providerState, prefix, Bundle_DP, _RC)

         ! Wet deposition
         call append_to_bundle('SUWT', providerState, prefix, Bundle_DP, _RC)

         ! Gravitational Settling
         call append_to_bundle('SUSD', providerState, prefix, Bundle_DP, _RC)
      end if

      self%instance = instance

      ! Create Radiation Mie Table
      call ESMF_ConfigGetAttribute(cfg, file_, label="aerosol_radBands_optics_file:", _RC)
      self%rad_Mie = GOCART2G_Mie(trim(file_), _RC)

      ! Trigger for photolysis calculations
      call ESMF_AttributeSet(aero, NAME="use_photolysis_table", VALUE=0, _RC)

      ! Create Photolysis Mie Table
      ! Get file names for the optical tables
      call ESMF_ConfigGetAttribute(cfg, file_, &
           label="aerosol_monochromatic_optics_file:", _RC)
      call ESMF_ConfigGetAttribute(universal_cfg, nmom_, label="n_phase_function_moments_photolysis:", default=0, &
           _RC)
      i = ESMF_ConfigGetLen(universal_cfg, label='aerosol_photolysis_wavelength_in_nm_from_LUT:', _RC)
      allocate(channels_(i), _STAT)
      call ESMF_ConfigGetAttribute(universal_cfg, channels_, &
           label="aerosol_photolysis_wavelength_in_nm_from_LUT:", _RC)
      self%phot_Mie = GOCART2G_Mie(trim(file_), channels_ * 1.e-9, nmom=nmom_, _RC)
      deallocate(channels_)

      ! Create Diagnostics Mie Table
      ! Get file names for the optical tables
      call ESMF_ConfigGetAttribute(cfg, file_, &
           label="aerosol_monochromatic_optics_file:", _RC)
      call ESMF_ConfigGetAttribute(cfg, nmom_, label="n_moments:", default=0, _RC)
      i = ESMF_ConfigGetLen(universal_cfg, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', _RC)
      allocate(channels_(i), _STAT)
      call ESMF_ConfigGetAttribute(universal_cfg, channels_, &
           label="aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:", _RC)
      self%diag_Mie = GOCART2G_Mie(trim(file_), channels_ * 1.e-9, nmom=nmom_, _RC)
      deallocate(channels_)

      ! Mie Table instance/index
      call ESMF_AttributeSet(aero, NAME='mie_table_instance', VALUE=instance, _RC)

      ! Add variables to SU instance's aero state. This is used in aerosol optics calculations
      call add_aero_named_alias(aero, label='air_pressure_for_aerosol_optics', label2='PLE', grid=grid, typekind=&
           MAPL_R4, _RC)
      call add_aero_named_alias(aero, label='relative_humidity_for_aerosol_optics', label2='RH', grid=grid, typekind=&
           MAPL_R4, _RC)
      ! call ESMF_StateGet (import, 'PLE', field, __RC__)
      ! call MAPL_StateAdd (aero, field, __RC__)
      ! call ESMF_StateGet (import, 'RH2', field, __RC__)
      ! call MAPL_StateAdd (aero, field, __RC__)

      ! Add variables to SU instance aero state for chemistry
      call add_aero_named_alias(aero, label='surface_area_density', label2='SAREA', grid=grid, typekind=MAPL_R4, _RC)
      call add_aero_named_alias(aero, label='effective_radius_in_microns', label2='REFF', grid=grid, typekind=MAPL_R4, &
           _RC)

      call add_aero_named_alias(aero, label='extinction_in_air_due_to_ambient_aerosol', label2='EXT', grid=grid, &
           typekind=MAPL_R8, _RC)
      call add_aero_named_alias(aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, &
           typekind=MAPL_R8, _RC)
      call add_aero_named_alias(aero, label='asymmetry_parameter_of_ambient_aerosol', label2='ASY', grid=grid, typekind&
           =MAPL_R8, _RC)
      call ESMF_ConfigGetAttribute(universal_cfg, nmom_, label='n_phase_function_moments_photolysis:', default=0, &
           _RC)
      if (nmom_ > 0) then
         call add_aero_named_alias(aero, label='legendre_coefficients_of_p11_for_photolysis', label2='MOM', &
              grid=grid, typekind=MAPL_R8, ungrid=nmom_, _RC)
      end if
      call add_aero_named_alias(aero, label='monochromatic_extinction_in_air_due_to_ambient_aerosol', &
           label2='monochromatic_EXT', grid=grid, typekind=MAPL_R4, _RC)
      call add_aero_named_alias(aero, label='sum_of_internalState_aerosol', label2='aerosolSum', grid=grid, typekind=&
           MAPL_R4, _RC)

      call ESMF_AttributeSet(aero, NAME='band_for_aerosol_optics', VALUE=0, _RC)
      call ESMF_AttributeSet(aero, NAME='wavelength_for_aerosol_optics', VALUE=0., _RC)

      mieTable_pointer = transfer(c_loc(self), [1])
      call ESMF_AttributeSet(aero, NAME='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer&
           ), _RC)

      call ESMF_AttributeSet(aero, NAME='internal_variable_name', VALUE='SO4', _RC)

      call ESMF_MethodAdd(aero, label='aerosol_optics', userRoutine=aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label='monochromatic_aerosol_optics', userRoutine=monochromatic_aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label='get_mixR', userRoutine=get_mixR, _RC)

      ! Deal with replenishment alarm (formerly the daily_alarm subroutine)
      self%alarm = createReplenishAlarm(gc, clock, 30000, _RC)
      _RETURN(ESMF_SUCCESS)

   end subroutine Initialize

   !IROUTINE: Run0

   !INTERFACE:
   subroutine Run0(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Clears klid to 0.0 for Sulfates

      ! Locals
      character(len=ESMF_MAXSTR) :: COMP_NAME
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal
      type(wrap_) :: wrap
      type(SU2G_GridComp), pointer :: self
      real, pointer, dimension(:, :, :) :: ple
      real, pointer, dimension(:, :, :) :: ptr3d_int

      integer :: status
      character(len=255) :: Iam = 'Run0'

      ! Begin...

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(gc, NAME=COMP_NAME, _RC)
      Iam = trim(COMP_NAME) // '::' // Iam

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      ! Get parameters from generic state.
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, _RC)

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, 'SU2G_GridComp', wrap, status)
      _VERIFY(status)
      self => wrap%PTR

      ! Set klid and Set internal values to 0 above klid
      call MAPL_GetPointer(import, ple, 'PLE', _RC)
      call findKlid(self%klid, self%plid, ple, _RC)
      call MAPL_GetPointer(internal, NAME='SO4', PTR=ptr3d_int, _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_GetPointer(internal, NAME='SO2', PTR=ptr3d_int, _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_GetPointer(internal, NAME='DMS', PTR=ptr3d_int, _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_GetPointer(internal, NAME='MSA', PTR=ptr3d_int, _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)

      _RETURN(ESMF_SUCCESS)

   end subroutine Run0

   !IROUTINE: Run

   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run method for the Sulfate Component. Determines whether to run
      ! data or computational run method.

      ! Locals
      character(len=ESMF_MAXSTR) :: COMP_NAME
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal

      logical :: data_driven

      integer :: status
      character(len=255) :: Iam = 'Run'

      ! Begin...

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(gc, NAME=COMP_NAME, _RC)
      Iam = trim(COMP_NAME) // '::' // Iam

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      ! Get parameters from generic state.
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, _RC)

      ! Is SU data driven?
      call determine_data_driven(COMP_NAME, data_driven, _RC)

      ! Update INTERNAL state variables with ExtData
      if (data_driven) then
         call Run_data(gc, import, export, internal, _RC)
      else
         call Run1(gc, import, export, clock, _RC)
      end if

      _RETURN(ESMF_SUCCESS)

   end subroutine Run

   !IROUTINE: Run1

   !INTERFACE:
   subroutine Run1(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Computes emissions/sources for Sea Salt

      ! Locals
      character(len=ESMF_MAXSTR) :: COMP_NAME
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal
      type(ESMF_Grid) :: grid
      type(wrap_) :: wrap
      type(SU2G_GridComp), pointer :: self
      type(ESMF_Time) :: time

      character(len=3) :: cdow
      integer :: idow
      integer :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
      real, pointer, dimension(:, :) :: lats
      real, pointer, dimension(:, :) :: lons
      real, dimension(:, :, :), allocatable :: aircraft_fuel_src, so2_ocs_src
      real, dimension(:, :), allocatable :: so2biomass_src, so2biomass_src_, so2anthro_l1_src, &
           so2anthro_l2_src, so2ship_src, so4ship_src, dmso_conc, &
           aviation_lto_src, aviation_cds_src, aviation_crs_src
      integer, dimension(:), allocatable :: iPoint, jPoint
      real, dimension(:, :, :), allocatable :: emissions_point
      character(len=ESMF_MAXSTR) :: fname ! file name for point source emissions
      logical :: fileExists

      real, pointer, dimension(:, :, :) :: dummyMSA !=> null() ! This is so the model can run without MSA enabled
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread

#include "SU2G_DeclarePointer___.h"

      integer :: status
      character(len=255) :: Iam = 'Run1'

      ! Begin...

      nullify(dummyMSA)

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(gc, NAME=COMP_NAME, _RC)
      Iam = trim(COMP_NAME) // '::' // Iam

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      call MAPL_Get(MAPL, grid=grid, _RC)

      ! Get parameters from generic state.
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, &
           lons=lons, &
           lats=lats, _RC)

#include "SU2G_GetPointer___.h"

      call MAPL_GetPointer(internal, dummyMSA, 'MSA', rc=status)

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, 'SU2G_GridComp', wrap, status)
      _VERIFY(status)
      self => wrap%PTR

      ! Extract nymd(yyyymmdd) from clock
      call ESMF_ClockGet(clock, currTime=time, _RC)
      call ESMF_TimeGet(time, YY=iyr, mm=imm, DD=idd, H=ihr, M=imn, S=isc, _RC)
      nymd = MAPL_PackedDateCreate(iyr, imm, idd)
      nhms = MAPL_PackedTimeCreate(ihr, imn, isc)

      ! Reset tracer to zero at 0Z on specific day of week
      idow = Chem_UtilIdow(nymd)
      if ((nhms == 0) .and. (idow == self%myDOW)) then
         cdow = Chem_UtilCdow(nymd)
         DMS = tiny(1.) ! avoid division by zero
         SO2 = tiny(1.) ! avoid division by zero
         SO4 = tiny(1.) ! avoid division by zero
         if (associated(dummyMSA)) dummyMSA = tiny(1.) ! avoid division by zero
         if (mapl_am_i_root()) then
            print *, '<> SU ' // cdow // ' tracer being set to zero on ', nymd, nhms
         end if
      end if

      ! Implicit allocation with Fortran 2003
      so2anthro_l1_src = SU_ANTHROL1
      so2anthro_l2_src = SU_ANTHROL2
      so2ship_src = SU_SHIPSO2
      so4ship_src = SU_SHIPSO4

      ! As a safety check, where value is undefined set to 0
      where (1.01 * so2anthro_l1_src > MAPL_UNDEF) so2anthro_l1_src = 0.
      where (1.01 * so2anthro_l2_src > MAPL_UNDEF) so2anthro_l2_src = 0.
      where (1.01 * so2ship_src > MAPL_UNDEF) so2ship_src = 0.
      where (1.01 * so4ship_src > MAPL_UNDEF) so4ship_src = 0.

      aircraft_fuel_src = SU_AIRCRAFT
      so2biomass_src = SU_BIOMASS
      dmso_conc = SU_DMSO
      aviation_lto_src = SU_AVIATION_LTO
      aviation_cds_src = SU_AVIATION_CDS
      aviation_crs_src = SU_AVIATION_CRS
      so2_ocs_src = pSO2_OCS

      ! PRC: turn off emissions
      if (self%disable_emissions) then
         so2anthro_l1_src = 0.
         so2anthro_l2_src = 0.
         so2ship_src = 0.
         so4ship_src = 0.
         aircraft_fuel_src = 0.
         so2biomass_src = 0.
         dmso_conc = 0.
         aviation_lto_src = 0.
         aviation_cds_src = 0.
         aviation_crs_src = 0.
         so2_ocs_src = 0.
      end if

      ! As a safety check, where value is undefined set to 0
      where (1.01 * so2biomass_src > MAPL_UNDEF) so2biomass_src = 0.
      where (1.01 * dmso_conc > MAPL_UNDEF) dmso_conc = 0.
      where (1.01 * aircraft_fuel_src > MAPL_UNDEF) aircraft_fuel_src = 0.
      where (1.01 * aviation_lto_src > MAPL_UNDEF) aviation_lto_src = 0.
      where (1.01 * aviation_cds_src > MAPL_UNDEF) aviation_cds_src = 0.
      where (1.01 * aviation_crs_src > MAPL_UNDEF) aviation_crs_src = 0.

      ! Start with a clean emission diagnostic
      if (associated(SUEM)) SUEM = 0.0

      ! Update emissions/production if necessary (daily)
      thread = MAPL_get_current_thread()
      workspace => self%workspaces(thread)
      ! Update Volcanic SO2 Emissions Daily
      if (workspace%nymd_last /= nymd) then
         workspace%nymd_last = nymd

         ! DEGASSING: Get pointwise SO2 and altitude of volcanoes from a daily file data base
         workspace%nVolc = 0 ! case of /dev/null (no volcanoes) or ill-formed filename
         if (index(self%volcano_srcfilen_degassing, 'volcanic_') /= 0) then
            call StrTemplate(fname, self%volcano_srcfilen_degassing, xid='unknown', &
                 nymd=nymd, nhms=120000)
            inquire(file=fname, exist=fileExists)
            if (fileExists) then
               call ReadPointEmissions(nymd, fname, workspace%nVolc, workspace%vLat, workspace%vLon, &
                    workspace%vElev, workspace%vCloud, workspace%vSO2, workspace%vStart, &
                    workspace%vEnd, label='volcano', _RC)
               workspace%vSO2 = workspace%vSO2 * fMassSO2 / fMassSulfur
            end if
         end if

         ! EXPLOSIVE: Get pointwise SO2 and altitude of volcanoes from a daily file data base
         workspace%nVolcE = 0 ! case of /dev/null (no volcanoes) or ill-formed filename
         if (index(self%volcano_srcfilen_explosive, 'volcanic_') /= 0) then
            call StrTemplate(fname, self%volcano_srcfilen_explosive, xid='unknown', &
                 nymd=nymd, nhms=120000)
            inquire(file=fname, exist=fileExists)
            if (fileExists) then
               call ReadPointEmissions(nymd, fname, workspace%nVolcE, workspace%vLatE, workspace%vLonE, &
                    workspace%vElevE, workspace%vCloudE, workspace%vSO2E, workspace%vStartE, &
                    workspace%vEndE, label='volcano', _RC)
               workspace%vSO2E = workspace%vSO2E * fMassSO2 / fMassSulfur
            end if
         end if

      end if

      ! DEGASSING: Apply volcanic emissions
      if (workspace%nVolc > 0) then
         if (associated(SO2EMVN)) SO2EMVN = 0.0
         allocate(iPoint(workspace%nVolc), jPoint(workspace%nVolc), _STAT)
         call MAPL_GetHorzIJIndex(workspace%nVolc, iPoint, jPoint, &
              grid=grid, &
              lon=workspace%vLon / real(MAPL_RADIANS_TO_DEGREES), &
              lat=workspace%vLat / real(MAPL_RADIANS_TO_DEGREES), &
              rc=status)
         if (status /= 0) then
            if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
            _VERIFY(status)
         end if

         call SUvolcanicEmissions(workspace%nVolc, workspace%vStart, workspace%vEnd, workspace%vSO2, workspace%vElev, &
              workspace%vCloud, iPoint, jPoint, nhms, SO2EMVN, SO2, nSO2, SUEM, &
              self%km, self%CDT, MAPL_GRAV, zle, delp, area, workspace%vLat, workspace%vLon, _RC)
         deallocate(iPoint, jPoint, _STAT)
      end if

      ! EXPLOSIVE: Apply volcanic emissions
      if (workspace%nVolcE > 0) then
         if (associated(SO2EMVE)) SO2EMVE = 0.0
         allocate(iPoint(workspace%nVolcE), jPoint(workspace%nVolcE), _STAT)
         call MAPL_GetHorzIJIndex(workspace%nVolcE, iPoint, jPoint, &
              grid=grid, &
              lon=workspace%vLonE / real(MAPL_RADIANS_TO_DEGREES), &
              lat=workspace%vLatE / real(MAPL_RADIANS_TO_DEGREES), &
              rc=status)
         if (status /= 0) then
            if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
            _VERIFY(status)
         end if

         call SUvolcanicEmissions(workspace%nVolcE, workspace%vStartE, workspace%vEndE, workspace%vSO2E, workspace%&
              vElevE, &
              workspace%vCloudE, iPoint, jPoint, nhms, SO2EMVE, SO2, nSO2, SUEM, &
              self%km, self%CDT, MAPL_GRAV, zle, delp, area, workspace%vLatE, workspace%vLonE, _RC)
         deallocate(iPoint, jPoint, _STAT)
      end if

      ! Apply diurnal cycle if so desired
      if (self%diurnal_bb) then
         so2biomass_src_ = so2biomass_src
         call Chem_BiomassDiurnal(so2biomass_src, so2biomass_src_, &
              lons(:, :) * real(MAPL_RADIANS_TO_DEGREES), &
              lats(:, :) * real(MAPL_RADIANS_TO_DEGREES), &
              nhms, self%CDT)
      end if

      ! Apply sulfate emissions
      call SulfateDistributeEmissions(self%km, self%nbins, self%CDT, MAPL_GRAV, nymd, nhms, &
           fMassSO4, fMassSO2, self%fSO4anth, self%eAircraftFuel, &
           nSO2, nSO4, &
           so2anthro_l1_src, so2anthro_l2_src, &
           so2biomass_src, dmso_conc, &
           so2ship_src, so4ship_src, &
           aircraft_fuel_src, &
           SO2, SO4, &
           lwi, u10m, v10m, zle, zpbl, &
           t, airdens, delp, workspace%nVolc, &
           SUEM, SO4EMAN, SO2EMAN, SO2EMBB, &
           self%aviation_layers, &
           aviation_lto_src, &
           aviation_cds_src, &
           aviation_crs_src, _RC)

      if (associated(DMS)) then
         call DMSemission(self%km, self%CDT, MAPL_GRAV, t, u10m, v10m, lwi, delp, &
              fMassDMS, dmso_conc, DMS, SUEM, nDMS, _RC)
      end if

      ! Add source of OCS-produced SO2
      SO2 = SO2 + so2_ocs_src * self%CDT

      ! Read any pointwise emissions, if requested
      if (self%doing_point_emissions) then
         call StrTemplate(fname, self%point_emissions_srcfilen, xid='unknown', &
              nymd=nymd, nhms=120000)
         inquire(file=fname, exist=fileExists)
         if (fileExists) then
            call ReadPointEmissions(nymd, fname, workspace%nPts, workspace%pLat, workspace%pLon, &
                 workspace%pBase, workspace%pTop, workspace%pEmis, workspace%pstart, &
                 workspace%pend, label='source', _RC)
         else if (.not. fileExists) then
            if (mapl_am_i_root()) print*, 'GOCART2G ', trim(COMP_NAME), ': ', trim(fname), ' not found; proceeding.'
            workspace%nPts = -1 ! set this back to -1 so the "if (workspace%nPts > 0)" conditional is not exercised.
         end if
      end if

      ! Get indices for point emissions
      if (workspace%nPts > 0) then
         allocate(emissions_point, mold=delp,  _STAT)
         emissions_point = 0.0
         allocate(iPoint(workspace%nPts), jPoint(workspace%nPts), _STAT)
         call MAPL_GetHorzIJIndex(workspace%nPts, iPoint, jPoint, &
              grid=grid, &
              lon=workspace%pLon / real(MAPL_RADIANS_TO_DEGREES), &
              lat=workspace%pLat / real(MAPL_RADIANS_TO_DEGREES), &
              rc=status)
         if (status /= 0) then
            if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
            _VERIFY(status)
         end if

         call updatePointwiseEmissions(self%km, workspace%pBase, workspace%pTop, workspace%pEmis, workspace%nPts, &
              workspace%pstart, workspace%pend, zle, &
              area, iPoint, jPoint, nhms, emissions_point, _RC)

         deallocate(iPoint, jPoint, _STAT)
         SO4 = SO4 + self%CDT * MAPL_GRAV / delp * emissions_point
      end if

      _RETURN(ESMF_SUCCESS)

   end subroutine Run1

   !IROUTINE: Run2

   !INTERFACE:

   subroutine Run2(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run2 method for the Sulfate Grid Component.

      ! Locals
      character(len=ESMF_MAXSTR) :: COMP_NAME
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(wrap_) :: wrap
      type(SU2G_GridComp), pointer :: self
      type(ESMF_Time) :: time
      type(MAPL_VarSpec), pointer :: InternalSpec(:)

      integer :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
      integer :: n
      logical :: KIN
      real, pointer, dimension(:, :) :: lats
      real, pointer, dimension(:, :) :: lons
      character(len=ESMF_MAXSTR) :: SHORT_NAME, fld_name
      real, pointer, dimension(:, :, :) :: int_ptr

      real, dimension(:, :, :), allocatable :: xoh, xno3, xh2o2

      real, dimension(:, :), allocatable :: drydepositionf
      real, pointer, dimension(:, :, :) :: susd_vel

      real, pointer, dimension(:, :, :) :: dummyMSA !=> null() ! this is so the model can run without MSA enabled
      logical :: alarm_is_ringing
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread
      integer :: i1, j1, i2, j2, km
      real, dimension(3) :: rainout_eff
      real, target, allocatable, dimension(:, :, :) :: RH20, RH80
      real, pointer, dimension(:, :) :: flux_ptr
      integer :: settling_opt
      type(ESMF_StateItem_Flag) :: item_type
#include "SU2G_DeclarePointer___.h"

      integer :: status
      character(len=255) :: Iam = 'Run2'

      ! Begin...

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(gc, NAME=COMP_NAME, _RC)
      Iam = trim(COMP_NAME) // '::' // Iam

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      ! Get parameters from generic state.
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, &
           InternalSpec=InternalSpec, &
           lons=lons, &
           lats=lats, _RC)

      ! Get the aero state
      call ESMF_StateGet(export, trim(COMP_NAME) // '_AERO', aero, _RC)

#include "SU2G_GetPointer___.h"

      call MAPL_GetPointer(internal, dummyMSA, 'MSA', rc=status)

      ! Extract nymd(yyyymmdd) from clock
      call ESMF_ClockGet(clock, currTime=time, _RC)
      call ESMF_TimeGet(time, YY=iyr, mm=imm, DD=idd, H=ihr, M=imn, S=isc, _RC)
      nymd = MAPL_PackedDateCreate(iyr, imm, idd)
      nhms = MAPL_PackedTimeCreate(ihr, imn, isc)

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, 'SU2G_GridComp', wrap, status)
      _VERIFY(status)
      self => wrap%PTR

      ! Set klid and Set internal values to 0 above klid
      call findKlid(self%klid, self%plid, ple, _RC)
      call setZeroKlid(self%km, self%klid, SO4)
      call setZeroKlid(self%km, self%klid, SO2)
      call setZeroKlid(self%km, self%klid, DMS)
      call setZeroKlid(self%km, self%klid, MSA)

      thread = MAPL_get_current_thread()
      workspace => self%workspaces(thread)

      alarm_is_ringing = ESMF_AlarmIsRinging(self%alarm, _RC)
      ! recycle H2O2 every 3 hours
      if (alarm_is_ringing) then
         workspace%recycle_H2O2 = .true.
#ifdef DEBUG
         if (mapl_am_i_root()) then
            print *, 'DEBUG:: SU replenish alarm is ringing'
         end if
#endif
      end if

      allocate(xoh, mold=airdens, _STAT)
      allocate(xno3, mold=airdens, _STAT)
      allocate(xh2o2, mold=airdens, _STAT)

      if (self%using_GMI) then
         xoh = GMI_OH
         xh2o2 = GMI_H2O2
         xno3 = GMI_NO3
         xoh = xoh * (MAPL_AVOGAD / 1000.) / MAPL_AIRMW * 1000. * airdens * 1.00E-06
         call MAPL_MaxMin('GMI:OH   ', xoh)
         call MAPL_MaxMin('GMI:H2O2 ', xh2o2)
         call MAPL_MaxMin('GMI:NO3  ', xno3)
         call MAPL_MaxMin('GMI:rhoa ', airdens)
      else

         xoh = 0.0
         xno3 = 0.0
         xh2o2 = h2o2_init

         call SulfateUpdateOxidants(nymd, nhms, lons, lats, airdens, self%km, self%CDT, &
              workspace%nymd_oxidants, MAPL_UNDEF, real(MAPL_RADIANS_TO_DEGREES), &
              MAPL_AVOGAD / 1000., MAPL_PI, MAPL_AIRMW, &
              SU_OH, SU_NO3, SU_H2O2, &
              xoh, xno3, xh2o2, workspace%recycle_H2O2, _RC)
      end if

      ! SU Settling
      select case (self%settling_scheme)
      case ('gocart')
         settling_opt = 1
      case ('ufs')
         settling_opt = 2
      case default
         _ASSERT_RC(.false., 'Unsupported settling scheme: ' // trim(self%settling_scheme), ESMF_RC_NOT_IMPL)
      end select

      ! Set default export value to 0.0 for all tracers
      if (associated(SUSD)) SUSD(:, :, :) = 0.0
      if (associated(SUSD_V)) SUSD_V(:, :, :, :) = 0.0
      ! Do settling only for sulfate aerosol tracer NSO4
      call MAPL_VarSpecGet(InternalSpec(nSO4), SHORT_NAME=SHORT_NAME, _RC)
      call MAPL_GetPointer(internal, NAME=SHORT_NAME, PTR=int_ptr, _RC)
      nullify(flux_ptr)
      if (associated(SUSD)) flux_ptr => SUSD(:, :, nSO4)
      nullify(susd_vel)
      if (associated(SUSD_V)) susd_vel => SUSD_V(:, :, :, nSO4)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 1, self%CDT, MAPL_GRAV, &
           int_ptr, t, airdens, &
           rh2, zle, delp, flux_ptr, susd_vel, settling_scheme=settling_opt, _RC)

      allocate(drydepositionf, mold=lwi, _STAT)
      call SulfateChemDriver(self%km, self%klid, self%CDT, MAPL_PI, real(MAPL_RADIANS_TO_DEGREES), MAPL_KARMAN, &
           MAPL_AIRMW, MAPL_AVOGAD / 1000., cpd, MAPL_GRAV, &
           fMassMSA, fMassDMS, fMassSO2, fMassSO4, &
           nymd, nhms, lons, lats, &
           DMS, SO2, SO4, dummyMSA, &
           nDMS, nSO2, nSO4, nMSA, &
           xoh, xno3, xh2o2, h2o2_init, &
           delp, t, fcld, airdens, zle, &
           ustar, sh, lwi, zpbl, z0h, &
           SUDP, SUPSO2, SUPMSA, &
           SUPSO4, SUPSO4g, SUPSO4aq, &
           pso2, pmsa, pso4, pso4g, pso4aq, drydepositionf, & ! 3d diagnostics
           _RC)

      KIN = .true.
      select case (self%wet_removal_scheme)
      case ('gocart')
         call SU_Wet_Removal(self%km, self%nbins, self%klid, self%CDT, KIN, MAPL_GRAV, MAPL_AIRMW, &
              delp, fMassSO4, fMassSO2, &
              h2o2_init, ple, airdens, cn_prcp, ncn_prcp, pfl_lsan, pfi_lsan, t, &
              nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, dummyMSA, &
              SUWT, SUPSO4, SUPSO4WT, pso4, PSO4WET, _RC)

      case ('ufs')

         do n = 1, self%nbins
            rainout_eff = 0.0
            rainout_eff(1) = self%fwet_ice(n)
            rainout_eff(2) = self%fwet_snow(n)
            rainout_eff(3) = self%fwet_rain(n)

            call MAPL_VarSpecGet(InternalSpec(n), SHORT_NAME=SHORT_NAME, _RC)
            call MAPL_GetPointer(internal, NAME=SHORT_NAME, PTR=int_ptr, _RC)
            call WetRemovalUFS(self%km, self%klid, n, self%CDT, 'sulfate', &
                 KIN, MAPL_GRAV, self%radius(n), rainout_eff, self%washout_tuning, &
                 self%wet_radius_thr, int_ptr, ple, t, airdens, pfl_lsan, pfi_lsan, SUWT, _RC)
         end do

      case default
         _ASSERT_RC(.false., 'Unsupported wet removal scheme: ' // trim(self%wet_removal_scheme), ESMF_RC_NOT_IMPL)
      end select

      ! Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
      call SU_Compute_Diags(self%km, self%klid, self%rhop(nSO4), &
           MAPL_GRAV, MAPL_PI, nSO4, self%diag_Mie, &
           self%wavelengths_profile * 1.0e-9, self%wavelengths_vertint * 1.0e-9, &
           t, airdens, delp, ple, tropp, rh2, u, v, DMS, SO2, SO4, dummyMSA, &
           DMSSMASS, DMSCMASS, &
           MSASMASS, MSACMASS, &
           SO2SMASS, SO2CMASS, &
           SO4SMASS, SO4CMASS, &
           SUEXTTAU, SUSTEXTTAU, SUSCATAU, SUSTSCATAU, SO4MASS, SUCONC, SUEXTCOEF, &
           SUSCACOEF, SUBCKCOEF, SUANGSTR, SUFLUXU, SUFLUXV, SO4SAREA, SO4SNUM, SO4REFF, _RC)

      if (associated(SO4SAREA)) then
         nullify(int_ptr)
         call ESMF_StateGet(aero, 'surface_area_density', item_type, _RC)
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_GetPointer(aero, int_ptr, 'surface_area_density', _RC)
            int_ptr = SO4SAREA
         end if
      end if

      if (associated(SO4REFF)) then ! Note unit conversion below to microns
         nullify(int_ptr)
         call ESMF_StateGet(aero, 'effective_radius_in_microns', item_type, _RC)
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_GetPointer(aero, int_ptr, 'effective_radius_in_microns', _RC)
            int_ptr = SO4REFF * 1.e6
         end if
      end if

      i1 = lbound(rh2, 1)
      i2 = ubound(rh2, 1)
      j1 = lbound(rh2, 2)
      j2 = ubound(rh2, 2)
      km = ubound(rh2, 3)

      allocate(RH20(i1:i2, j1:j2, km), _STAT)
      allocate(RH80(i1:i2, j1:j2, km), _STAT)

      RH20(:, :, :) = 0.20
      call SU_Compute_Diags(km=self%km, klid=self%klid, &
           rhop=self%rhop(nSO4), &
           grav=MAPL_GRAV, pi=MAPL_PI, nSO4=nSO4, mie=self%diag_Mie, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           tmpu=t, rhoa=airdens, delp=delp, ple=ple, tropp=tropp, rh=RH20, u=u, v=v, &
           DMS=DMS, SO2=SO2, SO4=SO4, MSA=dummyMSA, extcoef=SUEXTCOEFRH20, &
           scacoef=SUSCACOEFRH20, _RC)

      RH80(:, :, :) = 0.80
      call SU_Compute_Diags(km=self%km, klid=self%klid, &
           rhop=self%rhop(nSO4), &
           grav=MAPL_GRAV, pi=MAPL_PI, nSO4=nSO4, mie=self%diag_Mie, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           tmpu=t, rhoa=airdens, delp=delp, ple=ple, tropp=tropp, rh=RH80, u=u, v=v, &
           DMS=DMS, SO2=SO2, SO4=SO4, MSA=dummyMSA, extcoef=SUEXTCOEFRH80,&
           scacoef=SUSCACOEFRH80, _RC)

      deallocate(xoh, xh2o2, xno3, stat=status)

      _RETURN(ESMF_SUCCESS)

   end subroutine Run2

   !IROUTINE: Run_data -- ExtData Sulfate Grid Component

   !INTERFACE:

   subroutine Run_data(gc, import, export, internal, rc)

      !ARGUMENTS:

      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_State), intent(inout) :: internal ! Interal state
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION: Updates pointers in Internal state with fields from ExtData.

      ! Locals
      character(len=ESMF_MAXSTR) :: COMP_NAME
      type(wrap_) :: wrap
      type(SU2G_GridComp), pointer :: self

      real, pointer, dimension(:, :, :) :: ptr3d_int, ptr3d_imp

      integer :: status
      character(len=255) :: Iam = 'Run_data'

      ! Begin...

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet(gc, NAME=COMP_NAME, _RC)
      Iam = trim(COMP_NAME) // '::' // Iam

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, 'SU2G_GridComp', wrap, status)
      _VERIFY(status)
      self => wrap%PTR

      ! Update interal data pointers with ExtData
      call MAPL_GetPointer(internal, NAME='SO4', PTR=ptr3d_int, _RC)
      call MAPL_GetPointer(import, NAME='climSO4', PTR=ptr3d_imp, _RC)
      ptr3d_int = ptr3d_imp

      _RETURN(ESMF_SUCCESS)

   end subroutine Run_data

   subroutine aerosol_optics(state, rc)

      implicit none

      ! ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      ! Local
      integer, parameter :: DP = kind(1.0d0)
      real, dimension(:, :, :), pointer :: ple, rh
      real(kind=DP), dimension(:, :, :), pointer :: var
      real(kind=DP), dimension(:, :, :, :), pointer :: var4d
      real, dimension(:, :, :), pointer :: q
      real, dimension(:, :, :, :), pointer :: q_4d
      integer, allocatable :: opaque_self(:)
      type(c_ptr) :: address
      type(SU2G_GridComp), pointer :: self

      type(ESMF_Field) :: fld
      character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)

      real(kind=DP), dimension(:, :, :), allocatable :: ext_s, ssa_s, asy_s ! (lon:,lat:,lev:)
      real(kind=DP), dimension(:, :, :, :), allocatable :: pmom_s ! (lon:,lat:,lev:,nmom:)
      real :: x
      integer :: instance
      integer :: n, nbins
      integer :: i1, j1, i2, j2, km
      integer :: band
      integer :: usePhotTable
      real :: wavelength

      integer :: i, j, k
      type(ESMF_StateItem_Flag) :: item_type

      integer :: status
      character(len=255) :: Iam = 'SU2G::aerosol_optics'

      ! Begin...

      ! Mie Table instance/index
      call ESMF_AttributeGet(state, NAME='mie_table_instance', VALUE=instance, _RC)

      ! Get aerosol names
      call ESMF_AttributeGet(state, NAME='internal_variable_name', itemCount=nbins, _RC)
      allocate(aerosol_names(nbins), _STAT)
      call ESMF_AttributeGet(state, NAME='internal_variable_name', valueList=aerosol_names, _RC)

      ! Radiation band
      band = 0
      call ESMF_AttributeGet(state, NAME='band_for_aerosol_optics', VALUE=band, _RC)

      ! Are we doing a photolysis calculation?
      usePhotTable = 0
      call ESMF_AttributeGet(state, NAME='use_photolysis_table', VALUE=usePhotTable, _RC)

      ! Pressure at layer edges
      call MAPL_GetPointer(state, ple, 'air_pressure_for_aerosol_optics', _RC)

      ! call MAPL_GetPointer (state, ple, 'PLE', __RC__)

      i1 = lbound(ple, 1)
      i2 = ubound(ple, 1)
      j1 = lbound(ple, 2)
      j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      ! Relative humidity
      call MAPL_GetPointer(state, rh, 'relative_humidity_for_aerosol_optics', _RC)

      ! call MAPL_GetPointer (state, rh, 'RH2', __RC__)

      allocate(ext_s(i1:i2, j1:j2, km), &
           ssa_s(i1:i2, j1:j2, km), &
           asy_s(i1:i2, j1:j2, km), _STAT)

      allocate(q_4d(i1:i2, j1:j2, km, nbins), _STAT)

      do n = 1, nbins
         call ESMF_StateGet(state, trim(aerosol_names(n)), field=fld, _RC)
         call ESMF_FieldGet(fld, farrayPtr=q, _RC)

         do k = 1, km
            do j = j1, j2
               do i = i1, i2
                  x = ((ple(i, j, k) - ple(i, j, k - 1)) * 0.01) * (100. / MAPL_GRAV)
                  q_4d(i, j, k, n) = x * q(i, j, k)
               end do
            end do
         end do
      end do

      call ESMF_AttributeGet(state, NAME='mieTable_pointer', itemCount=n, _RC)
      allocate(opaque_self(n), _STAT)
      call ESMF_AttributeGet(state, NAME='mieTable_pointer', valueList=opaque_self, _RC)

      address = transfer(opaque_self, address)
      call c_f_pointer(address, self)

      if (usePhotTable /= 0) then
         wavelength = band * 1.e-9
         allocate(pmom_s(i1:i2, j1:j2, km, self%phot_Mie%nmom), _STAT)
         call miephot_(self%phot_Mie, nbins, wavelength, q_4d, rh, ext_s, ssa_s, pmom_s, _RC)
      else
         call mie_(self%rad_Mie, nbins, band, q_4d, rh, ext_s, ssa_s, asy_s, _RC)
      end if

      call ESMF_StateGet(state, 'extinction_in_air_due_to_ambient_aerosol', item_type, _RC)
      if (item_type == ESMF_STATEITEM_FIELD) then
         call MAPL_GetPointer(state, var, 'extinction_in_air_due_to_ambient_aerosol', _RC)
         var = ext_s(:, :, :)
      end if

      if (item_type == ESMF_STATEITEM_FIELD) then
         call MAPL_GetPointer(state, var, 'single_scattering_albedo_of_ambient_aerosol', _RC)
         var = ssa_s(:, :, :)
      end if

      if (usePhotTable /= 0) then
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_GetPointer(state, var4d, 'legendre_coefficients_of_p11_for_photolysis', _RC)
            var4d = pmom_s(:, :, :, :)
         end if
      else
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_GetPointer(state, var, 'asymmetry_parameter_of_ambient_aerosol', _RC)
            var = asy_s(:, :, :)
         end if
      end if

      deallocate(ext_s, ssa_s, asy_s, _STAT)
      if (usePhotTable /= 0) deallocate(pmom_s, _STAT)
      deallocate(q_4d, _STAT)

      _RETURN(ESMF_SUCCESS)

   contains

      subroutine mie_(mie, nbins, band, q, rh, bext_s, bssa_s, basym_s, rc)

         implicit none

         type(GOCART2G_Mie), intent(inout) :: mie ! mie table
         integer, intent(in) :: nbins ! number of bins
         integer, intent(in) :: band ! channel
         real, intent(in) :: q(:, :, :, :) ! aerosol mass mixing ratio, kg kg-1
         real, intent(in) :: rh(:, :, :) ! relative humidity
         real(kind=8), intent(out) :: bext_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         real(kind=8), intent(out) :: bssa_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         real(kind=8), intent(out) :: basym_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         integer, intent(out) :: rc
         ! local
         integer :: l
         real :: bext(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! extinction
         real :: bssa(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! SSA
         real :: gasym(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! asymmetry parameter

         integer :: status
         character(len=255) :: Iam = 'SU2G::aerosol_optics::mie_'

         bext_s = 0.0d0
         bssa_s = 0.0d0
         basym_s = 0.0d0

         do l = 1, nbins
            ! tau is converted to bext
            call mie%Query(band, l, q(:, :, :, l), rh, tau=bext, gasym=gasym, ssa=bssa, _RC)

            bext_s = bext_s + bext ! extinction
            bssa_s = bssa_s + (bssa * bext) ! scattering
            basym_s = basym_s + gasym * (bssa * bext) ! asymmetry parameter multiplied by scattering

         end do

         _RETURN(ESMF_SUCCESS)

      end subroutine mie_

      subroutine miephot_(mie, nbins, wavelength, q, rh, bext_s, bssa_s, bpmom_s, rc)

         implicit none

         type(GOCART2G_Mie), intent(inout) :: mie ! mie table
         integer, intent(in) :: nbins ! number of bins
         real, intent(in) :: wavelength ! wavelength in nm
         real, intent(in) :: q(:, :, :, :) ! aerosol mass mixing ratio, kg kg-1
         real, intent(in) :: rh(:, :, :) ! relative humidity
         real(kind=DP), intent(out) :: bext_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         real(kind=DP), intent(out) :: bssa_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         real(kind=DP), intent(out) :: bpmom_s(size(pmom_s, 1), size(pmom_s, 2), size(pmom_s, 3), size(pmom_s, 4))
         integer, intent(out) :: rc

         ! local
         integer :: l, M
         real :: bext(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! extinction
         real :: bssa(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! SSA
         real :: pmom(size(pmom_s, 1), size(pmom_s, 2), size(pmom_s, 3), size(pmom_s, 4), 6)

         integer :: status
         character(len=255) :: Iam = 'SU2G::aerosol_optics::miephot_'

         bext_s = 0.0d0
         bssa_s = 0.0d0
         bpmom_s = 0.0d0

         do l = 1, nbins
            ! tau is converted to bext
            call mie%Query(wavelength, l, q(:, :, :, l), rh, tau=bext, pmom=pmom, ssa=bssa, _RC)
            bext_s = bext_s + bext ! extinction
            bssa_s = bssa_s + (bssa * bext) ! scattering
            do M = 1, mie%nmom
               bpmom_s(:, :, :, M) = bpmom_s(:, :, :, M) + pmom(:, :, :, M, 1) * (bssa * bext) ! moments multiplied by scattering
            end do
         end do

         _RETURN(ESMF_SUCCESS)

      end subroutine miephot_

   end subroutine aerosol_optics

   subroutine monochromatic_aerosol_optics(state, rc)

      implicit none

      ! ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      ! Local
      real, dimension(:, :, :), pointer :: ple, rh
      real, dimension(:, :), pointer :: var
      real, dimension(:, :, :), pointer :: q
      real, dimension(:, :, :, :), pointer :: q_4d
      integer, allocatable :: opaque_self(:)
      type(c_ptr) :: address
      type(SU2G_GridComp), pointer :: self

      type(ESMF_Field) :: fld
      character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)

      real, dimension(:, :, :), allocatable :: tau_s, tau ! (lon:,lat:,lev:)
      real :: x
      integer :: instance
      integer :: n, nbins
      integer :: i1, j1, i2, j2, km
      real :: wavelength
      integer :: i, j, k
      type(ESMF_StateItem_Flag) :: item_type

      integer :: status
      character(len=255) :: Iam = 'SU2G::monochromatic_aerosol_optics'

      ! Begin...

      ! Mie Table instance/index
      call ESMF_AttributeGet(state, NAME='mie_table_instance', VALUE=instance, _RC)

      ! Get aerosol names
      call ESMF_AttributeGet(state, NAME='internal_variable_name', itemCount=nbins, _RC)
      allocate(aerosol_names(nbins), _STAT)
      call ESMF_AttributeGet(state, NAME='internal_variable_name', valueList=aerosol_names, _RC)

      ! Radiation band
      call ESMF_AttributeGet(state, NAME='wavelength_for_aerosol_optics', VALUE=wavelength, _RC)

      ! Pressure at layer edges
      call MAPL_GetPointer(state, ple, 'air_pressure_for_aerosol_optics', _RC)
      ! call MAPL_GetPointer (state, ple, 'PLE', __RC__)

      i1 = lbound(ple, 1)
      i2 = ubound(ple, 1)
      j1 = lbound(ple, 2)
      j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      ! Relative humidity
      call MAPL_GetPointer(state, rh, 'relative_humidity_for_aerosol_optics', _RC)
      ! call MAPL_GetPointer (state, rh, 'RH2', __RC__)

      allocate(tau_s(i1:i2, j1:j2, km), &
           tau(i1:i2, j1:j2, km), _STAT)
      tau_s = 0.0
      tau = 0.0

      allocate(q_4d(i1:i2, j1:j2, km, nbins), _STAT)

      do n = 1, nbins
         call ESMF_StateGet(state, trim(aerosol_names(n)), field=fld, _RC)
         call ESMF_FieldGet(fld, farrayPtr=q, _RC)

         do k = 1, km
            do j = j1, j2
               do i = i1, i2
                  x = (ple(i, j, k) - ple(i, j, k - 1)) / MAPL_GRAV
                  q_4d(i, j, k, n) = x * q(i, j, k)
               end do
            end do
         end do
      end do

      call ESMF_AttributeGet(state, NAME='mieTable_pointer', itemCount=n, _RC)
      allocate(opaque_self(n), _STAT)
      call ESMF_AttributeGet(state, NAME='mieTable_pointer', valueList=opaque_self, _RC)

      address = transfer(opaque_self, address)
      call c_f_pointer(address, self)

      do n = 1, nbins
         call self%diag_Mie%Query(wavelength, n, q_4d(:, :, :, n), rh, tau=tau, _RC)
         tau_s = tau_s + tau
      end do

      call ESMF_StateGet(state, 'monochromatic_extinction_in_air_due_to_ambient_aerosol', item_type, _RC)
      if (item_type == ESMF_STATEITEM_FIELD) then
         call MAPL_GetPointer(state, var, 'monochromatic_extinction_in_air_due_to_ambient_aerosol', _RC)
         var = sum(tau_s, dim=3)
      end if

      deallocate(q_4d, _STAT)

      _RETURN(ESMF_SUCCESS)

   end subroutine monochromatic_aerosol_optics

end module SU2G_GridCompMod
