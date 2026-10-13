#include "MAPL.h"

!MODULE: SU2G_GridCompMod - GOCART refactoring of the SU gridded component

!INTERFACE:
module SU2G_GridCompMod

   !USES:
   use ESMF
   use pflogger, only: logger_t => logger
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

   ! MAPL2's undefined value; MAPL3's MAPL_UNDEF is -huge, which breaks the 1.01*x > undef masks
   real, parameter :: MAPL2_UNDEF = huge(1.)

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

   character(*), parameter :: PRIVATE_STATE = "SU2G_GridComp"

contains

   !IROUTINE: SetServices

   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, intent(out) :: rc ! return code

      !DESCRIPTION: Declares SU's private state, entry points (Initialize and
      ! the Run0/Run1/Run2 phases) and state specs.

      !REVISION HISTORY:
      ! 08July2020   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

      ! Locals
      character(len=:), allocatable :: comp_name
      type(SU2G_GridComp), pointer :: self
      class(logger_t), pointer :: logger
      character(len=ESMF_MAXSTR) :: field_name
      logical :: data_driven = .true.
      integer :: num_threads
      type(MAPL_UngriddedDim) :: ungrd_nbins
      type(MAPL_UngriddedDim) :: ungrd_wavelengths_profile, ungrd_wavelengths_vertint
      integer :: i, status

      call MAPL_GridCompGet(gc, name=comp_name, logger=logger, _RC)

      ! Wrap gridcomp's private state and store it in gridcomp
      _SET_NAMED_PRIVATE_STATE(gc, SU2G_GridComp, PRIVATE_STATE)

      ! Retrieve the private state
      _GET_NAMED_PRIVATE_STATE(gc, SU2G_GridComp, PRIVATE_STATE, self)

      num_threads = MAPL_get_num_threads()
      allocate(self%workspaces(0:num_threads - 1), _STAT)

      ! process generic config items
      call self%GA_Environment%load_from_config(gc, _RC)

      ! Defined UngriddedDim items
      ungrd_nbins = MAPL_UngriddedDim(self%nbins, name="nbins", units="1")
      ungrd_wavelengths_profile = MAPL_UngriddedDim( &
           size(self%wavelengths_profile), &
           name="wavelengths_profile", &
           units="nm")
      ungrd_wavelengths_vertint = MAPL_UngriddedDim( &
           size(self%wavelengths_vertint), &
           name="wavelengths_vertint", &
           units="nm")

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
      call determine_data_driven(comp_name, data_driven, _RC)

      ! Set entry points
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name="Run1", _RC)
      if (data_driven .neqv. .true.) then
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run2, phase_name="Run2", _RC)
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run0, phase_name="Run0", _RC)
      end if

      ! Import and Internal states if data instance
      if (data_driven) then
         _FAIL("data driven section has not been activated yet")

         ! Pressure at layer edges
         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='PLE', &
              standard_name='air_pressure', &
              units='Pa', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_EDGE, &
              restart_mode=MAPL_RESTART_SKIP, _RC)

         ! RH: is between 0 and 1
         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='RH2', &
              standard_name='Rel_Hum_after_moist', &
              units='1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              restart_mode=MAPL_RESTART_SKIP, _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_INTERNAL, &
              short_name='SO4', &
              standard_name='Sulphate aerosol', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              restart_mode=MAPL_RESTART_OPTIONAL, &
              fill_value=0.0, &
              add_to_export=.true., _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='climSO4', &
              standard_name='Sulphate aerosol', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              restart_mode=MAPL_RESTART_OPTIONAL, _RC)

         do i = 1, self%nbins
            write(field_name, '(A, I0.3)') '', i
            ! dry deposition
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climSUDP' // trim(field_name), &
                 standard_name='Sulfate dry deposition (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)

            ! wet deposition
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climSUWT' // trim(field_name), &
                 standard_name='Sulfate wet deposition (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)

            ! gravitational settling
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climSUSD' // trim(field_name), &
                 standard_name='Sulfate settling (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)

            ! convective scavenging
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climSUSV' // trim(field_name), &
                 standard_name='Sulfate convective scavenging (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)
         end do
      end if ! (data_driven)

      ! Import, Export, Internal states for computational instance
      if (.not. data_driven) then
#include "SU2G_Export___.h"
         associate(using_GMI => self%using_GMI)
#include "SU2G_Import___.h"
         end associate
#include "SU2G_Internal___.h"
         call logger%info("Wet removal scheme is %a", trim(self%wet_removal_scheme))
         call logger%info("Settling scheme is %a", trim(self%settling_scheme))
      end if

      ! This state holds fields needed by radiation
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name=trim(comp_name) // "_AERO", &
           standard_name="aerosols_from_" // trim(comp_name), &
           units="kg kg-1", &
           dims="xyz", &
           vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
           itemtype=MAPL_STATEITEM_STATE, _RC)

      ! This bundle is needed by surface for snow albedo modification
      ! by aerosol settling and deposition
      ! DEVELOPMENT NOTE - Change to StateItem in future
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name=trim(comp_name) // "_AERO_DP", &
           standard_name="aerosol_deposition_from_" // trim(comp_name), &
           units="kg m-2 s-1", &
           dims="xy", &
           vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
           itemtype=MAPL_STATEITEM_FIELDBUNDLE, _RC)

      _RETURN(_SUCCESS)

   end subroutine SetServices

   !IROUTINE: Initialize

   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, intent(out) :: rc ! Error code

      !DESCRIPTION: This initializes SU Grid Component.

      !REVISION HISTORY:
      ! 08July2019   E.Sherman  First attempt at refactoring

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_Geom) :: geom
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(ESMF_State) :: provider_state
      type(ESMF_FieldBundle) :: bundle_dp
      type(ESMF_Field) :: field
      type(ESMF_Info) :: field_info, aero_info
      type(ESMF_StateItem_Flag) :: item_type
      type(SU2G_GridComp), pointer :: self

      integer, allocatable :: mieTable_pointer(:)
      integer :: km
      integer :: instance
      character(len=ESMF_MAXSTR) :: prefix
      character(len=:), allocatable :: diurnal_bb
      real :: CDT ! chemistry timestep (secs)
      real :: HDT ! model timestep (secs)
      logical :: data_driven

      real, dimension(4) :: Vect_Hcts
      integer, allocatable, dimension(:) :: channels_
      integer :: nmom_
      character(len=:), allocatable :: file_
      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, geom=geom, num_levels=km, _RC)

      ! Get my internal private state
      _GET_NAMED_PRIVATE_STATE(gc, SU2G_GridComp, PRIVATE_STATE, self)

      self%km = km

      ! Get DTs
      call MAPL_ClockGet(clock, dt=HDT, _RC)
      call MAPL_GridCompGetResource(gc, "GOCART2G_DT", CDT, default=real(HDT), _RC)
      self%CDT = CDT

      ! Check whether to de-activate diurnal biomass burning (default is *on*)
      call MAPL_GridCompGetResource(gc, "DIURNAL_BIOMASS_BURNING", diurnal_bb, default='YES', _RC)
      diurnal_bb = ESMF_UtilStringUpperCase(diurnal_bb, _RC)
      self%diurnal_bb = (trim(diurnal_bb) == 'YES')

      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is SU data driven?
      call determine_data_driven(comp_name, data_driven, _RC)

      ! If this is a data component, the data is provided in the import
      ! state via ExtData instead of the actual GOCART children
      if (data_driven) then
         provider_state = import
         prefix = 'clim'
      else
         provider_state = export
         prefix = ''
      end if

      ! Add attribute information
      if (.not. data_driven) then
         call ESMF_StateGet(internal, 'DMS', field, _RC)
         call ESMF_InfoGetFromHost(field, field_info, _RC)
         call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(1), _RC)
         ! Vect_Hcts=(/0.56, 3500.0, 0.0, 0.0/)
         call get_HenrysLawCts('DMS', Vect_Hcts(1), Vect_Hcts(2), Vect_Hcts(3), Vect_Hcts(4), _RC)
         call ESMF_InfoSet(field_info, key='SetofHenryLawCts', values=Vect_Hcts, _RC)

         call ESMF_StateGet(internal, 'SO2', field, _RC)
         call ESMF_InfoGetFromHost(field, field_info, _RC)
         call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(2), _RC)
         ! Vect_Hcts=(/1.4, 2900.0, 1.3E-02, 2000.0/)
         call get_HenrysLawCts('SO2', Vect_Hcts(1), Vect_Hcts(2), Vect_Hcts(3), Vect_Hcts(4), _RC)
         call ESMF_InfoSet(field_info, key='SetofHenryLawCts', values=Vect_Hcts, _RC)

         call ESMF_StateGet(internal, 'MSA', item_type, _RC)
         if (item_type == ESMF_STATEITEM_FIELD) then
            call ESMF_StateGet(internal, 'MSA', field, _RC)
            call ESMF_InfoGetFromHost(field, field_info, _RC)
            call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(4), _RC)
         end if
      end if

      ! Fill AERO State with SO4
      call ESMF_StateGet(export, trim(comp_name) // '_AERO', aero, _RC)
      call ESMF_StateGet(export, trim(comp_name) // '_AERO_DP', bundle_dp, _RC)

      call ESMF_StateGet(internal, 'SO4', field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(3), _RC)
      call ESMF_StateAdd(aero, [field], _RC)

      if (data_driven) then
         instance = instanceData

         ! (SO4; only aerosol component; bin 003)
         ! Dry deposition
         call append_to_bundle('SUDP003', provider_state, prefix, bundle_dp, _RC)

         ! Wet deposition (Convective scavenging)
         call append_to_bundle('SUSV003', provider_state, prefix, bundle_dp, _RC)

         ! Wet deposition
         call append_to_bundle('SUWT003', provider_state, prefix, bundle_dp, _RC)

         ! Gravitational Settling
         call append_to_bundle('SUSD003', provider_state, prefix, bundle_dp, _RC)
      else
         instance = instanceComputational

         ! Dry deposition
         call append_to_bundle('SUDP', provider_state, prefix, bundle_dp, _RC)

         ! Wet deposition (Convective scavenging)
         call append_to_bundle('SUSV', provider_state, prefix, bundle_dp, _RC)

         ! Wet deposition
         call append_to_bundle('SUWT', provider_state, prefix, bundle_dp, _RC)

         ! Gravitational Settling
         call append_to_bundle('SUSD', provider_state, prefix, bundle_dp, _RC)
      end if

      self%instance = instance

      ! Create Radiation Mie Table
      call MAPL_GridCompGetResource(gc, "aerosol_radBands_optics_file", file_, _RC)
      self%rad_Mie = GOCART2G_Mie(trim(file_), _RC)

      ! Trigger for photolysis calculations
      call ESMF_InfoGetFromHost(aero, aero_info, _RC)
      call ESMF_InfoSet(aero_info, key="use_photolysis_table", value=0, _RC)

      ! Create Photolysis Mie Table
      ! Get file names for the optical tables
      call MAPL_GridCompGetResource(gc, "aerosol_monochromatic_optics_file", file_, _RC)
      call MAPL_GridCompGetResource(gc, "n_phase_function_moments_photolysis", nmom_, default=0, _RC)
      call MAPL_GridCompGetResource(gc, "aerosol_photolysis_wavelength_in_nm_from_LUT", channels_, _RC)
      self%phot_Mie = GOCART2G_Mie(trim(file_), channels_ * 1.e-9, nmom=nmom_, _RC)

      ! Create Diagnostics Mie Table
      ! Get file names for the optical tables
      call MAPL_GridCompGetResource(gc, "aerosol_monochromatic_optics_file", file_, _RC)
      call MAPL_GridCompGetResource(gc, "n_moments", nmom_, default=0, _RC)
      call MAPL_GridCompGetResource(gc, "aerosol_monochromatic_optics_wavelength_in_nm_from_LUT", channels_, _RC)
      _ASSERT(allocated(channels_), "problem with reading wavelengths")
      self%diag_Mie = GOCART2G_Mie(trim(file_), channels_ * 1.e-9, nmom=nmom_, _RC)

      ! Mie Table instance/index
      call ESMF_InfoSet(aero_info, key="mie_table_instance", value=instance, _RC)

      ! Add variables to SU instance's aero state. This is used in aerosol optics calculations
      call add_aero(aero, label="air_pressure_for_aerosol_optics", label2="PLE", geom=geom, km=self%km, _RC)
      call add_aero(aero, label="relative_humidity_for_aerosol_optics", label2="RH", geom=geom, km=self%km, _RC)

      ! Add variables to SU instance aero state for chemistry
      call add_aero(aero, label="surface_area_density", label2="SAREA", geom=geom, km=self%km, _RC)
      call add_aero(aero, label="effective_radius_in_microns", label2="REFF", geom=geom, km=self%km, _RC)

      call add_aero( &
           aero, &
           label="extinction_in_air_due_to_ambient_aerosol", label2="EXT", &
           geom=geom, km=self%km, typekind=ESMF_TYPEKIND_R8, _RC)
      call add_aero( &
           aero, &
           label="single_scattering_albedo_of_ambient_aerosol", label2="SSA", &
           geom=geom, km=self%km, typekind=ESMF_TYPEKIND_R8, _RC)
      call add_aero(aero, &
           label="asymmetry_parameter_of_ambient_aerosol", label2="ASY", &
           geom=geom, km=self%km, typekind=ESMF_TYPEKIND_R8, _RC)
      call MAPL_GridCompGetResource(gc, "n_phase_function_moments_photolysis", nmom_, default=0, _RC)
      if (nmom_ > 0) then
         call add_aero(aero, &
              label="legendre_coefficients_of_p11_for_photolysis", label2="MOM", &
              geom=geom, km=self%km, typekind=ESMF_TYPEKIND_R8, ungrid=nmom_, _RC)
      end if
      call add_aero(aero, &
           label="monochromatic_extinction_in_air_due_to_ambient_aerosol", label2="monochromatic_EXT", &
           geom=geom, typekind=ESMF_TYPEKIND_R4, _RC)
      call add_aero(aero, label="sum_of_internalState_aerosol", label2="aerosolSum", geom=geom, km=self%km, _RC)

      call ESMF_InfoSet(aero_info, key="band_for_aerosol_optics", value=0, _RC)
      call ESMF_InfoSet(aero_info, key="wavelength_for_aerosol_optics", value=0., _RC)

      mieTable_pointer = transfer(c_loc(self), [1])
      call ESMF_InfoSet(aero_info, key="mieTable_pointer", values=mieTable_pointer, _RC)

      call ESMF_InfoSet(aero_info, key="internal_variable_name", value="SO4", _RC)

      call ESMF_MethodAdd(aero, label="aerosol_optics", userRoutine=aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label="monochromatic_aerosol_optics", userRoutine=monochromatic_aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label="get_mixR", userRoutine=get_mixR, _RC)

      ! Deal with replenishment alarm (formerly the daily_alarm subroutine)
      self%alarm = createReplenishAlarm(gc, clock, 30000, _RC)

      _RETURN(_SUCCESS)

   end subroutine Initialize

   !IROUTINE: Run0

   !INTERFACE:
   subroutine Run0(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Clears klid to 0.0 for Sulfates

      ! Locals
      type(ESMF_State) :: internal
      type(SU2G_GridComp), pointer :: self
      real, pointer, dimension(:, :, :) :: ple, ple0
      real, pointer, dimension(:, :, :) :: ptr3d_int
      integer :: status

      ! Begin...

      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, SU2G_GridComp, PRIVATE_STATE, self)

      ! Edge variable PLE is expected to be 0-based
      call MAPL_StateGetPointer(import, ple, 'PLE', _RC)
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple

      ! Set klid and Set internal values to 0 above klid
      call findKlid(self%klid, self%plid, ple0, _RC)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'SO4', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'SO2', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'DMS', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'MSA', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)

   end subroutine Run0

   !IROUTINE: Run

   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run method for the Sulfate Component. Determines whether to run
      ! data or computational run method.

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      logical :: data_driven
      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is SU data driven?
      call determine_data_driven(comp_name, data_driven, _RC)

      ! Update INTERNAL state variables with ExtData
      if (data_driven) then
         call Run_data(gc, import, export, internal, _RC)
      else
         call Run1(gc, import, export, clock, _RC)
      end if

      _RETURN(_SUCCESS)

   end subroutine Run

   !IROUTINE: Run1

   !INTERFACE:
   subroutine Run1(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Computes emissions/sources for Sulfate

      ! Locals
      type(ESMF_State) :: internal
      type(ESMF_Grid) :: grid
      type(ESMF_Geom) :: geom
      type(SU2G_GridComp), pointer :: self
      type(ESMF_Time) :: time
      class(logger_t), pointer :: logger

      character(len=3) :: cdow
      integer :: idow
      integer :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
      real, allocatable, dimension(:, :) :: lats
      real, allocatable, dimension(:, :) :: lons
      real, dimension(:, :, :), allocatable :: aircraft_fuel_src, so2_ocs_src
      real, dimension(:, :), allocatable :: so2biomass_src, so2biomass_src_, so2anthro_l1_src, &
           so2anthro_l2_src, so2ship_src, so4ship_src, dmso_conc, &
           aviation_lto_src, aviation_cds_src, aviation_crs_src
      integer, dimension(:), allocatable :: iPoint, jPoint
      real, dimension(:, :, :), allocatable :: emissions_point
      character(len=ESMF_MAXSTR) :: fname ! file name for point source emissions
      logical :: fileExists

      real, pointer, dimension(:, :, :) :: dummyMSA !=> null() ! This is so the model can run without MSA enabled
      real, pointer, dimension(:, :, :) :: zle0
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread

#include "SU2G_DeclarePointer___.h"

      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, logger=logger, geom=geom, grid=grid, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, SU2G_GridComp, PRIVATE_STATE, self)

      associate(using_GMI => self%using_GMI)
#include "SU2G_GetPointer___.h"
      end associate

      dummyMSA => MSA

      ! Edge variable ZLE is expected to be 0-based
      zle0(lbound(zle, 1):, lbound(zle, 2):, 0:) => zle

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
         call logger%info('<> SU ' // cdow // ' tracer being set to zero on %i0 %i0', nymd, nhms)
      end if

      ! Implicit allocation with Fortran 2003
      so2anthro_l1_src = SU_ANTHROL1
      so2anthro_l2_src = SU_ANTHROL2
      so2ship_src = SU_SHIPSO2
      so4ship_src = SU_SHIPSO4

      ! As a safety check, where value is undefined set to 0
      where (1.01 * so2anthro_l1_src > MAPL2_UNDEF) so2anthro_l1_src = 0.
      where (1.01 * so2anthro_l2_src > MAPL2_UNDEF) so2anthro_l2_src = 0.
      where (1.01 * so2ship_src > MAPL2_UNDEF) so2ship_src = 0.
      where (1.01 * so4ship_src > MAPL2_UNDEF) so4ship_src = 0.

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
      where (1.01 * so2biomass_src > MAPL2_UNDEF) so2biomass_src = 0.
      where (1.01 * dmso_conc > MAPL2_UNDEF) dmso_conc = 0.
      where (1.01 * aircraft_fuel_src > MAPL2_UNDEF) aircraft_fuel_src = 0.
      where (1.01 * aviation_lto_src > MAPL2_UNDEF) aviation_lto_src = 0.
      where (1.01 * aviation_cds_src > MAPL2_UNDEF) aviation_cds_src = 0.
      where (1.01 * aviation_crs_src > MAPL2_UNDEF) aviation_crs_src = 0.

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
            call MAPL_StrTemplate(fname, self%volcano_srcfilen_degassing, xid='unknown', &
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
            call MAPL_StrTemplate(fname, self%volcano_srcfilen_explosive, xid='unknown', &
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
         call MAPL_GeomGetHorzIJIndex(geom, &
              workspace%vLon / real(MAPL_RADIANS_TO_DEGREES), &
              workspace%vLat / real(MAPL_RADIANS_TO_DEGREES), &
              iPoint, jPoint, local_indices=.true., _RC)

         call SUvolcanicEmissions(workspace%nVolc, workspace%vStart, workspace%vEnd, workspace%vSO2, workspace%vElev, &
              workspace%vCloud, iPoint, jPoint, nhms, SO2EMVN, SO2, nSO2, SUEM, &
              self%km, self%CDT, MAPL_GRAV, zle0, delp, area, workspace%vLat, workspace%vLon, _RC)
         deallocate(iPoint, jPoint, _STAT)
      end if

      ! EXPLOSIVE: Apply volcanic emissions
      if (workspace%nVolcE > 0) then
         if (associated(SO2EMVE)) SO2EMVE = 0.0
         call MAPL_GeomGetHorzIJIndex(geom, &
              workspace%vLonE / real(MAPL_RADIANS_TO_DEGREES), &
              workspace%vLatE / real(MAPL_RADIANS_TO_DEGREES), &
              iPoint, jPoint, local_indices=.true., _RC)

         call SUvolcanicEmissions(workspace%nVolcE, workspace%vStartE, workspace%vEndE, workspace%vSO2E, workspace%&
              vElevE, &
              workspace%vCloudE, iPoint, jPoint, nhms, SO2EMVE, SO2, nSO2, SUEM, &
              self%km, self%CDT, MAPL_GRAV, zle0, delp, area, workspace%vLatE, workspace%vLonE, _RC)
         deallocate(iPoint, jPoint, _STAT)
      end if

      ! Apply diurnal cycle if so desired
      if (self%diurnal_bb) then
         call MAPL_GridGetCoordinates(grid, latitudes=lats, longitudes=lons, _RC)
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
           lwi, u10m, v10m, zle0, zpbl, &
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
         call MAPL_StrTemplate(fname, self%point_emissions_srcfilen, xid='unknown', &
              nymd=nymd, nhms=120000)
         inquire(file=fname, exist=fileExists)
         if (fileExists) then
            call ReadPointEmissions(nymd, fname, workspace%nPts, workspace%pLat, workspace%pLon, &
                 workspace%pBase, workspace%pTop, workspace%pEmis, workspace%pstart, &
                 workspace%pend, label='source', _RC)
         else if (.not. fileExists) then
            !$omp critical (SU2G_1)
            call logger%info("[" // trim(fname) // "] not found; proceeding...")
            !$omp end critical (SU2G_1)
            workspace%nPts = -1 ! set this back to -1 so the "if (workspace%nPts > 0)" conditional is not exercised.
         end if
      end if

      ! Get indices for point emissions
      if (workspace%nPts > 0) then
         allocate(emissions_point, mold=delp,  _STAT)
         emissions_point = 0.0
         call MAPL_GeomGetHorzIJIndex(geom, &
              workspace%pLon / real(MAPL_RADIANS_TO_DEGREES), &
              workspace%pLat / real(MAPL_RADIANS_TO_DEGREES), &
              iPoint, jPoint, local_indices=.true., _RC)

         call updatePointwiseEmissions(self%km, workspace%pBase, workspace%pTop, workspace%pEmis, workspace%nPts, &
              workspace%pstart, workspace%pend, zle0, &
              area, iPoint, jPoint, nhms, emissions_point, _RC)

         deallocate(iPoint, jPoint, _STAT)
         SO4 = SO4 + self%CDT * MAPL_GRAV / delp * emissions_point
      end if

      _RETURN(_SUCCESS)

   end subroutine Run1

   !IROUTINE: Run2

   !INTERFACE:

   subroutine Run2(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run2 method for the Sulfate Grid Component.

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(ESMF_Info) :: aero_info
      type(ESMF_Grid) :: grid
      type(SU2G_GridComp), pointer :: self
      type(ESMF_Time) :: time

      ! Internal state order, as in SU2G_StateSpecs.rc (nDMS, nSO2, nSO4, nMSA)
      character(len=3), parameter :: internal_names(4) = ['DMS', 'SO2', 'SO4', 'MSA']

      integer :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
      integer :: n
      logical :: KIN
      real, allocatable, dimension(:, :) :: lats
      real, allocatable, dimension(:, :) :: lons
      character(len=ESMF_MAXSTR) :: fld_name
      real, pointer, dimension(:, :, :) :: int_ptr

      real, dimension(:, :, :), allocatable :: xoh, xno3, xh2o2

      real, dimension(:, :), allocatable :: drydepositionf
      real, pointer, dimension(:, :, :) :: susd_vel

      real, pointer, dimension(:, :, :) :: dummyMSA !=> null() ! this is so the model can run without MSA enabled
      real, pointer, dimension(:, :, :) :: ple0, zle0, pfl_lsan0, pfi_lsan0
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

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, grid=grid, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, SU2G_GridComp, PRIVATE_STATE, self)

      call MAPL_GridGetCoordinates(grid, latitudes=lats, longitudes=lons, _RC)

      ! Get the aero state
      call ESMF_StateGet(export, trim(comp_name) // '_AERO', aero, _RC)
      call ESMF_InfoGetFromHost(aero, aero_info, _RC)

      associate(using_GMI => self%using_GMI)
#include "SU2G_GetPointer___.h"
      end associate

      dummyMSA => MSA

      ! Edge variables are expected to be 0-based
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple
      zle0(lbound(zle, 1):, lbound(zle, 2):, 0:) => zle
      pfl_lsan0(lbound(pfl_lsan, 1):, lbound(pfl_lsan, 2):, 0:) => pfl_lsan
      pfi_lsan0(lbound(pfi_lsan, 1):, lbound(pfi_lsan, 2):, 0:) => pfi_lsan

      ! Extract nymd(yyyymmdd) from clock
      call ESMF_ClockGet(clock, currTime=time, _RC)
      call ESMF_TimeGet(time, YY=iyr, mm=imm, DD=idd, H=ihr, M=imn, S=isc, _RC)
      nymd = MAPL_PackedDateCreate(iyr, imm, idd)
      nhms = MAPL_PackedTimeCreate(ihr, imn, isc)

      ! Set klid and Set internal values to 0 above klid
      call findKlid(self%klid, self%plid, ple0, _RC)
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
              workspace%nymd_oxidants, MAPL2_UNDEF, real(MAPL_RADIANS_TO_DEGREES), &
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
      call MAPL_StateGetPointer(internal, int_ptr, internal_names(nSO4), _RC)
      nullify(flux_ptr)
      if (associated(SUSD)) flux_ptr => SUSD(:, :, nSO4)
      nullify(susd_vel)
      if (associated(SUSD_V)) susd_vel => SUSD_V(:, :, :, nSO4)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 1, self%CDT, MAPL_GRAV, &
           int_ptr, t, airdens, &
           rh2, zle0, delp, flux_ptr, susd_vel, settling_scheme=settling_opt, _RC)

      allocate(drydepositionf, mold=lwi, _STAT)
      call SulfateChemDriver(self%km, self%klid, self%CDT, MAPL_PI, real(MAPL_RADIANS_TO_DEGREES), MAPL_KARMAN, &
           MAPL_AIRMW, MAPL_AVOGAD / 1000., cpd, MAPL_GRAV, &
           fMassMSA, fMassDMS, fMassSO2, fMassSO4, &
           nymd, nhms, lons, lats, &
           DMS, SO2, SO4, dummyMSA, &
           nDMS, nSO2, nSO4, nMSA, &
           xoh, xno3, xh2o2, h2o2_init, &
           delp, t, fcld, airdens, zle0, &
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
              h2o2_init, ple0, airdens, cn_prcp, ncn_prcp, pfl_lsan0, pfi_lsan0, t, &
              nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, dummyMSA, &
              SUWT, SUPSO4, SUPSO4WT, pso4, PSO4WET, _RC)

      case ('ufs')

         do n = 1, self%nbins
            rainout_eff = 0.0
            rainout_eff(1) = self%fwet_ice(n)
            rainout_eff(2) = self%fwet_snow(n)
            rainout_eff(3) = self%fwet_rain(n)

            call MAPL_StateGetPointer(internal, int_ptr, internal_names(n), _RC)
            call WetRemovalUFS(self%km, self%klid, n, self%CDT, 'sulfate', &
                 KIN, MAPL_GRAV, self%radius(n), rainout_eff, self%washout_tuning, &
                 self%wet_radius_thr, int_ptr, ple0, t, airdens, pfl_lsan0, pfi_lsan0, SUWT, _RC)
         end do

      case default
         _ASSERT_RC(.false., 'Unsupported wet removal scheme: ' // trim(self%wet_removal_scheme), ESMF_RC_NOT_IMPL)
      end select

      ! Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
      call SU_Compute_Diags(self%km, self%klid, self%rhop(nSO4), &
           MAPL_GRAV, MAPL_PI, nSO4, self%diag_Mie, &
           self%wavelengths_profile * 1.0e-9, self%wavelengths_vertint * 1.0e-9, &
           t, airdens, delp, ple0, tropp, rh2, u, v, DMS, SO2, SO4, dummyMSA, &
           DMSSMASS, DMSCMASS, &
           MSASMASS, MSACMASS, &
           SO2SMASS, SO2CMASS, &
           SO4SMASS, SO4CMASS, &
           SUEXTTAU, SUSTEXTTAU, SUSCATAU, SUSTSCATAU, SO4MASS, SUCONC, SUEXTCOEF, &
           SUSCACOEF, SUBCKCOEF, SUANGSTR, SUFLUXU, SUFLUXV, SO4SAREA, SO4SNUM, SO4REFF, _RC)

      if (associated(SO4SAREA)) then
         nullify(int_ptr)
         call ESMF_InfoGet(aero_info, key='surface_area_density', value=fld_name, _RC)
         call ESMF_StateGet(aero, trim(fld_name), item_type, _RC)
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_StateGetPointer(aero, int_ptr, trim(fld_name), _RC)
            int_ptr = SO4SAREA
         end if
      end if

      if (associated(SO4REFF)) then ! Note unit conversion below to microns
         nullify(int_ptr)
         call ESMF_InfoGet(aero_info, key='effective_radius_in_microns', value=fld_name, _RC)
         call ESMF_StateGet(aero, trim(fld_name), item_type, _RC)
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_StateGetPointer(aero, int_ptr, trim(fld_name), _RC)
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
           tmpu=t, rhoa=airdens, delp=delp, ple=ple0, tropp=tropp, rh=RH20, u=u, v=v, &
           DMS=DMS, SO2=SO2, SO4=SO4, MSA=dummyMSA, extcoef=SUEXTCOEFRH20, &
           scacoef=SUSCACOEFRH20, _RC)

      RH80(:, :, :) = 0.80
      call SU_Compute_Diags(km=self%km, klid=self%klid, &
           rhop=self%rhop(nSO4), &
           grav=MAPL_GRAV, pi=MAPL_PI, nSO4=nSO4, mie=self%diag_Mie, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           tmpu=t, rhoa=airdens, delp=delp, ple=ple0, tropp=tropp, rh=RH80, u=u, v=v, &
           DMS=DMS, SO2=SO2, SO4=SO4, MSA=dummyMSA, extcoef=SUEXTCOEFRH80,&
           scacoef=SUSCACOEFRH80, _RC)

      deallocate(xoh, xh2o2, xno3, stat=status)

      _RETURN(_SUCCESS)

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
      real, pointer, dimension(:, :, :) :: ptr3d_int, ptr3d_imp

      integer :: status

      ! Begin...

      ! Update interal data pointers with ExtData
      call MAPL_StateGetPointer(internal, ptr3d_int, 'SO4', _RC)
      call MAPL_StateGetPointer(import, ptr3d_imp, 'climSO4', _RC)
      ptr3d_int = ptr3d_imp

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(gc)
      _UNUSED_DUMMY(export)

   end subroutine Run_data

   subroutine aerosol_optics(state, rc)

      ! ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      ! Local
      type(ESMF_Info) :: info
      real, dimension(:, :, :), pointer :: ple, ple0, rh
      real(kind=DP), dimension(:, :, :), pointer :: var
      real(kind=DP), dimension(:, :, :, :), pointer :: var4d
      real, dimension(:, :, :), pointer :: q
      real, dimension(:, :, :, :), allocatable :: q_4d
      integer, allocatable :: opaque_self(:)
      type(c_ptr) :: address
      type(SU2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: fld_name

      real(kind=DP), dimension(:, :, :), allocatable :: ext_s, ssa_s, asy_s ! (lon:,lat:,lev:)
      real(kind=DP), dimension(:, :, :, :), allocatable :: pmom_s ! (lon:,lat:,lev:,nmom:)
      integer :: instance
      integer :: nbins
      integer :: i1, j1, i2, j2, km
      integer :: band
      integer :: usePhotTable
      real :: wavelength
      integer :: k
      integer :: status

      ! Begin...

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Mie Table instance/index
      call ESMF_InfoGet(info, key='mie_table_instance', value=instance, _RC)

      ! Radiation band
      call ESMF_InfoGet(info, key='band_for_aerosol_optics', value=band, default=0, _RC)

      ! Are we doing a photolysis calculation?
      call ESMF_InfoGet(info, key='use_photolysis_table', value=usePhotTable, default=0, _RC)

      ! Pressure at layer edges, expected to be 0-based
      call ESMF_InfoGet(info, key='air_pressure_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, ple, trim(fld_name), _RC)
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple

      i1 = lbound(ple0, 1)
      i2 = ubound(ple0, 1)
      j1 = lbound(ple0, 2)
      j2 = ubound(ple0, 2)
      km = ubound(ple0, 3)

      ! Relative humidity
      call ESMF_InfoGet(info, key='relative_humidity_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, rh, trim(fld_name), _RC)

      allocate(ext_s(i1:i2, j1:j2, km), &
           ssa_s(i1:i2, j1:j2, km), &
           asy_s(i1:i2, j1:j2, km), _STAT)

      ! SO4 is the only aerosol tracer
      call ESMF_InfoGet(info, key='internal_variable_name', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, q, trim(fld_name), _RC)
      nbins = 1

      allocate(q_4d(i1:i2, j1:j2, km, nbins), _STAT)
      do k = 1, km
         q_4d(:, :, k, 1) = ((ple0(:, :, k) - ple0(:, :, k - 1)) * 0.01) * (100. / MAPL_GRAV) * q(:, :, k)
      end do

      call ESMF_InfoGet(info, key='mieTable_pointer', values=opaque_self, _RC)

      address = transfer(opaque_self, address)
      call c_f_pointer(address, self)

      if (usePhotTable /= 0) then
         wavelength = band * 1.e-9
         allocate(pmom_s(i1:i2, j1:j2, km, self%phot_Mie%nmom), _STAT)
         call miephot_(self%phot_Mie, nbins, wavelength, q_4d, rh, ext_s, ssa_s, pmom_s, _RC)
      else
         call mie_(self%rad_Mie, nbins, band, q_4d, rh, ext_s, ssa_s, asy_s, _RC)
      end if

      call ESMF_InfoGet(info, key='extinction_in_air_due_to_ambient_aerosol', value=fld_name, _RC)
      if (fld_name /= '') then
         call MAPL_StateGetPointer(state, var, trim(fld_name), _RC)
         var = ext_s(:, :, :)
      end if

      call ESMF_InfoGet(info, key='single_scattering_albedo_of_ambient_aerosol', value=fld_name, _RC)
      if (fld_name /= '') then
         call MAPL_StateGetPointer(state, var, trim(fld_name), _RC)
         var = ssa_s(:, :, :)
      end if

      if (usePhotTable /= 0) then
         call ESMF_InfoGet(info, key='legendre_coefficients_of_p11_for_photolysis', value=fld_name, _RC)
         if (fld_name /= '') then
            call MAPL_StateGetPointer(state, var4d, trim(fld_name), _RC)
            var4d = pmom_s(:, :, :, :)
         end if
      else
         call ESMF_InfoGet(info, key='asymmetry_parameter_of_ambient_aerosol', value=fld_name, _RC)
         if (fld_name /= '') then
            call MAPL_StateGetPointer(state, var, trim(fld_name), _RC)
            var = asy_s(:, :, :)
         end if
      end if

      deallocate(ext_s, ssa_s, asy_s, _STAT)
      if (usePhotTable /= 0) deallocate(pmom_s, _STAT)
      deallocate(q_4d, _STAT)

      _RETURN(_SUCCESS)

   contains

      subroutine mie_(mie, nbins, band, q, rh, bext_s, bssa_s, basym_s, rc)

         type(GOCART2G_Mie), intent(inout) :: mie ! mie table
         integer, intent(in) :: nbins ! number of bins
         integer, intent(in) :: band ! channel
         real, intent(in) :: q(:, :, :, :) ! aerosol mass mixing ratio, kg kg-1
         real, intent(in) :: rh(:, :, :) ! relative humidity
         real(kind=DP), intent(out) :: bext_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         real(kind=DP), intent(out) :: bssa_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         real(kind=DP), intent(out) :: basym_s(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3))
         integer, intent(out) :: rc
         ! local
         integer :: l
         real :: bext(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! extinction
         real :: bssa(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! SSA
         real :: gasym(size(ext_s, 1), size(ext_s, 2), size(ext_s, 3)) ! asymmetry parameter

         integer :: status

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

         _RETURN(_SUCCESS)

      end subroutine mie_

      subroutine miephot_(mie, nbins, wavelength, q, rh, bext_s, bssa_s, bpmom_s, rc)

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

         _RETURN(_SUCCESS)

      end subroutine miephot_

   end subroutine aerosol_optics

   subroutine monochromatic_aerosol_optics(state, rc)

      ! ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      ! Local
      type(ESMF_Info) :: info
      real, dimension(:, :, :), pointer :: ple, ple0, rh
      real, dimension(:, :), pointer :: var
      real, dimension(:, :, :), pointer :: q
      real, dimension(:, :, :, :), allocatable :: q_4d
      integer, allocatable :: opaque_self(:)
      type(c_ptr) :: address
      type(SU2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: fld_name

      real, dimension(:, :, :), allocatable :: tau_s, tau ! (lon:,lat:,lev:)
      integer :: instance
      integer :: n, nbins
      integer :: i1, j1, i2, j2, km
      real :: wavelength
      integer :: k
      integer :: status

      ! Begin...

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Mie Table instance/index
      call ESMF_InfoGet(info, key='mie_table_instance', value=instance, _RC)

      ! Radiation band
      call ESMF_InfoGet(info, key='wavelength_for_aerosol_optics', value=wavelength, default=0., _RC)

      ! Pressure at layer edges, expected to be 0-based
      call ESMF_InfoGet(info, key='air_pressure_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, ple, trim(fld_name), _RC)
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple

      i1 = lbound(ple0, 1)
      i2 = ubound(ple0, 1)
      j1 = lbound(ple0, 2)
      j2 = ubound(ple0, 2)
      km = ubound(ple0, 3)

      ! Relative humidity
      call ESMF_InfoGet(info, key='relative_humidity_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, rh, trim(fld_name), _RC)

      allocate(tau_s(i1:i2, j1:j2, km), &
           tau(i1:i2, j1:j2, km), _STAT)
      tau_s = 0.0
      tau = 0.0

      ! SO4 is the only aerosol tracer
      call ESMF_InfoGet(info, key='internal_variable_name', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, q, trim(fld_name), _RC)
      nbins = 1

      allocate(q_4d(i1:i2, j1:j2, km, nbins), _STAT)
      do k = 1, km
         q_4d(:, :, k, 1) = (ple0(:, :, k) - ple0(:, :, k - 1)) / MAPL_GRAV * q(:, :, k)
      end do

      call ESMF_InfoGet(info, key='mieTable_pointer', values=opaque_self, _RC)

      address = transfer(opaque_self, address)
      call c_f_pointer(address, self)

      do n = 1, nbins
         call self%diag_Mie%Query(wavelength, n, q_4d(:, :, :, n), rh, tau=tau, _RC)
         tau_s = tau_s + tau
      end do

      call ESMF_InfoGet(info, key='monochromatic_extinction_in_air_due_to_ambient_aerosol', value=fld_name, _RC)
      if (fld_name /= '') then
         call MAPL_StateGetPointer(state, var, trim(fld_name), _RC)
         var = sum(tau_s, dim=3)
      end if

      deallocate(q_4d, _STAT)

      _RETURN(_SUCCESS)

   end subroutine monochromatic_aerosol_optics

end module SU2G_GridCompMod

subroutine SU2G_SetServices(gc, rc)
   use ESMF
   use SU2G_GridCompMod, only : mySetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine SU2G_SetServices
