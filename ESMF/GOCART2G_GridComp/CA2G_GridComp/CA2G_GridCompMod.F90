#include "MAPL.h"

!MODULE: CA2G_GridCompMod - GOCART refactoring of the Carbonaceous Aerosol (CA) gridded component

!INTERFACE:
module CA2G_GridCompMod

   !USES:
   use ESMF
   use MAPL
   use GOCART2G_MieMod
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process ! GOCART2G process library
   use GA_EnvironmentMod
   !$ use omp_lib

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData = 2

   integer, parameter :: DP = kind(1.0d0)

   ! Number of carbonaceous aerosol bins (phobic, philic)
   integer, parameter :: n_ca_bins = 2

   real, parameter :: cpd = 1004.16

   ! MAPL2's MAPL_UNDEF; MAPL3's is -huge, which breaks the 1.01*x > undef masks
   real, parameter :: MAPL2_UNDEF = 1.0e15

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices

   !DESCRIPTION: This module implements GOCART's Carbonaceous Aerosol (CA) Gridded Component.

   !REVISION HISTORY:
   ! 04January2024  Collow - Update to ChemSettling Call
   ! 15June2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

   ! Carbonaceous aerosol state
   type :: ThreadWorkspace
      ! Point emissions of carbonaceous aerosol
      integer :: nPts = -1
      integer, allocatable, dimension(:) :: pstart, pend
      real, allocatable, dimension(:) :: pLat, pLon, pBase, pTop, pEmis
   end type ThreadWorkspace

   type, extends(GA_Environment) :: CA2G_GridComp
      integer :: myDOW = -1 ! my Day of the week: Sun=1, Mon=2,...,Sat=7
      real :: ratPOM = 1.0 ! Ratio of POM to OC mass
      real :: fMonoterpenes = 0.0 ! Fraction of monoterpene emissions -> aerosol
      real :: fIsoprene = 0.0 ! Fraction of isoprene emissions -> aerosol
      real :: fHydrophobic ! Initially hydrophobic portion
      real :: tConvPhobicToPhilic ! e-folding time [days] hydrophobic to hydrophilic
      real, allocatable :: tChemLoss(:) ! e-folding time [days] for parameterized chemistry loss
      logical :: diurnal_bb ! diurnal biomass burning
      real :: eAircraftFuel ! Aircraft emission factor: go from kg fuel to kg C
      real, allocatable :: aviation_layers(:) ! heights of the LTO, CDS and CRS layers

      ! Workspace for point emissions
      logical :: doing_point_emissions = .false.
      character(len=:), allocatable :: point_emissions_srcfilen ! filename for pointwise emissions
      type(ThreadWorkspace), allocatable :: workspaces(:)
   end type CA2G_GridComp

   character(*), parameter :: PRIVATE_STATE = "CA2G_GridComp"

contains

   !IROUTINE: SetServices
   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, intent(out) :: rc ! return code

      !DESCRIPTION: Declares CA's private state, entry points (Initialize and
      ! the Run0/Run1/Run2 phases) and state specs.

      !REVISION HISTORY:
      ! 15June2020   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

      ! Locals
      character(len=:), allocatable :: comp_name
      type(CA2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: field_name
      logical :: data_driven = .true.
      integer :: num_threads
      type(MAPL_UngriddedDim) :: ungrd_nbins
      type(MAPL_UngriddedDim) :: ungrd_wavelengths_profile, ungrd_wavelengths_vertint
      integer :: i, status

      call MAPL_GridCompGet(gc, name=comp_name, _RC)

      ! Wrap gridcomp's private state and store it in gridcomp
      _SET_NAMED_PRIVATE_STATE(gc, CA2G_GridComp, PRIVATE_STATE)

      ! Retrieve the private state
      _GET_NAMED_PRIVATE_STATE(gc, CA2G_GridComp, PRIVATE_STATE, self)

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

      ! process CA-specific items
      call MAPL_GridCompGetResource(gc, "my_day_of_week", self%myDOW, default=-1, _RC)
      call MAPL_GridCompGetResource(gc, "hydrophobic_fraction", self%fHydrophobic, _RC)
      call MAPL_GridCompGetResource(gc, &
           "time_days_hydrophobic_to_hydrophilic", self%tConvPhobicToPhilic, default=2.5, _RC)
      call MAPL_GridCompGetResource(gc, "time_days_chemical_destruction", self%tChemLoss, _RC)
      _ASSERT(size(self%tChemLoss) == self%nbins, "time_days_chemical_destruction must have nbins values")
      call MAPL_GridCompGetResource(gc, "pom_ca_ratio", self%ratPOM, default=1.0, _RC)
      call MAPL_GridCompGetResource(gc, "monoterpenes_emission_fraction", self%fMonoterpenes, default=0.0, _RC)
      call MAPL_GridCompGetResource(gc, "isoprene_emission_fraction", self%fIsoprene, default=0.0, _RC)
      call MAPL_GridCompGetResource(gc, "aircraft_fuel_emission_factor", self%eAircraftFuel, _RC)
      call MAPL_GridCompGetResource(gc, "aviation_vertical_layers", self%aviation_layers, _RC)
      _ASSERT(size(self%aviation_layers) == 4, "aviation_vertical_layers must have 4 values")

      call MAPL_GridCompGetResource(gc, &
           "point_emissions_srcfilen", self%point_emissions_srcfilen, default='/dev/null', _RC)
      if ((index(self%point_emissions_srcfilen, '/dev/null') > 0)) then
         self%doing_point_emissions = .false. ! disable it if no file specified
      else
         self%doing_point_emissions = .true. ! we are good to go
      end if

      ! Is CA data driven?
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
              short_name=trim(comp_name) // 'phobic', &
              standard_name='Carbonaceous aerosol, hydrophobic mode', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              fill_value=0.0, &
              add_to_export=.true., _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_INTERNAL, &
              short_name=trim(comp_name) // 'philic', &
              standard_name='Carbonaceous aerosol, hydrophilic mode', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              fill_value=0.0, &
              add_to_export=.true., _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='clim' // trim(comp_name) // 'phobic', &
              standard_name='Carbonaceous aerosol, hydrophobic mode', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='clim' // trim(comp_name) // 'philic', &
              standard_name='Carbonaceous aerosol, hydrophilic mode', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, _RC)

         do i = 1, self%nbins
            write(field_name, '(A, I0.3)') '', i
            ! dry deposition
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climCADP' // trim(field_name), &
                 standard_name='Carbon dry deposition (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)

            ! wet deposition
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climCAWT' // trim(field_name), &
                 standard_name='Carbon wet deposition (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)

            ! gravitational settling
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climCASD' // trim(field_name), &
                 standard_name='Carbon settling (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)

            ! convective scavenging
            call MAPL_GridCompAddSpec(gc, &
                 state_intent=ESMF_STATEINTENT_IMPORT, &
                 short_name='climCASV' // trim(field_name), &
                 standard_name='Carbon convective scavenging (bin ' // trim(field_name) // ')', &
                 units='kg kg-1 s-1', &
                 dims='xy', &
                 vertical_stagger=MAPL_VERTICAL_STAGGER_NONE, &
                 restart_mode=MAPL_RESTART_SKIP, _RC)
         end do
      end if ! (data_driven)

      ! Import, Export, Internal states for computational instance
      if (.not. data_driven) then
#include "CA2G_Export___.h"
#include "CA2G_Import___.h"
#include "CA2G_Internal___.h"
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
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code

      !DESCRIPTION: This initializes CA's Grid Component. It primarily fills
      ! GOCART's AERO states with its carbonaceous aerosol fields.

      !REVISION HISTORY:
      ! 15June2020   E.Sherman  First attempt at refactoring

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_Geom) :: geom
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(ESMF_State) :: provider_state
      type(ESMF_FieldBundle) :: bundle_dp
      type(ESMF_Field) :: field
      type(ESMF_Info) :: field_info, aero_info
      type(CA2G_GridComp), pointer :: self

      integer, allocatable :: mieTable_pointer(:)
      integer :: km
      integer :: instance
      character(len=ESMF_MAXSTR) :: prefix
      real :: CDT ! chemistry timestep (secs)
      real :: HDT ! model timestep (secs)
      logical :: data_driven

      integer, allocatable, dimension(:) :: channels_
      integer :: nmom_
      character(len=:), allocatable :: file_
      integer :: i
      character(len=3) :: bin_index
      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, geom=geom, num_levels=km, _RC)

      ! Get my internal private state
      _GET_NAMED_PRIVATE_STATE(gc, CA2G_GridComp, PRIVATE_STATE, self)

      self%km = km

      ! Get DTs
      call MAPL_ClockGet(clock, dt=HDT, _RC)
      call MAPL_GridCompGetResource(gc, "GOCART2G_DT", CDT, default=real(HDT), _RC)
      self%CDT = CDT

      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is CA data driven?
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

      ! Fill AERO State with phobic/philic fields
      call ESMF_StateGet(export, trim(comp_name) // '_AERO', aero, _RC)
      call ESMF_StateGet(export, trim(comp_name) // '_AERO_DP', bundle_dp, _RC)

      call ESMF_StateGet(internal, trim(comp_name) // 'phobic', field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(1), _RC)
      call ESMF_StateAdd(aero, [field], _RC)

      call ESMF_StateGet(internal, trim(comp_name) // 'philic', field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(2), _RC)
      call ESMF_StateAdd(aero, [field], _RC)

      if (data_driven) then
         instance = instanceData

         do i = 1, self%nbins
            write(bin_index, '(A, I0.3)') '', i
            ! Dry deposition
            call append_to_bundle('CA' // 'DP' // trim(bin_index), provider_state, prefix, bundle_dp, _RC)

            ! Wet deposition (Convective scavenging)
            call append_to_bundle('CA' // 'SV' // trim(bin_index), provider_state, prefix, bundle_dp, _RC)

            ! Wet deposition
            call append_to_bundle('CA' // 'WT' // trim(bin_index), provider_state, prefix, bundle_dp, _RC)

            ! Gravitational Settling
            call append_to_bundle('CA' // 'SD' // trim(bin_index), provider_state, prefix, bundle_dp, _RC)
         end do
      else
         instance = instanceComputational

         ! Dry deposition
         call append_to_bundle(trim(comp_name) // 'DP', provider_state, prefix, bundle_dp, _RC)

         ! Wet deposition (Convective scavenging)
         call append_to_bundle(trim(comp_name) // 'SV', provider_state, prefix, bundle_dp, _RC)

         ! Wet deposition
         call append_to_bundle(trim(comp_name) // 'WT', provider_state, prefix, bundle_dp, _RC)

         ! Gravitational Settling
         call append_to_bundle(trim(comp_name) // 'SD', provider_state, prefix, bundle_dp, _RC)
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

      ! Add variables to CA instance's aero state. This is used in aerosol optics calculations
      call add_aero(aero, label="air_pressure_for_aerosol_optics", label2="PLE", geom=geom, km=self%km, _RC)
      call add_aero(aero, label="relative_humidity_for_aerosol_optics", label2="RH", geom=geom, km=self%km, _RC)

      ! Add variables to CA instance aero state for chemistry
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

      ! A single name identifying this species; the phobic/philic bin
      ! field names are reconstructed at runtime in aerosol_optics/get_mixR
      ! by appending 'phobic'/'philic' to this scalar (see Chem_AeroGeneric's
      ! get_mixR, which does the same).
      call ESMF_InfoSet(aero_info, key="internal_variable_name", value=trim(comp_name), _RC)

      call ESMF_MethodAdd(aero, label="aerosol_optics", userRoutine=aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label="monochromatic_aerosol_optics", userRoutine=monochromatic_aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label="get_mixR", userRoutine=get_mixR, _RC)

      _RETURN(_SUCCESS)

   end subroutine Initialize

   !IROUTINE: Run0
   !INTERFACE:
   subroutine Run0(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Clears klid to 0.0 for Carbon

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      type(CA2G_GridComp), pointer :: self
      real, pointer, dimension(:, :, :) :: ple, ple0
      real, pointer, dimension(:, :, :) :: ptr3d_int
      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, CA2G_GridComp, PRIVATE_STATE, self)

      ! Edge variable PLE is expected to be 0-based
      call MAPL_StateGetPointer(import, ple, 'PLE', _RC)
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple

      ! Set klid and Set internal values to 0 above klid
      call findKlid(self%klid, self%plid, ple0, _RC)
      call MAPL_StateGetPointer(internal, ptr3d_int, trim(comp_name) // 'phobic', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, trim(comp_name) // 'philic', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(export)
      _UNUSED_DUMMY(clock)

   end subroutine Run0

   !IROUTINE: Run
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run method for the Carbon Component. Determines whether to run
      ! data or computational run method.

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      logical :: data_driven
      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is CA data driven?
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
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Computes emissions/sources for Carbon

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      type(ESMF_Grid) :: grid
      type(ESMF_Geom) :: geom
      type(CA2G_GridComp), pointer :: self
      type(ESMF_Time) :: time

      character(len=2) :: GCsuffix
      character(len=3) :: cdow
      integer :: idow
      integer :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
      real, allocatable, dimension(:, :) :: lats
      real, allocatable, dimension(:, :) :: lons
      real, dimension(:, :, :), allocatable :: aircraft_fuel_src
      real, dimension(:, :), allocatable :: biomass_src, biofuel_src, biogvoc_src, &
           eocant1_src, eocant2_src, oc_ship_src, aviation_lto_src, aviation_cds_src, &
           aviation_crs_src, biomass_src_
      integer, dimension(:), allocatable :: iPoint, jPoint
      real, dimension(:, :, :), allocatable :: emissions_point
      character(len=ESMF_MAXSTR) :: fname ! file name for point source emissions
      logical :: fileExists

      real, pointer, dimension(:, :, :) :: intPtr_phobic, intPtr_philic
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread

#include "CA2G_DeclarePointer___.h"

      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, geom=geom, grid=grid, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, CA2G_GridComp, PRIVATE_STATE, self)

#include "CA2G_GetPointer___.h"

      if (trim(comp_name) == 'CA.oc') then
         GCsuffix = 'OC'
      else if (trim(comp_name) == 'CA.bc') then
         GCsuffix = 'BC'
      else if (trim(comp_name) == 'CA.br') then
         GCsuffix = 'BR'
      end if

      call MAPL_StateGetPointer(internal, intPtr_phobic, trim(comp_name) // 'phobic', _RC)
      call MAPL_StateGetPointer(internal, intPtr_philic, trim(comp_name) // 'philic', _RC)

      allocate(emissions_point, mold=delp, _STAT)
      emissions_point = 0.0

      ! Extract nymd(yyyymmdd) from clock
      call ESMF_ClockGet(clock, currTime=time, _RC)
      call ESMF_TimeGet(time, YY=iyr, mm=imm, DD=idd, H=ihr, M=imn, S=isc, _RC)
      nymd = MAPL_PackedDateCreate(iyr, imm, idd)
      nhms = MAPL_PackedTimeCreate(ihr, imn, isc)

      ! Reset tracer to zero at 0Z on specific day of week
      idow = Chem_UtilIdow(nymd)
      if ((nhms == 0) .and. (idow == self%myDOW)) then
         cdow = Chem_UtilCdow(nymd)
         intPtr_phobic = tiny(1.) ! avoid division by zero
         intPtr_philic = tiny(1.) ! avoid division by zero
         if (MAPL_Am_I_Root()) then
            !$omp critical (CA2G_1)
            print *, '<> CA ' // cdow // ' tracer being set to zero on ', nymd, nhms
            !$omp end critical (CA2G_1)
         end if
      end if

      ! Implicit allocation with Fortran 2003
      if (trim(comp_name) == 'CA.oc') then
         aircraft_fuel_src = OC_AIRCRAFT
         biomass_src = OC_BIOMASS
         biofuel_src = OC_BIOFUEL
         eocant1_src = OC_ANTEOC1
         eocant2_src = OC_ANTEOC2
         oc_ship_src = OC_SHIP
         aviation_lto_src = OC_AVIATION_LTO
         aviation_cds_src = OC_AVIATION_CDS
         aviation_crs_src = OC_AVIATION_CRS
         allocate(biogvoc_src, mold=OC_MTPA, _STAT)
         biogvoc_src = ((OC_MTPA + OC_MTPO + OC_LIMO) * self%fMonoterpenes) + (OC_ISOPRENE * self%fIsoprene)
      else if (trim(comp_name) == 'CA.bc') then
         aircraft_fuel_src = BC_AIRCRAFT
         biomass_src = BC_BIOMASS
         biofuel_src = BC_BIOFUEL
         eocant1_src = BC_ANTEBC1
         eocant2_src = BC_ANTEBC2
         oc_ship_src = BC_SHIP
         aviation_lto_src = BC_AVIATION_LTO
         aviation_cds_src = BC_AVIATION_CDS
         aviation_crs_src = BC_AVIATION_CRS
         ! Black carbon has no biogvoc_src, so we set it to zero.
         ! biogvoc_src is still needed for the call to CAEmission, however it
         ! effectively does nothing since we set all its values to zero.
         allocate(biogvoc_src, mold=BC_BIOMASS, _STAT)
         biogvoc_src = 0.0
      else if (trim(comp_name) == 'CA.br') then
         aircraft_fuel_src = BRC_AIRCRAFT
         biomass_src = BRC_BIOMASS
         biogvoc_src = BRC_TERPENE
         biofuel_src = BRC_BIOFUEL
         eocant1_src = BRC_ANTEBRC1
         eocant2_src = BRC_ANTEBRC2
         oc_ship_src = BRC_SHIP
         aviation_lto_src = BRC_AVIATION_LTO
         aviation_cds_src = BRC_AVIATION_CDS
         aviation_crs_src = BRC_AVIATION_CRS
      end if

      ! As a safety check, where value is undefined set to 0
      where (1.01 * biomass_src > MAPL2_UNDEF) biomass_src = 0.
      where (1.01 * biogvoc_src > MAPL2_UNDEF) biogvoc_src = 0.
      where (1.01 * biofuel_src > MAPL2_UNDEF) biofuel_src = 0.
      where (1.01 * eocant1_src > MAPL2_UNDEF) eocant1_src = 0.
      where (1.01 * eocant2_src > MAPL2_UNDEF) eocant2_src = 0.
      where (1.01 * oc_ship_src > MAPL2_UNDEF) oc_ship_src = 0.
      where (1.01 * aircraft_fuel_src > MAPL2_UNDEF) aircraft_fuel_src = 0.
      where (1.01 * aviation_lto_src > MAPL2_UNDEF) aviation_lto_src = 0.
      where (1.01 * aviation_cds_src > MAPL2_UNDEF) aviation_cds_src = 0.
      where (1.01 * aviation_crs_src > MAPL2_UNDEF) aviation_crs_src = 0.

      ! Apply diurnal cycle if so desired
      if (self%diurnal_bb) then
         call MAPL_GridGetCoordinates(grid, latitudes=lats, longitudes=lons, _RC)
         biomass_src_ = biomass_src
         call Chem_BiomassDiurnal(biomass_src, biomass_src_, &
              lons(:, :) * real(MAPL_RADIANS_TO_DEGREES), &
              lats(:, :) * real(MAPL_RADIANS_TO_DEGREES), &
              nhms, self%cdt)
      end if

      call CAEmission(self%diag_Mie, self%km, self%nbins, self%cdt, &
           MAPL_GRAV, GCsuffix, self%ratPOM, &
           self%eAircraftFuel, aircraft_fuel_src, &
           aviation_lto_src, aviation_cds_src, &
           aviation_crs_src, self%fHydrophobic, zpbl, t, airdens, rh2, &
           intPtr_philic, intPtr_phobic, delp, self%aviation_layers, biomass_src, &
           biogvoc_src, eocant1_src, eocant2_src, oc_ship_src, biofuel_src, &
           EM, EMAN, EMBB, EMBF, EMBG, _RC)

      ! Read any pointwise emissions, if requested
      thread = MAPL_get_current_thread()
      workspace => self%workspaces(thread)
      if (self%doing_point_emissions) then
         call MAPL_StrTemplate(fname, self%point_emissions_srcfilen, xid='unknown', &
              nymd=nymd, nhms=120000)
         inquire(file=fname, exist=fileExists)
         if (fileExists) then
            call ReadPointEmissions(nymd, fname, workspace%nPts, workspace%pLat, workspace%pLon, &
                 workspace%pBase, workspace%pTop, workspace%pEmis, workspace%pstart, &
                 workspace%pend, label='source', _RC)
         else if (.not. fileExists) then
            !$omp critical (CA2G_2)
            if (mapl_am_i_root()) print *, 'GOCART2G ', trim(comp_name), ': ', trim(fname), ' not found; proceeding.'
            !$omp end critical (CA2G_2)
            workspace%nPts = -1 ! set this back to -1 so the "if (workspace%nPts > 0)" conditional is not exercised.
         end if
      end if

      ! Get indices for point emissions
      if (workspace%nPts > 0) then
         call MAPL_GeomGetHorzIJIndex(geom, &
              workspace%pLon / real(MAPL_RADIANS_TO_DEGREES), &
              workspace%pLat / real(MAPL_RADIANS_TO_DEGREES), &
              iPoint, jPoint, local_indices=.true., _RC)

         call updatePointwiseEmissions( &
              self%km, workspace%pBase, workspace%pTop, workspace%pEmis, workspace%nPts, &
              workspace%pstart, workspace%pend, zle, &
              area, iPoint, jPoint, nhms, emissions_point, _RC)

         deallocate(iPoint, jPoint, _STAT)
         intPtr_phobic = intPtr_phobic + self%fHydrophobic * self%cdt * MAPL_GRAV / delp * emissions_point
         intPtr_philic = intPtr_philic + (1 - self%fHydrophobic) * self%cdt * MAPL_GRAV / delp * emissions_point
      end if

      _RETURN(_SUCCESS)

   end subroutine Run1

   !IROUTINE: Run2
   !INTERFACE:
   subroutine Run2(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run2 method for the Carbon Grid Component.

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(ESMF_Info) :: aero_info
      type(CA2G_GridComp), pointer :: self

      ! Internal state order, as in CA2G_StateSpecs.rc (phobic, philic)
      character(len=ESMF_MAXSTR) :: internal_names(n_ca_bins)

      character(len=2) :: GCsuffix
      integer :: n
      logical :: KIN
      real, allocatable, dimension(:, :) :: drydepositionfrequency, dqa
      real, pointer, dimension(:, :, :) :: casd_vel
      real, dimension(3) :: rainout_eff
      real, allocatable, dimension(:, :, :) :: pSOA_VOC
      real, pointer, dimension(:, :, :) :: int_ptr
      real, allocatable, dimension(:, :, :, :) :: int_arr
      character(len=ESMF_MAXSTR) :: fld_name
      real, pointer, dimension(:, :, :) :: intPtr_phobic, intPtr_philic
      real, pointer, dimension(:, :) :: flux_ptr

      integer :: i1, j1, i2, j2, km
      real, target, allocatable, dimension(:, :, :) :: RH20, RH80
      integer :: settling_opt
      type(ESMF_StateItem_Flag) :: item_type
      real, pointer, dimension(:, :, :) :: ple0, zle0, pfl_lsan0, pfi_lsan0

#include "CA2G_DeclarePointer___.h"

      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, CA2G_GridComp, PRIVATE_STATE, self)

      if (trim(comp_name) == 'CA.oc') then
         GCsuffix = 'OC'
      else if (trim(comp_name) == 'CA.bc') then
         GCsuffix = 'BC'
      else if (trim(comp_name) == 'CA.br') then
         GCsuffix = 'BR'
      end if

      internal_names(1) = trim(comp_name) // 'phobic'
      internal_names(2) = trim(comp_name) // 'philic'

      ! Get the aero state
      call ESMF_StateGet(export, trim(comp_name) // '_AERO', aero, _RC)
      call ESMF_InfoGetFromHost(aero, aero_info, _RC)

#include "CA2G_GetPointer___.h"

      call MAPL_StateGetPointer(internal, intPtr_phobic, trim(comp_name) // 'phobic', _RC)
      call MAPL_StateGetPointer(internal, intPtr_philic, trim(comp_name) // 'philic', _RC)

      ! Edge variables are expected to be 0-based
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple
      zle0(lbound(zle, 1):, lbound(zle, 2):, 0:) => zle
      pfl_lsan0(lbound(pfl_lsan, 1):, lbound(pfl_lsan, 2):, 0:) => pfl_lsan
      pfi_lsan0(lbound(pfi_lsan, 1):, lbound(pfi_lsan, 2):, 0:) => pfi_lsan

      ! Set klid and Set internal values to 0 above klid
      call findKlid(self%klid, self%plid, ple0, _RC)
      call setZeroKlid(self%km, self%klid, intPtr_phobic)
      call setZeroKlid(self%km, self%klid, intPtr_philic)

      ! Add on SOA from Anthropogenic/Biomass-burning VOC oxidation
      if (trim(comp_name) == 'CA.oc') then
         pSOA_VOC = pSOA_ANTHRO_VOC
         where (1.01 * pSOA_VOC > MAPL2_UNDEF) pSOA_VOC = 0.0

         intPtr_philic = intPtr_philic + self%cdt * pSOA_VOC / airdens
         if (associated(PSOA)) &
              PSOA = sum(pSOA_VOC * delp / airdens / MAPL_GRAV, 3)
      end if

      if (trim(comp_name) == 'CA.br') then
         pSOA_VOC = pSOA_BIOB_VOC
         where (1.01 * pSOA_VOC > MAPL2_UNDEF) pSOA_VOC = 0.0

         intPtr_philic = intPtr_philic + self%cdt * pSOA_VOC / airdens
         if (associated(PSOA)) &
              PSOA = sum(pSOA_VOC * delp / airdens / MAPL_GRAV, 3)
      end if

      ! Ad Hoc transfer of hydrophobic to hydrophilic aerosols
      ! Rate controlled in RC file; tConvPhobicToPhilic < 0 means no transfer
      call phobicToPhilic(intPtr_phobic, intPtr_philic, HYPHIL, &
           self%tConvPhobicToPhilic, self%km, self%cdt, MAPL_GRAV, delp, _RC)

      ! Ad Hoc chemical destruction of carbon
      ! This applies a simple exponential decay to both hydrophobic and
      ! hydrophilic modes with the time constant tChemLoss (e-folding
      ! time in days)
      do n = 1, self%nbins
         call MAPL_StateGetPointer(internal, int_ptr, internal_names(n), _RC)
         call carbonChemLoss(self%km, self%klid, n, self%cdt, MAPL_GRAV, delp, &
              self%tChemLoss(n), int_ptr, CH, _RC)
      end do

      ! CA Settling
      select case (self%settling_scheme)
      case ('gocart')
         settling_opt = 1
      case ('ufs')
         settling_opt = 2
      case default
         _ASSERT_RC(.false., 'Unsupported settling scheme: ' // trim(self%settling_scheme), ESMF_RC_NOT_IMPL)
      end select

      do n = 1, self%nbins
         call MAPL_StateGetPointer(internal, int_ptr, internal_names(n), _RC)
         nullify(flux_ptr)
         if (associated(SD)) flux_ptr => SD(:, :, n)
         nullify(casd_vel)
         if (associated(SD_V)) casd_vel => SD_V(:, :, :, n)
         call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, n, self%cdt, MAPL_GRAV, &
              int_ptr, t, airdens, &
              rh2, zle0, delp, flux_ptr, casd_vel, settling_scheme=settling_opt, _RC)
      end do

      ! CA Deposition
      allocate(dqa, mold=lwi, _STAT)
      allocate(drydepositionfrequency, mold=lwi, _STAT)

      drydepositionfrequency = 0.
      call DryDeposition(self%km, t, airdens, zle0, lwi, ustar, zpbl, sh, &
           MAPL_KARMAN, cpd, MAPL_GRAV, z0h, drydepositionfrequency, _RC)

      do n = 1, self%nbins
         call MAPL_StateGetPointer(internal, int_ptr, internal_names(n), _RC)
         dqa = 0.
         dqa = max(0.0, int_ptr(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%cdt)))
         int_ptr(:, :, self%km) = int_ptr(:, :, self%km) - dqa
         if (associated(DP)) then
            DP(:, :, n) = dqa * delp(:, :, self%km) / MAPL_GRAV / self%cdt
         end if
      end do

      ! Large-scale Wet Removal
      KIN = .true.

      select case (self%wet_removal_scheme)

      case ('gocart')
         ! Hydrophobic mode (first tracer) is not removed
         if (associated(WT)) WT(:, :, 1) = 0.0

         ! Hydrophilic mode (second tracer) is removed
         call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 2, self%cdt, GCsuffix, &
              KIN, MAPL_GRAV, self%fwet(2), philic, ple0, t, airdens, &
              pfl_lsan0, pfi_lsan0, cn_prcp, ncn_prcp, WT, _RC)

      case ('ufs')
         ! Both hydrophobic and hydrophilic modes can be removed
         do n = 1, self%nbins
            rainout_eff = 0.0
            rainout_eff(1) = self%fwet_ice(n) ! remove with ice
            rainout_eff(2) = self%fwet_snow(n) ! remove with snow
            rainout_eff(3) = self%fwet_rain(n) ! remove with rain

            call MAPL_StateGetPointer(internal, int_ptr, internal_names(n), _RC)
            call WetRemovalUFS(self%km, self%klid, n, self%cdt, GCsuffix, &
                 KIN, MAPL_GRAV, self%radius(n), rainout_eff, self%washout_tuning, &
                 self%wet_radius_thr, int_ptr, ple0, t, airdens, pfl_lsan0, pfi_lsan0, WT, _RC)
         end do
      case default
         _ASSERT_RC(.false., 'Unsupported wet removal scheme: ' // trim(self%wet_removal_scheme), ESMF_RC_NOT_IMPL)
      end select

      ! Compute diagnostics
      allocate(int_arr(ubound(intPtr_phobic, 1), ubound(intPtr_phobic, 2), self%km, self%nbins), _STAT)
      int_arr(:, :, :, 1) = intPtr_phobic
      int_arr(:, :, :, 2) = intPtr_philic

      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, nbins=2, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, aerosol=int_arr, grav=MAPL_GRAV, &
           tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, delp=delp, ple=ple0, tropp=tropp, &
           sfcmass=SMASS, colmass=CMASS, mass=MASS, &
           exttau=EXTTAU, stexttau=STEXTTAU, scatau=SCATAU, stscatau=STSCATAU, &
           fluxu=FLUXU, fluxv=FLUXV, &
           conc=CONC, extcoef=EXTCOEF, scacoef=SCACOEF, bckcoef=BCKCOEF, angstrom=ANGSTR, &
           aerindx=AERIDX, NO3nFlag=.false., SAREA=SAREA, REFF=REFF, _RC)

      if (associated(SAREA)) then
         nullify(int_ptr)
         call ESMF_InfoGet(aero_info, key='surface_area_density', value=fld_name, _RC)
         call ESMF_StateGet(aero, trim(fld_name), item_type, _RC)
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_StateGetPointer(aero, int_ptr, trim(fld_name), _RC)
            int_ptr = SAREA
         end if
      end if

      if (associated(REFF)) then ! Note unit conversion below to microns
         nullify(int_ptr)
         call ESMF_InfoGet(aero_info, key='effective_radius_in_microns', value=fld_name, _RC)
         call ESMF_StateGet(aero, trim(fld_name), item_type, _RC)
         if (item_type == ESMF_STATEITEM_FIELD) then
            call MAPL_StateGetPointer(aero, int_ptr, trim(fld_name), _RC)
            int_ptr = REFF * 1.e6
         end if
      end if

      i1 = lbound(rh2, 1); i2 = ubound(rh2, 1)
      j1 = lbound(rh2, 2); j2 = ubound(rh2, 2)
      km = ubound(rh2, 3)

      allocate(RH20(i1:i2, j1:j2, km), _STAT)
      allocate(RH80(i1:i2, j1:j2, km), _STAT)

      RH20(:, :, :) = 0.20
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, nbins=2, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, aerosol=int_arr, grav=MAPL_GRAV, &
           tmpu=t, rhoa=airdens, rh=RH20, u=u, v=v, delp=delp, ple=ple0, tropp=tropp, &
           extcoef=EXTCOEFRH20, scacoef=SCACOEFRH20, NO3nFlag=.false., _RC)

      RH80(:, :, :) = 0.80
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, nbins=2, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, aerosol=int_arr, grav=MAPL_GRAV, &
           tmpu=t, rhoa=airdens, rh=RH80, u=u, v=v, delp=delp, ple=ple0, tropp=tropp, &
           extcoef=EXTCOEFRH80, scacoef=SCACOEFRH80, NO3nFlag=.false., _RC)

      deallocate(RH20, RH80, _STAT)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(clock)

   end subroutine Run2

   !IROUTINE: Run_data -- ExtData Carbon Grid Component
   !INTERFACE:
   subroutine Run_data(gc, import, export, internal, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_State), intent(inout) :: internal ! Interal state
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION: Updates pointers in Internal state with fields from ExtData.

      ! Locals
      character(len=:), allocatable :: comp_name
      character(len=2) :: GCsuffix
      real, pointer, dimension(:, :, :) :: ptr3d_int, ptr3d_imp
      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, _RC)

      if (trim(comp_name) == 'CA.oc') then
         GCsuffix = 'OC'
      else if (trim(comp_name) == 'CA.bc') then
         GCsuffix = 'BC'
      else if (trim(comp_name) == 'CA.br') then
         GCsuffix = 'BR'
      end if

      ! Update internal data pointers with ExtData
      call MAPL_StateGetPointer(internal, ptr3d_int, trim(comp_name) // 'phobic', _RC)
      call MAPL_StateGetPointer(import, ptr3d_imp, 'clim' // trim(GCsuffix) // 'phobic', _RC)
      ptr3d_int = ptr3d_imp

      call MAPL_StateGetPointer(internal, ptr3d_int, trim(comp_name) // 'philic', _RC)
      call MAPL_StateGetPointer(import, ptr3d_imp, 'clim' // trim(GCsuffix) // 'philic', _RC)
      ptr3d_int = ptr3d_imp

      _RETURN(_SUCCESS)
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
      type(CA2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: fld_name

      ! CA's two bins (phobic, philic) are plain named fields derived at
      ! runtime from the scalar 'internal_variable_name' attribute (see
      ! Initialize), rather than an ungridded array.
      character(len=:), allocatable :: comp_name
      character(len=ESMF_MAXSTR) :: aerosol_names(n_ca_bins)

      real(kind=DP), dimension(:, :, :), allocatable :: ext_s, ssa_s, asy_s ! (lon:,lat:,lev:)
      real(kind=DP), dimension(:, :, :, :), allocatable :: pmom_s ! (lon:,lat:,lev:,nmom:)
      integer :: instance
      integer :: i1, j1, i2, j2, km
      integer :: band
      integer :: usePhotTable
      real :: wavelength
      integer :: n, k
      integer :: status

      ! Begin...

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Mie Table instance/index
      call ESMF_InfoGet(info, key='mie_table_instance', value=instance, _RC)

      ! Radiation band
      call ESMF_InfoGet(info, key='band_for_aerosol_optics', value=band, default=0, _RC)

      ! Are we doing a photolysis calculation?
      call ESMF_InfoGet(info, key='use_photolysis_table', value=usePhotTable, default=0, _RC)

      ! Build the two bin field names from the scalar component name
      call ESMF_InfoGetCharAlloc(info, key='internal_variable_name', value=comp_name, _RC)
      aerosol_names(1) = trim(comp_name) // 'phobic'
      aerosol_names(2) = trim(comp_name) // 'philic'

      ! Pressure at layer edges, expected to be 0-based
      call ESMF_InfoGet(info, key='air_pressure_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, ple, trim(fld_name), _RC)
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple

      i1 = lbound(ple0, 1); i2 = ubound(ple0, 1)
      j1 = lbound(ple0, 2); j2 = ubound(ple0, 2)
      km = ubound(ple0, 3)

      ! Relative humidity
      call ESMF_InfoGet(info, key='relative_humidity_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, rh, trim(fld_name), _RC)

      allocate(ext_s(i1:i2, j1:j2, km), &
           ssa_s(i1:i2, j1:j2, km), &
           asy_s(i1:i2, j1:j2, km), _STAT)

      allocate(q_4d(i1:i2, j1:j2, km, n_ca_bins), _STAT)

      do n = 1, n_ca_bins
         call MAPL_StateGetPointer(state, q, trim(aerosol_names(n)), _RC)
         do k = 1, km
            q_4d(:, :, k, n) = ((ple0(:, :, k) - ple0(:, :, k - 1)) * 0.01) * (100. / MAPL_GRAV) * q(:, :, k)
         end do
      end do

      call ESMF_InfoGetAlloc(info, key='mieTable_pointer', values=opaque_self, _RC)

      address = transfer(opaque_self, address)
      call c_f_pointer(address, self)

      if (usePhotTable /= 0) then
         wavelength = band * 1.e-9
         allocate(pmom_s(i1:i2, j1:j2, km, self%phot_Mie%nmom), _STAT)
         call miephot_(self%phot_Mie, n_ca_bins, wavelength, q_4d, rh, ext_s, ssa_s, pmom_s, _RC)
      else
         call mie_(self%rad_Mie, n_ca_bins, band, q_4d, rh, ext_s, ssa_s, asy_s, _RC)
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
         integer :: l, m
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
            do m = 1, mie%nmom
               bpmom_s(:, :, :, m) = bpmom_s(:, :, :, m) + pmom(:, :, :, m, 1) * (bssa * bext) ! moments multiplied by scattering
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
      type(CA2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: fld_name
      character(len=:), allocatable :: comp_name
      character(len=ESMF_MAXSTR) :: aerosol_names(n_ca_bins)

      real, dimension(:, :, :), allocatable :: tau_s, tau ! (lon:,lat:,lev:)
      integer :: instance
      integer :: n, i1, j1, i2, j2, km
      real :: wavelength
      integer :: k
      integer :: status

      ! Begin...

      call ESMF_InfoGetFromHost(state, info, _RC)

      ! Mie Table instance/index
      call ESMF_InfoGet(info, key='mie_table_instance', value=instance, _RC)

      ! Radiation band
      call ESMF_InfoGet(info, key='wavelength_for_aerosol_optics', value=wavelength, default=0., _RC)

      ! Build the two bin field names from the scalar component name
      call ESMF_InfoGetCharAlloc(info, key='internal_variable_name', value=comp_name, _RC)
      aerosol_names(1) = trim(comp_name) // 'phobic'
      aerosol_names(2) = trim(comp_name) // 'philic'

      ! Pressure at layer edges, expected to be 0-based
      call ESMF_InfoGet(info, key='air_pressure_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, ple, trim(fld_name), _RC)
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple

      i1 = lbound(ple0, 1); i2 = ubound(ple0, 1)
      j1 = lbound(ple0, 2); j2 = ubound(ple0, 2)
      km = ubound(ple0, 3)

      ! Relative humidity
      call ESMF_InfoGet(info, key='relative_humidity_for_aerosol_optics', value=fld_name, _RC)
      call MAPL_StateGetPointer(state, rh, trim(fld_name), _RC)

      allocate(tau_s(i1:i2, j1:j2, km), &
           tau(i1:i2, j1:j2, km), _STAT)
      tau_s = 0.0
      tau = 0.0

      allocate(q_4d(i1:i2, j1:j2, km, n_ca_bins), _STAT)

      do n = 1, n_ca_bins
         call MAPL_StateGetPointer(state, q, trim(aerosol_names(n)), _RC)
         do k = 1, km
            q_4d(:, :, k, n) = ((ple0(:, :, k) - ple0(:, :, k - 1)) * 0.01) * (100. / MAPL_GRAV) * q(:, :, k)
         end do
      end do

      call ESMF_InfoGetAlloc(info, key='mieTable_pointer', values=opaque_self, _RC)

      address = transfer(opaque_self, address)
      call c_f_pointer(address, self)

      do n = 1, n_ca_bins
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

end module CA2G_GridCompMod

subroutine CA2G_SetServices(gc, rc)
   use ESMF
   use CA2G_GridCompMod, only : mySetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine CA2G_SetServices
