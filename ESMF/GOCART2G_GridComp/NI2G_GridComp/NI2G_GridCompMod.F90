#include "MAPL.h"

!MODULE: NI2G_GridCompMod - GOCART refactoring of the Nitrate gridded component

!INTERFACE:
module NI2G_GridCompMod

   !USES:
   use ESMF
   use pflogger, only: logger_t => logger
   use gFTL2_StringVector, only: StringVector
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

   ! Number of dust/sea salt bins NI imports for heterogeneous chemistry
   integer, parameter :: n_du_ss_bins = 5
   ! Number of nitrate size bins
   integer, parameter :: n_no3_bins = 3

   real, parameter :: OCEAN = 0.0, LAND = 1.0, SEA_ICE = 2.0
   real, parameter :: fMassHNO3 = 63., fMassNO3 = 62.
   real, parameter :: cpd = 1004.16

   ! MAPL2's MAPL_UNDEF; MAPL3's is -huge, which breaks the 1.01*x > undef masks
   real, parameter :: MAPL2_UNDEF = 1.0e15

   ! relative position of nitrate tracers, matching NI2G_StateSpecs.rc order
   integer, parameter :: nNH3 = 1
   integer, parameter :: nNH4a = 2
   integer, parameter :: nNO3an1 = 3
   integer, parameter :: nNO3an2 = 4
   integer, parameter :: nNO3an3 = 5

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices

   !DESCRIPTION: This module implements GOCART's Nitrate (NI) Gridded Component.

   !REVISION HISTORY:
   ! 4January2024   Collow - Updated call for ChemSettling
   ! 01July2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

   ! Nitrate state
   type :: ThreadWorkspace
      logical :: first = .true.
   end type ThreadWorkspace

   type, extends(GA_Environment) :: NI2G_GridComp
      logical :: recycle_HNO3 = .false.
      real, allocatable :: rmedDU(:), rmedSS(:) ! DU and SS radius
      real, allocatable :: fnumDU(:), fnumSS(:) ! DU and SS particles per kg mass
      type(ThreadWorkspace), allocatable :: workspaces(:)
      type(ESMF_Alarm) :: alarm
      logical :: using_GMI ! logic for HNO3 GMI coupling
   end type NI2G_GridComp

   character(*), parameter :: PRIVATE_STATE = "NI2G_GridComp"

contains

   !IROUTINE: SetServices
   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! gridded component
      integer, intent(out) :: rc ! return code

      !DESCRIPTION: Declares NI's private state, entry points (Initialize and
      ! the Run0/Run1/Run2 phases) and state specs.

      !REVISION HISTORY:
      ! 24oct2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

      ! Locals
      character(len=:), allocatable :: comp_name
      type(NI2G_GridComp), pointer :: self
      logical :: data_driven = .true.
      integer :: num_threads
      type(MAPL_UngriddedDim) :: ungrd_wavelengths_profile, ungrd_wavelengths_vertint
      type(MAPL_UngriddedDim) :: ungrd_no3_bins
      type(MAPL_UngriddedDim) :: ungrd_du_ss_bins
      type(StringVector) :: advected_tracers
      integer :: status

      call MAPL_GridCompGet(gc, name=comp_name, _RC)

      ! Wrap gridcomp's private state and store it in gridcomp
      _SET_NAMED_PRIVATE_STATE(gc, NI2G_GridComp, PRIVATE_STATE)

      ! Retrieve the private state
      _GET_NAMED_PRIVATE_STATE(gc, NI2G_GridComp, PRIVATE_STATE, self)

      num_threads = MAPL_get_num_threads()
      allocate(self%workspaces(0:num_threads - 1), _STAT)

      ! process generic config items
      call self%GA_Environment%load_from_config(gc, _RC)

      ! Defined UngriddedDim items
      ungrd_wavelengths_profile = MAPL_UngriddedDim( &
           size(self%wavelengths_profile), &
           name="wavelengths_profile", &
           units="nm")
      ungrd_wavelengths_vertint = MAPL_UngriddedDim( &
           size(self%wavelengths_vertint), &
           name="wavelengths_vertint", &
           units="nm")
      ungrd_no3_bins = MAPL_UngriddedDim(n_no3_bins, name="no3_bins", units="1")
      ungrd_du_ss_bins = MAPL_UngriddedDim(n_du_ss_bins, name="du_ss_bins", units="1")

      ! process NI-specific items
      call MAPL_GridCompGetResource(gc, "using_GMI", self%using_GMI, _RC)

      ! Is NI data driven?
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
              short_name='NO3an1', &
              standard_name='Nitrate size bin 001', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              fill_value=0.0, &
              add_to_export=.true., _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_INTERNAL, &
              short_name='NO3an2', &
              standard_name='Nitrate size bin 002', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              fill_value=0.0, &
              add_to_export=.true., _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_INTERNAL, &
              short_name='NO3an3', &
              standard_name='Nitrate size bin 003', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, &
              fill_value=0.0, &
              add_to_export=.true., _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='climNO3an1', &
              standard_name='Nitrate size bin 001', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='climNO3an2', &
              standard_name='Nitrate size bin 002', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, _RC)

         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name='climNO3an3', &
              standard_name='Nitrate size bin 003', &
              units='kg kg-1', &
              dims='xyz', &
              vertical_stagger=MAPL_VERTICAL_STAGGER_CENTER, _RC)
      end if ! (data_driven)

      ! Import, Export, Internal states for computational instance
      if (.not. data_driven) then
#include "NI2G_Export___.h"
         associate(using_GMI => self%using_GMI)
#include "NI2G_Import___.h"
         end associate
#include "NI2G_Internal___.h"
         ! Replaces MAPL2's FRIENDLYTO=DYNAMICS:TURBULENCE:MOIST on NH3, NH4a,
         ! NO3an1, NO3an2, NO3an3
         call advected_tracers%push_back("NH3")
         call advected_tracers%push_back("NH4a")
         call advected_tracers%push_back("NO3an1")
         call advected_tracers%push_back("NO3an2")
         call advected_tracers%push_back("NO3an3")
         call MAPL_GridCompAddSpec(gc, &
              state_intent=ESMF_STATEINTENT_IMPORT, &
              short_name="TRADV", &
              standard_name="advected_quantities", &
              itemtype=MAPL_STATEITEM_SERVICE, &
              service_items=advected_tracers, _RC)
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

      !DESCRIPTION: This initializes the Nitrate gridded component. It
      ! primarily fills GOCART's AERO state with its nitrate fields.

      !REVISION HISTORY:
      ! 30June2020   E.Sherman  First attempt at refactoring

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_Geom) :: geom
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(ESMF_Field) :: field
      type(ESMF_Info) :: field_info, aero_info
      type(NI2G_GridComp), pointer :: self

      integer, allocatable :: mieTable_pointer(:)
      integer :: km
      integer :: instance
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
      _GET_NAMED_PRIVATE_STATE(gc, NI2G_GridComp, PRIVATE_STATE, self)

      self%km = km

      ! Get DTs
      call MAPL_ClockGet(clock, dt=HDT, _RC)
      call MAPL_GridCompGetResource(gc, "GOCART2G_DT", CDT, default=real(HDT), _RC)
      self%CDT = CDT

      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is NI data driven?
      call determine_data_driven(comp_name, data_driven, _RC)

      ! Get DU and SS attribute information for use in heterogeneous chemistry
      if (.not. data_driven) then
         call ESMF_StateGet(import, 'DU', field, _RC)
         call ESMF_InfoGetFromHost(field, field_info, _RC)
         call ESMF_InfoGetAlloc(field_info, key='radius', values=self%rmedDU, _RC)
         call ESMF_InfoGetAlloc(field_info, key='fnum', values=self%fnumDU, _RC)

         call ESMF_StateGet(import, 'SS', field, _RC)
         call ESMF_InfoGetFromHost(field, field_info, _RC)
         call ESMF_InfoGetAlloc(field_info, key='radius', values=self%rmedSS, _RC)
         call ESMF_InfoGetAlloc(field_info, key='fnum', values=self%fnumSS, _RC)
      end if

      ! Add attribute information to internal state variables
      if (.not. data_driven) then
         call ESMF_StateGet(internal, 'NH3', field, _RC)
         call ESMF_InfoGetFromHost(field, field_info, _RC)
         call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(nNH3), _RC)
         call get_HenrysLawCts('NH3', Vect_Hcts(1), Vect_Hcts(2), Vect_Hcts(3), Vect_Hcts(4), _RC)
         call ESMF_InfoSet(field_info, key='SetofHenryLawCts', values=Vect_Hcts, _RC)

         call ESMF_StateGet(internal, 'NH4a', field, _RC)
         call ESMF_InfoGetFromHost(field, field_info, _RC)
         call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(nNH4a), _RC)
      end if

      ! Fill AERO State with NO3an(1,2,3) fields
      call ESMF_StateGet(export, trim(comp_name) // '_AERO', aero, _RC)

      call ESMF_StateGet(internal, 'NO3an1', field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(nNO3an1), _RC)
      call ESMF_StateAdd(aero, [field], _RC)

      call ESMF_StateGet(internal, 'NO3an2', field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(nNO3an2), _RC)
      call ESMF_StateAdd(aero, [field], _RC)

      call ESMF_StateGet(internal, 'NO3an3', field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key='ScavengingFractionPerKm', value=self%fscav(nNO3an3), _RC)
      call ESMF_StateAdd(aero, [field], _RC)

      if (data_driven) then
         instance = instanceData
      else
         instance = instanceComputational
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

      ! Add variables to NI instance's aero state. This is used in aerosol optics calculations
      call add_aero(aero, label="air_pressure_for_aerosol_optics", label2="PLE", geom=geom, km=self%km, _RC)
      call add_aero(aero, label="relative_humidity_for_aerosol_optics", label2="RH", geom=geom, km=self%km, _RC)

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

      ! A single name identifying this species; the three bin fields
      ! (NO3an1, NO3an2, NO3an3) are hardcoded in aerosol_optics/get_mixR
      call ESMF_InfoSet(aero_info, key="internal_variable_name", value="NI", _RC)

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
      type(ESMF_GridComp) :: gc ! Gridded component
      type(ESMF_State) :: import ! Import state
      type(ESMF_State) :: export ! Export state
      type(ESMF_Clock) :: clock ! The clock
      integer, intent(out) :: rc ! Error code:

      !DESCRIPTION: Clears klid to 0.0 for Nitrate

      ! Locals
      type(ESMF_State) :: internal
      type(NI2G_GridComp), pointer :: self
      real, pointer, dimension(:, :, :) :: ple, ple0
      real, pointer, dimension(:, :, :) :: ptr3d_int
      integer :: status

      ! Begin...

      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, NI2G_GridComp, PRIVATE_STATE, self)

      ! Edge variable PLE is expected to be 0-based
      call MAPL_StateGetPointer(import, ple, 'PLE', _RC)
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple

      ! Set klid and Set internal values to 0 above klid
      call findKlid(self%klid, self%plid, ple0, _RC)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NO3an1', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NO3an2', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NO3an3', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NH4a', _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NH3', _RC)
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

      !DESCRIPTION: Run method for the Nitrate Grid Component. Determines
      ! whether to run data or computational run method.

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      logical :: data_driven
      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is NI data driven?
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

      !DESCRIPTION: Computes emissions/sources for Nitrate

      ! Locals
      type(ESMF_State) :: internal
      type(NI2G_GridComp), pointer :: self

#include "NI2G_DeclarePointer___.h"

      integer :: status

      ! Begin...

      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, NI2G_GridComp, PRIVATE_STATE, self)

      associate (using_GMI => self%using_GMI)
#include "NI2G_GetPointer___.h"
      end associate

      ! NH3 Emissions
      if (associated(NH3EM)) then
         NH3EM = 0.
         if (associated(EMI_NH3_BB)) NH3EM = NH3EM + EMI_NH3_BB
         if (associated(EMI_NH3_AG)) NH3EM = NH3EM + EMI_NH3_AG
         if (associated(EMI_NH3_EN)) NH3EM = NH3EM + EMI_NH3_EN
         if (associated(EMI_NH3_TR)) NH3EM = NH3EM + EMI_NH3_TR
         if (associated(EMI_NH3_RE)) NH3EM = NH3EM + EMI_NH3_RE
         if (associated(EMI_NH3_IN)) NH3EM = NH3EM + EMI_NH3_IN
         if (associated(EMI_NH3_OC)) NH3EM = NH3EM + EMI_NH3_OC
      end if

      if (associated(EMI_NH3_BB)) &
         NH3(:, :, self%km) = NH3(:, :, self%km) + self%cdt * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_BB
      if (associated(EMI_NH3_AG)) &
         NH3(:, :, self%km) = NH3(:, :, self%km) + self%cdt * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_AG
      if (associated(EMI_NH3_EN)) &
         NH3(:, :, self%km) = NH3(:, :, self%km) + self%cdt * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_EN
      if (associated(EMI_NH3_IN)) &
         NH3(:, :, self%km) = NH3(:, :, self%km) + self%cdt * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_IN
      if (associated(EMI_NH3_RE)) &
         NH3(:, :, self%km) = NH3(:, :, self%km) + self%cdt * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_RE
      if (associated(EMI_NH3_TR)) &
         NH3(:, :, self%km) = NH3(:, :, self%km) + self%cdt * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_TR
      if (associated(EMI_NH3_OC)) &
         NH3(:, :, self%km) = NH3(:, :, self%km) + self%cdt * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_OC

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(clock)

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

      !DESCRIPTION: Run2 method for the Nitrate Grid Component.

      ! Locals
      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      type(NI2G_GridComp), pointer :: self

      real, allocatable, dimension(:, :) :: drydepositionfrequency, dqa
      real, pointer, dimension(:, :, :) :: nisd_vel
      logical :: KIN
      real, allocatable, target, dimension(:, :, :) :: fluxoutWT
      real, allocatable, dimension(:, :, :, :) :: aerosol
      real, pointer, dimension(:, :) :: flux_ptr
      real, pointer, dimension(:, :, :) :: fluxWT_ptr

      logical :: alarm_is_ringing
      integer :: i1, j1, i2, j2, km
      real, target, allocatable, dimension(:, :, :) :: RH20, RH80
      integer :: i, j
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread
      integer :: settling_opt
      real, pointer, dimension(:, :, :) :: ple0, zle0, pfl_lsan0, pfi_lsan0

#include "NI2G_DeclarePointer___.h"

      integer :: status

      ! Begin...

      call MAPL_GridCompGet(gc, name=comp_name, _RC)
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      _GET_NAMED_PRIVATE_STATE(gc, NI2G_GridComp, PRIVATE_STATE, self)

      associate (using_GMI => self%using_GMI)
#include "NI2G_GetPointer___.h"
      end associate

      ! Edge variables are expected to be 0-based
      ple0(lbound(ple, 1):, lbound(ple, 2):, 0:) => ple
      zle0(lbound(zle, 1):, lbound(zle, 2):, 0:) => zle
      pfl_lsan0(lbound(pfl_lsan, 1):, lbound(pfl_lsan, 2):, 0:) => pfl_lsan
      pfi_lsan0(lbound(pfi_lsan, 1):, lbound(pfi_lsan, 2):, 0:) => pfi_lsan

      ! Set klid and Set internal values to 0 above klid
      call findKlid(self%klid, self%plid, ple0, _RC)
      call setZeroKlid(self%km, self%klid, NO3an1)
      call setZeroKlid(self%km, self%klid, NO3an2)
      call setZeroKlid(self%km, self%klid, NO3an3)
      call setZeroKlid(self%km, self%klid, NH4a)
      call setZeroKlid(self%km, self%klid, NH3)

      allocate(dqa, mold=lwi, _STAT)
      allocate(drydepositionfrequency, mold=lwi, _STAT)

      if (self%using_GMI) then
         xhno3 = GMI_HNO3
      else
         alarm_is_ringing = ESMF_AlarmIsRinging(self%alarm, _RC)
#ifdef DEBUG
         if (alarm_is_ringing) then
            if (MAPL_Am_I_Root()) then
               print *, 'DEBUG:: NI replenish alarm is ringing'
            end if
         end if
#endif

         ! Save local copy of HNO3 for first pass through run method regardless
         thread = MAPL_get_current_thread()
         workspace => self%workspaces(thread)

         ! Recycle HNO3 every 3 hours
         if (alarm_is_ringing) then
            xhno3 = NITRATE_HNO3
         end if
      end if

      if (associated(NIPNO3AQ)) NIPNO3AQ(:, :) = 0.
      if (associated(NIPNH4AQ)) NIPNH4AQ(:, :) = 0.
      if (associated(NIPNH3AQ)) NIPNH3AQ(:, :) = 0.
      if (associated(NIPHNO3AQ)) NIPHNO3AQ(:, :) = 0.
      ! 3D diag
      if (associated(NIPNO3AQV)) NIPNO3AQV(:, :, :) = 0.
      if (associated(NIPNH4AQV)) NIPNH4AQV(:, :, :) = 0.
      if (associated(NIPNH3AQV)) NIPNH3AQV(:, :, :) = 0.
      if (associated(NIPHNO3AQV)) NIPHNO3AQV(:, :, :) = 0.

      call NIthermo(self%km, self%klid, self%cdt, MAPL_GRAV, delp, airdens, &
           t, rh2, fMassHNO3, MAPL_AIRMW, SO4, NH3, NO3an1, NH4a, &
           xhno3, NIPNO3AQ, NIPNH4AQ, NIPNH3AQ, NIPHNO3AQ, &
           NIPNO3AQV, NIPNH4AQV, NIPNH3AQV, NIPHNO3AQV, _RC)

      call NIheterogenousChem(NIHT, NIHTV, xhno3, MAPL2_UNDEF, MAPL_AVOGAD, MAPL_AIRMW, &
           MAPL_PI, MAPL_RUNIV / 1000., airdens, t, rh2, delp, DU, &
           SS, self%rmedDU * 1.e-6, self%rmedSS * 1.e-6, &
           self%fnumDU, self%fnumSS, self%km, self%klid, &
           self%cdt, MAPL_GRAV, fMassHNO3, fMassNO3, NO3an1, NO3an2, &
           NO3an3, HNO3CONC, HNO3SMASS, HNO3CMASS, _RC)

      ! Update GMI_HNO3 after nitrate processes for 2-way coupling
      if (self%using_GMI) then
         GMI_HNO3 = xhno3
      end if

      ! NI Settling
      select case (self%settling_scheme)
      case ('gocart')
         settling_opt = 1
      case ('ufs')
         settling_opt = 2
      case default
         _ASSERT_RC(.false., 'Unsupported settling scheme: ' // trim(self%settling_scheme), ESMF_RC_NOT_IMPL)
      end select

      ! Ammonium - settles like bin 1 of nitrate
      nullify(nisd_vel)
      if (associated(NH4SD_V)) nisd_vel => NH4SD_V
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 1, self%cdt, MAPL_GRAV, &
           NH4a, t, airdens, rh2, zle0, delp, NH4SD, nisd_vel, &
           settling_scheme=settling_opt, _RC)

      ! Nitrate Bin 1
      nullify(flux_ptr)
      if (associated(NISD)) flux_ptr => NISD(:, :, 1)
      nullify(nisd_vel)
      if (associated(NISD_V)) nisd_vel => NISD_V(:, :, :, 1)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 1, self%cdt, MAPL_GRAV, &
           NO3an1, t, airdens, &
           rh2, zle0, delp, flux_ptr, nisd_vel, settling_scheme=settling_opt, _RC)
      ! Nitrate Bin 2
      nullify(flux_ptr)
      if (associated(NISD)) flux_ptr => NISD(:, :, 2)
      nullify(nisd_vel)
      if (associated(NISD_V)) nisd_vel => NISD_V(:, :, :, 2)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 2, self%cdt, MAPL_GRAV, &
           NO3an2, t, airdens, &
           rh2, zle0, delp, flux_ptr, nisd_vel, settling_scheme=settling_opt, _RC)
      ! Nitrate Bin 3
      nullify(flux_ptr)
      if (associated(NISD)) flux_ptr => NISD(:, :, 3)
      nullify(nisd_vel)
      if (associated(NISD_V)) nisd_vel => NISD_V(:, :, :, 3)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 3, self%cdt, MAPL_GRAV, &
           NO3an3, t, airdens, &
           rh2, zle0, delp, flux_ptr, nisd_vel, settling_scheme=settling_opt, _RC)

      ! NI Deposition
      drydepositionfrequency = 0.
      call DryDeposition(self%km, t, airdens, zle0, lwi, ustar, zpbl, sh, &
           MAPL_KARMAN, cpd, MAPL_GRAV, z0h, drydepositionfrequency, _RC)

      ! NH3
      dqa = 0.
      do i = 1, ubound(lwi, 1)
         do j = 1, ubound(lwi, 2)
            if (abs(lwi(i, j) - OCEAN) < 0.5) then
               dqa(i, j) = max(0.0, NH3(i, j, self%km) * (1. - exp(-10.0 * drydepositionfrequency(i, j) * self%cdt)))
            else
               dqa(i, j) = max(0.0, NH3(i, j, self%km) * (1. - exp(-3.0 * drydepositionfrequency(i, j) * self%cdt)))
            end if
         end do
      end do

      NH3(:, :, self%km) = NH3(:, :, self%km) - dqa
      if (associated(NH3DP)) NH3DP = dqa * delp(:, :, self%km) / MAPL_GRAV / self%cdt

      ! NH4a
      dqa = 0.
      dqa = max(0.0, NH4a(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%cdt)))
      NH4a(:, :, self%km) = NH4a(:, :, self%km) - dqa
      if (associated(NH4DP)) NH4DP = dqa * delp(:, :, self%km) / MAPL_GRAV / self%cdt

      ! NO3anx
      dqa = 0.
      dqa = max(0.0, NO3an1(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%cdt)))
      NO3an1(:, :, self%km) = NO3an1(:, :, self%km) - dqa
      if (associated(NIDP)) NIDP(:, :, 1) = dqa * delp(:, :, self%km) / MAPL_GRAV / self%cdt

      dqa = 0.
      dqa = max(0.0, NO3an2(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%cdt)))
      NO3an2(:, :, self%km) = NO3an2(:, :, self%km) - dqa
      if (associated(NIDP)) NIDP(:, :, 2) = dqa * delp(:, :, self%km) / MAPL_GRAV / self%cdt

      dqa = 0.
      dqa = max(0.0, NO3an3(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%cdt)))
      NO3an3(:, :, self%km) = NO3an3(:, :, self%km) - dqa
      if (associated(NIDP)) NIDP(:, :, 3) = dqa * delp(:, :, self%km) / MAPL_GRAV / self%cdt

      ! NI Large-scale Wet Removal
      if (associated(NH3WT) .or. associated(NH4WT)) then
         allocate(fluxoutWT(ubound(t, 1), ubound(t, 2), 1), _STAT)
      end if
      ! NH3
      KIN = .false.
      nullify(fluxWT_ptr)
      if (associated(NH3WT)) fluxWT_ptr => fluxoutWT
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'NH3', &
           KIN, MAPL_GRAV, self%fwet(nNH3), NH3, ple0, t, airdens, &
           pfl_lsan0, pfi_lsan0, cn_prcp, ncn_prcp, fluxWT_ptr, _RC)
      if (associated(NH3WT)) NH3WT = fluxWT_ptr(:, :, 1)

      ! NH4a
      KIN = .true.
      nullify(fluxWT_ptr)
      if (associated(NH4WT)) fluxWT_ptr => fluxoutWT
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'NH4a', &
           KIN, MAPL_GRAV, self%fwet(nNH4a), NH4a, ple0, t, airdens, &
           pfl_lsan0, pfi_lsan0, cn_prcp, ncn_prcp, fluxWT_ptr, _RC)
      if (associated(NH4WT)) NH4WT = fluxWT_ptr(:, :, 1)

      if (allocated(fluxoutWT)) then
         deallocate(fluxoutWT, _STAT)
      end if

      KIN = .true.
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'nitrate', &
           KIN, MAPL_GRAV, self%fwet(nNO3an1), NO3an1, ple0, t, airdens, &
           pfl_lsan0, pfi_lsan0, cn_prcp, ncn_prcp, NIWT, _RC)

      KIN = .true.
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 2, self%cdt, 'nitrate', &
           KIN, MAPL_GRAV, self%fwet(nNO3an2), NO3an2, ple0, t, airdens, &
           pfl_lsan0, pfi_lsan0, cn_prcp, ncn_prcp, NIWT, _RC)

      KIN = .true.
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 3, self%cdt, 'nitrate', &
           KIN, MAPL_GRAV, self%fwet(nNO3an3), NO3an3, ple0, t, airdens, &
           pfl_lsan0, pfi_lsan0, cn_prcp, ncn_prcp, NIWT, _RC)

      ! Compute desired output diagnostics
      ! Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
      allocate(aerosol(ubound(NH4a, 1), ubound(NH4a, 2), ubound(NH4a, 3), n_no3_bins), _STAT)
      aerosol(:, :, :, :) = 0.0
      aerosol(:, :, :, 1) = NH4a
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=1, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple0, tropp=tropp, &
           sfcmass=NH4SMASS, colmass=NH4CMASS, mass=NH4MASS, conc=NH4CONC, _RC)

      aerosol(:, :, :, 1) = NH3
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=1, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple0, tropp=tropp, &
           sfcmass=NH3SMASS, colmass=NH3CMASS, mass=NH3MASS, conc=NH3CONC, _RC)

      aerosol(:, :, :, 1) = NO3an1
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=1, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple0, tropp=tropp, &
           sfcmass=NISMASS25, colmass=NICMASS25, mass=NIMASS25, conc=NICONC25, &
           exttau25=NIEXTT25, scatau25=NISCAT25, exttaufm=NIEXTTFM, scataufm=NISCATFM, &
           NO3nFlag=.true., _RC)

      aerosol(:, :, :, 1) = NO3an1
      aerosol(:, :, :, 2) = NO3an2
      aerosol(:, :, :, 3) = NO3an3
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=n_no3_bins, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple0, tropp=tropp, sfcmass=NISMASS, colmass=NICMASS, mass=NIMASS, conc=NICONC, &
           exttau=NIEXTTAU, stexttau=NISTEXTTAU, scatau=NISCATAU, stscatau=NISTSCATAU, &
           fluxu=NIFLUXU, fluxv=NIFLUXV, extcoef=NIEXTCOEF, scacoef=NISCACOEF, &
           bckcoef=NIBCKCOEF, angstrom=NIANGSTR, _RC)

      i1 = lbound(rh2, 1); i2 = ubound(rh2, 1)
      j1 = lbound(rh2, 2); j2 = ubound(rh2, 2)
      km = ubound(rh2, 3)

      allocate(RH20(i1:i2, j1:j2, km), _STAT)
      allocate(RH80(i1:i2, j1:j2, km), _STAT)

      RH20(:, :, :) = 0.20

      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=n_no3_bins, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, aerosol=aerosol, &
           grav=MAPL_GRAV, tmpu=t, rhoa=airdens, &
           rh=RH20, u=u, v=v, delp=delp, ple=ple0, tropp=tropp, &
           extcoef=NIEXTCOEFRH20, scacoef=NISCACOEFRH20, _RC)

      RH80(:, :, :) = 0.80
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=n_no3_bins, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, aerosol=aerosol, &
           grav=MAPL_GRAV, tmpu=t, rhoa=airdens, &
           rh=RH80, u=u, v=v, delp=delp, ple=ple0, tropp=tropp, &
           extcoef=NIEXTCOEFRH80, scacoef=NISCACOEFRH80, _RC)

      deallocate(RH20, RH80)

      _RETURN(_SUCCESS)
      _UNUSED_DUMMY(clock)

   end subroutine Run2

   !IROUTINE: Run_data -- ExtData Nitrate Grid Component
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
      real, pointer, dimension(:, :, :) :: ptr3d_int, ptr3d_imp
      integer :: status

      ! Begin...

      ! Update internal data pointers with ExtData
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NO3an1', _RC)
      call MAPL_StateGetPointer(import, ptr3d_imp, 'climNO3an1', _RC)
      ptr3d_int = ptr3d_imp
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NO3an2', _RC)
      call MAPL_StateGetPointer(import, ptr3d_imp, 'climNO3an2', _RC)
      ptr3d_int = ptr3d_imp
      call MAPL_StateGetPointer(internal, ptr3d_int, 'NO3an3', _RC)
      call MAPL_StateGetPointer(import, ptr3d_imp, 'climNO3an3', _RC)
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
      type(NI2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: fld_name

      ! Nitrate's three size bins are plain named fields, not an ungridded
      ! array; hardcode the names (as in NI2G_StateSpecs.rc / the generated
      ! NI2G_GetPointer___.h) rather than reading a list from Info.
      character(len=6), parameter :: aerosol_names(n_no3_bins) = ['NO3an1', 'NO3an2', 'NO3an3']

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

      allocate(q_4d(i1:i2, j1:j2, km, n_no3_bins), _STAT)

      do n = 1, n_no3_bins
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
         call miephot_(self%phot_Mie, n_no3_bins, wavelength, q_4d, rh, ext_s, ssa_s, pmom_s, _RC)
      else
         call mie_(self%rad_Mie, n_no3_bins, band, q_4d, rh, ext_s, ssa_s, asy_s, _RC)
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
      type(NI2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: fld_name
      character(len=6), parameter :: aerosol_names(n_no3_bins) = ['NO3an1', 'NO3an2', 'NO3an3']

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

      allocate(q_4d(i1:i2, j1:j2, km, n_no3_bins), _STAT)

      do n = 1, n_no3_bins
         call MAPL_StateGetPointer(state, q, trim(aerosol_names(n)), _RC)
         do k = 1, km
            q_4d(:, :, k, n) = ((ple0(:, :, k) - ple0(:, :, k - 1)) * 0.01) * (100. / MAPL_GRAV) * q(:, :, k)
         end do
      end do

      call ESMF_InfoGetAlloc(info, key='mieTable_pointer', values=opaque_self, _RC)

      address = transfer(opaque_self, address)
      call c_f_pointer(address, self)

      do n = 1, n_no3_bins
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

end module NI2G_GridCompMod

subroutine NI2G_SetServices(gc, rc)
   use ESMF
   use NI2G_GridCompMod, only : mySetServices => SetServices
   type(ESMF_GridComp) :: gc
   integer, intent(out) :: rc
   call mySetServices(gc, rc=rc)
end subroutine NI2G_SetServices
