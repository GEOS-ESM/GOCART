#include "MAPL_Generic.h"

!BOP
!MODULE: NI2G_GridCompMod - GOCART Nitrate gridded component

!INTERFACE:
module NI2G_GridCompMod

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

   real, parameter :: OCEAN = 0.0, LAND = 1.0, SEA_ICE = 2.0
   real, parameter :: fMassHNO3 = 63., fMassNO3 = 62.
   real, parameter :: cpd = 1004.16

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData = 2

   integer, parameter :: nNH3 = 1
   integer, parameter :: nNH4a = 2
   integer, parameter :: nNO3an1 = 3
   integer, parameter :: nNO3an2 = 4
   integer, parameter :: nNO3an3 = 5
   integer, parameter :: DP = kind(1.0d0)

   !PUBLIC MEMBER FUNCTIONS:
   public SetServices

   !DESCRIPTION: This module implements GOCART's Nitrate (NI) Gridded Component.

   !REVISION HISTORY:
   ! 4January2024   Collow - Updated call for ChemSettling
   ! 01July2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

   !EOP

   !Nitrate state
   type :: ThreadWorkspace
      logical :: first = .true.
   end type ThreadWorkspace

   type, extends(GA_Environment) :: NI2G_GridComp
      ! logical :: first
      logical :: recycle_HNO3 = .false.
      real, allocatable :: rmedDU(:), rmedSS(:) ! DU and SS radius
      real, allocatable :: fnumDU(:), fnumSS(:) ! DU and SS particles per kg mass
      type(ThreadWorkspace), allocatable :: workspaces(:)
      type(ESMF_Alarm) :: alarm
   end type NI2G_GridComp

   type wrap_
      type(NI2G_GridComp), pointer :: PTR !=> null()
   end type wrap_

contains

   !BOP
   !IROUTINE: SetServices
   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      integer, intent(out) :: rc

      !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
      !     the Initialize and Finalize services to generic versions. It also
      !     allocates our instance of a generic state and puts it in the
      !     gridded component (GC). Here we only set the two-stage run method
      !     and declare the data services.

      !REVISION HISTORY:
      !   24oct2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(ESMF_Config) :: cfg
      type(ESMF_Config) :: universal_cfg
      type(wrap_) :: wrap
      type(NI2G_GridComp), pointer :: self

      real :: defval
      logical :: data_driven = .true.
      logical :: file_exists
      integer :: num_threads

      __Iam__('SetServices')

      !   Get my name and set-up traceback handle
      !   ---------------------------------------
      call ESMF_GridCompGet(gc, NAME=comp_name, config=universal_cfg, _RC)
      Iam = trim(comp_name) // '::' // Iam

      !   Wrap internal state for storing in GC
      !   -------------------------------------
      allocate(self, _STAT)
      wrap%PTR => self

      num_threads = MAPL_get_num_threads()
      allocate(self%workspaces(0:num_threads - 1), _STAT)

      !   Load resource file
      !   -------------------
      cfg = ESMF_ConfigCreate(_RC)
      inquire(file='NI2G_instance_' // trim(comp_name) // '.rc', exist=file_exists)
      if (file_exists) then
         call ESMF_ConfigLoadFile(cfg, 'NI2G_instance_' // trim(comp_name) // '.rc', _RC)
      else
         if (mapl_am_i_root()) print*, 'NI2G_instance_' // trim(comp_name) // '.rc does not exist! Loading' // &
              ' NI2G_instance_NI.rc instead'
         call ESMF_ConfigLoadFile(cfg, 'NI2G_instance_NI.rc', _RC)
      end if

      ! process generic config items
      call self%GA_Environment%load_from_config(cfg, universal_cfg, _RC)

      !   Is NI data driven?
      !   ------------------
      call determine_data_driven(comp_name, data_driven, _RC)

      !   Set entry points
      !   ------------------------
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, _RC)
      if (data_driven .neqv. .true.) then
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run2, _RC)
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run0, _RC)
      end if

      defval = 0.0

      !   Import and Internal states if data instance
      !   -------------------------------------------
      if (data_driven) then

         !      Pressure at layer edges
         !      -----------------------
         call MAPL_AddImportSpec(gc, &
              short_name='PLE', &
              long_name='air_pressure', &
              units='Pa', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationEdge, &
              restart=MAPL_RestartSkip, _RC)

         !      RH: is between 0 and 1
         !      ----------------------
         call MAPL_AddImportSpec(gc, &
              short_name='RH2', &
              long_name='Rel_Hum_after_moist', &
              units='1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, &
              restart=MAPL_RestartSkip, _RC)

         call MAPL_AddInternalSpec(gc,&
              short_name='NO3an1', &
              long_name='Nitrate size bin 001', &
              units='kg kg-1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, &
              restart=MAPL_RestartOptional, &
              ! friendlyto='DYNAMICS:TURBULENCE:MOIST', &
              & add2export=.true., _RC)

         call MAPL_AddInternalSpec(gc,&
              short_name='NO3an2', &
              long_name='Nitrate size bin 002', &
              units='kg kg-1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, &
              restart=MAPL_RestartOptional, &
              ! friendlyto='DYNAMICS:TURBULENCE:MOIST', &
              add2export=.true., _RC)

         call MAPL_AddInternalSpec(gc,&
              short_name='NO3an3', &
              long_name='Nitrate size bin 003', &
              units='kg kg-1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, &
              restart=MAPL_RestartOptional, &
              ! friendlyto='DYNAMICS:TURBULENCE:MOIST', &
              add2export=.true., _RC)

         call MAPL_AddImportSpec(gc,&
              short_name='climNO3an1', &
              long_name='Nitrate size bin 001', &
              units='kg kg-1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, &
              restart=MAPL_RestartOptional, _RC)

         call MAPL_AddImportSpec(gc,&
              short_name='climNO3an2', &
              long_name='Nitrate size bin 002', &
              units='kg kg-1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, &
              restart=MAPL_RestartOptional, _RC)

         call MAPL_AddImportSpec(gc,&
              short_name='climNO3an3', &
              long_name='Nitrate size bin 003', &
              units='kg kg-1', &
              dims=MAPL_DimsHorzVert, &
              vlocation=MAPL_VLocationCenter, &
              restart=MAPL_RestartOptional, _RC)
      end if ! (data_driven)

      !   Import, Export, Internal states for computational instance
      !   ----------------------------------------------------------
      if (.not.data_driven) then
#include "NI2G_Export___.h"
#include "NI2G_Import___.h"
#include "NI2G_Internal___.h"
      end if

      !   This state holds fields needed by radiation
      !   ---------------------------------------------
      call MAPL_AddExportSpec(gc, &
           short_name=trim(comp_name) // '_AERO', &
           long_name='aerosols_from_' // trim(comp_name), &
           units='kg kg-1', &
           dims=MAPL_DimsHorzVert, &
           vlocation=MAPL_VLocationCenter, &
           DATATYPE=MAPL_StateItem, _RC)

      !   Store internal state in GC
      !   --------------------------
      call ESMF_UserCompSetInternalState(gc, 'NI2G_GridComp', wrap, _RC)

      !   Set generic services
      !   ----------------------------------
      call MAPL_GenericSetServices(gc, _RC)

      RETURN_(ESMF_SUCCESS)

   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_State), intent(inout) :: import
      type(ESMF_State), intent(inout) :: export
      type(ESMF_Clock), intent(inout) :: clock
      integer, optional, intent(out) :: rc

      !DESCRIPTION: This initializes the Nitrate gridded component. It primaryily
      !               fills GOCART's AERO states with its nitrate fields.

      !REVISION HISTORY:
      ! 30June2020   E.Sherman  First attempt at refactoring

      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_Config) :: universal_cfg
      type(ESMF_Grid) :: grid
      type(ESMF_State) :: internal
      type(ESMF_State) :: aero
      type(ESMF_State) :: providerState
      type(ESMF_Config) :: cfg
      type(ESMF_FieldBundle) :: bundle_dp
      type(wrap_) :: wrap
      type(NI2G_GridComp), pointer :: self

      integer, allocatable :: mieTable_pointer(:)
      integer :: i, dims(3), km
      integer :: instance
      type(ESMF_Field) :: field, fld
      character(len=ESMF_MAXSTR) :: prefix
      real :: CDT ! chemistry timestep (secs)
      integer :: HDT ! model     timestep (secs)
      logical :: data_driven
      logical :: bands_are_present
      integer :: NUM_BANDS
      character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)

      type(ESMF_Calendar) :: calendar
      type(ESMF_Time) :: currentTime
      type(ESMF_Time) :: ringTime
      type(ESMF_TimeInterval) :: ringInterval
      integer :: year, month, day, hh, mm, ss

      real, dimension(4) :: Vect_Hcts
      !    real, allocatable, dimension(:) :: rmedDU, rmedSS, fnumDU, fnumSS
      integer :: itemCount
      integer, allocatable, dimension(:) :: channels_
      integer :: nmom_
      character(len=ESMF_MAXSTR) :: file_
      logical :: file_exists

      __Iam__('Initialize')

      !   Get the target components name and set-up traceback handle.
      !   -----------------------------------------------------------
      call ESMF_GridCompGet(gc, grid=grid, NAME=comp_name, config=universal_cfg, _RC)
      Iam = trim(comp_name) // '::' // trim(Iam)

      !   Get my internal MAPL_Generic state
      !   -----------------------------------
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      !   Get my internal private state
      !   -----------------------------
      call ESMF_UserCompGetInternalState(gc, 'NI2G_GridComp', wrap, _RC)
      self => wrap%PTR

      !   Get dimensions
      !   ---------------
      call MAPL_GridGet(grid, localCellCountPerDim=dims, _RC)
      km = dims(3)
      self%km = km

      !   Get DTs
      !   -------
      call MAPL_GetResource(MAPL, HDT, Label='RUN_DT:', _RC)
      call MAPL_GetResource(MAPL, CDT, Label='GOCART2G_DT:', default=real(HDT), _RC)
      self%CDT = CDT

      !  Load resource file and get number of bins
      !  -------------------------------------------
      cfg = ESMF_ConfigCreate(_RC)
      inquire(file='NI2G_instance_' // trim(comp_name) // '.rc', exist=file_exists)
      if (file_exists) then
         call ESMF_ConfigLoadFile(cfg, 'NI2G_instance_' // trim(comp_name) // '.rc', _RC)
      else
         if (mapl_am_i_root()) print*, 'NI2G_instance_' // trim(comp_name) // '.rc does not exist! Loading' // &
              ' NI2G_instance_NI.rc instead'
         call ESMF_ConfigLoadFile(cfg, 'NI2G_instance_NI.rc', _RC)
      end if

      !self%first = .true.

      !   Call Generic Initialize
      !   ----------------------------------------
      call MAPL_GenericInitialize(gc, import, export, clock, _RC)

      !   Get parameters from generic state.
      !   -----------------------------------
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, _RC)

      !   Is NI data driven?
      !   ------------------
      call determine_data_driven(comp_name, data_driven, _RC)

      !   Get DU and SS attribute information for use in heterogenous chemistry
      if (.not.data_driven) then
         call ESMF_StateGet(import, 'DU', field, _RC)
         call ESMF_AttributeGet(field, NAME='radius', itemCount=itemCount, _RC)
         allocate(self%rmedDU(itemCount), _STAT)
         allocate(self%fnumDU(itemCount), _STAT)
         call ESMF_AttributeGet(field, NAME='radius', valueList=self%rmedDU, _RC)
         call ESMF_AttributeGet(field, NAME='fnum', valueList=self%fnumDU, _RC)

         call ESMF_StateGet(import, 'SS', field, _RC)
         call ESMF_AttributeGet(field, NAME='radius', itemCount=itemCount, _RC)
         allocate(self%rmedSS(itemCount), _STAT)
         allocate(self%fnumSS(itemCount), _STAT)
         call ESMF_AttributeGet(field, NAME='radius', valueList=self%rmedSS, _RC)
         call ESMF_AttributeGet(field, NAME='fnum', valueList=self%fnumSS, _RC)
      end if

      !   If this is a data component, the data is provided in the import
      !   state via ExtData instead of the actual GOCART children
      !   ----------------------------------------------------------------
      if (data_driven) then
         providerState = import
         prefix = 'clim'
      else
         providerState = export
         prefix = ''
      end if

      if (.not.data_driven) then
         call ESMF_StateGet(internal, 'NH3', field, _RC)
         call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), _RC)
         call get_HenrysLawCts('NH3', Vect_Hcts(1), Vect_Hcts(2), Vect_Hcts(3), Vect_Hcts(4), _RC)
         call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, _RC)

         call ESMF_StateGet(internal, 'NH4a', field, _RC)
         call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(2), _RC)
      end if

      !   Fill AERO State with N03an(1,2,3) fields
      !   ----------------------------------------
      call ESMF_StateGet(export, trim(comp_name) // '_AERO', aero, _RC)

      call ESMF_StateGet(internal, 'NO3an1', field, _RC)
      call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(3), _RC)
      fld = MAPL_FieldCreate(field, 'NO3an1', _RC)
      call MAPL_StateAdd(aero, fld, _RC)

      call ESMF_StateGet(internal, 'NO3an2', field, _RC)
      call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(4), _RC)
      fld = MAPL_FieldCreate(field, 'NO3an2', _RC)
      call MAPL_StateAdd(aero, fld, _RC)

      call ESMF_StateGet(internal, 'NO3an3', field, _RC)
      call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(5), _RC)
      fld = MAPL_FieldCreate(field, 'NO3an3', _RC)
      call MAPL_StateAdd(aero, fld, _RC)

      if (data_driven) then
         instance = instanceData
      else
         instance = instanceComputational
      end if

      self%instance = instance

      !   Create Radiation Mie Table
      !   --------------------------
      call ESMF_ConfigGetAttribute(cfg, file_, Label="aerosol_radBands_optics_file:", _RC)
      self%rad_Mie = GOCART2G_Mie(trim(file_), _RC)

      !   Create Diagnostics Mie Table
      !   -----------------------------
      !   Get file names for the optical tables
      call ESMF_ConfigGetAttribute(cfg, file_, &
           Label="aerosol_monochromatic_optics_file:", _RC)
      call ESMF_ConfigGetAttribute(cfg, nmom_, Label="n_moments:", default=0, _RC)
      i = ESMF_ConfigGetLen(universal_cfg, Label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', _RC)
      allocate(channels_(i), _STAT)
      call ESMF_ConfigGetAttribute(universal_cfg, channels_, &
           Label="aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:", _RC)
      self%diag_Mie = GOCART2G_Mie(trim(file_), channels_ * 1.e-9, nmom=nmom_, _RC)
      deallocate(channels_)

      ! Mie Table instance/index
      call ESMF_AttributeSet(aero, NAME='mie_table_instance', VALUE=instance, _RC)

      ! Add variables to NI instance's aero state. This is used in aerosol optics calculations
      call add_aero(aero, Label='air_pressure_for_aerosol_optics', label2='PLE', grid=grid, typekind=MAPL_R4, _RC)
      call add_aero(aero, Label='relative_humidity_for_aerosol_optics', label2='RH', grid=grid, typekind=MAPL_R4, &
           _RC)
      !   call ESMF_StateGet (import, 'PLE', field, _RC)
      !   call MAPL_StateAdd (aero, field, _RC)
      !   call ESMF_StateGet (import, 'RH2', field, _RC)
      !   call MAPL_StateAdd (aero, field, _RC)

      call add_aero(aero, Label='extinction_in_air_due_to_ambient_aerosol', label2='EXT', grid=grid, typekind=MAPL_R8, &
           _RC)
      call add_aero(aero, Label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=&
           MAPL_R8, _RC)
      call add_aero(aero, Label='asymmetry_parameter_of_ambient_aerosol', label2='ASY', grid=grid, typekind=MAPL_R8, &
           _RC)
      call add_aero(aero, Label='monochromatic_extinction_in_air_due_to_ambient_aerosol', &
           label2='monochromatic_EXT', grid=grid, typekind=MAPL_R4, _RC)
      call add_aero(aero, Label='sum_of_internalState_aerosol', label2='aerosolSum', grid=grid, typekind=MAPL_R4, &
           _RC)

      call ESMF_AttributeSet(aero, NAME='band_for_aerosol_optics', VALUE=0, _RC)
      call ESMF_AttributeSet(aero, NAME='wavelength_for_aerosol_optics', VALUE=0., _RC)

      mieTable_pointer = transfer(c_loc(self), [1])
      call ESMF_AttributeSet(aero, NAME='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer&
           ), _RC)

      allocate(aerosol_names(3), _STAT)
      aerosol_names(1) = 'NO3an1'
      aerosol_names(2) = 'NO3an2'
      aerosol_names(3) = 'NO3an3'
      call ESMF_AttributeSet(aero, NAME='internal_variable_name', valueList=aerosol_names, &
           itemCount=size(aerosol_names), _RC)

      call ESMF_MethodAdd(aero, Label='aerosol_optics', userRoutine=aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, Label='monochromatic_aerosol_optics', userRoutine=monochromatic_aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, Label='get_mixR', userRoutine=get_mixR, _RC)

      ! Deal with replenishment alarm (formerly the daily_alarm subroutine)
      ! ===================================================================
      self%alarm = createReplenishAlarm(gc, clock, 30000, _RC)
      _RETURN(ESMF_SUCCESS)

   end subroutine Initialize

   !BOP
   !IROUTINE: Run0
   !INTERFACE:
   subroutine Run0(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Clears klid to 0.0 for Nitrate
      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal
      type(wrap_) :: wrap
      type(NI2G_GridComp), pointer :: self
      real, pointer, dimension(:, :, :) :: ple
      real, pointer, dimension(:, :, :) :: ptr3d_int

      __Iam__('Run0')

      !   Get my name and set-up traceback handle
      !   ---------------------------------------
      call ESMF_GridCompGet(gc, NAME=comp_name, _RC)
      Iam = trim(comp_name) // '::' // Iam

      !   Get my internal MAPL_Generic state
      !   -----------------------------------
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      !   Get parameters from generic state.
      !   -----------------------------------
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, _RC)

      !   Get my private internal state
      !   ------------------------------
      call ESMF_UserCompGetInternalState(gc, 'NI2G_GridComp', wrap, _RC)
      self => wrap%PTR

      !   Set klid and Set internal values to 0 above klid
      !   ---------------------------------------------------
      call MAPL_GetPointer(import, ple, 'PLE', _RC)
      call findKlid(self%klid, self%plid, ple, _RC)
      call MAPL_GetPointer(internal, NAME='NO3an1', PTR=ptr3d_int, _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_GetPointer(internal, NAME='NO3an2', PTR=ptr3d_int, _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)
      call MAPL_GetPointer(internal, NAME='NO3an3', PTR=ptr3d_int, _RC)
      call setZeroKlid(self%km, self%klid, ptr3d_int)

      RETURN_(ESMF_SUCCESS)

   end subroutine Run0

   !BOP
   !IROUTINE: Run
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run method for the Nitrate Grid Component. Determines whether to run
      !               data or computational run method.
      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal
      logical :: data_driven

      __Iam__('Run')

      !   Get my name and set-up traceback handle
      !   ---------------------------------------
      call ESMF_GridCompGet(gc, NAME=comp_name, _RC)
      Iam = trim(comp_name) // '::' // Iam

      !   Get my internal MAPL_Generic state
      !   -----------------------------------
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      !   Get parameters from generic state.
      !   -----------------------------------
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, _RC)

      !   Is NI data driven?
      !   ------------------
      call determine_data_driven(comp_name, data_driven, _RC)

      !   Update INTERNAL state variables with ExtData
      !   ---------------------------------------------
      if (data_driven) then
         call Run_data(gc, import, export, internal, _RC)
      else
         call Run1(gc, import, export, clock, _RC)
      end if

      RETURN_(ESMF_SUCCESS)

   end subroutine Run

   !BOP
   !IROUTINE: Run1
   !INTERFACE:
   subroutine Run1(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION:  Computes emissions/sources for Nitrate
      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal
      type(ESMF_Grid) :: grid
      type(wrap_) :: wrap
      type(NI2G_GridComp), pointer :: self

#include "NI2G_DeclarePointer___.h"

      __Iam__('Run1')

      !   Get my name and set-up traceback handle
      !   ---------------------------------------
      call ESMF_GridCompGet(gc, NAME=comp_name, _RC)
      Iam = trim(comp_name) // '::' // Iam

      !   Get my internal MAPL_Generic state
      !   -----------------------------------
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      call MAPL_Get(MAPL, grid=grid, _RC)

      !   Get parameters from generic state.
      !   -----------------------------------
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, _RC)

#include "NI2G_GetPointer___.h"

      !   Get my private internal state
      !   ------------------------------
      call ESMF_UserCompGetInternalState(gc, 'NI2G_GridComp', wrap, _RC)
      self => wrap%PTR

      !   NH3 Emissions
      !   -------------
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
           NH3(:, :, self%km) = NH3(:, :, self%km) + self%CDT * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_BB
      if (associated(EMI_NH3_AG)) &
           NH3(:, :, self%km) = NH3(:, :, self%km) + self%CDT * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_AG
      if (associated(EMI_NH3_EN)) &
           NH3(:, :, self%km) = NH3(:, :, self%km) + self%CDT * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_EN
      if (associated(EMI_NH3_IN)) &
           NH3(:, :, self%km) = NH3(:, :, self%km) + self%CDT * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_IN
      if (associated(EMI_NH3_RE)) &
           NH3(:, :, self%km) = NH3(:, :, self%km) + self%CDT * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_RE
      if (associated(EMI_NH3_TR)) &
           NH3(:, :, self%km) = NH3(:, :, self%km) + self%CDT * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_TR
      if (associated(EMI_NH3_OC)) &
           NH3(:, :, self%km) = NH3(:, :, self%km) + self%CDT * MAPL_GRAV / delp(:, :, self%km) * EMI_NH3_OC

      RETURN_(ESMF_SUCCESS)

   end subroutine Run1

   !BOP
   !IROUTINE: Run2
   !INTERFACE:
   subroutine Run2(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_Clock), intent(inout) :: clock ! The clock
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION: Run2 method for the Nitrate Grid Component.
      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(MAPL_MetaComp), pointer :: MAPL
      type(ESMF_State) :: internal
      type(wrap_) :: wrap
      type(NI2G_GridComp), pointer :: self
      real, allocatable, dimension(:, :) :: drydepositionfrequency, dqa
      real :: fwet
      logical :: KIN
      real, allocatable, target, dimension(:, :, :) :: fluxoutWT
      real, allocatable, dimension(:, :, :, :) :: aerosol
      real, pointer, dimension(:, :) :: flux_ptr
      real, pointer, dimension(:, :, :) :: fluxWT_ptr
      logical :: alarm_is_ringing
      integer :: i1, j1, i2, j2, km
      real, target, allocatable, dimension(:, :, :) :: RH20, RH80
      integer :: rhFlag
      integer :: i, j
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread
      integer :: settling_opt

#include "NI2G_DeclarePointer___.h"

      __Iam__('Run2')

      !   Get my name and set-up traceback handle
      !   ---------------------------------------
      call ESMF_GridCompGet(gc, NAME=comp_name, _RC)
      Iam = trim(comp_name) // '::' // Iam

      !   Get my internal MAPL_Generic state
      !   -----------------------------------
      call MAPL_GetObjectFromGC(gc, MAPL, _RC)

      !   Get parameters from generic state.
      !   -----------------------------------
      call MAPL_Get(MAPL, INTERNAL_ESMF_STATE=internal, _RC)

#include "NI2G_GetPointer___.h"

      !   Get my private internal state
      !   ------------------------------
      call ESMF_UserCompGetInternalState(gc, 'NI2G_GridComp', wrap, _RC)
      self => wrap%PTR

      !   Set klid and Set internal values to 0 above klid
      !   ---------------------------------------------------
      call findKlid(self%klid, self%plid, ple, _RC)
      call setZeroKlid(self%km, self%klid, NO3an1)
      call setZeroKlid(self%km, self%klid, NO3an2)
      call setZeroKlid(self%km, self%klid, NO3an3)

      allocate(dqa, mold=lwi) !, _STAT)
      allocate(drydepositionfrequency, mold=lwi) !, _STAT)

      alarm_is_ringing = ESMF_AlarmIsRinging(self%alarm, _RC)
#ifdef DEBUG
      if (alarm_is_ringing) then
         if (mapl_am_i_root()) then
            print *, 'DEBUG:: NI replenish alarm is ringing'
         end if
      end if
#endif

      !   Save local copy of HNO3 for first pass through run method regardless
      thread = MAPL_get_current_thread()
      workspace => self%workspaces(thread)

      !if (workspace%first) then
      !xhno3 = MAPL_UNDEF
      !workspace%first = .false.
      !end if

      !   Recycle HNO3 every 3 hours
      if (alarm_is_ringing) then
         xhno3 = NITRATE_HNO3
      end if

      if (associated(NIPNO3AQ)) NIPNO3AQ(:, :) = 0.
      if (associated(NIPNH4AQ)) NIPNH4AQ(:, :) = 0.
      if (associated(NIPNH3AQ)) NIPNH3AQ(:, :) = 0.

      call NIthermo(self%km, self%klid, self%CDT, MAPL_GRAV, delp, airdens, &
           t, rh2, fMassHNO3, MAPL_AIRMW, SO4, NH3, NO3an1, NH4a, &
           xhno3, NIPNO3AQ, NIPNH4AQ, NIPNH3AQ, _RC)

      call NIheterogenousChem(NIHT, xhno3, MAPL_UNDEF, MAPL_AVOGAD, MAPL_AIRMW, &
           MAPL_PI, MAPL_RUNIV / 1000., airdens, t, rh2, delp, DU, &
           ss, self%rmedDU * 1.e-6, self%rmedSS * 1.e-6, &
           self%fnumDU, self%fnumSS, self%km, self%klid, &
           self%CDT, MAPL_GRAV, fMassHNO3, fMassNO3, NO3an1, NO3an2, &
           NO3an3, HNO3CONC, HNO3SMASS, HNO3CMASS, _RC)
      !   Save local copy of HNO3 for first pass through run method regardless

      !   NI Settling
      !   -----------
      select case (self%settling_scheme)
      case ('gocart')
         settling_opt = 1
      case ('ufs')
         settling_opt = 2
      case default
         _ASSERT_RC(.false., 'Unsupported settling scheme: ' // trim(self%settling_scheme), ESMF_RC_NOT_IMPL)
      end select

      !   Ammonium - settles like bin 1 of nitrate
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 1, self%CDT, MAPL_GRAV, &
           NH4a, t, airdens, rh2, zle, delp, NH4SD, settling_scheme=settling_opt, _RC)
      !   Nitrate Bin 1
      nullify(flux_ptr)
      if (associated(NISD)) flux_ptr => NISD(:, :, 1)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 1, self%CDT, MAPL_GRAV, &
           NO3an1, t, airdens, &
           rh2, zle, delp, flux_ptr, settling_scheme=settling_opt, _RC)
      !   Nitrate Bin 2
      nullify(flux_ptr)
      if (associated(NISD)) flux_ptr => NISD(:, :, 2)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 2, self%CDT, MAPL_GRAV, &
           NO3an2, t, airdens, &
           rh2, zle, delp, flux_ptr, settling_scheme=settling_opt, _RC)
      !   Nitrate Bin 3
      nullify(flux_ptr)
      if (associated(NISD)) flux_ptr => NISD(:, :, 3)
      call Chem_SettlingSimple(self%km, self%klid, self%diag_Mie, 3, self%CDT, MAPL_GRAV, &
           NO3an3, t, airdens, &
           rh2, zle, delp, flux_ptr, settling_scheme=settling_opt, _RC)

      !  NI Deposition
      !  -----------
      drydepositionfrequency = 0.
      call DryDeposition(self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
           MAPL_KARMAN, cpd, MAPL_GRAV, z0h, drydepositionfrequency, _RC)

      !   Save local copy of HNO3 for first pass through run method regardless
      !  NH3
      dqa = 0.
      do i = 1, ubound(lwi, 1)
         do j = 1, ubound(lwi, 2)
            if (abs(lwi(i, j) - OCEAN) < 0.5) then
               dqa(i, j) = max(0.0, NH3(i, j, self%km) * (1. - exp(-10.0 * drydepositionfrequency(i, j) * self%CDT)))
            else
               dqa(i, j) = max(0.0, NH3(i, j, self%km) * (1. - exp(-3.0 * drydepositionfrequency(i, j) * self%CDT)))
            end if
         end do
      end do
      !   Save local copy of HNO3 for first pass through run method regardless

      NH3(:, :, self%km) = NH3(:, :, self%km) - dqa
      if (associated(NH3DP)) NH3DP = dqa * delp(:, :, self%km) / MAPL_GRAV / self%CDT

      !  NH4a
      dqa = 0.
      dqa = max(0.0, NH4a(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%CDT)))
      NH4a(:, :, self%km) = NH4a(:, :, self%km) - dqa
      if (associated(NH4DP)) NH4DP = dqa * delp(:, :, self%km) / MAPL_GRAV / self%CDT

      !  NO3anx
      dqa = 0.
      dqa = max(0.0, NO3an1(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%CDT)))
      NO3an1(:, :, self%km) = NO3an1(:, :, self%km) - dqa
      if (associated(NIDP)) NIDP(:, :, 1) = dqa * delp(:, :, self%km) / MAPL_GRAV / self%CDT

      dqa = 0.
      dqa = max(0.0, NO3an2(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%CDT)))
      NO3an2(:, :, self%km) = NO3an2(:, :, self%km) - dqa
      if (associated(NIDP)) NIDP(:, :, 2) = dqa * delp(:, :, self%km) / MAPL_GRAV / self%CDT

      dqa = 0.
      dqa = max(0.0, NO3an3(:, :, self%km) * (1. - exp(-drydepositionfrequency * self%CDT)))
      NO3an3(:, :, self%km) = NO3an3(:, :, self%km) - dqa
      if (associated(NIDP)) NIDP(:, :, 3) = dqa * delp(:, :, self%km) / MAPL_GRAV / self%CDT

      !  NI Large-scale Wet Removal
      !  --------------------------
      if (associated(NH3WT) .or. associated(NH4WT)) then
         allocate(fluxoutWT(ubound(t, 1), ubound(t, 2), 1), _STAT)
      end if
      !  NH3
      KIN = .false.
      fwet = 1.
      nullify(fluxWT_ptr)
      if (associated(NH3WT)) fluxWT_ptr => fluxoutWT
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%CDT, 'NH3', &
           KIN, MAPL_GRAV, fwet, NH3, ple, t, airdens, &
           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, fluxWT_ptr, _RC)
      !   Save local copy of HNO3 for first pass through run method regardless
      if (associated(NH3WT)) NH3WT = fluxWT_ptr(:, :, 1)

      !  NH4a
      KIN = .true.
      fwet = 1.
      nullify(fluxWT_ptr)
      if (associated(NH4WT)) fluxWT_ptr => fluxoutWT
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%CDT, 'NH4a', &
           KIN, MAPL_GRAV, fwet, NH4a, ple, t, airdens, &
           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, fluxWT_ptr, _RC)
      if (associated(NH4WT)) NH4WT = fluxWT_ptr(:, :, 1)

      if (allocated(fluxoutWT)) then
         deallocate(fluxoutWT, _STAT)
      end if

      KIN = .true.
      fwet = 1.
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%CDT, 'nitrate', &
           KIN, MAPL_GRAV, fwet, NO3an1, ple, t, airdens, &
           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, _RC)

      KIN = .true.
      fwet = 1.
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 2, self%CDT, 'nitrate', &
           KIN, MAPL_GRAV, fwet, NO3an2, ple, t, airdens, &
           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, _RC)

      KIN = .true.
      fwet = 0.3
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 3, self%CDT, 'nitrate', &
           KIN, MAPL_GRAV, fwet, NO3an3, ple, t, airdens, &
           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, _RC)

      !   Save local copy of HNO3 for first pass through run method regardless
      !  Compute desired output diagnostics
      !  ----------------------------------
      !  Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
      allocate(aerosol(ubound(NH4a, 1), ubound(NH4a, 2), ubound(NH4a, 3), 3), _STAT)
      aerosol(:, :, :, :) = 0.0
      aerosol(:, :, :, 1) = NH4a
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=1, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple, tropp=tropp,&
           sfcmass=NH4SMASS, colmass=NH4CMASS, mass=NH4MASS, conc=NH4CONC, _RC)
      !   Save local copy of HNO3 for first pass through run method regardless

      aerosol(:, :, :, 1) = NH3
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=1, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple, tropp=tropp,&
           sfcmass=NH3SMASS, colmass=NH3CMASS, mass=NH3MASS, conc=NH3CONC, _RC)
      !   Save local copy of HNO3 for first pass through run method regardless

      aerosol(:, :, :, 1) = NO3an1
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=1, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple, tropp=tropp,&
           sfcmass=NISMASS25, colmass=NICMASS25, mass=NIMASS25, conc=NICONC25, &
           exttau25=NIEXTT25, scatau25=NISCAT25, exttaufm=NIEXTTFM, scataufm=NISCATFM, &
           NO3nFlag=.true., _RC)
      !   Save local copy of HNO3 for first pass through run method regardless

      aerosol(:, :, :, 1) = NO3an1
      aerosol(:, :, :, 2) = NO3an2
      aerosol(:, :, :, 3) = NO3an3
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=3, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, &
           aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
           delp=delp, ple=ple, tropp=tropp, sfcmass=NISMASS, colmass=NICMASS, mass=NIMASS, conc=NICONC, &
           exttau=NIEXTTAU, stexttau=NISTEXTTAU, scatau=NISCATAU, stscatau=NISTSCATAU,&
           fluxu=NIFLUXU, fluxv=NIFLUXV, extcoef=NIEXTCOEF, scacoef=NISCACOEF, &
           bckcoef=NIBCKCOEF, angstrom=NIANGSTR, _RC)

      i1 = lbound(rh2, 1)
      i2 = ubound(rh2, 1)
      j1 = lbound(rh2, 2)
      j2 = ubound(rh2, 2)
      km = ubound(rh2, 3)

      allocate(RH20(i1:i2, j1:j2, km), _STAT)
      allocate(RH80(i1:i2, j1:j2, km), _STAT)

      RH20(:, :, :) = 0.20

      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=3, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, aerosol=aerosol, &
           grav=MAPL_GRAV, tmpu=t, rhoa=airdens, &
           rh=RH20, u=u, v=v, delp=delp, ple=ple, tropp=tropp, &
           extcoef=NIEXTCOEFRH20, scacoef=NISCACOEFRH20, _RC)
      !   Save local copy of HNO3 for first pass through run method regardless

      RH80(:, :, :) = 0.80
      call Aero_Compute_Diags(mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=3, &
           wavelengths_profile=self%wavelengths_profile * 1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint * 1.0e-9, aerosol=aerosol, &
           grav=MAPL_GRAV, tmpu=t, rhoa=airdens, &
           rh=RH80, u=u, v=v, delp=delp, ple=ple, tropp=tropp, &

           extcoef=NIEXTCOEFRH80, scacoef=NISCACOEFRH80, _RC)
      !   Save local copy of HNO3 for first pass through run method regardless

      deallocate(RH20, RH80)

      RETURN_(ESMF_SUCCESS)

   end subroutine Run2

   !BOP
   !IROUTINE: Run_data -- ExtData Nitrate Grid Component
   !INTERFACE:
   subroutine Run_data(gc, import, export, internal, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc ! Gridded component
      type(ESMF_State), intent(inout) :: import ! Import state
      type(ESMF_State), intent(inout) :: export ! Export state
      type(ESMF_State), intent(inout) :: internal ! Interal state
      integer, optional, intent(out) :: rc ! Error code:

      !DESCRIPTION: Updates pointers in Internal state with fields from ExtData.
      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(wrap_) :: wrap
      type(NI2G_GridComp), pointer :: self
      real, pointer, dimension(:, :, :) :: ptr3d_int, ptr3d_imp

      __Iam__('Run_data')

      !   Get my name and set-up traceback handle
      !   ---------------------------------------
      call ESMF_GridCompGet(gc, NAME=comp_name, _RC)
      Iam = trim(comp_name) // '::' // Iam

      !   Get my private internal state
      !   ------------------------------
      call ESMF_UserCompGetInternalState(gc, 'NI2G_GridComp', wrap, _RC)
      self => wrap%PTR

      !   Update interal data pointers with ExtData
      !   -----------------------------------------
      call MAPL_GetPointer(internal, NAME='NO3an1', PTR=ptr3d_int, _RC)
      call MAPL_GetPointer(import, NAME='climNO3an1', PTR=ptr3d_imp, _RC)
      ptr3d_int = ptr3d_imp
      call MAPL_GetPointer(internal, NAME='NO3an2', PTR=ptr3d_int, _RC)
      call MAPL_GetPointer(import, NAME='climNO3an2', PTR=ptr3d_imp, _RC)
      ptr3d_int = ptr3d_imp
      call MAPL_GetPointer(internal, NAME='NO3an3', PTR=ptr3d_int, _RC)
      call MAPL_GetPointer(import, NAME='climNO3an3', PTR=ptr3d_imp, _RC)
      ptr3d_int = ptr3d_imp

      RETURN_(ESMF_SUCCESS)

   end subroutine Run_data

   subroutine aerosol_optics(state, rc)

      !ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      !Local
      integer, parameter :: DP = kind(1.0d0)
      real, dimension(:, :, :), pointer :: ple, rh
      real(kind=DP), dimension(:, :, :), pointer :: var
      real, dimension(:, :, :), pointer :: q
      real, dimension(:, :, :, :), pointer :: q_4d
      integer, allocatable :: opaque_self(:)
      type(c_ptr) :: address
      type(NI2G_GridComp), pointer :: self
      character(len=ESMF_MAXSTR) :: fld_name
      type(ESMF_Field) :: fld
      character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)
      real(kind=DP), dimension(:, :, :), allocatable :: ext_s, ssa_s, asy_s ! (lon:,lat:,lev:)
      real :: x
      integer :: instance
      integer :: n, nbins
      integer :: i1, j1, i2, j2, km
      integer :: band

      integer :: i, j, k

      __Iam__('NI2G::aerosol_optics')

      !   Mie Table instance/index
      !   ------------------------
      call ESMF_AttributeGet(state, NAME='mie_table_instance', VALUE=instance, _RC)

      !   Get aerosol names
      !   -----------------
      call ESMF_AttributeGet(state, NAME='internal_variable_name', itemCount=nbins, _RC)
      allocate(aerosol_names(nbins), _STAT)
      call ESMF_AttributeGet(state, NAME='internal_variable_name', valueList=aerosol_names, _RC)

      !   Radiation band
      !   --------------
      band = 0
      call ESMF_AttributeGet(state, NAME='band_for_aerosol_optics', VALUE=band, _RC)

      !   Pressure at layer edges
      !   ------------------------
      call ESMF_AttributeGet(state, NAME='air_pressure_for_aerosol_optics', VALUE=fld_name, _RC)
      call MAPL_GetPointer(state, ple, trim(fld_name), _RC)

      !    call MAPL_GetPointer (state, ple, 'PLE', _RC)

      i1 = lbound(ple, 1)
      i2 = ubound(ple, 1)
      j1 = lbound(ple, 2)
      j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      !   Relative humidity
      !   -----------------
      call ESMF_AttributeGet(state, NAME='relative_humidity_for_aerosol_optics', VALUE=fld_name, _RC)
      call MAPL_GetPointer(state, rh, trim(fld_name), _RC)

      !    call MAPL_GetPointer (state, rh, 'RH2', _RC)

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

      call mie_(self%rad_Mie, nbins, band, q_4d, rh, ext_s, ssa_s, asy_s, _RC)

      call ESMF_AttributeGet(state, NAME='extinction_in_air_due_to_ambient_aerosol', VALUE=fld_name, _RC)
      if (fld_name /= '') then
         call MAPL_GetPointer(state, var, trim(fld_name), _RC)
         var = ext_s(:, :, :)
      end if

      call ESMF_AttributeGet(state, NAME='single_scattering_albedo_of_ambient_aerosol', VALUE=fld_name, _RC)
      if (fld_name /= '') then
         call MAPL_GetPointer(state, var, trim(fld_name), _RC)
         var = ssa_s(:, :, :)
      end if

      call ESMF_AttributeGet(state, NAME='asymmetry_parameter_of_ambient_aerosol', VALUE=fld_name, _RC)
      if (fld_name /= '') then
         call MAPL_GetPointer(state, var, trim(fld_name), _RC)
         var = asy_s(:, :, :)
      end if

      deallocate(ext_s, ssa_s, asy_s, _STAT)
      deallocate(q_4d, _STAT)

      RETURN_(ESMF_SUCCESS)

   contains

      !    subroutine mie_(mie_table, aerosol_names, nb, offset, q, rh, bext_s, bssa_s, basym_s, rc)
      subroutine mie_(mie, nbins, band, q, rh, bext_s, bssa_s, basym_s, rc)

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

         __Iam__('NI2G::aerosol_optics::mie_')

         bext_s = 0.0d0
         bssa_s = 0.0d0
         basym_s = 0.0d0

         do l = 1, nbins
            !tau is converted to bext
            call mie%Query(band, l, q(:, :, :, l), rh, tau=bext, gasym=gasym, ssa=bssa, _RC)

            bext_s = bext_s + bext ! extinction
            bssa_s = bssa_s + (bssa * bext) ! scattering extinction
            basym_s = basym_s + gasym * (bssa * bext) ! asymetry parameter multiplied by scatering extiction

         end do

         RETURN_(ESMF_SUCCESS)

      end subroutine mie_

   end subroutine aerosol_optics

   subroutine monochromatic_aerosol_optics(state, rc)

      !ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      real, dimension(:, :, :), pointer :: ple, rh
      real, dimension(:, :), pointer :: var
      real, dimension(:, :, :), pointer :: q
      real, dimension(:, :, :, :), pointer :: q_4d
      integer, allocatable :: opaque_self(:)
      type(c_ptr) :: address
      type(NI2G_GridComp), pointer :: self

      character(len=ESMF_MAXSTR) :: fld_name
      type(ESMF_Field) :: fld
      character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)

      real, dimension(:, :, :), allocatable :: tau_s, tau ! (lon:,lat:,lev:)
      real :: x
      integer :: instance
      integer :: n, nbins
      integer :: i1, j1, i2, j2, km
      real :: wavelength
      integer :: i, j, k

      __Iam__('NI2G:: monochromatic_aerosol_optics')

      !   Begin...

      !   Mie Table instance/index
      !   ------------------------
      call ESMF_AttributeGet(state, NAME='mie_table_instance', VALUE=instance, _RC)

      !   Get aerosol names
      !   -----------------
      call ESMF_AttributeGet(state, NAME='internal_variable_name', itemCount=nbins, _RC)
      allocate(aerosol_names(nbins), _STAT)
      call ESMF_AttributeGet(state, NAME='internal_variable_name', valueList=aerosol_names, _RC)

      !   Radiation band
      !   --------------
      call ESMF_AttributeGet(state, NAME='wavelength_for_aerosol_optics', VALUE=wavelength, _RC)

      !   Pressure at layer edges
      !   ------------------------
      call ESMF_AttributeGet(state, NAME='air_pressure_for_aerosol_optics', VALUE=fld_name, _RC)
      call MAPL_GetPointer(state, ple, trim(fld_name), _RC)
      !    call MAPL_GetPointer (state, ple, 'PLE', _RC)

      i1 = lbound(ple, 1)
      i2 = ubound(ple, 1)
      j1 = lbound(ple, 2)
      j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      !   Relative humidity
      !   -----------------
      call ESMF_AttributeGet(state, NAME='relative_humidity_for_aerosol_optics', VALUE=fld_name, _RC)
      call MAPL_GetPointer(state, rh, trim(fld_name), _RC)
      !    call MAPL_GetPointer (state, rh, 'RH2', _RC)

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

      do n = 1, nbins
         call self%diag_Mie%Query(wavelength, n, q_4d(:, :, :, n), rh, tau=tau, _RC)
         tau_s = tau_s + tau
      end do

      call ESMF_AttributeGet(state, NAME='monochromatic_extinction_in_air_due_to_ambient_aerosol', VALUE=fld_name, &
           _RC)
      if (fld_name /= '') then
         call MAPL_GetPointer(state, var, trim(fld_name), _RC)
         var = sum(tau_s, dim=3)
      end if

      deallocate(q_4d, _STAT)

      RETURN_(ESMF_SUCCESS)

   end subroutine monochromatic_aerosol_optics

end module NI2G_GridCompMod
