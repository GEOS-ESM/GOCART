#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: NI2G_GridCompMod - GOCART Nitrate gridded component 

! !INTERFACE:
module NI2G_GridCompMod

! !USES:
   use ESMF
   use MAPL
   use GOCART2G_MieMod 
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process       ! GOCART2G process library
   use GA_EnvironmentMod

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   real, parameter :: fMassHNO3 = 63., fMassNO3 = 62.
   integer, parameter :: nNH3 = 1
   integer, parameter :: nNH4a = 2
   integer, parameter :: nNO3an1 = 3
   integer, parameter :: nNO3an2 = 4
   integer, parameter :: nNO3an3 = 5

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices

real, parameter ::  cpd    = 1004.16
integer, parameter     :: DP = kind(1.0d0)

! !DESCRIPTION: This module implements GOCART's Nitrate (NI) Gridded Component.

! !REVISION HISTORY:
! 01July2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================

!  !Nitrate state
   type, extends(GA_Environment) :: NI2G_GridComp
       logical           :: first
       logical           :: recycle_HNO3 = .false.
       real, allocatable :: xhno3(:,:,:)   ! buffer for NITRATE_HNO3 [kg/(m^2 sec)]
       real, allocatable :: rmedDU(:), rmedSS(:) ! DU and SS radius
       real, allocatable :: fnumDU(:), fnumSS(:) ! DU and SS particles per kg mass
   end type NI2G_GridComp

   type wrap_
      type (NI2G_GridComp), pointer     :: PTR => null()
   end type wrap_

contains

!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices ( GC, RC )

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code

!    DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!     the Initialize and Finalize services to generic versions. It also
!     allocates our instance of a generic state and puts it in the 
!     gridded component (GC). Here we only set the two-stage run method
!     and declare the data services.

!   !REVISION HISTORY:
!   24oct2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

!EOP
!============================================================================
!

!   !Locals
    character (len=ESMF_MAXSTR)                 :: COMP_NAME
    type (ESMF_Config)                          :: cfg
    type (ESMF_Config)                          :: universal_cfg
    type (wrap_)                                :: wrap
    type (NI2G_GridComp), pointer               :: self

    real                                        :: DEFVAL
    logical                                     :: data_driven=.true.

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, config=universal_cfg, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file 
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'NI2G_instance_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
       if (mapl_am_i_root()) print*,'NI2G_instance_'//trim(COMP_NAME)//'.rc does not exist! Loading NI2G_instance_NI.rc instead'
       call ESMF_ConfigLoadFile (cfg, 'NI2G_instance_NI.rc', __RC__)
    end if

    ! process generic config items
    call self%GA_Environment%load_from_config( cfg, universal_cfg, __RC__)

!   Is NI data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, __RC__)
    if (data_driven .neqv. .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if

    DEFVAL = 0.0

!   Import and Internal states if data instance 
!   -------------------------------------------
    if (data_driven) then

!      Pressure at layer edges
!      -----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'PLE',                                 &
          LONG_NAME  = 'air_pressure',                        &
          UNITS      = 'Pa',                                  &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationEdge,                    &
          RESTART    = MAPL_RestartSkip,     __RC__)

!      RH: is between 0 and 1
!      ----------------------
       call MAPL_AddImportSpec(GC,                            &
          SHORT_NAME = 'RH2',                                 &
          LONG_NAME  = 'Rel_Hum_after_moist',                 &
          UNITS      = '1',                                   &
          DIMS       = MAPL_DimsHorzVert,                     &
          VLOCATION  = MAPL_VLocationCenter,                  &
          RESTART    = MAPL_RestartSkip,     __RC__)

       call MAPL_AddInternalSpec(gc,&
          & short_name='NO3an1', &
          & long_name='Nitrate size bin 001', &
          & units='kg kg-1', &
          & dims=MAPL_DimsHorzVert, &
          & vlocation=MAPL_VlocationCenter, &
          & restart=MAPL_RestartOptional, &
          & friendlyto='DYNAMICS:TURBULENCE:MOIST', &
          & add2export=.true., __RC__)

       call MAPL_AddInternalSpec(gc,&
          & short_name='NO3an2', &
          & long_name='Nitrate size bin 002', &
          & units='kg kg-1', &
          & dims=MAPL_DimsHorzVert, &
          & vlocation=MAPL_VlocationCenter, &
          & restart=MAPL_RestartOptional, &
          & friendlyto='DYNAMICS:TURBULENCE:MOIST', &
          & add2export=.true., __RC__)

       call MAPL_AddInternalSpec(gc,&
          & short_name='NO3an3', &
          & long_name='Nitrate size bin 003', &
          & units='kg kg-1', &
          & dims=MAPL_DimsHorzVert, &
          & vlocation=MAPL_VlocationCenter, &
          & restart=MAPL_RestartOptional, &
          & friendlyto='DYNAMICS:TURBULENCE:MOIST', &
          & add2export=.true., __RC__)

       call MAPL_AddImportSpec(gc,&
          short_name='climNO3an1', &
          long_name='Nitrate size bin 001', &
          units='kg kg-1', &
          dims=MAPL_DimsHorzVert, &
          vlocation=MAPL_VlocationCenter, &
          restart=MAPL_RestartOptional, __RC__)

       call MAPL_AddImportSpec(gc,&
          short_name='climNO3an2', &
          long_name='Nitrate size bin 002', &
          units='kg kg-1', &
          dims=MAPL_DimsHorzVert, &
          vlocation=MAPL_VlocationCenter, &
          restart=MAPL_RestartOptional, __RC__)

       call MAPL_AddImportSpec(gc,&
          short_name='climNO3an3', &
          long_name='Nitrate size bin 003', &
          units='kg kg-1', &
          dims=MAPL_DimsHorzVert, &
          vlocation=MAPL_VlocationCenter, &
          restart=MAPL_RestartOptional, __RC__)
    end if ! (data_driven)


!   Import, Export, Internal states for computational instance 
!   ----------------------------------------------------------
    if (.not. data_driven) then
#include "NI2G_Export___.h"
#include "NI2G_Import___.h"
#include "NI2G_Internal___.h"
    end if


!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec(GC,                                 &
      SHORT_NAME = trim(COMP_NAME)//'_AERO',                   &
       LONG_NAME  = 'aerosols_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                  &
       DIMS       = MAPL_DimsHorzVert,                          &
       VLOCATION  = MAPL_VLocationCenter,                       &
       DATATYPE   = MAPL_StateItem, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'NI2G_GridComp', wrap, STATUS )
    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine SetServices


!============================================================================
!BOP

! !IROUTINE: Initialize 

! !INTERFACE:
  subroutine Initialize (GC, IMPORT, EXPORT, CLOCK, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT ! Export state
    type (ESMF_Clock),    intent(inout) :: CLOCK  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes the Nitrate gridded component. It primaryily 
!               fills GOCART's AERO states with its nitrate fields. 

! !REVISION HISTORY: 
! 30June2020   E.Sherman  First attempt at refactoring

!EOP
!============================================================================
!   !Locals 
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),      pointer   :: MAPL
    type (ESMF_Config)                   :: universal_cfg
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (NI2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: prefix
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)
    real, pointer, dimension(:,:,:)      :: int_ptr
    logical                              :: data_driven
    logical                              :: bands_are_present
    integer                              :: NUM_BANDS
    character (len=ESMF_MAXSTR), allocatable    :: aerosol_names(:)
    real, pointer, dimension(:,:,:)      :: ple

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Alarm)        :: alarm_HNO3
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss

    real, dimension(4)   :: Vect_Hcts
!    real, allocatable, dimension(:) :: rmedDU, rmedSS, fnumDU, fnumSS
    integer :: itemCount
    integer, allocatable, dimension(:)   :: channels_
    integer                              :: nmom_
    character(len=ESMF_MAXSTR)           :: file_
    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, config=universal_cfg, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)

!   Get my internal MAPL_Generic state
!   ----------------------------------- 
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)
    
!   Get my internal private state  
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, globalCellCountPerDim=dims, __RC__ )
    km = dims(3)
    self%km = km

    allocate(self%xhno3(dims(1),dims(2),dims(3)), __STAT__)

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!  Load resource file and get number of bins 
!  -------------------------------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'NI2G_instance_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'NI2G_instance_'//trim(COMP_NAME)//'.rc does not exist! &
                                    loading NI2G_instance_NI.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'NI2G_instance_NI.rc', __RC__)
    end if

    self%first = .true.

!   Call Generic Initialize 
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( mapl, INTERNAL_ESMF_STATE = internal, __RC__)

!   Is NI data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Get DU and SS attribute information for use in heterogenous chemistry
    if (.not. data_driven) then
       call ESMF_StateGet(import, 'DU', field, __RC__)
       call ESMF_AttributeGet(field, name='radius', itemCount=itemCount, __RC__)
       allocate(self%rmedDU(itemCount), __STAT__)
       allocate(self%fnumDU(itemCount), __STAT__)
       call ESMF_AttributeGet(field, name='radius', valueList=self%rmedDU, __RC__)
       call ESMF_AttributeGet(field, name='fnum', valueList=self%fnumDU, __RC__)

       call ESMF_StateGet(import, 'SS', field, __RC__)
       call ESMF_AttributeGet(field, name='radius', itemCount=itemCount, __RC__)
       allocate(self%rmedSS(itemCount), __STAT__)
       allocate(self%fnumSS(itemCount), __STAT__)
       call ESMF_AttributeGet(field, name='radius', valueList=self%rmedSS, __RC__)
       call ESMF_AttributeGet(field, name='fnum', valueList=self%fnumSS, __RC__)
    end if

!   Se HNO3 recycle alarm
    if (.not. data_driven) then
        call ESMF_ClockGet(clock, calendar=calendar, currTime=currentTime, __RC__)
        call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, __RC__)
        call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, __RC__)
        call ESMF_TimeIntervalSet(ringInterval, H=3, calendar=calendar, __RC__)

        do while (ringTime < currentTime)! DO WE NEED THIS?
            ringTime = currentTime + ringInterval
        end do

        alarm_HNO3 = ESMF_AlarmCreate(Clock        = clock,        &
                                      Name         = 'HNO3_RECYCLE_ALARM', &
                                      RingInterval = ringInterval, &
                                      RingTime     = currentTime,  &
                                      Enabled      = .true.   ,    &
                                      Sticky       = .false.  , __RC__)
    end if

!   If this is a data component, the data is provided in the import
!   state via ExtData instead of the actual GOCART children
!   ----------------------------------------------------------------
    if ( data_driven ) then
       providerState = import
       prefix = 'clim'
    else
       providerState = export
       prefix = ''
    end if

    if (.not. data_driven) then
       call ESMF_StateGet (internal, 'NH3', field, __RC__)
       call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), __RC__)
       call get_HenrysLawCts('NH3',Vect_Hcts(1),Vect_Hcts(2),Vect_Hcts(3),Vect_Hcts(4),__RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, __RC__)

       call ESMF_StateGet (internal, 'NH4a', field, __RC__)
       call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(2), __RC__)
!    end if

!      Set klid
       call MAPL_GetPointer(import, ple, 'PLE', __RC__)
       call findKlid (self%klid, self%plid, ple, __RC__)
    end if

!   Fill AERO State with N03an(1,2,3) fields
!   ----------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)

    call ESMF_StateGet (internal, 'NO3an1', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(3), __RC__)
    fld = MAPL_FieldCreate (field, 'NO3an1', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    if (.not. data_driven) then
!      Set internal NO3an1 values to 0 where above klid
       call MAPL_GetPointer (internal, int_ptr, 'NO3an1', __RC__)
       call setZeroKlid(self%km, self%klid, int_ptr)
    end if

    call ESMF_StateGet (internal, 'NO3an2', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(4), __RC__)
    fld = MAPL_FieldCreate (field, 'NO3an2', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    if (.not. data_driven) then
!      Set internal NO3an2 values to 0 where above klid
       call MAPL_GetPointer (internal, int_ptr, 'NO3an2', __RC__)
       call setZeroKlid(self%km, self%klid, int_ptr)
    end if

    call ESMF_StateGet (internal, 'NO3an3', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(5), __RC__)
    fld = MAPL_FieldCreate (field, 'NO3an3', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    if (.not. data_driven) then
!      Set internal NO3an3 values to 0 where above klid
       call MAPL_GetPointer (internal, int_ptr, 'NO3an3', __RC__)
       call setZeroKlid(self%km, self%klid, int_ptr)
    end if

    if (data_driven) then
       instance = instanceData
    else
       instance = instanceComputational
    end if

    self%instance = instance
    
!   Create Radiation Mie Table
!   --------------------------
    call ESMF_ConfigGetAttribute (cfg, file_, label="aerosol_radBands_optics_file:", __RC__ )
    self%rad_Mie = GOCART2G_Mie(trim(file_), __RC__)

!   Create Diagnostics Mie Table
!   -----------------------------
!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, file_, &
                                  label="aerosol_monochromatic_optics_file:", __RC__ )
    call ESMF_ConfigGetAttribute (cfg, nmom_, label="n_moments:", default=0,  __RC__)
    i = ESMF_ConfigGetLen (universal_cfg, label='aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:', __RC__)
    allocate (channels_(i), __STAT__ )
    call ESMF_ConfigGetAttribute (universal_cfg, channels_, &
                                  label= "aerosol_monochromatic_optics_wavelength_in_nm_from_LUT:", __RC__)
    self%diag_Mie = GOCART2G_Mie(trim(file_), channels_*1.e-9, nmom=nmom_, __RC__)
    deallocate(channels_)

    ! Mie Table instance/index
    call ESMF_AttributeSet(aero, name='mie_table_instance', value=instance, __RC__)

    ! Add variables to NI instance's aero state. This is used in aerosol optics calculations
    call add_aero (aero, label='air_pressure_for_aerosol_optics',             label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='relative_humidity_for_aerosol_optics',        label2='RH',  grid=grid, typekind=MAPL_R4,__RC__)
!   call ESMF_StateGet (import, 'PLE', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)
!   call ESMF_StateGet (import, 'RH2', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)

    call add_aero (aero, label='extinction_in_air_due_to_ambient_aerosol',    label2='EXT', grid=grid, typekind=MAPL_R8,__RC__)
    call add_aero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R8,__RC__)
    call add_aero (aero, label='asymmetry_parameter_of_ambient_aerosol',      label2='ASY', grid=grid, typekind=MAPL_R8,__RC__)
    call add_aero (aero, label='monochromatic_extinction_in_air_due_to_ambient_aerosol', &
                   label2='monochromatic_EXT', grid=grid, typekind=MAPL_R4,__RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol', label2='aerosolSum', grid=grid, typekind=MAPL_R4, __RC__)

    call ESMF_AttributeSet (aero, name='band_for_aerosol_optics', value=0, __RC__)
    call ESMF_AttributeSet (aero, name='wavelength_for_aerosol_optics', value=0., __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet (aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

    allocate(aerosol_names(3), __STAT__)
    aerosol_names(1) = 'NO3an1'
    aerosol_names(2) = 'NO3an2'
    aerosol_names(3) = 'NO3an3'
    call ESMF_AttributeSet (aero, name='internal_variable_name', valueList=aerosol_names, &
                            itemCount=size(aerosol_names), __RC__)

    call ESMF_MethodAdd (aero, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)
    call ESMF_MethodAdd (aero, label='monochromatic_aerosol_optics', userRoutine=monochromatic_aerosol_optics, __RC__)
    call ESMF_MethodAdd (aero, label='get_mixR', userRoutine=get_mixR, __RC__)

    RETURN_(ESMF_SUCCESS)

  end subroutine Initialize

!============================================================================

!BOP
! !IROUTINE: Run 

! !INTERFACE:
  subroutine Run (GC, import, export, clock, rc)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: rc     ! Error code:

! !DESCRIPTION: Run method for the Sea Salt Grid Component. Determines whether to run
!               data or computational run method.

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal

    logical                           :: data_driven

    __Iam__('Run')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is NI data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Update INTERNAL state variables with ExtData
!   ---------------------------------------------
    if (data_driven) then
       call Run_data (GC, import, export, internal, __RC__)
    else
       call Run1 (GC, import, export, clock, __RC__)
    end if

    RETURN_(ESMF_SUCCESS)

  end subroutine Run

!============================================================================
!BOP
! !IROUTINE: Run1 

! !INTERFACE:
  subroutine Run1 (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION:  Computes emissions/sources for Nitrate

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (NI2G_GridComp), pointer     :: self

#include "NI2G_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=comp_name, __RC__)
    Iam = trim(comp_name) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, __RC__)

#include "NI2G_GetPointer___.h"

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

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
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*MAPL_GRAV/delp(:,:,self%km)*EMI_NH3_BB
    if (associated(EMI_NH3_AG)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*MAPL_GRAV/delp(:,:,self%km)*EMI_NH3_AG
    if (associated(EMI_NH3_EN)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*MAPL_GRAV/delp(:,:,self%km)*EMI_NH3_EN
    if (associated(EMI_NH3_IN)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*MAPL_GRAV/delp(:,:,self%km)*EMI_NH3_IN
    if (associated(EMI_NH3_RE)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*MAPL_GRAV/delp(:,:,self%km)*EMI_NH3_RE
    if (associated(EMI_NH3_TR)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*MAPL_GRAV/delp(:,:,self%km)*EMI_NH3_TR
    if (associated(EMI_NH3_OC)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*MAPL_GRAV/delp(:,:,self%km)*EMI_NH3_OC

    RETURN_(ESMF_SUCCESS)

  end subroutine Run1


!============================================================================
!BOP
! !IROUTINE: Run2 

! !INTERFACE:

  subroutine Run2 (GC, import, export, clock, RC)

    ! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code:

! !DESCRIPTION: Run2 method for the Dust Grid Component.

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: MAPL
    type (ESMF_State)                 :: internal
    type (wrap_)                      :: wrap
    type (NI2G_GridComp), pointer     :: self

    real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
    real                              :: fwet
    logical                           :: KIN
    real, allocatable, target, dimension(:,:,:) :: fluxoutWT
    real, allocatable, dimension(:,:,:,:) :: aerosol
    real, pointer, dimension(:,:)     :: flux_ptr
    real, pointer, dimension(:,:,:)   :: fluxWT_ptr

    type (ESMF_ALARM)               :: alarm
    logical                         :: alarm_is_ringing

    integer :: rhFlag
    integer :: i, j

#include "NI2G_DeclarePointer___.h"

    __Iam__('Run2')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=INTERNAL, __RC__)

#include "NI2G_GetPointer___.h"

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    allocate(dqa, mold=lwi, __STAT__)
    allocate(drydepositionfrequency, mold=lwi, __STAT__)

!   check hno3 alarm
    call ESMF_ClockGetAlarm(clock, 'HNO3_RECYCLE_ALARM', alarm, __RC__)
    alarm_is_ringing = ESMF_AlarmIsRinging(alarm, __RC__)

!   Save local copy of HNO3 for first pass through run method regardless
    if (self%first) then
       self%xhno3 = MAPL_UNDEF
       self%first = .false.
    end if

!   Recycle HNO3 every 3 hours
    if (alarm_is_ringing) then
       self%xhno3 = NITRATE_HNO3
       call ESMF_AlarmRingerOff(alarm, __RC__)
    end if

    if (associated(NIPNO3AQ)) NIPNO3AQ(:,:) = 0.
    if (associated(NIPNH4AQ)) NIPNH4AQ(:,:) = 0.
    if (associated(NIPNH3AQ)) NIPNH3AQ(:,:) = 0.

    call NIthermo (self%km, self%klid, self%cdt, MAPL_GRAV, delp, airdens, &
                   t, rh2, fMassHNO3, MAPL_AIRMW, SO4, NH3, NO3an1, NH4a, &
                   self%xhno3, NIPNO3AQ, NIPNH4AQ, NIPNH3AQ, __RC__)


    call NIheterogenousChem (NIHT, self%xhno3, MAPL_UNDEF, MAPL_AVOGAD, MAPL_AIRMW, &
                             MAPL_PI, MAPL_RUNIV/1000., airdens, t, rh2, delp, DU, &
                             SS, self%rmedDU*1.e-6, self%rmedSS*1.e-6, &
                             self%fnumDU, self%fnumSS, 5, 5, self%km, self%klid, &
                             self%cdt, MAPL_GRAV, fMassHNO3, fMassNO3, NO3an1, NO3an2, & 
                             NO3an3, HNO3CONC, HNO3SMASS,  HNO3CMASS, __RC__)


!   NI Settling
!   -----------
!   Because different bins having different swelling coefficients I need to
!   handle the call to settling differently.

!   Ammonium - settles like ammonium sulfate (rhflag = 3)
    rhflag = 3
!    call Chem_SettlingSimpleOrig (self%km, self%klid, rhflag, MAPL_GRAV, self%cdt, &
!                                  1.e-6*self%radius(nNH4a), self%rhop(nNH4a), &
!                                  NH4a, t, airdens, rh2, delp, zle, NH4SD, __RC__)
    call Chem_SettlingSimple (self%km, self%klid, rhFlag, self%cdt, MAPL_GRAV, &
                              self%radius(nNH4a)*1.e-6, self%rhop(nNH4a), NH4a, t, &
                              airdens, rh2, zle, delp, NH4SD, __RC__)

!  Nitrate bin 1 - settles like ammonium sulfate (rhflag = 3)
    rhflag = 3
    nullify(flux_ptr)
    if (associated(NISD)) flux_ptr => NISD(:,:,1)
!    call Chem_SettlingSimpleOrig (self%km, self%klid, rhFlag, MAPL_GRAV, self%cdt, &
!                                  1.e-6*self%radius(nNO3an1), self%rhop(nNO3an1), &
!                                  NO3an1, t, airdens, rh2, delp, zle, flux_ptr, __RC__)
    call Chem_SettlingSimple (self%km, self%klid, rhFlag, self%cdt, MAPL_GRAV, &
                              self%radius(nNO3an1)*1.e-6, self%rhop(nNO3an1), NO3an1, &
                              t, airdens, rh2, zle, delp, flux_ptr, __RC__)

!  Nitrate bin 2 - settles like sea salt (rhflag = 2)
    rhflag = 2
    nullify(flux_ptr)
    if (associated(NISD)) flux_ptr => NISD(:,:,2)
!    call Chem_SettlingSimpleOrig (self%km, self%klid, rhFlag, MAPL_GRAV, self%cdt, &
!                                  1.e-6*self%radius(nNO3an2), self%rhop(nNO3an2), &
!                                  NO3an2, t, airdens, rh2, delp, zle, flux_ptr, __RC__)
    call Chem_SettlingSimple (self%km, self%klid, rhFlag, self%cdt, MAPL_GRAV, &
                              self%radius(nNO3an2)*1.e-6, self%rhop(nNO3an2), NO3an2, &
                              t, airdens, rh2, zle, delp, flux_ptr, __RC__)

!  Nitrate bin 1 - settles like dust (rhflag = 0)
    rhflag = 0
    nullify(flux_ptr)
    if (associated(NISD)) flux_ptr => NISD(:,:,3)
!    call Chem_SettlingSimpleOrig (self%km, self%klid, rhFlag, MAPL_GRAV, self%cdt, &
!                                  1.e-6*self%radius(nNO3an3), self%rhop(nNO3an3), &
!                                  NO3an3, t, airdens, rh2, delp, zle, flux_ptr, __RC__)
    call Chem_SettlingSimple (self%km, self%klid, rhFlag, self%cdt, MAPL_GRAV, &
                              self%radius(nNO3an3)*1.e-6, self%rhop(nNO3an3), NO3an3, &
                              t, airdens, rh2, zle, delp, flux_ptr, __RC__)

!  NI Deposition
!  -----------
    drydepositionfrequency = 0.
    call DryDeposition(self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
                       MAPL_KARMAN, cpd, MAPL_GRAV, z0h, drydepositionfrequency, __RC__ )

!  NH3
   dqa = 0.
   do i=1,ubound(lwi,1)
      do j =1,ubound(lwi,2)
         if (abs(lwi(i,j) - OCEAN) < 0.5) then
            dqa(i,j) = max(0.0, NH3(i,j,self%km)*(1.-exp(-10.0*drydepositionfrequency(i,j)*self%cdt)))
         else
            dqa(i,j) = max(0.0, NH3(i,j,self%km)*(1.-exp( -3.0*drydepositionfrequency(i,j)*self%cdt)))
         end if
      end do
   end do

   NH3(:,:,self%km) = NH3(:,:,self%km) - dqa
   if( associated(NH3DP) ) NH3DP = dqa*delp(:,:,self%km)/MAPL_GRAV/self%cdt

!  NH4a
   dqa = 0.
   dqa = max(0.0, NH4a(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NH4a(:,:,self%km) = NH4a(:,:,self%km) - dqa
   if( associated(NH4DP) ) NH4DP = dqa*delp(:,:,self%km)/MAPL_GRAV/self%cdt

!  NO3anx
   dqa = 0.
   dqa = max(0.0, NO3an1(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NO3an1(:,:,self%km) = NO3an1(:,:,self%km) - dqa
   if( associated(NIDP) ) NIDP(:,:,1) = dqa*delp(:,:,self%km)/MAPL_GRAV/self%cdt

   dqa = 0.
   dqa = max(0.0, NO3an2(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NO3an2(:,:,self%km) = NO3an2(:,:,self%km) - dqa
   if( associated(NIDP) ) NIDP(:,:,2) = dqa*delp(:,:,self%km)/MAPL_GRAV/self%cdt

   dqa = 0.
   dqa = max(0.0, NO3an3(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NO3an3(:,:,self%km) = NO3an3(:,:,self%km) - dqa
   if( associated(NIDP) ) NIDP(:,:,3) = dqa*delp(:,:,self%km)/MAPL_GRAV/self%cdt

!  NI Large-scale Wet Removal
!  --------------------------
   if (associated(NH3WT) .or. associated(NH4WT)) then
      allocate(fluxoutWT(ubound(t,1), ubound(t,2), 1), __STAT__)
   end if
!  NH3
   KIN = .false.
   fwet = 1.
   nullify(fluxWT_ptr)
   if (associated(NH3WT)) fluxWT_ptr => fluxoutWT
   call WetRemovalGOCART2G (self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'NH3', &
                            KIN, MAPL_GRAV, fwet, NH3, ple, t, airdens, &
                            pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, fluxWT_ptr, __RC__)
   if (associated(NH3WT)) NH3WT = fluxWT_ptr(:,:,1)

!  NH4a
   KIN = .true.
   fwet = 1.
   nullify(fluxWT_ptr)
   if (associated(NH4WT)) fluxWT_ptr => fluxoutWT
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'NH4a', &
                           KIN, MAPL_GRAV, fwet, NH4a, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, fluxWT_ptr, __RC__)
   if (associated(NH4WT)) NH4WT = fluxWT_ptr(:,:,1)

   if (allocated(fluxoutWT)) then
      deallocate(fluxoutWT, __STAT__)
   end if

   KIN = .true.
   fwet = 1.
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'nitrate', &
                           KIN, MAPL_GRAV, fwet, NO3an1, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, __RC__)

   KIN = .true.
   fwet = 1.
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 2, self%cdt, 'nitrate', &
                           KIN, MAPL_GRAV, fwet, NO3an2, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, __RC__)

   KIN = .true.
   fwet = 0.3
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 3, self%cdt, 'nitrate', &
                           KIN, MAPL_GRAV, fwet, NO3an3, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, __RC__)

!  Compute desired output diagnostics
!  ----------------------------------
!  Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
   allocate(aerosol(ubound(NH4a,1), ubound(NH4a,2), ubound(NH4a,3), 3), __STAT__)
   aerosol(:,:,:,:) = 0.0
   aerosol(:,:,:,1) = NH4a
   call Aero_Compute_Diags (mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
                            nbins=1, &
                            wavelengths_profile=self%wavelengths_profile*1.0e-9, &
                            wavelengths_vertint=self%wavelengths_vertint*1.0e-9, &
                            aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, ple=ple, tropp=tropp,&
                            sfcmass=NH4SMASS, colmass=NH4CMASS, mass=NH4MASS, conc=NH4CONC, __RC__)

   aerosol(:,:,:,1) = NH3
   call Aero_Compute_Diags (mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
                            nbins=1, &
                            wavelengths_profile=self%wavelengths_profile*1.0e-9, &
                            wavelengths_vertint=self%wavelengths_vertint*1.0e-9, &
                            aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, ple=ple, tropp=tropp,&
                            sfcmass=NH3SMASS, colmass=NH3CMASS, mass=NH3MASS, conc=NH3CONC, __RC__)

   aerosol(:,:,:,1) = NO3an1
   call Aero_Compute_Diags (mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
                            nbins=1, &
                            wavelengths_profile=self%wavelengths_profile*1.0e-9, &
                            wavelengths_vertint=self%wavelengths_vertint*1.0e-9, &
                            aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, ple=ple, tropp=tropp,&
                            sfcmass=NISMASS25, colmass=NICMASS25, mass=NIMASS25, conc=NICONC25, &
                            exttau25=NIEXTT25, scatau25=NISCAT25, exttaufm=NIEXTTFM, scataufm=NISCATFM, &
                            NO3nFlag=.true., __RC__)

   aerosol(:,:,:,1) = NO3an1
   aerosol(:,:,:,2) = NO3an2
   aerosol(:,:,:,3) = NO3an3
   call Aero_Compute_Diags (mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
                            nbins=3,  &
                            wavelengths_profile=self%wavelengths_profile*1.0e-9, &
                            wavelengths_vertint=self%wavelengths_vertint*1.0e-9, &
                            aerosol=aerosol, grav=MAPL_GRAV, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, ple=ple, tropp=tropp,sfcmass=NISMASS, colmass=NICMASS, mass=NIMASS, conc=NICONC, &
                            exttau=NIEXTTAU, stexttau=NISTEXTTAU,scatau=NISCATAU, stscatau=NISTSCATAU,&
                            fluxu=NIFLUXU, fluxv=NIFLUXV, extcoef=NIEXTCOEF, scacoef=NISCACOEF, &
                            angstrom=NIANGSTR, __RC__ )

    RETURN_(ESMF_SUCCESS)
  
  end subroutine Run2


!============================================================================
!BOP
! !IROUTINE: Run_data -- ExtData Sea Salt Grid Component

! !INTERFACE:

  subroutine Run_data (GC, IMPORT, EXPORT, INTERNAL, RC)

    ! !ARGUMENTS:

    type (ESMF_GridComp), intent(inout) :: GC       ! Gridded component 
    type (ESMF_State),    intent(inout) :: IMPORT   ! Import state
    type (ESMF_State),    intent(inout) :: EXPORT   ! Export state
    type (ESMF_State),    intent(inout) :: INTERNAL ! Interal state
    integer, optional,    intent(  out) :: RC       ! Error code:

! !DESCRIPTION: Updates pointers in Internal state with fields from ExtData. 

!EOP
!============================================================================
! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (wrap_)                       :: wrap
    type (NI2G_GridComp), pointer      :: self

    real, pointer, dimension(:,:,:)  :: ptr3d_int, ptr3d_imp

    __Iam__('Run_data')

!*****************************************************************************
! Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'//Iam

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Update interal data pointers with ExtData
!   -----------------------------------------
    call MAPL_GetPointer (internal, name='NO3an1', ptr=ptr3d_int, __RC__)
    call MAPL_GetPointer (import, name='climNO3an1', ptr=ptr3d_imp, __RC__)
    ptr3d_int = ptr3d_imp
    call MAPL_GetPointer (internal, name='NO3an2', ptr=ptr3d_int, __RC__)
    call MAPL_GetPointer (import, name='climNO3an2', ptr=ptr3d_imp, __RC__)
    ptr3d_int = ptr3d_imp
    call MAPL_GetPointer (internal, name='NO3an3', ptr=ptr3d_int, __RC__)
    call MAPL_GetPointer (import, name='climNO3an3', ptr=ptr3d_imp, __RC__)
    ptr3d_int = ptr3d_imp

    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data

!-------------------------------------------------------------------------------------

 subroutine aerosol_optics(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    integer, parameter                               :: DP=kind(1.0d0)
    real, dimension(:,:,:), pointer                  :: ple, rh
    real(kind=DP), dimension(:,:,:), pointer         :: var
    real, dimension(:,:,:), pointer                  :: q
    real, dimension(:,:,:,:), pointer                :: q_4d
    integer, allocatable                             :: opaque_self(:)
    type(C_PTR)                                      :: address
    type(NI2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR),allocatable          :: aerosol_names(:)

    real(kind=DP), dimension(:,:,:), allocatable     :: ext_s, ssa_s, asy_s  ! (lon:,lat:,lev:)
    real                                             :: x
    integer                                          :: instance
    integer                                          :: n, nbins
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band

    integer :: i, j, k

    __Iam__('NI2G::aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

!   Get aerosol names
!   -----------------
    call ESMF_AttributeGet (state, name='internal_variable_name', itemCount=nbins, __RC__)
    allocate (aerosol_names(nbins), __STAT__)
    call ESMF_AttributeGet (state, name='internal_variable_name', valueList=aerosol_names, __RC__)

!   Radiation band
!   --------------
    band = 0
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)

!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, ple, 'PLE', __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, rh, 'RH2', __RC__)

    allocate(ext_s(i1:i2, j1:j2, km), &
             ssa_s(i1:i2, j1:j2, km), &
             asy_s(i1:i2, j1:j2, km), __STAT__)

    allocate(q_4d(i1:i2, j1:j2, km, nbins), __STAT__)

    do n = 1, nbins
       call ESMF_StateGet (state, trim(aerosol_names(n)), field=fld, __RC__)
       call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

        do k = 1, km
           do j = j1, j2
              do i = i1, i2
                 x = ((ple(i,j,k) - ple(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                 q_4d(i,j,k,n) = x * q(i,j,k)
              end do
           end do
        end do
    end do

    call ESMF_AttributeGet(state, name='mieTable_pointer', itemCount=n, __RC__)
    allocate (opaque_self(n), __STAT__)
    call ESMF_AttributeGet(state, name='mieTable_pointer', valueList=opaque_self, __RC__)

    address = transfer(opaque_self, address)
    call c_f_pointer(address, self)

    call mie_ (self%rad_Mie, nbins, band, q_4d, rh, ext_s, ssa_s, asy_s, __RC__)

    call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ext_s(:,:,:)
    end if

    call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = ssa_s(:,:,:)
    end if

   call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = asy_s(:,:,:)
    end if

    deallocate(ext_s, ssa_s, asy_s, __STAT__)
    deallocate(q_4d, __STAT__)

    RETURN_(ESMF_SUCCESS)

  contains

!    subroutine mie_(mie_table, aerosol_names, nb, offset, q, rh, bext_s, bssa_s, basym_s, rc)
    subroutine mie_(mie, nbins, band, q, rh, bext_s, bssa_s, basym_s, rc)

    implicit none

    type(GOCART2G_Mie) ,           intent(inout) :: mie              ! mie table
    integer,                       intent(in   ) :: nbins            ! number of bins
    integer,                       intent(in )   :: band             ! channel
    real,                          intent(in )   :: q(:,:,:,:)       ! aerosol mass mixing ratio, kg kg-1
    real,                          intent(in )   :: rh(:,:,:)        ! relative humidity
    real(kind=8), intent(  out) :: bext_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=8), intent(  out) :: bssa_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=8), intent(  out) :: basym_s(size(ext_s,1),size(ext_s,2),size(ext_s,3))
    integer,       intent(out)  :: rc
    ! local
    integer                     :: l
    real                        :: bext (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! extinction
    real                        :: bssa (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! SSA
    real                        :: gasym(size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! asymmetry parameter

    __Iam__('NI2G::aerosol_optics::mie_')

    bext_s  = 0.0d0
    bssa_s  = 0.0d0
    basym_s = 0.0d0

    do l = 1, nbins
       !tau is converted to bext
       call mie%Query(band, l, q(:,:,:,l), rh, tau=bext, gasym=gasym, ssa=bssa, __RC__)

       bext_s  = bext_s  +             bext     ! extinction
       bssa_s  = bssa_s  +       (bssa*bext)    ! scattering extinction
       basym_s = basym_s + gasym*(bssa*bext)    ! asymetry parameter multiplied by scatering extiction 

    end do

    RETURN_(ESMF_SUCCESS)

    end subroutine mie_

  end subroutine aerosol_optics

!-------------------------------------------------------------------------------------

 subroutine  monochromatic_aerosol_optics(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    real, dimension(:,:,:), pointer                  :: ple, rh
    real, dimension(:,:), pointer                    :: var
    real, dimension(:,:,:), pointer                  :: q
    real, dimension(:,:,:,:), pointer                :: q_4d
    integer, allocatable                             :: opaque_self(:)
    type(C_PTR)                                      :: address
    type(NI2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR),allocatable          :: aerosol_names(:)

    real, dimension(:,:,:), allocatable              :: tau_s, tau  ! (lon:,lat:,lev:)
    real                                             :: x
    integer                                          :: instance
    integer                                          :: n, nbins
    integer                                          :: i1, j1, i2, j2, km
    real                                             :: wavelength
    integer :: i, j, k

    __Iam__('NI2G:: monochromatic_aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

!   Get aerosol names
!   -----------------
    call ESMF_AttributeGet (state, name='internal_variable_name', itemCount=nbins, __RC__)
    allocate (aerosol_names(nbins), __STAT__)
    call ESMF_AttributeGet (state, name='internal_variable_name', valueList=aerosol_names, __RC__)

!   Radiation band
!   --------------
    call ESMF_AttributeGet(state, name='wavelength_for_aerosol_optics', value=wavelength, __RC__)


!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)
!    call MAPL_GetPointer (state, ple, 'PLE', __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)
!    call MAPL_GetPointer (state, rh, 'RH2', __RC__)

    allocate(tau_s(i1:i2, j1:j2, km), &
               tau(i1:i2, j1:j2, km), __STAT__)
    tau_s = 0.0
      tau = 0.0

    allocate(q_4d(i1:i2, j1:j2, km, nbins), __STAT__)

    do n = 1, nbins
       call ESMF_StateGet (state, trim(aerosol_names(n)), field=fld, __RC__)
       call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

       do k = 1, km
          do j = j1, j2
             do i = i1, i2
                x = ((ple(i,j,k) - ple(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                q_4d(i,j,k,n) = x * q(i,j,k)
             end do
          end do
       end do
    end do

    call ESMF_AttributeGet(state, name='mieTable_pointer', itemCount=n, __RC__)
    allocate (opaque_self(n), __STAT__)
    call ESMF_AttributeGet(state, name='mieTable_pointer', valueList=opaque_self, __RC__)

    address = transfer(opaque_self, address)
    call c_f_pointer(address, self)

    do n = 1, nbins
       call self%diag_Mie%Query(wavelength, n, q_4d(:,:,:,n), rh, tau=tau, __RC__)
       tau_s = tau_s + tau
    end do

    call ESMF_AttributeGet(state, name='monochromatic_extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
        var = sum(tau_s, dim=3)
    end if

    deallocate(q_4d, __STAT__)

    RETURN_(ESMF_SUCCESS)

  end subroutine monochromatic_aerosol_optics


end module NI2G_GridCompMod
