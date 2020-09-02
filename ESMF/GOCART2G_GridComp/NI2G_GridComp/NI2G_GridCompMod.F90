#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: NI2G_GridCompMod - GOCART Nitrate gridded component 

! !INTERFACE:
module NI2G_GridCompMod

! !USES:
   use ESMF
   use MAPL
   use Chem_MieTableMod2G
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process       ! GOCART2G process library
   use GA_GridCompMod

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   real, parameter :: fMassHNO3 = 63., fMassNO3 = 62., fMassAir = 29.
   integer, parameter :: nNH3 = 1
   integer, parameter :: nNH4a = 2
   integer, parameter :: nNO3an1 = 3
   integer, parameter :: nNO3an2 = 4
   integer, parameter :: nNO3an3 = 5

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices

real, parameter ::  chemgrav   = 9.80616
real, parameter ::  cpd    = 1004.16

! !DESCRIPTION: This module implements GOCART's Nitrate (NI) Gridded Component.

! !REVISION HISTORY:
! 01July2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================

!  !Nitrate state
   type, extends(GA_GridComp) :: NI2G_GridComp
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
    type (wrap_)                                :: wrap
    type (NI2G_GridComp), pointer               :: self

    real                                        :: DEFVAL
    logical                                     :: data_driven=.true.

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

if(mapl_am_i_root()) print*,trim(comp_name),' SetServices BEGIN'

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file 
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'NI2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
       if (mapl_am_i_root()) print*,'NI2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! Loading NI2G_GridComp_NI.rc instead'
       call ESMF_ConfigLoadFile (cfg, 'NI2G_GridComp_NI.rc', __RC__)
    end if

    ! process generic config items
    call self%GA_GridComp%load_from_config( cfg, __RC__)

!   Is NI data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, __RC__)
    if (data_driven /= .true.) then
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

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
!    call MAPL_AddExportSpec(GC,                                                  &
!       SHORT_NAME = trim(COMP_NAME)//'_AERO_ACI',                                &
!       LONG_NAME  = 'aerosol_cloud_interaction_aerosols_from_'//trim(COMP_NAME),  &
!       UNITS      = 'kg kg-1',                                                   &
!       DIMS       = MAPL_DimsHorzVert,                                           &
!       VLOCATION  = MAPL_VLocationCenter,                                        &
!       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   DEVELOPMENT NOTE - Change to StateItem in future
!   ---------------------------------------------------------------
!    call MAPL_AddExportSpec(GC,                                   &
!       SHORT_NAME = trim(COMP_NAME)//'_AERO_DP',                  &
!       LONG_NAME  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
!       UNITS      = 'kg m-2 s-1',                                 &
!       DIMS       = MAPL_DimsHorzOnly,                            &
!       DATATYPE   = MAPL_BundleItem, __RC__)


!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'NI2G_GridComp', wrap, STATUS )
    VERIFY_(STATUS)

!   Set generic services
!   ----------------------------------
    call MAPL_GenericSetServices (GC, __RC__)


if(mapl_am_i_root()) print*,trim(comp_name),' SetServices END'

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
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero, aero_aci
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
    integer                              :: NUM_BANDS
    character (len=ESMF_MAXSTR), allocatable    :: aerosol_names(:)
    real, pointer, dimension(:,:,:)      :: ple
    real, pointer, dimension(:,:)        :: area

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Alarm)        :: alarm_HNO3
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss

    real, dimension(4)   :: Vect_Hcts
!    real, allocatable, dimension(:) :: rmedDU, rmedSS, fnumDU, fnumSS
    integer :: itemCount

    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)

if(mapl_am_i_root()) print*,trim(comp_name),' Init BEGIN'

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
    call ESMF_ConfigLoadFile (cfg, 'NI2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'NI2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                    loading NI2G_GridComp_NI.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'NI2G_GridComp_NI.rc', __RC__)
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
    call MAPL_GetResource (MAPL, NUM_BANDS, 'NUM_BANDS:', __RC__)

!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, self%rad_MieTable(instance)%optics_file, &
                                  label="aerosol_radBands_optics_file:", __RC__ )

    allocate (self%rad_MieTable(instance)%channels(NUM_BANDS), __STAT__ )

    call ESMF_ConfigGetAttribute (cfg, self%rad_MieTable(instance)%channels, label= "BANDS:", &
                                 count=self%rad_MieTable(instance)%nch, rc=status)

    if (rc /= 0) then
       do i = 1, NUM_BANDS
          self%rad_MieTable(instance)%channels(i) = i
       end do
    end if

    allocate (self%rad_MieTable(instance)%mie_aerosol, __STAT__)
    self%rad_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%rad_MieTable(instance)%optics_file, rc)
    call Chem_MieTableRead (self%rad_MieTable(instance)%mie_aerosol, NUM_BANDS, self%rad_MieTable(instance)%channels, rc)

!   Create Diagnostics Mie Table
!   -----------------------------
!   Get file names for the optical tables
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%optics_file, &
                                  label="aerosol_monochromatic_optics_file:", __RC__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nch, label="n_channels:", __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%nmom, label="n_moments:", default=0,  __RC__)
    allocate (self%diag_MieTable(instance)%channels(self%diag_MieTable(instance)%nch), __STAT__ )
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%channels, &
                                  label="aerosol_monochromatic_optics_wavelength:", __RC__)
    allocate (self%diag_MieTable(instance)%mie_aerosol, __STAT__)
    self%diag_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%diag_MieTable(instance)%optics_file, __RC__ )
    call Chem_MieTableRead (self%diag_MieTable(instance)%mie_aerosol, self%diag_MieTable(instance)%nch, &
                            self%diag_MieTable(instance)%channels, rc, nmom=self%diag_MieTable(instance)%nmom)

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

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',             value=0,     __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet(aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

!    call ESMF_AttributeSet(aero, name='internal_varaible_name', value='SS', __RC__)
    allocate(aerosol_names(3), __STAT__)
    aerosol_names(1) = 'NO3an1'
    aerosol_names(2) = 'NO3an2'
    aerosol_names(3) = 'NO3an3'
    call ESMF_AttributeSet(aero, name='internal_varaible_name', valueList=aerosol_names, &
                           itemCount=size(aerosol_names), __RC__)

    call ESMF_MethodAdd(AERO, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

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

!if(mapl_am_i_root()) print*,trim(comp_name),'2G Run1 BEGIN'

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, __RC__)

#include "NI2G_GetPointer___.h"

!if(mapl_am_i_root()) print*,'NI2G Run1 BEGIN sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G Run1 BEGIN sum(NH4a) = ',sum(NH4a)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'NI2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!if(mapl_am_i_root()) print*,'NI2G Run1 BEGIN sum(self%xhno3) = ',sum(self%xhno3)

!   NH3 Emissions
!   -------------
    if (associated(NH3EM)) then
       NH3EM = 0.
       if (associated(EMI_NH3_BB)) NH3EM = NH3EM + EMI_NH3_BB
       if (associated(EMI_NH3_AG)) NH3EM = NH3EM + EMI_NH3_AG
       if (associated(EMI_NH3_EN)) NH3EM = NH3EM + EMI_NH3_EN
       if (associated(EMI_NH3_RE)) NH3EM = NH3EM + EMI_NH3_RE
       if (associated(EMI_NH3_TR)) NH3EM = NH3EM + EMI_NH3_TR
       if (associated(EMI_NH3_IN)) NH3EM = NH3EM + EMI_NH3_IN
       if (associated(EMI_NH3_OC)) NH3EM = NH3EM + EMI_NH3_OC
    end if

    if (associated(EMI_NH3_BB)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_BB
    if (associated(EMI_NH3_AG)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_AG
    if (associated(EMI_NH3_EN)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_EN
    if (associated(EMI_NH3_IN)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_IN
    if (associated(EMI_NH3_RE)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_RE
    if (associated(EMI_NH3_TR)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_TR
    if (associated(EMI_NH3_OC)) &
       NH3(:,:,self%km) = NH3(:,:,self%km)+self%cdt*chemgrav/delp(:,:,self%km)*EMI_NH3_OC


!if(mapl_am_i_root()) print*,'NI2G sum(DU)',sum(DU)

!if(mapl_am_i_root()) print*,'NI2G Run1 END sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G Run1 END sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G Run1 END sum(self%xhno3) = ',sum(self%xhno3)
!if(mapl_am_i_root()) print*,trim(comp_name),'2G Run1 END'


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
    real, pointer, dimension(:,:)     :: fluxout
    real, pointer, dimension(:,:,:)   :: fluxoutWT
    real, allocatable, dimension(:,:,:,:) :: aerosol

!    character(len=15)                 :: ind

    type (ESMF_ALARM)               :: alarm
    logical                         :: alarm_is_ringing

    integer :: rhFlag

integer :: i,j
!real :: rmedDU(5), rmedSS(5), fnumDU(5), fnumSS(5)

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

!if(mapl_am_i_root()) print*,trim(comp_name),'2G Run2 BEGIN'

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
!if(mapl_am_i_root()) print*,'NI2G TEST1'
    end if

!   Recycle HNO3 every 3 hours
    if (alarm_is_ringing) then
       self%xhno3 = NITRATE_HNO3
       call ESMF_AlarmRingerOff(alarm, __RC__)
!if(mapl_am_i_root()) print*,'NI2G recycle alarm TRUE'
!if(mapl_am_i_root()) print*,'NI recycle alarm sum(self%xhno3)',sum(self%xhno3)
    end if

!if(mapl_am_i_root()) print*,'NI2G Run2 BEGIN sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G Run2 BEGIN sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G Run2 BEGIN sum(hno3) = ',sum(nitrate_hno3)
!if(mapl_am_i_root()) print*,'NI2G Run2 BEGIN before sum(self%xhno3) = ',sum(self%xhno3)
!if(mapl_am_i_root()) print*,'NI2G sum(DU) = ',sum(DU)
!if(mapl_am_i_root()) print*,'NI2G sum(SS) = ',sum(SS)

! This could be in incorrect alarm. This alarm is currently for HNO3_RECYCLE_ALARM, but a
! new alarm might need to be created just for this GC, or all of gocart2g?
!RUN_ALARM: if (alarm_is_ringing) then

    if (associated(NIPNO3AQ)) NIPNO3AQ(:,:) = 0.
    if (associated(NIPNH4AQ)) NIPNH4AQ(:,:) = 0.
    if (associated(NIPNH3AQ)) NIPNH3AQ(:,:) = 0.

!if(mapl_am_i_root()) print*,'NI2G before thermo sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G before thermo sum(NO3an1) = ',sum(NO3an1)
!if(mapl_am_i_root()) print*,'NI2G before thermo sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G before sum(SO4) = ',sum(SO4)

    call NIthermo (self%km, self%klid, self%cdt, chemgrav, delp, airdens, t, rh2, fMassHNO3, fMassAir, &
                   SO4, NH3, NO3an1, NH4a, self%xhno3, NIPNO3AQ, NIPNH4AQ, NIPNH3AQ, rc)

!if(mapl_am_i_root()) print*,'NI2G after thermo sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G after thermo sum(NO3an1) = ',sum(NO3an1)
!if(mapl_am_i_root()) print*,'NI2G after thermo sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G after thermo sum(xhno3) = ',sum(self%xhno3)
!if(mapl_am_i_root()) print*,'NI2G sum(NIPNO3AQ) = ',sum(NIPNO3AQ)
!if(mapl_am_i_root()) print*,'NI2G sum(NIPNH4AQ) = ',sum(NIPNH4AQ)
!if(mapl_am_i_root()) print*,'NI2G sum(NIPNH3AQ) = ',sum(NIPNH3AQ)

!    call NIheterogenousChem (NIHT, self%xhno3, MAPL_AVOGAD, MAPL_AIRMW, MAPL_PI, MAPL_RUNIV, &
    call NIheterogenousChem (NIHT, self%xhno3, MAPL_AVOGAD, MAPL_AIRMW, MAPL_PI, MAPL_RUNIV/1000., &
                             airdens, t, rh2, delp, DU, SS, self%rmedDU*1.e-6, self%rmedSS*1.e-6, &
                             self%fnumDU, self%fnumSS, 5, 5, self%km, self%klid, self%cdt, chemgrav, fMassHNO3, &
                             fMassNO3, fmassair, NO3an1, NO3an2, NO3an3, HNO3CONC, HNO3SMASS, &
                             HNO3CMASS, rc)

!if(mapl_am_i_root()) print*,'NI2G sum(NIHT(:,:,1)) = ',sum(NIHT(:,:,1))
!if(mapl_am_i_root()) print*,'NI2G sum(NIHT(:,:,2)) = ',sum(NIHT(:,:,2))
!if(mapl_am_i_root()) print*,'NI2G sum(NIHT(:,:,3)) = ',sum(NIHT(:,:,3))
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an1) = ',sum(NO3an1)
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an2) = ',sum(NO3an2)
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an3) = ',sum(NO3an3)
!if(mapl_am_i_root()) print*,'NI2G sum(HNO3CONC) = ',sum(HNO3CONC)
!if(mapl_am_i_root()) print*,'NI2G sum(HNO3SMASS) = ',sum(HNO3SMASS)
!if(mapl_am_i_root()) print*,'NI2G sum(HNO3CMASS) = ',sum(HNO3CMASS)

!if(mapl_am_i_root()) print*,'NI2G after hetchem sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G after hetchem sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G after hetchem sum(xhno3) = ',sum(self%xhno3)

!if(mapl_am_i_root()) print*,'NI2G NH4a array = ',NH4a

!   NI Settling
!   -----------
!   Because different bins having different swelling coefficients I need to
!   handle the call to settling differently.

!   Ammonium - settles like ammonium sulfate (rhflag = 3)
    rhflag = 3
    call Chem_SettlingSimpleOrig (self%km, self%klid, rhflag, chemgrav, self%cdt, &
                                  1.e-6*self%radius(nNH4a), self%rhop(nNH4a), &
                                  NH4a, t, airdens, rh2, delp, zle, NH4SD, rc)

!if(mapl_am_i_root()) print*,'NI2G sum(NH4SD) = ',sum(NH4SD)
!if(mapl_am_i_root()) print*,'NI2G sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G NH4a array = ',NH4a
!if(mapl_am_i_root()) print*,'NI2G NH4SD array = ',NH4SD


    allocate(fluxout, mold=lwi, __STAT__)
!  Nitrate bin 1 - settles like ammonium sulfate (rhflag = 3)
    rhflag = 3
    fluxout = 0.
    call Chem_SettlingSimpleOrig (self%km, self%klid, rhFlag, chemgrav, self%cdt, &
                                  1.e-6*self%radius(nNO3an1), self%rhop(nNO3an1), &
                                  NO3an1, t, airdens, rh2, delp, zle, fluxout, rc)
    if (associated(NISD)) NISD(:,:,1) = fluxout
!if(mapl_am_i_root()) print*,'NI2G sum(NISD(:,:,1)) = ',sum(NISD(:,:,1))
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an1) = ',sum(NO3an1)

!  Nitrate bin 2 - settles like sea salt (rhflag = 2)
    rhflag = 2
    fluxout = 0.
    call Chem_SettlingSimpleOrig (self%km, self%klid, rhFlag, chemgrav, self%cdt, &
                                  1.e-6*self%radius(nNO3an2), self%rhop(nNO3an2), &
                                  NO3an2, t, airdens, rh2, delp, zle, fluxout, rc)
    if (associated(NISD)) NISD(:,:,2) = fluxout
!if(mapl_am_i_root()) print*,'NI2G sum(NISD(:,:,2)) = ',sum(NISD(:,:,2))
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an2) = ',sum(NO3an2)

!  Nitrate bin 1 - settles like dust (rhflag = 0)
    rhflag = 0
    fluxout = 0.
    call Chem_SettlingSimpleOrig (self%km, self%klid, rhFlag, chemgrav, self%cdt, &
                                  1.e-6*self%radius(nNO3an3), self%rhop(nNO3an3), &
                                  NO3an3, t, airdens, rh2, delp, zle, fluxout, rc)
    if (associated(NISD)) NISD(:,:,3) = fluxout
!if(mapl_am_i_root()) print*,'NI2G sum(NISD(:,:,3)) = ',sum(NISD(:,:,3))
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an3) = ',sum(NO3an3)


!if(mapl_am_i_root()) print*,'NI2G after chemset sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G after chemset sum(NH4a) = ',sum(NH4a)

!  NI Deposition
!  -----------
    drydepositionfrequency = 0.
    call DryDeposition(self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
                       MAPL_KARMAN, cpd, chemGRAV, z0h, drydepositionfrequency, __RC__ )
!if(mapl_am_i_root()) print"(g25.17)",'NI2G drydep = ',drydepositionfrequency
!if(mapl_am_i_root()) print*,'NI2G NH3 array = ',NH3
!if(mapl_am_i_root()) print*,'NI2G lwi array = ',lwi

!  NH3
   dqa = 0.
!   where (abs(lwi - OCEAN) < 0.5)
!       dqa = max(0.0, NH3(:,:,self%km)*(1.-exp(-10.0*drydepositionfrequency*self%cdt)))
!   elsewhere
!       dqa = max(0.0, NH3(:,:,self%km)*(1.-exp( -3.0*drydepositionfrequency*self%cdt)))
!   end where

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
   if( associated(NH3DP) ) NH3DP = dqa*delp(:,:,self%km)/chemgrav/self%cdt
!if(mapl_am_i_root()) print*,'NI2G sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G sum(NH3DP) = ',sum(NH3DP)
!if(mapl_am_i_root()) print*,'NI2G dqa array = ',dqa
!if(mapl_am_i_root()) print"(g25.17)",'NI2G NH3 array = ',NH3


!  NH4a
   dqa = 0.
   dqa = max(0.0, NH4a(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NH4a(:,:,self%km) = NH4a(:,:,self%km) - dqa
   if( associated(NH4DP) ) NH4DP = dqa*delp(:,:,self%km)/chemgrav/self%cdt
!if(mapl_am_i_root()) print*,'NI2G sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G sum(NH4DP) = ',sum(NH4DP)

!  NO3anx
   dqa = 0.
   dqa = max(0.0, NO3an1(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NO3an1(:,:,self%km) = NO3an1(:,:,self%km) - dqa
   if( associated(NIDP) ) NIDP(:,:,1) = dqa*delp(:,:,self%km)/chemgrav/self%cdt
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an1) = ',sum(NO3an1)
!if(mapl_am_i_root()) print*,'NI2G sum(NIDP(:,:,1)) = ',sum(NIDP(:,:,1))

   dqa = 0.
   dqa = max(0.0, NO3an2(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NO3an2(:,:,self%km) = NO3an2(:,:,self%km) - dqa
   if( associated(NIDP) ) NIDP(:,:,2) = dqa*delp(:,:,self%km)/chemgrav/self%cdt
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an2) = ',sum(NO3an2)
!if(mapl_am_i_root()) print*,'NI2G sum(NIDP(:,:,2)) = ',sum(NIDP(:,:,2))

   dqa = 0.
   dqa = max(0.0, NO3an3(:,:,self%km)*(1.-exp(-drydepositionfrequency*self%cdt)))
   NO3an3(:,:,self%km) = NO3an3(:,:,self%km) - dqa
   if( associated(NIDP) ) NIDP(:,:,3) = dqa*delp(:,:,self%km)/chemgrav/self%cdt
!if(mapl_am_i_root()) print*,'NI2G sum(NO3an3) = ',sum(NO3an3)
!if(mapl_am_i_root()) print*,'NI2G sum(NIDP(:,:,3)) = ',sum(NIDP(:,:,3))

!  NI Large-scale Wet Removal
!  --------------------------
   allocate(fluxoutWT(ubound(t,1), ubound(t,2), 1), __STAT__)
   fluxoutWT = 0.
!  NH3
   KIN = .false.
   fwet = 1.
   call WetRemovalGOCART2G (self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'NH3', &
                            KIN, chemGRAV, fwet, NH3, ple, t, airdens, &
                            pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, fluxoutWT, rc)
   if (associated(NH3WT)) NH3WT = fluxoutWT(:,:,1)
!if(mapl_am_i_root()) print*,'NI2G sum(NH3WT) = ',sum(NH3WT)
!if(mapl_am_i_root()) print*,'NI2G sum(NH3) = ',sum(NH3)

!  NH4a
   fluxoutWT = 0.
   KIN = .true.
   fwet = 1.
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'NH4a', &
                           KIN, chemGRAV, fwet, NH4a, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, fluxoutWT, rc)
   if (associated(NH4WT)) NH4WT = fluxoutWT(:,:,1)
!if(mapl_am_i_root()) print*,'NI2G sum(NH4WT) = ',sum(NH4WT)
!if(mapl_am_i_root()) print*,'NI2G sum(NH4) = ',sum(NH4a)

   KIN = .true.
   fwet = 1.
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 1, self%cdt, 'nitrate', &
                           KIN, chemGRAV, fwet, NO3an1, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, rc)
!if(mapl_am_i_root()) print*,'NI2G sum(NIWT(:,:,1)) = ',sum(NIWT(:,:,1))
   KIN = .true.
   fwet = 1.
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 2, self%cdt, 'nitrate', &
                           KIN, chemGRAV, fwet, NO3an2, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, rc)
!if(mapl_am_i_root()) print*,'NI2G sum(NIWT(:,:,2)) = ',sum(NIWT(:,:,2))

   KIN = .true.
   fwet = 1.
   call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, 3, self%cdt, 'nitrate', &
                           KIN, chemGRAV, fwet, NO3an3, ple, t, airdens, &
                           pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, NIWT, rc)
!if(mapl_am_i_root()) print*,'NI2G sum(NIWT(:,:,3)) = ',sum(NIWT(:,:,3))

!  Compute desired output diagnostics
!  ----------------------------------
   allocate(aerosol(ubound(NH4a,1), ubound(NH4a,2), ubound(NH4a,3), 3), __STAT__)
   aerosol(:,:,:,:) = 0.0
   aerosol(:,:,:,1) = NH4a
   call Aero_Compute_Diags (mie_table=self%diag_MieTable(self%instance), km=self%km, klid=self%klid, nbegin=1, &
                            nbins=1, channels=self%diag_MieTable(self%instance)%channels, &
                            aerosol=aerosol, grav=chemgrav, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, sfcmass=NH4SMASS, colmass=NH4CMASS, mass=NH4MASS, conc=NH4CONC, __RC__)
!if(mapl_am_i_root()) print*,'NI2G sum(NH4SMASS) = ',sum(NH4SMASS) 
!if(mapl_am_i_root()) print*,'NI2G sum(NH4CMASS) = ',sum(NH4CMASS)
!if(mapl_am_i_root()) print*,'NI2G sum(NH4MASS) = ',sum(NH4MASS)
!if(mapl_am_i_root()) print*,'NI2G sum(NH4CONC) = ',sum(NH4CONC)

   aerosol(:,:,:,1) = NH3
   call Aero_Compute_Diags (mie_table=self%diag_MieTable(self%instance), km=self%km, klid=self%klid, nbegin=1, &
                            nbins=1, channels=self%diag_MieTable(self%instance)%channels, &
                            aerosol=aerosol, grav=chemgrav, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, sfcmass=NH3SMASS, colmass=NH3CMASS, mass=NH3MASS, conc=NH3CONC, __RC__)
!if(mapl_am_i_root()) print*,'NI2G sum(NH3SMASS) = ',sum(NH3SMASS)
!if(mapl_am_i_root()) print*,'NI2G sum(NH3CMASS) = ',sum(NH3CMASS)
!if(mapl_am_i_root()) print*,'NI2G sum(NH3MASS) = ',sum(NH3MASS)
!if(mapl_am_i_root()) print*,'NI2G sum(NH3CONC) = ',sum(NH3CONC)

   aerosol(:,:,:,1) = NO3an1
   call Aero_Compute_Diags (mie_table=self%diag_MieTable(self%instance), km=self%km, klid=self%klid, nbegin=1, &
                            nbins=1, channels=self%diag_MieTable(self%instance)%channels, &
                            aerosol=aerosol, grav=chemgrav, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, sfcmass=NISMASS25, colmass=NICMASS25, mass=NIMASS25, conc=NICONC25, __RC__)
!if(mapl_am_i_root()) print*,'NI2G sum(NISMASS25) = ',sum(NISMASS25)
!if(mapl_am_i_root()) print*,'NI2G sum(NICMASS25) = ',sum(NICMASS25)
!if(mapl_am_i_root()) print*,'NI2G sum(NIMASS25) = ',sum(NIMASS25)
!if(mapl_am_i_root()) print*,'NI2G sum(NICONC25) = ',sum(NICONC25)

!#if 0
   aerosol(:,:,:,1) = NO3an1
   aerosol(:,:,:,2) = NO3an2
   aerosol(:,:,:,3) = NO3an3
   call Aero_Compute_Diags (mie_table=self%diag_MieTable(self%instance), km=self%km, klid=self%klid, nbegin=1, &
                            nbins=3, channels=self%diag_MieTable(self%instance)%channels, &
                            aerosol=aerosol, grav=chemgrav, tmpu=t, rhoa=airdens, rh=rh2, u=u, v=v, &
                            delp=delp, sfcmass=NISMASS, colmass=NICMASS, mass=NIMASS, conc=NICONC, &
                            exttau=NIEXTTAU, scatau=NISCATAU, exttau25=NIEXTT25, scatau25=NISCAT25, &
                            fluxu=NIFLUXU, fluxv=NIFLUXV, extcoef=NIEXTCOEF, scacoef=NISCACOEF, &
                            exttaufm=NIEXTTFM, scataufm=NISCATFM, angstrom=NIANGSTR, __RC__ )
!#endif
!if(mapl_am_i_root()) print*,'NI2G sum(NIEXTTAU) = ',sum(NIEXTTAU)
!if(mapl_am_i_root()) print*,'NI2G sum(NISCATAU) = ',sum(NISCATAU)
!if(mapl_am_i_root()) print*,'NI2G sum(NIMASS) = ',sum(NIMASS)
!if(mapl_am_i_root()) print*,'NI2G sum(NIFLUXU) = ',sum(NIFLUXU)

!if(mapl_am_i_root()) print*,'NI2G Run2 END sum(NIANGSTR) = ',sum(NIANGSTR)
!if(mapl_am_i_root()) print*,'NI2G Run2 END sum(NH3) = ',sum(NH3)
!if(mapl_am_i_root()) print*,'NI2G Run2 END sum(NH4a) = ',sum(NH4a)
!if(mapl_am_i_root()) print*,'NI2G Run2 END sum(self%xhno3) = ',sum(self%xhno3)

!if(mapl_am_i_root()) print*,'NI2G Run2 END array NH3 = ',NH3
!if(mapl_am_i_root()) print*,'NI2G Run2 END array NH4a = ',NH4a
!if(mapl_am_i_root()) print*,'NI2G Run2 END sum(NO3an1) = ',sum(NO3an1)
!if(mapl_am_i_root()) print*,'NI2G Run2 END sum(NO3an2) = ',sum(NO3an2)
!if(mapl_am_i_root()) print*,'NI2G Run2 END sum(NO3an3) = ',sum(NO3an3)


if(mapl_am_i_root()) print*,trim(comp_name),'2G Run2 END'

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

if(mapl_am_i_root())print*,'NI2G Run_data END'

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
    integer                                          :: band, offset
    integer, parameter                               :: n_bands = 1

    integer :: i, j, k

    __Iam__('NI2G::aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

!   Get aerosol names
!   -----------------
    call ESMF_AttributeGet (state, name='internal_varaible_name', itemCount=nbins, __RC__)
    allocate (aerosol_names(nbins), __STAT__)
    call ESMF_AttributeGet (state, name='internal_varaible_name', valueList=aerosol_names, __RC__)

!   Radiation band
!   --------------
    band = 0
    call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)
    offset = band - n_bands

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

    call mie_ (self%rad_MieTable(instance), nbins, n_bands, offset, q_4d, rh, ext_s, ssa_s, asy_s, __RC__)

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
    subroutine mie_(mie_table, nbins, nb, offset, q, rh, bext_s, bssa_s, basym_s, rc)

    implicit none

    type(Chem_Mie),                intent(inout) :: mie_table        ! mie table
    integer,                       intent(in   ) :: nbins            ! number of bins
    integer,                       intent(in )   :: nb               ! number of bands
    integer,                       intent(in )   :: offset           ! bands offset 
    real,                          intent(in )   :: q(:,:,:,:)       ! aerosol mass mixing ratio, kg kg-1
    real,                          intent(in )   :: rh(:,:,:)        ! relative humidity
    real(kind=8), intent(  out) :: bext_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=8), intent(  out) :: bssa_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=8), intent(  out) :: basym_s(size(ext_s,1),size(ext_s,2),size(ext_s,3))
    integer,           intent(out)  :: rc
    ! local
    integer                           :: l
    real                              :: bext (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! extinction
    real                              :: bssa (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! SSA
    real                              :: gasym(size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! asymmetry parameter

    __Iam__('NI2G::aerosol_optics::mie_')

     bext_s  = 0.0d0
     bssa_s  = 0.0d0
     basym_s = 0.0d0

    do l = 1, nbins
       call Chem_MieQuery(mie_table, l, real(offset+1.), q(:,:,:,l), rh, bext, gasym=gasym, ssa=bssa)

       bext_s  = bext_s  +             bext     ! extinction
       bssa_s  = bssa_s  +       (bssa*bext)    ! scattering extinction
       basym_s = basym_s + gasym*(bssa*bext)    ! asymetry parameter multiplied by scatering extiction 

    end do

    RETURN_(ESMF_SUCCESS)

    end subroutine mie_

  end subroutine aerosol_optics




end module NI2G_GridCompMod
