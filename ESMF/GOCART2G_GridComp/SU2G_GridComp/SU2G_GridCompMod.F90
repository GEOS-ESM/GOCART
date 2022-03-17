#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: SU2G_GridCompMod - GOCART refactoring of the SU gridded component 

! !INTERFACE:
module SU2G_GridCompMod

! !USES:
   use ESMF
   use MAPL
   use GOCART2G_MieMod 
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process       ! GOCART2G process library
   use GA_EnvironmentMod
   use MAPL_StringTemplate, only: StrTemplate

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

   integer, parameter     :: DP=kind(1.0d0)

!  gram molecular weights of species
   real, parameter :: fMassSulfur = 32., fMassSO2 = 64., fMassSO4 = 96., &
                      fMassDMS = 62., fMassMSA = 96.

   real, parameter ::  cpd    = 1004.16

!  relative position of sulfate tracers
   integer, parameter :: nDMS = 1, &
                         nSO2 = 2, &
                         nSO4 = 3, &
                         nMSA = 4
real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices


! !DESCRIPTION: This module implements GOCART's Sulfer (SU) Gridded Component.

! !REVISION HISTORY:
! 08July2020  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring.

!EOP
!===========================================================================
!  !Sulfer state
   type, extends(GA_Environment) :: SU2G_GridComp
      integer :: myDOW = -1     ! my Day of the week: Sun=1, Mon=2,...,Sat=7
      logical :: diurnal_bb     ! diurnal biomass burning
      integer :: nymd_last = -1 ! Previous nymd. Updated daily
      integer :: nymd_oxidants = -1 ! Update the oxidant files?
      real    :: eAircraftFuel  ! Aircraft emission factor: go from kg fuel to kg SO2
      real    :: aviation_layers(4)  ! heights of the LTO, CDS and CRS layers
      real    :: fSO4anth  ! Fraction of anthropogenic emissions that are SO4
      logical :: recycle_H2O2 = .false.
      logical :: firstRun = .true.
      real, allocatable  :: sigma(:) ! Sigma of lognormal number distribution
      real, pointer :: h2o2_init(:,:,:)

!     Special handling for volcanic emissions
      character(len=255) :: volcano_srcfilen
      integer :: nVolc = 0
      real, allocatable, dimension(:)  :: vLat, &
                                          vLon, &
                                          vSO2, &
                                          vElev, &
                                          vCloud
      integer, allocatable, dimension(:) :: vStart, &
                                            vEnd
!     !Workspae for point emissions
      logical                :: doing_point_emissions = .false.
      character(len=255)     :: point_emissions_srcfilen   ! filename for pointwise emissions
      integer                         :: nPts = -1
      integer, allocatable, dimension(:)  :: pstart, pend
      real, allocatable, dimension(:)     :: pLat, &
                                             pLon, &
                                             pBase, &
                                             pTop, &
                                             pEmis
   end type SU2G_GridComp

   type wrap_
      type (SU2G_GridComp), pointer     :: PTR => null()
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
!   08July2020   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

!EOP
!============================================================================
!

!   !Locals
    character (len=ESMF_MAXSTR)                 :: COMP_NAME
    type (ESMF_Config)                          :: cfg
    type (ESMF_Config)                          :: universal_cfg
    type (wrap_)                                :: wrap
    type (SU2G_GridComp), pointer               :: self

    character (len=ESMF_MAXSTR)                 :: field_name

    integer                                     :: i
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
    call ESMF_ConfigLoadFile (cfg, 'SU2G_instance_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
       if (mapl_am_i_root()) print*,'SU2G_instance_'//trim(COMP_NAME)//'.rc does not exist! loading SU2G_instance_SU.rc instead'
       call ESMF_ConfigLoadFile (cfg, 'SU2G_instance_SU.rc', __RC__)
    end if

!   process generic config items
    call self%GA_Environment%load_from_config( cfg, universal_cfg, __RC__)

    allocate(self%sigma(self%nbins), __STAT__)

!   process SU-specific items
    call ESMF_ConfigGetAttribute(cfg, self%volcano_srcfilen, label='volcano_srcfilen:', __RC__)
    call ESMF_ConfigGetAttribute(cfg, self%eAircraftFuel, label='aircraft_fuel_emission_factor:', __RC__)
    call ESMF_ConfigGetAttribute(cfg, self%fSO4anth, label='so4_anthropogenic_fraction:', __RC__)
    call ESMF_ConfigGetAttribute(cfg, self%sigma, label='sigma:', __RC__)
    call ESMF_ConfigFindLabel (cfg, 'aviation_vertical_layers:', __RC__)
    do i=1,size(self%aviation_layers)
       call ESMF_ConfigGetAttribute (cfg, self%aviation_layers(i), __RC__)
    end do

    call ESMF_ConfigGetAttribute (cfg, self%point_emissions_srcfilen, &
                                  label='point_emissions_srcfilen:', default='/dev/null', __RC__)
    if ( (index(self%point_emissions_srcfilen,'/dev/null')>0) ) then
       self%doing_point_emissions = .false. ! disable it if no file specified
    else
       self%doing_point_emissions = .true.  ! we are good to go
    end if

!   Is SU data driven?
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
         short_name='SO4', &
         long_name='Sulphate aerosol', &
         units='kg kg-1', &
         dims=MAPL_DimsHorzVert, &
         vlocation=MAPL_VlocationCenter, &
         restart=MAPL_RestartOptional, &
!         friendlyto='DYNAMICS:TURBULENCE:MOIST', &
         add2export=.true., __RC__)

       call MAPL_AddImportSpec(gc,&
         short_name='climSO4', &
         long_name='Sulphate aerosol', &
         units='kg kg-1', &
         dims=MAPL_DimsHorzVert, &
         vlocation=MAPL_VlocationCenter, &
         restart=MAPL_RestartOptional, __RC__) 

        do i = 1, self%nbins
            write(field_name, '(A, I0.3)') '', i
!           ! dry deposition
            call MAPL_AddImportSpec(GC,                                           &
              SHORT_NAME = 'climSUDP'//trim(field_name),                          &
              LONG_NAME  = 'Sulfate dry deposition (bin '//trim(field_name)//')',  &
              UNITS      = 'kg kg-1 s-1',                                             &
              DIMS       = MAPL_DimsHorzOnly,                                     &
              VLOCATION  = MAPL_VLocationCenter,                                  &
              RESTART    = MAPL_RestartSkip, __RC__)

!           ! wet deposition    
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSUWT'//trim(field_name),                         &
               LONG_NAME  = 'Sulfate wet deposition (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1 s-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!           ! gravitational settling
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSUSD'//trim(field_name),                         &
               LONG_NAME  = 'Sulfate settling (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1 s-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!        ! convective scavenging
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSUSV'//trim(field_name),                         &
               LONG_NAME  = 'Sulfate convective scavenging (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1 s-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if ! (data_driven)


    if (.not. data_driven) then
       call MAPL_AddImportSpec(GC,             &
          SHORT_NAME = 'SU_H2O2', &
          LONG_NAME  = 'source species'  ,     &
          UNITS      = '1',                    &
          DIMS       = MAPL_DimsHorzVert,      &
          VLOCATION  = MAPL_VLocationCenter,   &
          RESTART    = MAPL_RestartSkip, __RC__) 

       call MAPL_AddImportSpec(GC,           &
          SHORT_NAME = 'SU_OH', &
          LONG_NAME  = 'source species'  ,   &
          UNITS      = '1',                  &
          DIMS       = MAPL_DimsHorzVert,    &
          VLOCATION  = MAPL_VLocationCenter, &
          RESTART    = MAPL_RestartSkip, __RC__) 

       call MAPL_AddImportSpec(GC,            &
          SHORT_NAME = 'SU_NO3', &
          LONG_NAME  = 'source species'  ,    &
          UNITS      = '1',                   &
          DIMS       = MAPL_DimsHorzVert,     &
          VLOCATION  = MAPL_VLocationCenter,  &
          RESTART    = MAPL_RestartSkip, __RC__)
    end if

!   Import, Export, Internal states for computational instance 
!   ----------------------------------------------------------
    if (.not. data_driven) then
#include "SU2G_Export___.h"
#include "SU2G_Import___.h"
#include "SU2G_Internal___.h"
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

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   DEVELOPMENT NOTE - Change to StateItem in future
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                   &
       SHORT_NAME = trim(COMP_NAME)//'_AERO_DP',                  &
       LONG_NAME  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg m-2 s-1',                                 &
       DIMS       = MAPL_DimsHorzOnly,                            &
       DATATYPE   = MAPL_BundleItem, __RC__)


!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'SU2G_GridComp', wrap, STATUS )
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

! !DESCRIPTION: This initializes SU Grid Component.  

! !REVISION HISTORY: 
! 08July2019   E.Sherman  First attempt at refactoring

!EOP
!============================================================================
!   !Locals 
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),      pointer   :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg, universal_cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (SU2G_GridComp), pointer        :: self
    type (ESMF_Alarm)                    :: alarm_H2O2

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: prefix, diurnal_bb
    real, pointer, dimension(:,:)        :: lats
    real, pointer, dimension(:,:)        :: lons
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)
    real, pointer, dimension(:,:,:)      :: int_ptr
    logical                              :: data_driven
    integer                              :: NUM_BANDS
    logical                              :: bands_are_present
    real, pointer, dimension(:,:,:)      :: ple

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss

    real, dimension(4)   :: Vect_Hcts
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
    call ESMF_UserCompGetInternalState(GC, 'SU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, globalCellCountPerDim=dims, __RC__ )
    km = dims(3)
    self%km = km

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!   Check whether to de-activate diurnal biomass burning (default is *on*)
!   ----------------------------------------------------------------------
    call ESMF_ConfigGetAttribute(universal_cfg, diurnal_bb, label='DIURNAL_BIOMASS_BURNING:', &
                                 default='YES', __RC__)
    diurnal_bb = ESMF_UtilStringUpperCase(diurnal_bb, __RC__)
    if (trim(diurnal_bb) == 'YES') then
       self%diurnal_bb = .true.
    else
       self%diurnal_bb = .false.
    end if

!  Load resource file and get number of bins 
!  -------------------------------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'SU2G_instance_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'SU2G_instance_'//trim(COMP_NAME)//'.rc does not exist! &
                                    loading SU2G_instance_SU.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'SU2G_instance_SU.rc', __RC__)
    end if

!   Call Generic Initialize 
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( mapl, INTERNAL_ESMF_STATE = internal, &
                         LONS = LONS, &
                         LATS = LATS, __RC__ )

    allocate(self%h2o2_init(size(lats,1),size(lats,2),self%km), __STAT__)

!   Is SU data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set H2O2 recycle alarm
!   ----------------------
    if (.not. data_driven) then
        call ESMF_ClockGet(clock, calendar=calendar, currTime=currentTime, __RC__)
        call ESMF_TimeGet(currentTime, YY=year, MM=month, DD=day, H=hh, M=mm, S=ss, __RC__)
        call ESMF_TimeSet(ringTime, YY=year, MM=month, DD=day, H=0, M=0, S=0, __RC__)
        call ESMF_TimeIntervalSet(ringInterval, H=3, calendar=calendar, __RC__)

        do while (ringTime < currentTime)! DO WE NEED THIS?
            ringTime = currentTime + ringInterval
        end do

        alarm_H2O2 = ESMF_AlarmCreate(Clock        = clock,        &
                                      Name         = 'H2O2_RECYCLE_ALARM', &
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

!   Add attribute information
    if (.not. data_driven) then
       call ESMF_StateGet (internal, 'DMS', field, __RC__)
       call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), __RC__)
!       Vect_Hcts=(/0.56, 3500.0, 0.0, 0.0/)
       call get_HenrysLawCts('DMS',Vect_Hcts(1),Vect_Hcts(2),Vect_Hcts(3),Vect_Hcts(4),__RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, __RC__)

       call ESMF_StateGet (internal, 'SO2', field, __RC__)
       call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(2), __RC__)
!       Vect_Hcts=(/1.4, 2900.0, 1.3E-02, 2000.0/)
       call get_HenrysLawCts('SO2',Vect_Hcts(1),Vect_Hcts(2),Vect_Hcts(3),Vect_Hcts(4),__RC__)
       call ESMF_AttributeSet(field, 'SetofHenryLawCts', Vect_Hcts, __RC__)

       call ESMF_StateGet (internal, 'MSA', field, rc=status)
       if (status == 0) then
          call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(4), __RC__)
       end if
    end if

!   Fill AERO State with SO4
!   ----------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

    call ESMF_StateGet (internal, 'SO4', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(3), __RC__)
    fld = MAPL_FieldCreate (field, 'SO4', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    if (.not. data_driven) then
!      Set klid
       call MAPL_GetPointer(import, ple, 'PLE', __RC__)
       call findKlid (self%klid, self%plid, ple, __RC__)
!      Set internal SO4 values to 0 where above klid
       call MAPL_GetPointer (internal, int_ptr, 'SO4', __RC__)
       call setZeroKlid (self%km, self%klid, int_ptr)
    end if

    if (data_driven) then
       instance = instanceData

!      (SO4; only aerosol component; bin 003)
!      Dry deposition
       call append_to_bundle('SUDP003', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition (Convective scavenging)
       call append_to_bundle('SUSV003', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition
       call append_to_bundle('SUWT003', providerState, prefix, Bundle_DP, __RC__)

!      Gravitational Settling
       call append_to_bundle('SUSD003', providerState, prefix, Bundle_DP, __RC__)
    else
       instance = instanceComputational

!      Dry deposition
       call append_to_bundle('SUDP', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition (Convective scavenging)
       call append_to_bundle('SUSV', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition
       call append_to_bundle('SUWT', providerState, prefix, Bundle_DP, __RC__)

!      Gravitational Settling
       call append_to_bundle('SUSD', providerState, prefix, Bundle_DP, __RC__)
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

    ! Add variables to SU instance's aero state. This is used in aerosol optics calculations
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
                   label2='monochromatic_EXT', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol', label2='aerosolSum', grid=grid, typekind=MAPL_R4, __RC__)

    call ESMF_AttributeSet (aero, name='band_for_aerosol_optics', value=0, __RC__)
    call ESMF_AttributeSet (aero, name='wavelength_for_aerosol_optics', value=0., __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet (aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

    call ESMF_AttributeSet (aero, name='internal_variable_name', value='SO4', __RC__)

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

!   Is SU data driven?
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

! !DESCRIPTION:  Computes emissions/sources for Sea Salt

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (SU2G_GridComp), pointer     :: self
    type(ESMF_Time)                   :: time

    character(len=3) :: cdow
    integer          :: idow
    integer          :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
    real, pointer, dimension(:,:)        :: lats
    real, pointer, dimension(:,:)        :: lons
    real, dimension(:,:,:), allocatable  :: aircraft_fuel_src
    real, dimension(:,:), allocatable :: so2biomass_src, so2biomass_src_, so2anthro_l1_src, &
                                         so2anthro_l2_src, so2ship_src, so4ship_src, dmso_conc, &
                                         aviation_lto_src, aviation_cds_src, aviation_crs_src
    integer, dimension(:), allocatable  :: iPointVolc, jPointVolc, iPoint, jPoint
    real, dimension(:,:,:), allocatable :: emissions_point
    character (len=ESMF_MAXSTR)  :: fname ! file name for point source emissions
    logical :: fileExists

    real, pointer, dimension(:,:,:) :: dummyMSA => null() ! This is so the model can run without MSA enabled

#include "SU2G_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    Iam = trim(comp_name) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, &
                        LONS = LONS, &
                        LATS = LATS, __RC__ )

#include "SU2G_GetPointer___.h"

    call MAPL_GetPointer(internal, dummyMSA, 'MSA', rc=status)


!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'SU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

!   Reset tracer to zero at 0Z on specific day of week
!   --------------------------------------------------
    idow = Chem_UtilIdow(nymd)
    if ( (nhms==0) .and. (idow == self%myDOW) ) then
       cdow = Chem_UtilCdow(nymd)
       DMS = tiny(1.) ! avoid division by zero
       SO2 = tiny(1.) ! avoid division by zero
       SO4 = tiny(1.) ! avoid division by zero
       if ( associated(dummyMSA)) dummyMSA = tiny(1.) ! avoid division by zero
       if ( MAPL_AM_I_ROOT() ) then
          print *, '<> SU '//cdow//' tracer being set to zero on ', nymd, nhms
       end if
    end if

!   Implicit allocation with Fortran 2003
    so2anthro_l1_src = SU_ANTHROL1
    so2anthro_l2_src = SU_ANTHROL2
    so2ship_src = SU_SHIPSO2
    so4ship_src = SU_SHIPSO4

!   As a safety check, where value is undefined set to 0
    where(1.01*so2anthro_l1_src > MAPL_UNDEF)  so2anthro_l1_src = 0.
    where(1.01*so2anthro_l2_src > MAPL_UNDEF)  so2anthro_l2_src = 0.
    where(1.01*so2ship_src > MAPL_UNDEF)       so2ship_src = 0.
    where(1.01*so4ship_src > MAPL_UNDEF)       so4ship_src = 0.

    aircraft_fuel_src = SU_AIRCRAFT
    so2biomass_src = SU_BIOMASS
    dmso_conc = SU_DMSO
    aviation_lto_src = SU_AVIATION_LTO
    aviation_cds_src = SU_AVIATION_CDS
    aviation_crs_src = SU_AVIATION_CRS

!   As a safety check, where value is undefined set to 0
    where(1.01*so2biomass_src > MAPL_UNDEF)    so2biomass_src = 0.
    where(1.01*dmso_conc > MAPL_UNDEF)         dmso_conc = 0.
    where(1.01*aircraft_fuel_src > MAPL_UNDEF) aircraft_fuel_src = 0.
    where(1.01*aviation_lto_src > MAPL_UNDEF ) aviation_lto_src = 0.
    where(1.01*aviation_cds_src > MAPL_UNDEF ) aviation_cds_src = 0.
    where(1.01*aviation_crs_src > MAPL_UNDEF ) aviation_crs_src = 0.

!   Update emissions/production if necessary (daily)
!   -----------------------------------------------
    if(self%nymd_last /= nymd) then
       self%nymd_last = nymd

!      Get pointwise SO2 and altitude of volcanoes from a daily file data base
       if(index(self%volcano_srcfilen,'volcanic_') /= 0) then
          call StrTemplate(fname, self%volcano_srcfilen, xid='unknown', &
                            nymd=nymd, nhms=120000 )
          call ReadPointEmissions (nymd, fname, self%nVolc, self%vLat, self%vLon, &
                                   self%vElev, self%vCloud, self%vSO2, self%vStart, &
                                   self%vEnd, label='volcano', __RC__)
          self%vSO2 = self%vSO2 * fMassSO2 / fMassSulfur
!         Special possible case
          if(self%volcano_srcfilen(1:9) == '/dev/null') self%nVolc = 0
       end if
    end if

!   Apply volcanic emissions
!   ------------------------
    if (self%nVolc > 0) then
       if (associated(SO2EMVE)) SO2EMVE=0.0
       if (associated(SO2EMVN)) SO2EMVN=0.0
       allocate(iPointVolc(self%nVolc), jPointVolc(self%nVolc),  __STAT__)
       call MAPL_GetHorzIJIndex(self%nVolc, iPointVolc, jPointVolc, &
                                grid = grid,               &
                                lon  = self%vLon/real(MAPL_RADIANS_TO_DEGREES), &
                                lat  = self%vLat/real(MAPL_RADIANS_TO_DEGREES), &
                                rc   = status)
           if ( status /= 0 ) then
              if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
              VERIFY_(status)
           end if

       call SUvolcanicEmissions (self%nVolc, self%vStart, self%vEnd, self%vSO2, self%vElev, &
                                 self%vCloud, iPointVolc, jPointVolc, nhms, SO2EMVN, SO2EMVE, SO2, nSO2, SUEM, &
                                 self%km, self%cdt, MAPL_GRAV, zle, delp, area, self%vLat, self%vLon, __RC__)
    end if

!   Apply diurnal cycle if so desired
    if ( self%diurnal_bb ) then
       so2biomass_src_ = so2biomass_src
       call Chem_BiomassDiurnal (so2biomass_src, so2biomass_src_, &
                                 lons(:,:)*real(MAPL_RADIANS_TO_DEGREES), &
                                 lats(:,:)*real(MAPL_RADIANS_TO_DEGREES), &
                                 nhms, self%cdt)
    end if

!   Apply sulfate emissions
!   -----------------------
    call SulfateDistributeEmissions ( self%km, self%nbins, self%cdt, MAPL_GRAV, nymd, nhms, &
                                      fMassSO4, fMassSO2, self%fSO4anth, self%eAircraftFuel, &
                                      nSO2, nSO4, &
                                      so2anthro_l1_src, so2anthro_l2_src, &
                                      so2biomass_src, dmso_conc, &
                                      so2ship_src, so4ship_src, &
                                      aircraft_fuel_src, &
                                      SO2, SO4, &
                                      lwi, u10m, v10m, zle, zpbl, &
                                      t, airdens, delp, self%nVolc, &
                                      SUEM, SO4EMAN, SO2EMAN, SO2EMBB, &
                                      self%aviation_layers,   &
                                      aviation_lto_src, &
                                      aviation_cds_src, &
                                      aviation_crs_src, __RC__) 

    if (associated(dms)) then 
       call DMSemission (self%km, self%cdt, MAPL_GRAV, t, u10m, v10m, lwi, delp, &
                         fMassDMS, SU_DMSO, dms, SUEM, nDMS, __RC__)
    end if

!   Add source of OCS-produced SO2
!   ------------------------------
    SO2 = SO2 + pSO2_OCS*self%cdt

!   Read any pointwise emissions, if requested
!   ------------------------------------------
    if(self%doing_point_emissions) then
       call StrTemplate(fname, self%point_emissions_srcfilen, xid='unknown', &
                         nymd=nymd, nhms=120000 )
       inquire( file=fname, exist=fileExists)
       if (fileExists) then
          call ReadPointEmissions (nymd, fname, self%nPts, self%pLat, self%pLon, &
                                   self%pBase, self%pTop, self%pEmis, self%pStart, &
                                   self%pEnd, label='source', __RC__)
       else if (.not. fileExists) then
         if(mapl_am_i_root()) print*,'GOCART2G ',trim(comp_name),': ',trim(fname),' not found; proceeding.'
         self%nPts = -1 ! set this back to -1 so the "if (self%nPts > 0)" conditional is not exercised.
       end if
    endif

!   Get indices for point emissions
!   -------------------------------
    if (self%nPts > 0) then
       allocate(emissions_point, mold=delp,  __STAT__)
       emissions_point = 0.0
       allocate(iPoint(self%nPts), jPoint(self%nPts),  __STAT__)
       call MAPL_GetHorzIJIndex(self%nPts, iPoint, jPoint, &
                                grid = grid,               &
                                lon  = self%pLon/real(MAPL_RADIANS_TO_DEGREES), &
                                lat  = self%pLat/real(MAPL_RADIANS_TO_DEGREES), &
                                rc   = status)
          if ( status /= 0 ) then
             if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
             VERIFY_(status)
          end if

        call updatePointwiseEmissions (self%km, self%pBase, self%pTop, self%pEmis, self%nPts, &
                                       self%pStart, self%pEnd, zle, &
                                       area, iPoint, jPoint, nhms, emissions_point, __RC__)
   
        SO4 = SO4 + self%cdt * MAPL_GRAV / delp * emissions_point
     end if

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
    type (SU2G_GridComp), pointer     :: self
    type (ESMF_Time)                  :: time
    type (ESMF_Alarm)                 :: ALARM
    type(MAPL_VarSpec), pointer       :: InternalSpec(:)

    integer         :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
    integer                           :: n
    logical                           :: KIN
    real, pointer, dimension(:,:)     :: lats
    real, pointer, dimension(:,:)     :: lons
    character(len=ESMF_MAXSTR)        :: short_name
    real, pointer, dimension(:,:,:)   :: int_ptr

    real, dimension(:,:,:), allocatable :: xoh, xno3, xh2o2

    real, dimension(:,:), allocatable   :: drydepositionf
    real, pointer, dimension(:,:,:)     :: dummyMSA => null() ! this is so the model can run without MSA enabled
    logical :: alarm_is_ringing  

#include "SU2G_DeclarePointer___.h"

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
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=INTERNAL, &
                         INTERNALSPEC = InternalSpec, &
                         LONS = LONS, &
                         LATS = LATS, __RC__ )

#include "SU2G_GetPointer___.h"
    
    call MAPL_GetPointer(internal, dummyMSA, 'MSA', rc=status)
    
!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'SU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    call ESMF_ClockGetAlarm(clock, 'H2O2_RECYCLE_ALARM', alarm, __RC__)
    alarm_is_ringing = ESMF_AlarmIsRinging(alarm, __RC__)
!   recycle H2O2 every 3 hours
    if (alarm_is_ringing) then
       self%recycle_h2o2 = ESMF_AlarmIsRinging(alarm, __RC__)
       call ESMF_AlarmRingerOff(alarm, __RC__)
    end if

    allocate(xoh, mold=airdens, __STAT__)
    allocate(xno3, mold=airdens, __STAT__)
    allocate(xh2o2, mold=airdens, __STAT__)
    xoh = 0.0
    xno3 = 0.0

    if (self%firstRun) then
       xh2o2          = MAPL_UNDEF
       self%h2o2_init = MAPL_UNDEF
       self%firstRun  = .false.
    end if

    xh2o2 = self%h2o2_init 

    call SulfateUpdateOxidants (nymd, nhms, LONS, LATS, airdens, self%km, self%cdt, &
                                self%nymd_oxidants, MAPL_UNDEF, real(MAPL_RADIANS_TO_DEGREES), &
                                MAPL_AVOGAD/1000., MAPL_PI, MAPL_AIRMW, &
                                SU_OH, SU_NO3, SU_H2O2, &
                                xoh, xno3, xh2o2, self%recycle_h2o2, __RC__)

!   SU Settling
!   -----------
    do n = 1, self%nbins
       ! if radius == 0 then we're dealing with a gas which has no settling losses
       if (self%radius(n) == 0.0) then
          if (associated(SUSD)) SUSD(:,:,n) = 0.0
          cycle
       end if

       call MAPL_VarSpecGet(InternalSpec(n), SHORT_NAME=short_name, __RC__)
       call MAPL_GetPointer(internal, NAME=short_name, ptr=int_ptr, __RC__)

       call Chem_Settling (self%km, self%klid, n, self%rhFlag, self%cdt, MAPL_GRAV, &
                           self%radius(n)*1.e-6, self%rhop(n), int_ptr, t, airdens, &
                           rh2, zle, delp, SUSD, __RC__)
    end do

    allocate(drydepositionf, mold=lwi, __STAT__)
    call SulfateChemDriver (self%km, self%klid, self%cdt, MAPL_PI, real(MAPL_RADIANS_TO_DEGREES), MAPL_KARMAN, &
                            MAPL_AIRMW, MAPL_AVOGAD/1000., cpd, MAPL_GRAV, &
                            fMassMSA, fMassDMS, fMassSO2, fMassSO4, &
                            nymd, nhms, lons, lats, &
                            dms, so2, so4, dummyMSA, &
                            nDMS, nSO2, nSO4, nMSA, &
                            xoh, xno3, xh2o2, self%h2o2_init, &
                            delp, t, fcld, airdens, zle, &
                            ustar, sh, lwi, zpbl, z0h, &
                            SUDP, SUPSO2, SUPMSA, &
                            SUPSO4, SUPSO4g, SUPSO4aq, &
                            pso2, pmsa, pso4, pso4g, pso4aq, drydepositionf, & ! 3d diagnostics
                            __RC__)

    KIN = .true.
    call SU_Wet_Removal ( self%km, self%nbins, self%klid, self%cdt, kin, MAPL_GRAV, MAPL_AIRMW, &
                          delp, fMassSO4, fMassSO2, &
                          self%h2o2_init, ple, airdens, cn_prcp, ncn_prcp, pfl_lsan, pfi_lsan, t, &
                          nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, dummyMSA, &
                          SUWT, SUPSO4, SUPSO4WT, PSO4, PSO4WET, __RC__ )

!   Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
    call SU_Compute_Diags ( self%km, self%klid, self%radius(nSO4), self%sigma(nSO4), self%rhop(nSO4), &
                            MAPL_GRAV, MAPL_PI, nSO4, self%diag_Mie, &
                            self%wavelengths_profile*1.0e-9, self%wavelengths_vertint*1.0e-9, &
                            t, airdens, delp, ple,tropp, rh2, u, v, DMS, SO2, SO4, dummyMSA, &
                            DMSSMASS, DMSCMASS, &
                            MSASMASS, MSACMASS, &
                            SO2SMASS, SO2CMASS, &
                            SO4SMASS, SO4CMASS, &
                            SUEXTTAU, SUSTEXTTAU,SUSCATAU,SUSTSCATAU, SO4MASS, SUCONC, SUEXTCOEF, &
                            SUSCACOEF, SUANGSTR, SUFLUXU, SUFLUXV, SO4SAREA, SO4SNUM, __RC__)

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
    type (SU2G_GridComp), pointer      :: self

    real, pointer, dimension(:,:,:)  :: ptr3d_int, ptr3d_imp

    __Iam__('Run_data')

!*****************************************************************************
! Begin... 

! Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'//Iam

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'SU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Update interal data pointers with ExtData
!   -----------------------------------------
    call MAPL_GetPointer (internal, NAME='SO4', ptr=ptr3d_int, __RC__)
    call MAPL_GetPointer (import, NAME='climSO4', ptr=ptr3d_imp, __RC__)
    ptr3d_int = ptr3d_imp

    RETURN_(ESMF_SUCCESS)

  end subroutine Run_data


!-----------------------------------------------------------------------------------
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
    type(SU2G_GridComp), pointer                     :: self

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

    __Iam__('SU2G::aerosol_optics')

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

    type(GOCART2G_Mie),            intent(inout) :: mie              ! mie table
    integer,                       intent(in   ) :: nbins            ! number of bins
    integer,                       intent(in )   :: band             ! channel
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

    __Iam__('SU2G::aerosol_optics::mie_')

    bext_s  = 0.0d0
    bssa_s  = 0.0d0
    basym_s = 0.0d0

    do l = 1, nbins
       ! tau is converted to bext
       call mie%Query(band, l, q(:,:,:,l), rh, tau=bext, gasym=gasym, ssa=bssa, __RC__)

       bext_s  = bext_s  +             bext     ! extinction
       bssa_s  = bssa_s  +       (bssa*bext)    ! scattering extinction
       basym_s = basym_s + gasym*(bssa*bext)    ! asymetry parameter multiplied by scatering extiction 

    end do

    RETURN_(ESMF_SUCCESS)

    end subroutine mie_

  end subroutine aerosol_optics

!-----------------------------------------------------------------------------------
  subroutine monochromatic_aerosol_optics(state, rc)

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
    type(SU2G_GridComp), pointer                     :: self

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

    __Iam__('SU2G::monochromatic_aerosol_optics')

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
                 x = (ple(i,j,k) - ple(i,j,k-1))/MAPL_GRAV
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


end module SU2G_GridCompMod

