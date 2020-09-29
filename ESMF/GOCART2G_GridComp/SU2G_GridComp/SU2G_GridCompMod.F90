#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: SU2G_GridCompMod - GOCART refactoring of the SU gridded component 

! !INTERFACE:
module SU2G_GridCompMod

! !USES:
   use ESMF
   use MAPL
   use Chem_MieTableMod2G
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process       ! GOCART2G process library
   use GA_GridCompMod
   use m_StrTemplate          ! string templates

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

   integer, parameter     :: DP=kind(1.0d0)

!  gram molecular weights of species
   real, parameter :: fMassSulfur = 32., fMassSO2 = 64., fMassSO4 = 96., &
                      fMassDMS = 62., fMassMSA = 96.

   real, parameter ::  airmw  = 28.97   ! molecular weight of air
   real, parameter :: nAvogadro  = 6.022e23 ! molecules per mole of air 

   real, parameter ::  cpd    = 1004.16
   real, parameter ::  undefval  = 1.e15   ! missing value
   real, parameter ::  chemgrav   = 9.80616

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
   type, extends(GA_GridComp) :: SU2G_GridComp
      integer :: myDOW = -1     ! my Day of the week: Sun=1, Mon=2,...,Sat=7
      logical :: using_GMI_OH
      logical :: using_GMI_NO3
      logical :: using_GMI_H2O2
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
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

if(mapl_am_i_root()) print*,trim(comp_name),'2G SetServices BEGIN'

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file 
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'SU2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
       if (mapl_am_i_root()) print*,'SU2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! loading SU2G_GridComp_SU.data.rc instead'
       call ESMF_ConfigLoadFile (cfg, 'SU2G_GridComp_SU.rc', __RC__)
    end if

!   process generic config items
    call self%GA_GridComp%load_from_config( cfg, __RC__)

    allocate(self%sigma(self%nbins), __STAT__)

!   process SU-specific items
    call ESMF_ConfigGetAttribute(cfg, self%using_GMI_H2O2, label='using_GMI_H2O2:', __RC__)
    call ESMF_ConfigGetAttribute(cfg, self%using_GMI_OH,   label='using_GMI_OH:',   __RC__)
    call ESMF_ConfigGetAttribute(cfg, self%using_GMI_NO3,  label='using_GMI_NO3:',  __RC__)
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

!   This state is needed by MOIST - It will contain aerosols
!   ----------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                                  &
       SHORT_NAME = trim(COMP_NAME)//'_AERO_ACI',                                &
       LONG_NAME  = 'aerosol_cloud_interaction_aerosols_from_'//trim(COMP_NAME),  &
       UNITS      = 'kg kg-1',                                                   &
       DIMS       = MAPL_DimsHorzVert,                                           &
       VLOCATION  = MAPL_VLocationCenter,                                        &
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
    type (ESMF_State)                    :: aero, aero_aci
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg, cf
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
    real, pointer, dimension(:,:,:)      :: ple

    type(ESMF_Calendar)     :: calendar
    type(ESMF_Time)         :: currentTime
    type(ESMF_Time)         :: ringTime
    type(ESMF_TimeInterval) :: ringInterval
    integer                 :: year, month, day, hh, mm, ss

    real, dimension(4)   :: Vect_Hcts

    __Iam__('Initialize')

!****************************************************************************
!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, config=cf, __RC__)
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
    call ESMF_ConfigGetAttribute(cf, diurnal_bb, label='DIURNAL_BIOMASS_BURNING:', &
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
    call ESMF_ConfigLoadFile (cfg, 'SU2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'SU2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                    loading SU2G_GridComp_SU.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'SU2G_GridComp_SU.rc', __RC__)
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
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_ACI', aero_aci, __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

    call ESMF_StateGet (internal, 'SO4', field, __RC__)
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(3), __RC__)
    fld = MAPL_FieldCreate (field, 'SO4', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)
    call MAPL_StateAdd (aero_aci, fld, __RC__)

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
                                  label= "aerosol_monochromatic_optics_wavelength:", __RC__)

    allocate (self%diag_MieTable(instance)%mie_aerosol, __STAT__)
    self%diag_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%diag_MieTable(instance)%optics_file, __RC__ )
    call Chem_MieTableRead (self%diag_MieTable(instance)%mie_aerosol, self%diag_MieTable(instance)%nch, &
                            self%diag_MieTable(instance)%channels, rc, nmom=self%diag_MieTable(instance)%nmom)

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

    call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',             value=0,     __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet(aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

    call ESMF_AttributeSet(aero, name='internal_varaible_name', value='SO4', __RC__)

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
    where(1.01*so2anthro_l1_src > undefval)  so2anthro_l1_src = 0.
    where(1.01*so2anthro_l2_src > undefval)  so2anthro_l2_src = 0.
    where(1.01*so2ship_src > undefval)       so2ship_src = 0.
    where(1.01*so4ship_src > undefval)       so4ship_src = 0.

    aircraft_fuel_src = SU_AIRCRAFT
    so2biomass_src = SU_BIOMASS
    dmso_conc = SU_DMSO
    aviation_lto_src = SU_AVIATION_LTO
    aviation_cds_src = SU_AVIATION_CDS
    aviation_crs_src = SU_AVIATION_CRS

!   As a safety check, where value is undefined set to 0
    where(1.01*so2biomass_src > undefval)    so2biomass_src = 0.
    where(1.01*dmso_conc > undefval)         dmso_conc = 0.
    where(1.01*aircraft_fuel_src > undefval) aircraft_fuel_src = 0.
    where(1.01*aviation_lto_src > undefval ) aviation_lto_src = 0.
    where(1.01*aviation_cds_src > undefval ) aviation_cds_src = 0.
    where(1.01*aviation_crs_src > undefval ) aviation_crs_src = 0.

!   Update emissions/production if necessary (daily)
!   -----------------------------------------------
    if(self%nymd_last /= nymd) then
       self%nymd_last = nymd

!      Get pointwise SO2 and altitude of volcanoes from a daily file data base
       if(index(self%volcano_srcfilen,'volcanic_') /= 0) then
          call StrTemplate (fname, self%volcano_srcfilen, xid='unknown', &
                            nymd=nymd, nhms=120000 )
          call ReadPointEmissions (nymd, fname, self%nVolc, self%vLat, self%vLon, &
                                   self%vElev, self%vCloud, self%vSO2, self%vStart, &
                                   self%vEnd, label='volcano')
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
                                 self%km, self%cdt, chemgrav, zle, delp, area, self%vLat, self%vLon, __RC__)
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
    call SulfateDistributeEmissions ( self%km, self%nbins, self%cdt, chemgrav, nymd, nhms, &
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
                                      aviation_crs_src, rc) 

    if (associated(dms)) then 
       call DMSemission (self%km, self%cdt, chemgrav, t, u10m, v10m, lwi, delp, &
                         fMassDMS, SU_DMSO, dms, SUEM, nDMS, rc)
    end if

!   Read any pointwise emissions, if requested
!   ------------------------------------------
    if(self%doing_point_emissions) then
       call StrTemplate (fname, self%point_emissions_srcfilen, xid='unknown', &
                         nymd=nymd, nhms=120000 )
       call ReadPointEmissions (nymd, fname, self%nPts, self%pLat, self%pLon, &
                                 self%pBase, self%pTop, self%pEmis, self%pStart, &
                                 self%pEnd, label='source')

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

    real, pointer, dimension(:,:,:)  ::  oh, no3, h2o2
    real, dimension(:,:,:), allocatable ::  xoh, xno3, xh2o2

    real, dimension(:,:), allocatable :: drydepositionf
    real, pointer, dimension(:,:,:) :: dummyMSA => null() ! this is so the model can run without MSA enabled
    

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
    self%recycle_h2o2 = ESMF_AlarmIsRinging(alarm, __RC__)

!   Get oxidant pointers from specified provider
!   ----------------------------------------------
    call GetOxidant (self%using_GMI_OH, oh, 'OH', __RC__)
    call GetOxidant (self%using_GMI_OH, no3, 'NO3', __RC__)
    call GetOxidant (self%using_GMI_OH, h2o2, 'H2O2', __RC__)

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

if(mapl_am_i_root()) print*,'SU2G Run2 sum(xh2o2) = ',sum(xh2o2)
if(mapl_am_i_root()) print*,'SU2G Run2 sum(self%h2o2_init) = ',sum(self%h2o2_init)
!if(mapl_am_i_root()) print*,'SU2G Run2 sum(xoh) = ',sum(xoh)

    call SulfateUpdateOxidants (nymd, nhms, LONS, LATS, &
                                airdens, self%km, self%cdt, &
                                self%nymd_oxidants, undefval, &
                                oh, no3, h2o2, &
                                xoh, xno3, xh2o2, self%recycle_h2o2, rc)

if(mapl_am_i_root()) print*,'SU2G Run2 UpdateOxidants sum(xh2o2) = ',sum(xh2o2)
if(mapl_am_i_root()) print*,'SU2G Run2 UpdateOxidants sum(self%h2o2_init) = ',sum(self%h2o2_init)
!if(mapl_am_i_root()) print*,'SU2G Run2 UpdateOxidants sum(xoh) = ',sum(xoh)

!   SU Settling
!   -----------
    do n = 1, self%nbins
       call MAPL_VarSpecGet(InternalSpec(n), SHORT_NAME=short_name, __RC__)
       call MAPL_GetPointer(internal, NAME=short_name, ptr=int_ptr, __RC__)

       call Chem_Settling2Gorig (self%km, self%klid, self%rhFlag, n, int_ptr, CHEMgrav, delp, &
                                 self%radius(n)*1.e-6, self%rhop(n), self%cdt, t, airdens, &
                                 rh2, zle, SUSD, rc=rc)
    end do

    allocate(drydepositionf, mold=lwi, __STAT__)
    call SulfateChemDriver (self%km, self%klid, self%cdt, MAPL_PI, real(MAPL_RADIANS_TO_DEGREES), MAPL_karman, &
                            airmw, nAvogadro, cpd, chemgrav, &
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
                            rc)

if(mapl_am_i_root()) print*,'SU2G Run2 ChemDriver sum(xh2o2) = ',sum(xh2o2)
if(mapl_am_i_root()) print*,'SU2G Run2 ChemDriver sum(self%h2o2_init) = ',sum(self%h2o2_init)
!if(mapl_am_i_root()) print*,'SU2G Run2 ChemDriver sum(xoh) = ',sum(xoh)
if(mapl_am_i_root()) print*,'SU2G Run2 sum(SUPSO4aq) = ',sum(SUPSO4aq)
if(mapl_am_i_root()) print*,'SU2G Run2 sum(SUPSO4g) = ',sum(SUPSO4g)

    KIN = .true.
    call SU_Wet_Removal ( self%km, self%nbins, self%klid, self%cdt, kin, chemgrav, airMW, &
                          delp, fMassSO4, fMassSO2, &
                          self%h2o2_init, ple, airdens, cn_prcp, ncn_prcp, pfl_lsan, pfi_lsan, t, &
                          nDMS, nSO2, nSO4, nMSA, DMS, SO2, SO4, dummyMSA, &
                          SUWT, SUPSO4, SUPSO4WT, PSO4, PSO4WET, rc )

if(mapl_am_i_root()) print*,'SU2G Run2 WetRemoval sum(xh2o2) = ',sum(xh2o2)
if(mapl_am_i_root()) print*,'SU2G Run2 WetRemoval sum(self%h2o2_init) = ',sum(self%h2o2_init)

    call SU_Compute_Diags ( self%km, self%klid, self%radius(nSO4), self%sigma(nSO4), self%rhop(nSO4), &
                            chemgrav, MAPL_pi, nSO4, self%diag_MieTable(self%instance), &
                            self%diag_MieTable(self%instance)%channels, &
                            t, airdens, delp, rh2, u, v, DMS, SO2, SO4, dummyMSA, &
                            DMSSMASS, DMSCMASS, &
                            MSASMASS, MSACMASS, &
                            SO2SMASS, SO2CMASS, &
                            SO4SMASS, SO4CMASS, &
                            SUEXTTAU, SUSCATAU, SO4MASS, SUCONC, SUEXTCOEF, &
                            SUSCACOEF, SUANGSTR, SUFLUXU, SUFLUXV, SO4SAREA, SO4SNUM, rc)

if(mapl_am_i_root()) print*,'SU2G Run2 E sum(SO2) = ',sum(SO2)
if(mapl_am_i_root()) print*,'SU2G Run2 E sum(SO4) = ',sum(SO4)
if(mapl_am_i_root()) print*,'SU2G Run2 E sum(DMS) = ',sum(DMS)
if (associated(dummyMSA)) then
   if(mapl_am_i_root()) print*,'SU2G Run2 E sum(MSA) = ',sum(dummyMSA)
end if

    RETURN_(ESMF_SUCCESS)

  contains
!..................................................................
  subroutine GetOxidant (using_GMI, ptr, ptr_name, rc)

     logical, intent(in) :: using_GMI
     real, pointer, dimension(:,:,:), intent(inout) :: ptr
     character(len=*), intent(in) :: ptr_name
     integer, optional, intent(out) :: rc

!      Begin...
       rc = 0
       if (using_GMI) then
          call MAPL_GetPointer(import, ptr, trim(ptr_name), __RC__)
       else
          call MAPL_GetPointer(import, ptr, 'SU_'//trim(ptr_name), __RC__)
       end if

  end subroutine GetOxidant
!..................................................................


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

if (mapl_am_I_root()) print*,trim(comp_name),' Run_data BEGIN'

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

if (mapl_am_I_root()) print*,trim(comp_name),' Run_data END'

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
    integer                                          :: band, offset
    integer, parameter                               :: n_bands = 1

    integer :: i, j, k

    __Iam__('SU2G::aerosol_optics')

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


    __Iam__('SU2G::aerosol_optics::mie_')

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



!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver

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
   real, pointer, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]  
   real, pointer, dimension(:,:,:), intent(in) :: tmpu   ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: cloud  ! cloud fraction for radiation [1]
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
   if( associated(su_pSO4) )  su_pSO4g  = 0.
   if( associated(su_pSO4) )  su_pSO4aq = 0.
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
                        von_karman, cpd, grav, z0h, drydepositionfrequency, rc )

if(mapl_am_i_root()) print*,'SU2G sum(drydepositionfrequency) = ',sum(drydepositionfrequency)


!  Now call the chemistry packages...
!  ----------------------------------
!  DMS source and oxidation to SO2 and MSA
   call SulfateChemDriver_DMS (km, klid, cdt, airMolWght, nAvogadro, cpd,&
                               fMassMSA, fMassDMS, fMassSO2, &
                               dms, nDMS, xoh, xno3, &
                               cossza, tmpu, rhoa, &
                               pSO2_DMS, pMSA_DMS, SU_dep, &
                               rc)
if(mapl_am_i_root()) print*,'SU2G sum(pSO2_DMS) = ',sum(pSO2_DMS)


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

if(mapl_am_i_root()) print*,'SU2G klid = ',klid


if(mapl_am_i_root()) then
   print*,' SU2G_ChemDrv_SO2 CHECK'
   print*,'=========================='
   print*,'xoh = ',sum(xoh)
   print*,'xh2o2 = ',sum(xh2o2)
   print*,'so2 = ',sum(so2)
   print*,'pSO4g_SO2 = ',sum(pSO4g_SO2)
   print*,'pSO4aq_SO2 = ',sum(pSO4aq_SO2)
   print *,'i1 = ' 
   print *,'i2 = ',ubound(tmpu,1)
   print *,'j1 = '
   print *,'j2 = ',ubound(tmpu,2)
   print *,'km = ',km
   print *,'cdt = ',cdt
   print *,'rhoa = ',sum(rhoa)
   print *,'delp = ',sum(delp)
   print *,'tmpu = ',sum(tmpu)
   print *,'cloud = ',sum(cloud)
   print *,'oro = ',sum(oro)
   print*,'======== END ============='
end if

!  SO2 source and oxidation to SO4
   call SulfateChemDriver_SO2 (km, klid, cdt, airMolWght, nAvogadro, cpd, grav, &
                               fMassSO4, fMassSO2, &
                               so2, nSO2, xoh, xh2o2, &
                               tmpu, rhoa, delp, oro, cloud, drydepositionfrequency, &
                               pSO2_DMS, pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                               rc)

   if( associated(pSO4g) )  pSO4g = pSO4g_SO2
   if( associated(su_pSO4g)) then
     do k = klid, km
      su_pSO4g(:,:) = su_pSO4g(:,:) + pSO4g_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

!if(mapl_am_i_root()) print*,'SU2G pSO4g_SO2 = ',sum(pSO4g_SO2)

   if( associated(pSO4aq) )  pSO4aq = pSO4aq_SO2
   if( associated(su_pSO4aq)) then
     do k = klid, km
      su_pSO4aq(:,:) = su_pSO4aq(:,:) + pSO4aq_SO2(:,:,k)*delp(:,:,k)/grav
     enddo
   endif

if(mapl_am_i_root()) print*,'SU2G pSO4aq_SO2 = ',sum(pSO4aq_SO2)

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
                               rc)

!  MSA source and loss
   if( associated(msa)) then
      call SulfateChemDriver_MSA (km, klid, cdt, grav, msa, nMSA, delp, &
                                  drydepositionfrequency, pMSA_DMS, SU_dep, &
                                  rc)
   end if

!  Save the h2o2 value after chemistry
   h2o2_init = xh2o2

   rc = 0

   end subroutine SulfateChemDriver

!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_DMS

   subroutine SulfateChemDriver_DMS (km, klid, cdt, airMolWght, nAvogadro, cpd, &
                                     fMassMSA, fMassDMS, fMassSO2, &
                                     qa, nDMS, xoh, xno3, &
                                     cossza, tmpu, rhoa, &
                                     pSO2_DMS, pMSA_DMS, SU_dep, &
                                     rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:  
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: nAvogadro  ! Avogadro's number [1/kmol]
   real, intent(in)    :: airMolWght ! molecular weight of air [kg/kmol]
   real, intent(in)    :: cpd
   real, intent(in)    :: fMassMSA, fMassDMS, fMassSO2 ! gram molecular weights of species 
   integer, intent(in) :: nDMS       !index position of sulfates
   real, dimension(:,:,:), intent(in) :: xoh, xno3  ! OH, NO3 respectievly [kg/kg]
   real, dimension(:,:), intent(in)   :: cossza
   real, dimension(:,:,:), intent(in) :: tmpu   ! temperature [K]
   real, dimension(:,:,:), intent(in) :: rhoa   ! layer air density [kg/m^3]


! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg] 
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   real, dimension(:,:,:), allocatable,  intent(out) :: pSO2_DMS ! SO2 production from DMS oxidation [kg kg-1 s-1]
   real, dimension(:,:,:), allocatable,  intent(out) :: pMSA_DMS ! MSA production from DMS oxidation [kg kg-1 s-1]
   integer, optional, intent(out)   :: rc

! !DESCRIPTION: Computes the production of SO2 and MSA due to DMS oxidation
!
!   R1:    DMS + OH  -> a*SO2 + b*MSA                OH addition channel
!          k1 = { 1.7d-42*exp(7810/T)*[O2] / (1+5.5e-31*exp(7460/T)*[O2] }
!          a = 0.75, b = 0.25
!
!   R2:    DMS + OH  ->   SO2 + ...                  OH abstraction channel
!          k2 = 1.2e-11*exp(-260/T)
!
!      DMS_OH = DMS0 * exp(-(r1+r2)*NDT1)
!          where DMS0 is the DMS concentration at the beginning,
!          r1 = k1*[OH], r2 = k2*[OH]
!
!   R3:    DMS + NO3 ->   SO2 + ...
!          k3 = 1.9e-13*exp(500/T)
!
!      DMS = DMS_OH * exp(-r3*NDT1)
!          where r3 = k3*[NO3]
!
!   R4:    DMS + X   ->   SO2 + ...
!          assume to be at the rate of DMS+OH and DMS+NO3 combined.
!
!   The production of SO2 and MSA here, PSO2_DMS and PMSA_DMS, are saved
!   for use in CHEM_SO2 and CHEM_MSA subroutines as a source term.  They
!   are in unit of MixingRatio/second.
!
! !REVISION HISTORY:
!
!  06Nov2003, Colarco
!  Based on Ginoux
!
!  03Aug2020 E.Sherman - ported to process library

! !Local Variables
   integer :: i, j, k, i1=1, j1=1, i2, j2
   real*8  :: Fx, b, eff
   real*8  :: rk1, rk2, rk3, rk4
   real*8  :: tk, o2, oh, no3, air
   real*8  :: dms, dms0, dms_oh

   data Fx  / 1.0 /
   data b   / 0.25 /
   data eff / 1. /

!EOP
!-------------------------------------------------------------------------
! Begin...

   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

   allocate(pSO2_DMS, mold=tmpu)
   allocate(pMSA_DMS, mold=tmpu)

!  spatial loop 
   do k = klid, km
    do j = j1, j2
     do i = i1, i2

      rk1 = 0.
      rk2 = 0.
      rk3 = 0.
      rk4 = 0.

      tk  = tmpu(i,j,k)
      oh  = xoh(i,j,k)
!     air molecules in # cm-3
      air = 1000.*rhoa(i,j,k) / airMolWght * nAvogadro * 1.e-6
!     oxygen molecules in # cm-3
      o2 = 0.21 * air
!     no3 -> go from volume mixing ratio to # cm-3
      no3 = xno3(i,j,k) * air

!     initial DMS concentration (kg kg-1)
      dms0 = qa(i,j,k)
      dms0 = max(dms0,tiny(dms0))

!     1 & 2) DMS + OH: RK1 = addition, RK2 = abstraction
      if(oh .gt. 0.) then
       rk1 = (1.7d-42 * exp(7810./tk) * o2) / &
             (1. + 5.5e-31 * exp(7460./tk) * o2) * oh
       rk2 = (1.2e-11 * exp(-260./tk)) * oh
      endif

!     3) DMS + NO3: only happens at night
      if(cossza(i,j) .le. 0.) then
       rk3 = (1.9e-13 * exp(500./tk)) * no3
      endif

!     Now do the DMS loss
      dms_oh = dms0   * exp( -(rk1+rk2)* Fx * cdt)
      dms    = dms_oh * exp( -(rk3)    * Fx * cdt)

!     SO2 and MSA production terms
!     MSA is formed from the DMS+OH addition step
!     Production should go as mass mixing ratio change in MSA
      if( (rk1+rk2) .eq. 0.) then
       pMSA_DMS(i,j,k) = 0.
      else
       pMSA_DMS(i,j,k) =  (dms0 - dms_oh) * b*rk1/((rk1+rk2)*Fx) * eff &
                         * (fMassMSA/fMassDMS) / cdt
      endif

!     Everything else goes into SO2 formation step
      pSO2_DMS(i,j,k) = ( dms0 - dms - &
                          pMSA_DMS(i,j,k)*cdt*(fMassDMS/fMassMSA) &
                        ) * (fMassSO2/fMassDMS) / cdt

!     4) Dry deposition of DMS (not in GOCART?)
!      if(k .eq. km) rk4 = drydepf(i,j)
!      dms0 = dms
!      dms  = dms0 * exp(-rk4*cdt)
!      dms    = max(dms,1.e-32)

!     Update the mass mixing ratio and the dry deposition flux out of DMS
      dms    = max(dms,tiny(dms))
      qa(i,j,k) = dms

     end do ! i
    end do  ! j
    if(k .eq. km .and. associated(SU_dep) ) SU_dep(:,:,nDMS) = 0.
   end do   ! k


   rc = 0

   end subroutine SulfateChemDriver_DMS


!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_SO2

   subroutine SulfateChemDriver_SO2 (km, klid, cdt, airMolWght, nAvogadro, cpd, grav, &
                                     fMassSO4, fMassSO2, &
                                     qa, nSO2, xoh, xh2o2, &
                                     tmpu, rhoa, delp, oro, cloud, drydepf, &
                                     pSO2_DMS, pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                                     rc)

! !USES:
   implicit NONE

! !INPUT PARAMETERS:  
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: nAvogadro  ! Avogadro's number [1/kmol]
   real, intent(in)    :: airMolWght ! molecular weight of air [kg/kmol]
   real, intent(in)    :: cpd
   real, intent(in)    :: grav       ! gravity [m/sec]
   real, intent(in)    :: fMassSO4, fMassSO2 ! gram molecular weights of species 
   integer, intent(in) :: nSO2       !index position of sulfates
   real, pointer, dimension(:,:,:), intent(in) :: tmpu   ! temperature [K]
   real, pointer, dimension(:,:,:), intent(in) :: rhoa   ! layer air density [kg/m^3]
   real, pointer, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]  
   real, pointer, dimension(:,:,:), intent(in) :: cloud  ! cloud fraction for radiation [1]
   real, dimension(:,:), intent(in)   :: drydepf  ! dry deposition frequency [s-1]
   real, pointer, dimension(:,:), intent(in) :: oro  ! land-ocean-ice mask
   real, dimension(:,:,:), intent(in) :: pSO2_DMS ! SO2 production from DMS oxidation [kg m-2 s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg] 
   real, dimension(:,:,:), intent(inout) :: xoh, xh2o2  ! OH, H2O2 respectievly [kg/kg]
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
!   real, dimension(:,:,:), allocatable, intent(out) :: pSO4g_SO2 ! SO4 production - gas phase [kg kg-1 s-1]
!   real, dimension(:,:,:), allocatable, intent(out) :: pSO4aq_SO2 ! SO4 production - aqueous [kg kg-1 s-1]
   real, dimension(:,:,:), intent(inout) :: pSO4g_SO2 ! SO4 production - gas phase [kg kg-1 s-1]
   real, dimension(:,:,:), intent(inout) :: pSO4aq_SO2 ! SO4 production - aqueous [kg kg-1 s-1]
   integer, optional, intent(out)   :: rc

! !DESCRIPTION: Computes the concentration of SO2 and production of SO4
!
!  SO2 production:

!    DMS + OH, DMS + NO3 (saved in SU_ChemDrv_DMS)
!
!  SO2 loss:
!    SO2 + OH  -> SO4
!    SO2       -> drydep
!    SO2 + H2O2 or O3 (aq) -> SO4
!
!  SO2 = SO2_0 * exp(-bt)
!      + PSO2_DMS*dt/bt * [1-exp(-bt)]
!    where b is the sum of the reaction rate of SO2 + OH and the dry
!    deposition rate of SO2, PSO2_DMS is SO2 production from DMS in
!    MixingRatio/timestep.
!
!  If there is cloud in the gridbox (fraction = fc), then the aqueous
!  phase chemistry also takes place in cloud. The amount of SO2 oxidized
!  by H2O2 in cloud is limited by the available H2O2; the rest may be
!  oxidized due to additional chemistry, e.g, reaction with O3 or O2
!  (catalyzed by trace metal).
!
! !REVISION HISTORY:
!  06Nov2003, Colarco - Based on Ginoux!
!  15Jul2010, Colarco - modularized
!  03Aug2020 E.Sherman - ported to process library


! !Local Variables
   integer :: i, j, k, j2, i2
   real*8  :: rk1, rk2, rk, rkt, f1
   real*8  :: L1, L2, Ld, SO2, SO2_cd, fc, fMR
   real*8  :: oh, h2o2, SO20, tk, air, k0, ki, kk
   real, dimension(:,:), allocatable :: fout

   data ki / 1.5e-12 /

!EOP
!-------------------------------------------------------------------------
! Begin...

   j2 = ubound(tmpu, 2)
   i2 = ubound(tmpu, 1)

!   allocate(pSO4g_SO2, mold=tmpu)
!   allocate(pSO4aq_SO2, mold=tmpu)
   allocate(fout(i2,j2))

!  Conversion of SO2 mmr to SO2 vmr
   fMR = airMolWght / fMassSO2

!  Initialize flux variable   
   fout = 0.

!  spatial loop 
   do k = klid, km
    do j = 1, j2
     do i = 1, i2

      rk1 = 0.
      rk2 = 0.
      L1  = 0.
      L2  = 0.
      Ld  = 0.

      tk   = tmpu(i,j,k)
      oh   = xoh(i,j,k)
      h2o2 = max(xh2o2(i,j,k),tiny(xh2o2(i,j,k)))

!     air molecules in # cm-3
      air  = 1000.*rhoa(i,j,k) / airMolWght * nAvogadro * 1.e-6
!     1) SO2 + OH(g) in s-1
      k0 = 3.0e-31 * (300./tk)**4.3
      kk = k0 * air / ki
      f1 = (1. + (log10(kk))**2.)**(-1.)
      rk1 = ( (k0*air/(1.+kk)) * 0.6**f1) * oh

!     2) rk2 is the loss of SO2 due to dry deposition.
      if(k .eq. km) then
!      drydepf calculated for aerosol
!      follow Walcek: ocean drydepf_so2 = 10*drydepf_aer
!      or if land drydepf_so2 = 3*drydepf_aer
       if(oro(i,j) .eq. OCEAN) then
        rk2 = 10.*drydepf(i,j)
       else
        rk2 = 3.*drydepf(i,j)
       endif
!       rk2 = drydepf(i,j)
      else
       rk2 = 0.
      endif

      rk = (rk1 + rk2)
      rkt = rk*cdt

!     Update the SO2 concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.

!     initial SO2 concentration (kg kg-1) after adding source
      SO20 = qa(i,j,k) + pSO2_DMS(i,j,k)*cdt
      SO20 = max(SO20,tiny(SO20))

      if(rk .gt. 0.) then
       SO2_cd =  SO20 * exp(-rkt)
       L1     = (SO20 - SO2_cd) * rk1/rk
       if(k .eq. km) then
        Ld    = (SO20 - SO2_cd) * rk2/rk
        fout(i,j) = Ld * delp(i,j,km)/grav/cdt
       else
        Ld    = 0.
       endif
      else
       SO2_cd = SO20
       L1     = 0.
      endif

!     Update SO2 concentration after cloud chemistry, if it occurs
      fc = cloud(i,j,k)
      if(fc .gt. 0. .and. SO2_cd .gt. 0. .and. tk .gt. 258.) then
!      Check on H2O2 vmr -> is SO2 vmr greater?
       if(fMr * SO2_cd .gt. h2o2) then
        fc = fc*(h2o2/(fMR*SO2_cd))
        h2o2 = h2o2*(1.-cloud(i,j,k))
       else
        h2o2 = h2o2*(1. - cloud(i,j,k)*(fMR*SO2_cd)/h2o2)
       endif
       SO2 = SO2_cd*(1.-fc)
!      aqueous loss rate (mixing ratio/timestep)
       L2 = SO2_cd * fc
      else
       SO2 = SO2_cd
       L2 = 0.
      endif

!     Ideally you would update the H2O2 mixing ratio at this point
!     and then reset it periodically
      xh2o2(i,j,k) = max(h2o2,tiny(h2o2))

      SO2 = max(SO2,tiny(SO2))
      qa(i,j,k) = SO2
      pSO4g_SO2(i,j,k) = L1 * (fMassSO4/fMassSO2) / cdt
      pSO4aq_SO2(i,j,k) = L2 * (fMassSO4/fMassSO2) / cdt

     end do
    end do
   end do

   if( associated(SU_dep) ) SU_dep(:,:,nSO2) = fout

   rc = 0

   end subroutine SulfateChemDriver_SO2

!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_SO4

   subroutine SulfateChemDriver_SO4 (km, klid, cdt, grav, qa, nSO4, delp, drydepf, &
                                     pSO4g_SO2, pSO4aq_SO2, SU_dep, &
                                     rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:  
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: grav   ! gravity [m/sec]
   integer, intent(in) :: nSO4   ! index position of sulfate
   real, pointer, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]  
   real, dimension(:,:), intent(in)   :: drydepf    ! dry deposition frequency [s-1]
   real, dimension(:,:,:), intent(in) :: pSO4g_SO2  ! SO4 production - gas phase [kg kg-1 s-1]
   real, dimension(:,:,:), intent(in) :: pSO4aq_SO2 ! SO4 production - aqueous [kg kg-1 s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg] 
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc

! !DESCRIPTION:
!
!  SO4 production:
!    The only production term is due to SO2 oxidation.
!    SO4 = SO4_0 * exp(-kt) + pSO4_SO2/kt * (1.-exp(-kt))
!     where k is the dry deposition
!
! !REVISION HISTORY:
!
!  15Jul2010, Colarco - Modularized
!  06Nov2003, Colarco
!  Based on Ginoux
!
!  03Aug2020 E.Sherman - ported to process library

!
! !Local Variables
   integer :: i, j, k, i2, j2
   real*8  :: rk, rkt, Ld
   real*8  :: SO4, SO40, pSO4
   real, dimension(:,:), allocatable :: fout

!EOP
!-------------------------------------------------------------------------

! Begin...

   j2 = ubound(qa, 2)
   i2 = ubound(qa, 1)

   allocate(fout(i2,j2))

!  Initialize flux variable
   fout = 0.

!  spatial loop 
   do k = klid, km
    do j = 1, j2
     do i = 1, i2

      pSO4 = pSO4g_SO2(i,j,k)+pSO4aq_SO2(i,j,k)

!     initial SO4 concentration (kg kg-1)
      SO40 = qa(i,j,k)
      SO40 = max(SO40,tiny(SO40))

!     Update the SO4 concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.
!     RK is the dry deposition frequency
      if(k .eq. km) then
       RK = drydepf(i,j)
       RKT = RK*cdt
       SO4 = (SO40 + pSO4*cdt) * exp(-rkt)
       Ld  = (SO40 - SO4 + pSO4*cdt)
       fout(i,j) = Ld * delp(i,j,km)/grav/cdt
      else
       SO4 = SO40 + pSO4*cdt
       Ld = 0.
      endif

      SO4 = max(SO4,tiny(SO4))
      qa(i,j,k) = SO4

     end do
    end do
   end do

   if( associated(SU_dep) ) SU_dep(:,:,nSO4) = fout

   rc = 0

   end subroutine SulfateChemDriver_SO4

!==================================================================================
!BOP
! !IROUTINE: SulfateChemDriver_MSA

   subroutine SulfateChemDriver_MSA (km, klid, cdt, grav, qa, nMSA, delp, drydepf, &
                                     pMSA_DMS, SU_dep, &
                                     rc)
! !USES:
   implicit NONE

! !INPUT PARAMETERS:  
   integer, intent(in) :: km     ! number of model levels
   integer, intent(in) :: klid   ! index for pressure lid
   real, intent(in)    :: cdt    ! chemisty model timestep [sec]
   real, intent(in)    :: grav   ! gravity [m/sec]
   integer, intent(in) :: nMSA   ! index position of sulfate
   real, pointer, dimension(:,:,:), intent(in) :: delp   ! pressure thickness [Pa]  
   real, dimension(:,:), intent(in)   :: drydepf   ! dry deposition frequency [s-1]
   real, dimension(:,:,:), intent(in) :: pMSA_DMS  ! MSA production - gas phase [kg kg-1 s-1]

! !INOUTPUT PARAMETERS:
   real, dimension(:,:,:), intent(inout) :: qa  ! dimethyl sulfide [kg/kg] 
   real, pointer, dimension(:,:,:), intent(inout) :: SU_dep ! Sulfate Dry Deposition All Bins [kg m-2 s-1]

! !OUTPUT PARAMETERS:
   integer, optional, intent(out)   :: rc

! !DESCRIPTION: 
!
!  MSA production:
!    The only production term is due to DMS oxidation.
!    MSA = MSA_0 * exp(-kt) + pSO4_SO2/kt * (1.-exp(-kt))
!     where k is the dry deposition
!
! !REVISION HISTORY:
!
!  15Jul2010, Colarco -- modularized
!  06Nov2003, Colarco
!  Based on Ginoux
!
!  03Aug2020 E.Sherman - ported to process library

! !Local Variables
   integer :: i, j, k, i2, j2
   real*8  :: rk, rkt, Ld
   real*8  :: MSA, MSA0
   real, dimension(:,:), allocatable :: fout

!EOP
!-------------------------------------------------------------------------
! Begin...

   j2 = ubound(qa, 2)
   i2 = ubound(qa, 1)

   allocate(fout(i2,j2))

!  spatial loop 
   do k = klid, km
    do j = 1, j2
     do i = 1, i2

!     initial MSA concentration (kg kg-1)
      MSA0 = qa(i,j,k)
      MSA0 = max(MSA0,tiny(MSA0))

!     Update the MSA concentration
!     Originally this was solved like a simple exponential solution
!     after Jacobson eq. 13.38, which is more accurate but not mass
!     conserving.  We've already timesplit everything, so accuracy is
!     out to lunch, and I'd prefer to conserve mass.
!     RK is the dry deposition frequency
      if(k .eq. km) then
       RK = drydepf(i,j)
       RKT = RK*cdt
       MSA = (MSA0 + pMSA_DMS(i,j,k)*cdt) * exp(-rkt)
       Ld  = (MSA0 + pMSA_DMS(i,j,k)*cdt - MSA)
       fout(i,j) = Ld * delp(i,j,km)/grav/cdt
      else
       MSA = MSA0 + pMSA_DMS(i,j,k)*cdt
       Ld = 0.
      endif

      MSA = max(MSA,tiny(MSA))
      qa(i,j,k) = MSA

     end do
    end do
   end do

   if( associated(SU_dep) ) SU_dep(:,:,nMSA) = fout

   rc = 0


   end subroutine SulfateChemDriver_MSA



















end module SU2G_GridCompMod

