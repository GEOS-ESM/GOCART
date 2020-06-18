#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: SS2G_GridCompMod - GOCART refactoring of the SS gridded component 

! !INTERFACE:
module SS2G_GridCompMod

! !USES:
   use ESMF
   use MAPL
   use Chem_MieTableMod2G
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use Chem_UtilMod
   use GOCART2G_Process       ! GOCART2G process library

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2
   real, parameter :: OCEAN=0.0, LAND = 1.0, SEA_ICE = 2.0
   integer, parameter     :: DP=kind(1.0d0)

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices

real, parameter ::  chemgrav   = 9.80616

! !DESCRIPTION: This module implements GOCARTS' Sea Salt (SS) Gridded Component.

! !REVISION HISTORY:
! 24Oct2019  E.Sherman  First attempt at refactoring.

!EOP
!===========================================================================
   integer, parameter         :: NHRES = 6  ! DEV NOTE!!! should this be allocatable, and not a parameter?

!  !Sea Salt state
   type SS2G_GridComp
       type(Chem_Mie), dimension(2)    :: rad_MieTable, diag_MieTable
       real, allocatable      :: radius(:)      ! particle effective radius [um]
       real, allocatable      :: rlow(:)        ! particle effective radius lower bound [um]
       real, allocatable      :: rup(:)         ! particle effective radius upper bound [um]
       real, allocatable      :: fscav(:)       ! scavenging efficiency
       real, allocatable      :: molwght(:)     ! molecular weight
       real, allocatable      :: rhop(:)        ! dry particle density
       real, allocatable      :: fnum(:)        ! number of particles per kg mass
       real                   :: maringFlag     ! maring settling velocity correction
       integer                :: rhFlag         ! RH swelling of Seasalt (1 for Fitzgerald 1975; 2 for Gerber 1985)
       integer                :: sstEmisFlag    ! Choice of SST correction to emissions: 0 - none; 1 - Jaegle et al. 2011; 2 - GEOS5
       logical                :: hoppelFlag     ! Apply the Hoppel correction to emissions (Fan and Toon, 2011)
       logical                :: weibullFlag    ! Apply the Weibull distribution to wind speed for emissions (Fan and Toon, 2011)
       integer                :: nbins
       integer                :: km             ! vertical grid dimension
       real                   :: CDT            ! chemistry timestep (secs)
       integer                :: instance       ! data or computational instance
       real, allocatable      :: deep_lakes_mask(:,:)
       integer                :: emission_scheme
       real                   :: emission_scale ! global scaling factor
       real                   :: emission_scale_res(NHRES) ! global scaling factor
   end type SS2G_GridComp

   type wrap_
      type (SS2G_GridComp), pointer     :: PTR => null()
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
    type (SS2G_GridComp), pointer               :: self
    type (Chem_Mie)                             :: this

    character (len=ESMF_MAXSTR)                 :: field_name
    character (len=ESMF_MAXSTR), allocatable    :: aerosol_names(:)

    integer                                     :: n, i, nCols, nbins
    real                                        :: DEFVAL
    logical                                     :: data_driven=.true.

    !development testing variables - to be deleted
    real, dimension(:,:), pointer       :: ptr_test

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' // Iam

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file 
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'SS2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
       if (mapl_am_i_root()) print*,'SS2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! loading SS2G_GridComp_SS.data.rc instead'
       call ESMF_ConfigLoadFile (cfg, 'SS2G_GridComp_SS.rc', __RC__)
    end if

!   Get nbins from cfg
    call ESMF_ConfigGetAttribute (cfg, self%nbins, label='nbins:', __RC__)
    nbins = self%nbins

!   Parse config file into private internal state
!   ----------------------------------------------
    allocate(self%radius(nbins), self%rlow(nbins), self%rup(nbins), self%fscav(nbins), &
             self%molwght(nbins), self%fnum(nbins), self%rhop(nbins), __STAT__)

    call ESMF_ConfigGetAttribute (cfg, self%radius,     label='particle_radius_microns:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rlow,       label='radius_lower:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rup,        label='radius_upper:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%fscav,      label='fscav:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%molwght,    label='molecular_weight:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rhop,       label='SS_density:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%fnum,       label='fnum:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%sstEmisFlag, label='sstEmisFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%weibullFlag,  label='weibullFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%hoppelFlag, label='hoppelFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rhFlag, label='rhFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%emission_scheme, label='emission_scheme:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%emission_scale_res, label='emission_scale:', __RC__)


!   Is SS data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_INITIALIZE,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_METHOD_RUN, Run, __RC__)
    if (data_driven /= .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if

!   INTERNAL STATE
!   ---------------
!   Default INTERNAL state values
!   -----------------------------
    DEFVAL = 0.0

#include "SS2G_Internal___.h"

!   IMPORT STATE
!   -------------
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

        do i = 1, nbins
            write(field_name, '(A, I0.3)') '', i
            call MAPL_AddImportSpec(GC,                                           &
              SHORT_NAME = 'climss'//trim(field_name),                            &
              LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')',  &
              UNITS      = 'kg kg-1',                                             &
              RESTART    = MAPL_RestartSkip,                                      &
              DIMS       = MAPL_DimsHorzVert,                                     &
              VLOCATION  = MAPL_VLocationCenter, __RC__)

!           ! dry deposition
            call MAPL_AddImportSpec(GC,                                           &
              SHORT_NAME = 'climSSDP'//trim(field_name),                          &
              LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')',  &
              UNITS      = 'kg kg-1',                                             &
              DIMS       = MAPL_DimsHorzOnly,                                     &
              VLOCATION  = MAPL_VLocationCenter,                                  &
              RESTART    = MAPL_RestartSkip, __RC__)

!           ! wet deposition    
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSWT'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!           ! gravitational settling
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSSD'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!        ! convective scavenging
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSSV'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if ! (data_driven)


!   EXPORT STATE
!   -------------
    if (.not. data_driven) then
#include "SS2G_Export___.h"
#include "SS2G_Import___.h"
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
    call ESMF_UserCompSetInternalState ( GC, 'SS2G_GridComp', wrap, STATUS )
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

! !DESCRIPTION: This initializes SS' Grid Component. It primaryily fills 
!               GOCART's AERO states with its sea salt fields. 

! !REVISION HISTORY: 
! 24oct2019   E.Sherman  First attempt at refactoring

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
    type (SS2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, j, nbins, nCols, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: field_name, prefix
    real, pointer, dimension(:,:)        :: lats 
    real, pointer, dimension(:,:)        :: lons
    real                                 :: dummylon
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)

    logical                              :: data_driven
    integer                              :: NUM_BANDS

real, pointer :: ssptr(:,:,:,:)

    __Iam__('Initialize')

!****************************************************************************

!   Begin... 

!   Get the target components name and set-up traceback handle.
!   -----------------------------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, name=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) // '::' //trim(Iam)

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, MAPL, __RC__)

!   Get my internal private state
!   -----------------------------
    call ESMF_UserCompGetInternalState(GC, 'SS2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Get dimensions
!   ---------------
    call MAPL_GridGet (grid, globalCellCountPerDim=dims, __RC__ )
    km = dims(3)
    self%km = km

!   Scaling factor to multiply calculated
!   emissions by.  Applies to all size bins.
!   ----------------------------------------
    self%emission_scale = Chem_UtilResVal(dims(1), dims(2), self%emission_scale_res(:), __RC__)

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!  Load resource file and get number of bins 
!  -------------------------------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'SS2G_GridComp_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'SS2G_GridComp_'//trim(COMP_NAME)//'.rc does not exist! &
                                    loading SS2G_GridComp_SS.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'SS2G_GridComp_SS.rc', __RC__)
    end if

!   Call Generic Initialize 
!   ----------------------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get ( mapl, INTERNAL_ESMF_STATE = internal, &
                         LONS = LONS, &
                         LATS = LATS, __RC__ )

!   Is SS data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

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

!   Fill AERO State with sea salt fields
!   ----------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_ACI', aero_aci, __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

    call ESMF_StateGet (internal, 'SS', field, __RC__)
    fld = MAPL_FieldCreate (field, 'SS', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)
    call MAPL_StateAdd (aero_aci, fld, __RC__)

    ! ADD OTHER ATTRIBUTE, HENTRY COEFFICIENTS?
    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', VALUE=self%fscav(1), __RC__)

!   Dry deposition
!   ---------------
    call append_to_bundle('SSDP', providerState, prefix, Bundle_DP, __RC__)

!   Wet deposition (Convective scavenging)
!   --------------------------------------
    call append_to_bundle('SSSV', providerState, prefix, Bundle_DP, __RC__)

!   Wet deposition
!   ---------------
    call append_to_bundle('SSWT', providerState, prefix, Bundle_DP, __RC__)

!   Gravitational Settling
!   ----------------------
    call append_to_bundle('SSSD', providerState, prefix, Bundle_DP, __RC__)

!   Set AERO States' attributes
!   ----------------------------
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
    call ESMF_ConfigGetAttribute (cfg, self%diag_MieTable(instance)%channels, label= "r_channels:", __RC__)

    allocate (self%diag_MieTable(instance)%mie_aerosol, __STAT__)
    self%diag_MieTable(instance)%mie_aerosol = Chem_MieTableCreate (self%diag_MieTable(instance)%optics_file, __RC__ )
    call Chem_MieTableRead (self%diag_MieTable(instance)%mie_aerosol, self%diag_MieTable(instance)%nch, &
                            self%diag_MieTable(instance)%channels, rc, nmom=self%diag_MieTable(instance)%nmom)

    ! Mie Table instance/index
    call ESMF_AttributeSet(aero, name='mie_table_instance', value=instance, __RC__)

    ! Add variables to SS instance's aero state. This is used in aerosol optics calculations
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

    call ESMF_AttributeSet(aero, name='internal_varaible_name', value='SS', __RC__)

    call ESMF_MethodAdd(AERO, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)

!   Mask to prevent emissions from the Great Lakes and the Caspian Sea
!   ------------------------------------------------------------------
    allocate(self%deep_lakes_mask(ubound(lons, 1),ubound(lons, 2)), __STAT__)
    call deepLakesMask (lons, lats, real(MAPL_RADIANS_TO_DEGREES), self%deep_lakes_mask, rc)

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

!   Is SS data driven?
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
    type (SS2G_GridComp), pointer     :: self

    real, allocatable, dimension(:,:) :: fgridefficiency
    real, allocatable, dimension(:,:) :: fsstemis
    real, allocatable, dimension(:,:) :: tskin_c
    real, allocatable, dimension(:,:) :: fhoppel, vsettle
    real, allocatable, dimension(:,:) :: memissions, nemissions, dqa

    real    :: radius_wet, rhop_wet, diff_coef

    real(kind=DP), allocatable, dimension(:,:) :: gweibull

    integer :: n, i, j 

#include "SS2G_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

!if(mapl_am_i_root()) print*,'SS2G Run1 BEGIN'

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, grid=grid, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC (GC, mapl, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, __RC__)

#include "SS2G_GetPointer___.h"

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'SS2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Sea Salt Source (and modifications)
!   -----------------------------------
!   Grid box efficiency to emission (fraction of sea water)
    allocate(fgridefficiency, mold=frocean, __STAT__ )
    fgridefficiency = min(max(0.,(frocean-fraci)*self%deep_lakes_mask),1.)

!   Apply SST correction following Jaegle et al. 2011 if needed
!   ------------------------------------------------------------
    allocate(fsstemis, mold=frocean, __STAT__ )
    call jeagleSSTcorrection(self%sstEmisFlag, fsstemis, ts, __RC__)

!   Apply a Weibull distribution to emissions wind speeds
!   -----------------------------------------------------    
    allocate(gweibull(ubound(u10m,1), ubound(u10m,2)), __STAT__ )
    call weibullDistribution (gweibull, self%weibullFlag, u10m, v10m, __RC__)

!   Loop over bins and do emission calculation
!   Possibly apply the Hoppel correction based on fall speed (Fan and Toon, 2011)
!   -----------------------------------------------
    allocate(fhoppel, mold=frocean, __STAT__ )
    allocate(memissions, mold=frocean, __STAT__ )
    allocate(nemissions, mold=frocean, __STAT__ )
    allocate(dqa, mold=frocean, __STAT__ )

    fhoppel = 1.0

    do n = 1, self%nbins
       memissions = 0.
       nemissions = 0.
       dqa = 0.

       call SeasaltEmission (self%rlow(n), self%rup(n), self%emission_scheme, u10m, &
                             v10m, ustar, memissions, nemissions, __RC__ )

!   For the Hoppel correction need to compute the wet radius and settling velocity
!   in the surface
       if (self%hoppelFlag) then
          call hoppelCorrection (self%radius(n)*1.e-6, self%rhop(n), rh2(:,:,self%km), &
                                 dz, ustar, self%rhFlag, airdens(:,:,self%km), t(:,:,self%km), &
                                 chemGRAV, MAPL_KARMAN, fhoppel, __RC__)
       end if
 
       memissions = self%emission_scale * fgridefficiency * fsstemis * fhoppel * gweibull * memissions
       dqa = memissions * self%cdt * chemgrav / delp(:,:,self%km)
       SS(:,:,self%km,n) = SS(:,:,self%km,n) + dqa

       if (associated(SSEM)) then
          SSEM(:,:,n) = memissions
       end if
    end do !n = 1

    deallocate(fhoppel, memissions, nemissions, dqa, gweibull, &
               fsstemis, fgridefficiency, __STAT__)

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
    type (SS2G_GridComp), pointer     :: self

    integer                           :: n
    real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
    real                              :: fwet
    logical                           :: KIN

real, parameter ::  cpd    = 1004.16


#include "SS2G_DeclarePointer___.h"

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

#include "SS2G_GetPointer___.h"

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'SS2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    allocate(dqa, mold=lwi, __STAT__)
    allocate(drydepositionfrequency, mold=lwi, __STAT__)

!   Sea Salt Settling
!   -----------------
    call Chem_Settling2Gorig (self%km, self%rhFlag, SS, CHEMgrav, delp, &
                              self%radius*1.e-6, self%rhop, self%cdt, t, airdens, &
                              rh2, zle, SSSD, __RC__)

!   Sea Salt Deposition
!   -------------------
    drydepositionfrequency = 0.
    call DryDeposition(self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
                       MAPL_KARMAN, cpd, chemGRAV, z0h, drydepositionfrequency, __RC__ )

    ! increase deposition velocity over land
    where (abs(lwi - LAND) < 0.5)
        drydepositionfrequency = 5.0 * drydepositionfrequency
    end where

    do n = 1, self%nbins
       dqa = 0.
       dqa = max(0.0, SS(:,:,self%km,n)*(1.-exp(-drydepositionfrequency*self%cdt)))
       SS(:,:,self%km,n) = SS(:,:,self%km,n) - dqa
       if (associated(SSDP)) then
          SSDP(:,:,n) = dqa * delp(:,:,self%km) / chemGRAV / self%cdt
       end if
    end do

!   Sea Salt Large-scale Wet Removal
!   -------------------------------
    KIN = .TRUE.
    do n = 1, self%nbins
       fwet = 1.

       call WetRemovalGOCART2G(self%km, self%nbins, self%nbins, n, self%cdt, 'sea_salt', &
                               KIN, chemGRAV, fwet, SS(:,:,:,n), ple, t, airdens, &
                               pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, SSWT, rc)
    end do

    call Aero_Compute_Diags (self%diag_MieTable(self%instance), self%km, self%nbins, self%rlow, self%rup, &
                           self%diag_MieTable(self%instance)%channels, SS, chemGRAV, t, airdens, &
                           rh2, u, v, delp, SSSMASS, SSCMASS, SSMASS, SSEXTTAU, SSSCATAU,     &
                           SSSMASS25, SSCMASS25, SSMASS25, SSEXTT25, SSSCAT25, &
                           SSFLUXU, SSFLUXV, SSCONC, SSEXTCOEF, SSSCACOEF,    &
                           SSEXTTFM, SSSCATFM ,SSANGSTR, SSAERIDX, __RC__)


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
    type (MAPL_MetaComp), pointer      :: MAPL
    type (ESMF_State)                  :: aero

    integer                            :: i, n
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:)    :: ptr3d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp

    __Iam__('Run_data')

!*****************************************************************************
! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'//Iam

if (mapl_am_I_root()) print*,trim(comp_name),' Run_data BEGIN'

    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO', aero, __RC__)
    call ESMF_AttributeGet (aero, name='aerosol_names', itemCount=n, __RC__)

!   Update interal data pointers with ExtData
!   -----------------------------------------
    do i = 1, n
    write(field_name, '(A, I0.3)') 'ss', i
        call MAPL_GetPointer (internal, NAME=trim(field_name), ptr=ptr3d_int, __RC__)
        call MAPL_GetPointer (import,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

        ptr3d_int = ptr3d_imp
    end do

if (mapl_am_I_root()) print*,trim(comp_name),' Run_data END'

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
    real, dimension(:,:,:,:), pointer                :: q, q_4d
    integer, allocatable                             :: opaque_self(:)
    type(C_PTR)                                      :: address
    type(SS2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name, int_fld_name
    type(ESMF_Field)                                 :: fld
    character (len=ESMF_MAXSTR)                      :: COMP_NAME

    real(kind=DP), dimension(:,:,:), allocatable     :: ext_s, ssa_s, asy_s  ! (lon:,lat:,lev:)
    real, dimension(:,:,:), allocatable              :: x
    integer                                          :: instance
    integer                                          :: n, nbins, dims(4)
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band, offset
    integer, parameter                               :: n_bands = 1

    integer :: i, j, k

    __Iam__('SS2G::aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

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
             asy_s(i1:i2, j1:j2, km), &
                 x(i1:i2, j1:j2, km), __STAT__)

    call ESMF_AttributeGet(state, name='internal_varaible_name', value=int_fld_name, __RC__)
    call ESMF_StateGet (state, trim(int_fld_name), field=fld, __RC__) !add as attribute - dont hard code?
    call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

    nbins = size(q,4)

    allocate(q_4d(i1:i2, j1:j2, km, nbins), __STAT__)

    do n = 1, nbins
       do k = 1, km
          x(:,:,k) = ((PLE(:,:,k) - PLE(:,:,k-1))*0.01)*(100./MAPL_GRAV)
          q_4d(:,:,k,n) = x(:,:,k) * q(:,:,k,n)
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

    __Iam__('SS2G::aerosol_optics::mie_')

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



end module SS2G_GridCompMod

