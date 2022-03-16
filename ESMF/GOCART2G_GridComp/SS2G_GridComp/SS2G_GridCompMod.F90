#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: SS2G_GridCompMod - GOCART refactoring of the SS gridded component 

! !INTERFACE:
module SS2G_GridCompMod

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
   integer, parameter     :: DP=kind(1.0d0)

! !PUBLIC MEMBER FUNCTIONS:
   PUBLIC  SetServices

real, parameter ::  cpd    = 1004.16

! !DESCRIPTION: This module implements GOCART's Sea Salt (SS) Gridded Component.

! !REVISION HISTORY:
! 24Oct2019  E.Sherman  First attempt at refactoring.

!EOP
!===========================================================================
   integer, parameter         :: NHRES = 6

!  !Sea Salt state
   type, extends(GA_Environment) :: SS2G_GridComp
       real, allocatable      :: rlow(:)        ! particle effective radius lower bound [um]
       real, allocatable      :: rup(:)         ! particle effective radius upper bound [um]
       real, allocatable      :: rmed(:)        ! number median radius [um]
       integer                :: sstEmisFlag    ! Choice of SST correction to emissions: 
!                                                 0 - none; 1 - Jaegle et al. 2011; 2 - GEOS5
       logical                :: hoppelFlag     ! Apply the Hoppel correction to emissions (Fan and Toon, 2011)
       logical                :: weibullFlag    ! Apply the Weibull distribution to wind speed for emissions (Fan and Toon, 2011)
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

!   !Locals
    character (len=ESMF_MAXSTR)                 :: COMP_NAME
    type (ESMF_Config)                          :: cfg
    type (ESMF_Config)                          :: universal_cfg
    type (wrap_)                                :: wrap
    type (SS2G_GridComp), pointer               :: self

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
    call ESMF_ConfigLoadFile (cfg, 'SS2G_instance_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
       if (mapl_am_i_root()) print*,'SS2G_instance_'//trim(COMP_NAME)//'.rc does not exist! loading SS2G_instance_SS.rc instead'
       call ESMF_ConfigLoadFile (cfg, 'SS2G_instance_SS.rc', __RC__)
    end if

    ! process generic config items
    call self%GA_Environment%load_from_config( cfg, universal_cfg, __RC__)

    allocate(self%rlow(self%nbins), self%rup(self%nbins), self%rmed(self%nbins), __STAT__)

    ! process SS-specific items
!    call ESMF_ConfigGetAttribute (cfg, self%fscav,      label='fscav:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%sstEmisFlag, label='sstEmisFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%weibullFlag,  label='weibullFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%hoppelFlag, label='hoppelFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%emission_scheme, label='emission_scheme:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%emission_scale_res, label='emission_scale:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rlow, label='radius_lower:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rup, label='radius_upper:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rmed, label='particle_radius_number:', __RC__)

!   Is SS data driven?
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

   call MAPL_AddInternalSpec(gc,&
         short_name='SS', &
         long_name='Sea Salt Mixing Ratio all bins', &
         units='kg kg-1', &
         dims=MAPL_DimsHorzVert, &
         vlocation=MAPL_VlocationCenter, &
         restart=MAPL_RestartOptional, &
         ungridded_dims=[self%nbins], &
!         friendlyto='DYNAMICS:TURBULENCE:MOIST', &
         add2export=.true., __RC__)


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

        do i = 1, self%nbins
            write(field_name, '(A, I0.3)') '', i
            call MAPL_AddImportSpec(GC,                                           &
              SHORT_NAME = 'climss'//trim(field_name),                            &
              LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')',  &
              UNITS      = 'kg kg-1 s-1',                                             &
              RESTART    = MAPL_RestartSkip,                                      &
              DIMS       = MAPL_DimsHorzVert,                                     &
              VLOCATION  = MAPL_VLocationCenter, __RC__)

!           ! dry deposition
            call MAPL_AddImportSpec(GC,                                           &
              SHORT_NAME = 'climSSDP'//trim(field_name),                          &
              LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')',  &
              UNITS      = 'kg kg-1 s-1',                                             &
              DIMS       = MAPL_DimsHorzOnly,                                     &
              VLOCATION  = MAPL_VLocationCenter,                                  &
              RESTART    = MAPL_RestartSkip, __RC__)

!           ! wet deposition    
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSWT'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt wet removal (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1 s-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!           ! gravitational settling
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSSD'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1 s-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)

!        ! convective scavenging
            call MAPL_AddImportSpec(GC,                                           &
               SHORT_NAME = 'climSSSV'//trim(field_name),                         &
               LONG_NAME  = 'Sea Salt Mixing Ratio (bin '//trim(field_name)//')', &
               UNITS      = 'kg kg-1 s-1',                                            &
               DIMS       = MAPL_DimsHorzOnly,                                    &
               VLOCATION  = MAPL_VLocationCenter,                                 &
               RESTART    = MAPL_RestartSkip, __RC__)
        end do
    end if ! (data_driven)


!   Import, Export, Internal states for computational instance 
!   ----------------------------------------------------------
    if (.not. data_driven) then
#include "SS2G_Export___.h"
#include "SS2G_Import___.h"
#include "SS2G_Internal___.h"
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
    type (ESMF_State)                    :: aero
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg, universal_cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (SS2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: prefix, bin_index
    real, pointer, dimension(:,:)        :: lats 
    real, pointer, dimension(:,:)        :: lons
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)
    real, pointer, dimension(:,:,:,:)    :: int_ptr
    logical                              :: data_driven
    integer                              :: NUM_BANDS
    logical                              :: bands_are_present
    real, pointer, dimension(:,:,:)      :: ple
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
    call ESMF_ConfigLoadFile (cfg, 'SS2G_instance_'//trim(COMP_NAME)//'.rc', rc=status)
    if (status /= 0) then
      if (mapl_am_i_root()) print*,'SS2G_instance_'//trim(COMP_NAME)//'.rc does not exist! &
                                    loading SS2G_instance_SS.rc instead'
      call ESMF_ConfigLoadFile( cfg, 'SS2G_instance_SS.rc', __RC__)
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

!   Add attribute information for SS export. Used in NI hetergenous chemistry.
    call ESMF_StateGet (export, 'SS', field, __RC__)
    call ESMF_AttributeSet(field, NAME='radius', valueList=self%rmed, itemCount=self%nbins, __RC__)
    call ESMF_AttributeSet(field, NAME='fnum', valueList=self%fnum, itemCount=self%nbins, __RC__)

!   Fill AERO State with sea salt fields
!   ----------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

    call ESMF_StateGet (internal, 'SS', field, __RC__)
    call ESMF_AttributeSet(field, NAME='klid', value=self%klid, __RC__)
    fld = MAPL_FieldCreate (field, 'SS', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    if (.not. data_driven) then
!      Set klid
       call MAPL_GetPointer(import, ple, 'PLE', __RC__)
       call findKlid (self%klid, self%plid, ple, __RC__)
!      Set SS values to 0 where above klid
       call MAPL_GetPointer (internal, int_ptr, 'SS', __RC__)
       call setZeroKlid4d (self%km, self%klid, int_ptr)
    end if

    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', value=self%fscav(1), __RC__)

    if (data_driven) then
       instance = instanceData

       do i = 1, self%nbins
          write (bin_index, '(A, I0.3)') '', i
!         Dry deposition
          call append_to_bundle('SSDP'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Wet deposition (Convective scavenging)
          call append_to_bundle('SSSV'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Wet deposition
          call append_to_bundle('SSWT'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Gravitational Settling
          call append_to_bundle('SSSD'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)
       end do
    else
       instance = instanceComputational

!      Dry deposition
       call append_to_bundle('SSDP', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition (Convective scavenging)
       call append_to_bundle('SSSV', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition
       call append_to_bundle('SSWT', providerState, prefix, Bundle_DP, __RC__)

!      Gravitational Settling
       call append_to_bundle('SSSD', providerState, prefix, Bundle_DP, __RC__)
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
    call add_aero (aero, label='monochromatic_extinction_in_air_due_to_ambient_aerosol', &
                   label2='monochromatic_EXT', grid=grid, typekind=MAPL_R4,__RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol', label2='aerosolSum', grid=grid, typekind=MAPL_R4, __RC__)

    call ESMF_AttributeSet (aero, name='band_for_aerosol_optics', value=0, __RC__)
    call ESMF_AttributeSet (aero, name='wavelength_for_aerosol_optics', value=0., __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet (aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

    call ESMF_AttributeSet (aero, name='internal_variable_name', value='SS', __RC__)

    call ESMF_MethodAdd (aero, label='aerosol_optics', userRoutine=aerosol_optics, __RC__)
    call ESMF_MethodAdd (aero, label='monochromatic_aerosol_optics', userRoutine=monochromatic_aerosol_optics, __RC__)
    call ESMF_MethodAdd (aero, label='get_mixR', userRoutine=get_mixR, __RC__)

!   Mask to prevent emissions from the Great Lakes and the Caspian Sea
!   ------------------------------------------------------------------
    allocate(self%deep_lakes_mask(ubound(lons, 1),ubound(lons, 2)), __STAT__)
    call deepLakesMask (lons, lats, real(MAPL_RADIANS_TO_DEGREES), self%deep_lakes_mask, __RC__)

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
    real, allocatable, dimension(:,:) :: fhoppel
    real, allocatable, dimension(:,:) :: memissions, nemissions, dqa

    real(kind=DP), allocatable, dimension(:,:) :: gweibull

    integer :: n 

#include "SS2G_DeclarePointer___.h"

   __Iam__('Run1')

!*****************************************************************************
!   Begin... 

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
                             v10m, ustar, MAPL_PI, memissions, nemissions, __RC__ )

!      For the Hoppel correction need to compute the wet radius and settling velocity
!      in the surface
       if (self%hoppelFlag) then
          call hoppelCorrection (self%radius(n)*1.e-6, self%rhop(n), rh2(:,:,self%km), &
                                 dz, ustar, self%rhFlag, airdens(:,:,self%km), t(:,:,self%km), &
                                 MAPL_GRAV, MAPL_KARMAN, fhoppel, __RC__)
       end if
 
       memissions = self%emission_scale * fgridefficiency * fsstemis * fhoppel * gweibull * memissions
       dqa = memissions * self%cdt * MAPL_GRAV / delp(:,:,self%km)
       SS(:,:,self%km,n) = SS(:,:,self%km,n) + dqa

       if (associated(SSEM)) then
          SSEM(:,:,n) = memissions
       end if
    end do !n = 1

    deallocate(fhoppel, memissions, nemissions, dqa, gweibull, &
               fsstemis, fgridefficiency, __STAT__)

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
    type (SS2G_GridComp), pointer     :: self

    integer                           :: n
    real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
    real                              :: fwet
    logical                           :: KIN


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
    do n = 1, self%nbins
       call Chem_Settling (self%km, self%klid, n, self%rhFlag, self%cdt, MAPL_GRAV, &
                           self%radius(n)*1.e-6, self%rhop(n), SS(:,:,:,n), t, airdens, &
                           rh2, zle, delp, SSSD, __RC__)
    end do

!   Deposition
!   -----------
    drydepositionfrequency = 0.
    call DryDeposition(self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
                       MAPL_KARMAN, cpd, MAPL_GRAV, z0h, drydepositionfrequency, __RC__ )

    ! increase deposition velocity over land
    where (abs(lwi - LAND) < 0.5)
        drydepositionfrequency = 5.0 * drydepositionfrequency
    end where

    do n = 1, self%nbins
       dqa = 0.
       dqa = max(0.0, SS(:,:,self%km,n)*(1.-exp(-drydepositionfrequency*self%cdt)))
       SS(:,:,self%km,n) = SS(:,:,self%km,n) - dqa
       if (associated(SSDP)) then
          SSDP(:,:,n) = dqa * delp(:,:,self%km) / MAPL_GRAV / self%cdt
       end if
    end do

!   Large-scale Wet Removal
!   ------------------------
    KIN = .TRUE.
    do n = 1, self%nbins
       fwet = 1.
       call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, n, self%cdt, 'sea_salt', &
                               KIN, MAPL_GRAV, fwet, SS(:,:,:,n), ple, t, airdens, &
                               pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, SSWT, __RC__)
    end do

!   Compute diagnostics
!   -------------------
!   Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
    call Aero_Compute_Diags (self%diag_Mie, self%km, self%klid, 1, self%nbins, self%rlow, &
                             self%rup, self%wavelengths_profile*1.0e-9, &
                             self%wavelengths_vertint*1.0e-9, SS, MAPL_GRAV, t, airdens,rh2, u, v, &
                             delp, ple, tropp,SSSMASS, SSCMASS, SSMASS, SSEXTTAU,SSSTEXTTAU, SSSCATAU,SSSTSCATAU, &
                             SSSMASS25, SSCMASS25, SSMASS25, SSEXTT25, SSSCAT25, &
                             SSFLUXU, SSFLUXV, SSCONC, SSEXTCOEF, SSSCACOEF,    &
                             SSEXTTFM, SSSCATFM ,SSANGSTR, SSAERIDX, NO3nFlag=.false.,__RC__)
 
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
    type (SS2G_GridComp), pointer      :: self

    integer                            :: i
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:,:)  :: ptr4d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp

    __Iam__('Run_data')

!*****************************************************************************
! Begin... 

! Get my name and set-up traceback handle
! ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'//Iam

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'SS2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Update interal data pointers with ExtData
!   -----------------------------------------
    call MAPL_GetPointer (internal, NAME='SS', ptr=ptr4d_int, __RC__)

    do i = 1, self%nbins
    write(field_name, '(A, I0.3)') 'ss', i
        call MAPL_GetPointer (import,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

        ptr4d_int(:,:,:,i) = ptr3d_imp
    end do

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

    real(kind=DP), dimension(:,:,:), allocatable     :: ext_s, ssa_s, asy_s  ! (lon:,lat:,lev:)
    real, dimension(:,:,:), allocatable              :: x
    integer                                          :: instance
    integer                                          :: n, nbins
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band
    integer ::  k

    __Iam__('SS2G::aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

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
             asy_s(i1:i2, j1:j2, km), &
                 x(i1:i2, j1:j2, km), __STAT__)

    call ESMF_AttributeGet(state, name='internal_variable_name', value=int_fld_name, __RC__)
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

    __Iam__('SS2G::aerosol_optics::mie_')

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

!-------------------------------------------------------------------------------------
  subroutine monochromatic_aerosol_optics(state, rc)

    implicit none

!   !ARGUMENTS:
    type (ESMF_State)                                :: state
    integer,            intent(out)                  :: rc

!   !Local
    real, dimension(:,:,:), pointer                  :: ple, rh
    real, dimension(:,:), pointer                    :: var
    real, dimension(:,:,:,:), pointer                :: q, q_4d
    integer, allocatable                             :: opaque_self(:)
    type(C_PTR)                                      :: address
    type(SS2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld

    real, dimension(:,:,:), allocatable              :: tau_s, tau, x ! (lon:,lat:,lev:)
    integer                                          :: instance
    integer                                          :: n, nbins, k
    integer                                          :: i1, j1, i2, j2, km, i, j
    real                                             :: wavelength

    __Iam__('SS2G::monochromatic_aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet (state, name='mie_table_instance', value=instance, __RC__)

!   Radiation band
!   --------------
    wavelength = 0.
    call ESMF_AttributeGet (state, name='wavelength_for_aerosol_optics', value=wavelength, __RC__)

!   Pressure at layer edges 
!   ------------------------
    call ESMF_AttributeGet (state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer (state, ple, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, ple, 'PLE', __RC__)

    i1 = lbound(ple, 1); i2 = ubound(ple, 1)
    j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                         km = ubound(ple, 3)

!   Relative humidity
!   -----------------
    call ESMF_AttributeGet (state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer (state, rh, trim(fld_name), __RC__)

!    call MAPL_GetPointer (state, rh, 'RH2', __RC__)

    allocate(tau_s(i1:i2, j1:j2, km), &
               tau(i1:i2, j1:j2, km), &
                 x(i1:i2, j1:j2, km), __STAT__)
    tau_s = 0.
    tau   = 0.

    call ESMF_StateGet (state, 'SS', field=fld, __RC__)
    call ESMF_FieldGet (fld, farrayPtr=q, __RC__)

    nbins = size(q,4)

    allocate(q_4d(i1:i2, j1:j2, km, nbins), __STAT__)
    q_4d = 0.

    do n = 1, nbins
       do k = 1, km
          x(:,:,k) = (PLE(:,:,k) - PLE(:,:,k-1)) / MAPL_GRAV
          q_4d(:,:,k,n) = x(:,:,k) * q(:,:,k,n)
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

    call ESMF_AttributeGet (state, name='monochromatic_extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = sum(tau_s, dim=3)
    end if

    deallocate(q_4d, __STAT__)

    RETURN_(ESMF_SUCCESS)

  end subroutine monochromatic_aerosol_optics


end module SS2G_GridCompMod

