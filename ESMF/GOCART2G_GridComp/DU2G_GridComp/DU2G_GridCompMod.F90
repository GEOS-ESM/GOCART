#include "MAPL_Generic.h"

!=============================================================================
!BOP

! !MODULE: DU2G_GridCompMod - GOCART refactoring of the DU gridded component 

! !INTERFACE:
module DU2G_GridCompMod

!  !USES:
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

! !PUBLIC MEMBER FUNCTIONS:
   public  SetServices

! !DESCRIPTION: This module implements GOCART's Dust (DU) Gridded Component.

! !REVISION HISTORY:
! 16Oct2019  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring

!EOP
!===========================================================================

   integer, parameter         :: NHRES = 6

!  !Dust state
   type, extends(GA_Environment) :: DU2G_GridComp
       real, allocatable      :: rlow(:)        ! particle effective radius lower bound [um]
       real, allocatable      :: rup(:)         ! particle effective radius upper bound [um]
       real, allocatable      :: sfrac(:)       ! fraction of total source
       real, allocatable      :: sdist(:)       ! FENGSHA aerosol fractional size distribution [1]
       real                   :: alpha          ! FENGSHA scaling factor
       real                   :: gamma          ! FENGSHA tuning exponent
       real                   :: kvhmax         ! FENGSHA max. vertical/horizontal mass flux ratio [1]
       real                   :: Ch_DU_res(NHRES) ! resolutions used for Ch_DU
       real                   :: Ch_DU          ! dust emission tuning coefficient [kg s2 m-5].
       logical                :: maringFlag=.false.  ! maring settling velocity correction
       integer                :: day_save = -1      
       character(len=:), allocatable :: emission_scheme     ! emission scheme selector
       integer       :: clayFlag       ! clay and silt term in K14
       real          :: f_swc          ! soil mosture scaling factor
       real          :: f_scl          ! clay content scaling factor
       real          :: uts_gamma      ! threshold friction velocity parameter 'gamma'
!      !Workspae for point emissions
       logical                :: doing_point_emissions = .false.
       character(len=255)     :: point_emissions_srcfilen   ! filename for pointwise emissions
       integer                         :: nPts = -1
       integer, allocatable, dimension(:)  :: pstart, pend
       real, allocatable, dimension(:)     :: pLat, &
                                              pLon, &
                                              pBase, &
                                              pTop, &
                                              pEmis
   end type DU2G_GridComp

   type wrap_
      type (DU2G_GridComp), pointer     :: PTR => null()
   end type wrap_

contains

!============================================================================
!BOP

! !IROUTINE: SetServices 

! !INTERFACE:
  subroutine SetServices (GC, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(INOUT)   :: GC  ! gridded component
    integer,              intent(  OUT)   :: RC  ! return code

!   !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
!     the Initialize and Finalize services to generic versions. It also
!     allocates our instance of a generic state and puts it in the 
!     gridded component (GC). Here we only set the two-stage run method
!     and declare the data services.

!   !REVISION HISTORY: 
!   16oct2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring

!EOP
!============================================================================

!
!   !Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (ESMF_Config)                 :: cfg
    type (ESMF_Config)                 :: universal_cfg
    type (wrap_)                       :: wrap
    type (DU2G_GridComp), pointer      :: self

    character (len=ESMF_MAXSTR)        :: field_name
    character (len=ESMF_MAXSTR)        :: emission_scheme
    integer                            :: i
    real                               :: DEFVAL
    logical                            :: data_driven = .true.

    __Iam__('SetServices')

!****************************************************************************
!   Begin...

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, config=universal_cfg, __RC__)
    Iam = trim(COMP_NAME) //'::'// Iam

!   Wrap internal state for storing in GC
!   -------------------------------------
    allocate (self, __STAT__)
    wrap%ptr => self

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'DU2G_instance_'//trim(COMP_NAME)//'.rc', rc=status)

    if (status /= 0) then
        if (mapl_am_i_root()) print*,'DU2G_instance_'//trim(COMP_NAME)//'.rc does not exist! Loading DU2G_GridComp_DU.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'DU2G_instance_DU.rc', __RC__)
    end if

    ! process generic config items
    call self%GA_Environment%load_from_config(cfg, universal_cfg, __RC__)

    allocate(self%sfrac(self%nbins), self%rlow(self%nbins), self%rup(self%nbins), __STAT__)
    ! process DU-specific items
    call ESMF_ConfigGetAttribute (cfg, self%maringFlag, label='maringFlag:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%sfrac,      label='source_fraction:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%Ch_DU_res,  label='Ch_DU:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rlow,       label='radius_lower:', __RC__)
    call ESMF_ConfigGetAttribute (cfg, self%rup,        label='radius_upper:', __RC__)

    call ESMF_ConfigGetAttribute (cfg, emission_scheme, label='emission_scheme:', default='ginoux', __RC__)
    self%emission_scheme = ESMF_UtilStringLowerCase(trim(emission_scheme), __RC__)

    ! Test if our scheme is allowed, if so, print it out
    _ASSERT(any(self%emission_scheme == [character(len=7) :: 'ginoux','k14','fengsha']), "Error. Unallowed emission scheme: "//trim(self%emission_scheme)//". Allowed: ginoux, k14, fengsha")
    if (MAPL_AM_I_ROOT()) then
       write (*,*) trim(Iam)//": Dust emission scheme is "//trim(self%emission_scheme)
    end if

    call ESMF_ConfigGetAttribute (cfg, self%point_emissions_srcfilen, &
                                  label='point_emissions_srcfilen:', default='/dev/null', __RC__)
    if ( (index(self%point_emissions_srcfilen,'/dev/null')>0) ) then
       self%doing_point_emissions = .false. ! disable it if no file specified
    else
       self%doing_point_emissions = .true.  ! we are good to go
    end if

!   read scheme-specific parameters
!   --------------------------------
    select case (self%emission_scheme)
    case ('fengsha')
       call ESMF_ConfigGetAttribute (cfg, self%alpha,      label='alpha:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%gamma,      label='gamma:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%kvhmax,     label='vertical_to_horizontal_flux_ratio_limit:', __RC__)
    case ('k14')
       call ESMF_ConfigGetAttribute (cfg, self%clayFlag,   label='clayFlag:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%f_swc,      label='soil_moisture_factor:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%f_scl,      label='soil_clay_factor:', __RC__)
       call ESMF_ConfigGetAttribute (cfg, self%uts_gamma,  label='uts_gamma:', __RC__)
    case ('ginoux')
       ! nothing to do
    case default
       _ASSERT_RC(.false., "Unallowed emission scheme: "//trim(self%emission_scheme)//". Allowed: ginoux, k14, fengsha", ESMF_RC_NOT_IMPL)
    end select

!   Is DU data driven?
!   ------------------
    call determine_data_driven (COMP_NAME, data_driven, __RC__)

!   Set entry points
!   ------------------------
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Initialize,  Initialize, __RC__)
    call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run, __RC__)
    if (data_driven .neqv. .true.) then
       call MAPL_GridCompSetEntryPoint (GC, ESMF_Method_Run, Run2, __RC__)
    end if

    DEFVAL = 0.0


!   IMPORT STATE
!   -------------
    if (data_driven) then

       call MAPL_AddInternalSpec(gc,&
          short_name='DU', &
          long_name='Dust Mixing Ratio all bins', &
          units='kg kg-1', &
          dims=MAPL_DimsHorzVert, &
          vlocation=MAPL_VlocationCenter, &
          restart=MAPL_RestartOptional, &
          ungridded_dims=[self%nbins], &
          friendlyto='DYNAMICS:TURBULENCE:MOIST', &
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
          write (field_name, '(A, I0.3)') '', i
          call MAPL_AddImportSpec(GC,                                     &
             short_name = 'climdu'//trim(field_name),                        &
             long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')',  &
             units      = 'kg kg-1',                                         &
             restart    = MAPL_RestartSkip,                                  &
             dims       = MAPL_DimsHorzVert,                                 &
             vlocation  = MAPL_VLocationCenter, __RC__)

!        !dry deposition
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climDUDP'//trim(field_name),                      &
             long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')',  &
             units      = 'kg kg-1',                                         &
             dims       = MAPL_DimsHorzOnly,                                 &
             vlocation  = MAPL_VLocationCenter,                              &
             restart    = MAPL_RestartSkip, __RC__)

!        !wet deposition    
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climDUWT'//trim(field_name),                     &
             long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)

!        !gravitational settling
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climDUSD'//trim(field_name),                     &
             long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)

!        !convective scavenging
          call MAPL_AddImportSpec(GC,                                       &
             short_name = 'climDUSV'//trim(field_name),                     &
             long_name  = 'Dust Mixing Ratio (bin '//trim(field_name)//')', &
             units      = 'kg kg-1',                                        &
             dims       = MAPL_DimsHorzOnly,                                &
             vlocation  = MAPL_VLocationCenter,                             &
             restart    = MAPL_RestartSkip, __RC__)
        end do
    end if ! (data_driven)

!   Computational Instance
!   ----------------------
    if (.not. data_driven) then
#include "DU2G_Export___.h"
      associate (scheme => self%emission_scheme)
#include "DU2G_Import___.h"
      end associate
#include "DU2G_Internal___.h"
    end if

!   This state holds fields needed by radiation
!   ---------------------------------------------
    call MAPL_AddExportSpec (GC,                             &
       short_name = trim(COMP_NAME)//'_AERO',                &
       long_name  = 'aerosols_from_'//trim(COMP_NAME),       &
       units      = 'kg kg-1',                               &
       dims       = MAPL_DimsHorzVert,                       &
       vlocation  = MAPL_VLocationCenter,                    &
       datatype   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo modification
!   by aerosol settling and deposition
!   ~~~DEVELOPERS NOTE~~~ Change to StateItem when possible
!   ---------------------------------------------------------------
    call MAPL_AddExportSpec (GC,                                  &
       short_name = trim(COMP_NAME)//'_AERO_DP',                                    &
       long_name  = 'aerosol_deposition_from_'//trim(COMP_NAME),  &
       units      = 'kg m-2 s-1',                                 &
       dims       = MAPL_DimsHorzOnly,                            &
       datatype   = MAPL_BundleItem, __RC__)

!   Store internal state in GC
!   --------------------------
    call ESMF_UserCompSetInternalState ( GC, 'DU2G_GridComp', wrap, STATUS )
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
  subroutine Initialize (GC, import, export, clock, RC)

!   !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC     ! Gridded component 
    type (ESMF_State),    intent(inout) :: import ! Import state
    type (ESMF_State),    intent(inout) :: export ! Export state
    type (ESMF_Clock),    intent(inout) :: clock  ! The clock
    integer, optional,    intent(  out) :: RC     ! Error code

! !DESCRIPTION: This initializes DU's Grid Component. It primaryily fills 
!               GOCART's AERO states with its dust fields. 

! !REVISION HISTORY: 
! 16oct2019  E.Sherman, A.da Silva, T.Clune, A.Darmenov - First attempt at refactoring

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)          :: COMP_NAME
    type (MAPL_MetaComp),       pointer  :: MAPL
    type (ESMF_Grid)                     :: grid
    type (ESMF_State)                    :: internal
    type (ESMF_State)                    :: aero
    type (ESMF_State)                    :: providerState
    type (ESMF_Config)                   :: cfg
    type (ESMF_Config)                   :: universal_cfg
    type (ESMF_FieldBundle)              :: Bundle_DP
    type (wrap_)                         :: wrap
    type (DU2G_GridComp), pointer        :: self

    integer, allocatable                 :: mieTable_pointer(:)
    integer                              :: i, dims(3), km
    integer                              :: instance
    type (ESMF_Field)                    :: field, fld
    character (len=ESMF_MAXSTR)          :: bin_index, prefix
    real                                 :: CDT         ! chemistry timestep (secs)
    integer                              :: HDT         ! model     timestep (secs)
    real, pointer, dimension(:,:,:,:)    :: int_ptr
    real, pointer, dimension(:,:,:)      :: ple
    logical                              :: data_driven
    integer                              :: NUM_BANDS
    logical                              :: bands_are_present
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
    call ESMF_UserCompGetInternalState(GC, 'DU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    call MAPL_GridGet ( grid, globalCellCountPerDim=dims, __RC__ )

!   Dust emission tuning coefficient [kg s2 m-5]. NOT bin specific.
!   ---------------------------------------------------------------
    self%Ch_DU = Chem_UtilResVal(dims(1), dims(2), self%Ch_DU_res(:), __RC__)
    self%Ch_DU = self%Ch_DU * 1.0e-9

!   Dust emission size distribution for FENGSHA
!   ---------------------------------------------------------------
    if (self%emission_scheme == 'fengsha') then
      allocate(self%sdist(self%nbins), __STAT__)
      call DustAerosolDistributionKok(self%radius, self%rup, self%rlow, self%sdist)
    end if

!   Get dimensions
!   ---------------
    km = dims(3)
    self%km = km

!   Get DTs
!   -------
    call MAPL_GetResource(mapl, HDT, Label='RUN_DT:', __RC__)                        
    call MAPL_GetResource(mapl, CDT, Label='GOCART_DT:', default=real(HDT), __RC__)
    self%CDT = CDT

!   Load resource file  
!   -------------------
    cfg = ESMF_ConfigCreate (__RC__)
    call ESMF_ConfigLoadFile (cfg, 'DU2G_instance_'//trim(COMP_NAME)//'.rc', RC=STATUS)
    if (status /= 0) then
        if (mapl_am_i_root()) print*,'DU2G_instance_'//trim(COMP_NAME)//'.rc does not exist! &
                                      loading DU2G_instance_DU.rc instead'
        call ESMF_ConfigLoadFile (cfg, 'DU2G_instance_DU.rc', __RC__)
    end if

!   Call Generic Initialize 
!   ------------------------
    call MAPL_GenericInitialize (GC, import, export, clock, __RC__)

!   Get parameters from generic state.
!   -----------------------------------
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Is DU data driven?
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

!   Add attribute information for DU export. Used in NI hetergenous chemistry.
    call ESMF_StateGet (export, 'DU', field, __RC__)
    call ESMF_AttributeSet(field, NAME='radius', valueList=self%radius, itemCount=self%nbins, __RC__)
    call ESMF_AttributeSet(field, NAME='fnum', valueList=self%fnum, itemCount=self%nbins, __RC__)

!   Add attribute information to internal state variables
!   -----------------------------------------------------
!   Fill AERO States with dust fields
!   ------------------------------------
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO'    , aero    , __RC__)
    call ESMF_StateGet (export, trim(COMP_NAME)//'_AERO_DP' , Bundle_DP, __RC__)

    call ESMF_StateGet (internal, 'DU', field, __RC__)
    fld = MAPL_FieldCreate (field, 'DU', __RC__)
    call MAPL_StateAdd (aero, fld, __RC__)

    if (.not. data_driven) then
!      Set klid
       call MAPL_GetPointer(import, ple, 'PLE', __RC__)
       call findKlid (self%klid, self%plid, ple, __RC__)
!      Set internal DU values to 0 where above klid
       call MAPL_GetPointer (internal, int_ptr, 'DU', __RC__)
       call setZeroKlid4d (self%km, self%klid, int_ptr)
    end if

    call ESMF_AttributeSet(field, NAME='ScavengingFractionPerKm', value=self%fscav(1), __RC__)

    if (data_driven) then
       instance = instanceData

       do i = 1, self%nbins
          write (bin_index, '(A, I0.3)') '', i
!         Dry deposition
          call append_to_bundle('DUDP'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Wet deposition (Convective scavenging)
          call append_to_bundle('DUSV'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Wet deposition
          call append_to_bundle('DUWT'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)

!         Gravitational Settling
          call append_to_bundle('DUSD'//trim(bin_index), providerState, prefix, Bundle_DP, __RC__)
       end do
    else
       instance = instanceComputational

!      Dry deposition
       call append_to_bundle('DUDP', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition (Convective scavenging)
       call append_to_bundle('DUSV', providerState, prefix, Bundle_DP, __RC__)

!      Wet deposition
       call append_to_bundle('DUWT', providerState, prefix, Bundle_DP, __RC__)

!      Gravitational Settling
       call append_to_bundle('DUSD', providerState, prefix, Bundle_DP, __RC__)
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

!   Mie Table instance/index
    call ESMF_AttributeSet (aero, name='mie_table_instance', value=instance, __RC__)

!   Add variables to DU instance's aero state. This is used in aerosol optics calculations
!   --------------------------------------------------------------------------------------
    call add_aero (aero, label='air_pressure_for_aerosol_optics',      label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='relative_humidity_for_aerosol_optics', label2='RH',  grid=grid, typekind=MAPL_R4, __RC__)
!   call ESMF_StateGet (import, 'PLE', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)
!   call ESMF_StateGet (import, 'RH2', field, __RC__)
!   call MAPL_StateAdd (aero, field, __RC__)

    call add_aero (aero, label='extinction_in_air_due_to_ambient_aerosol',    label2='EXT', grid=grid, typekind=MAPL_R8, __RC__)
    call add_aero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R8, __RC__)
    call add_aero (aero, label='asymmetry_parameter_of_ambient_aerosol',      label2='ASY', grid=grid, typekind=MAPL_R8, __RC__)
    call add_aero (aero, label='monochromatic_extinction_in_air_due_to_ambient_aerosol', label2='monochromatic_EXT', &
                   grid=grid, typekind=MAPL_R4, __RC__)
    call add_aero (aero, label='sum_of_internalState_aerosol', label2='aerosolSum', grid=grid, typekind=MAPL_R4, __RC__)

    call ESMF_AttributeSet (aero, name='band_for_aerosol_optics', value=0, __RC__)
    call ESMF_AttributeSet (aero, name='wavelength_for_aerosol_optics', value=0., __RC__)

    mieTable_pointer = transfer(c_loc(self), [1])
    call ESMF_AttributeSet (aero, name='mieTable_pointer', valueList=mieTable_pointer, itemCount=size(mieTable_pointer), __RC__)

    call ESMF_AttributeSet (aero, name='internal_variable_name', value='DU', __RC__)

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

!   !DESCRIPTION: Run method for the Dust Grid Component. Determines whether to
!                 run data or computational run method.

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

!   Is DU data driven?
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

! !DESCRIPTION:  Computes emissions/sources for Dust

!EOP
!============================================================================
!   !Locals
    character (len=ESMF_MAXSTR)       :: COMP_NAME
    type (MAPL_MetaComp), pointer     :: mapl
    type (ESMF_State)                 :: internal
    type (ESMF_Grid)                  :: grid
    type (wrap_)                      :: wrap
    type (DU2G_GridComp), pointer     :: self
    type(ESMF_Time)                   :: time

    integer  :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
    integer  :: import_shape(2), i2, j2
!    real, dimension(:,:), pointer     :: emissions_surface
    real, dimension(:,:,:), allocatable     :: emissions_surface
    real, dimension(:,:,:,:), allocatable :: emissions
    real, dimension(:,:,:), allocatable   :: emissions_point
    character (len=ESMF_MAXSTR)  :: fname ! file name for point source emissions
    integer, pointer, dimension(:)  :: iPoint, jPoint
    logical :: fileExists
    real :: qmax, qmin
    integer :: n, ijl
    real, dimension(:,:), allocatable   :: z_
    real, dimension(:,:), allocatable   :: ustar_
    real, dimension(:,:), allocatable   :: ustar_t_
    real, dimension(:,:), allocatable   :: ustar_ts_
    real, dimension(:,:), allocatable   :: R_
    real, dimension(:,:), allocatable   :: H_w_
    real, dimension(:,:), allocatable   :: f_erod_

#include "DU2G_DeclarePointer___.h"

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

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'DU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Extract nymd(yyyymmdd) from clock
!   ---------------------------------
    call ESMF_ClockGet (clock, currTime=time, __RC__)
    call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
    call MAPL_PackTime (nymd, iyr, imm , idd)
    call MAPL_PackTime (nhms, ihr, imn, isc)

    associate (scheme => self%emission_scheme)
#include "DU2G_GetPointer___.h"
    end associate

!   Set du_src to 0 where undefined
!   --------------------------------
    if (associated(du_src)) then
       where (1.01*du_src > MAPL_UNDEF) du_src = 0.
    endif
!   Get dimensions
!   ---------------
    import_shape = shape(wet1)
    i2 = import_shape(1)
    j2 = import_shape(2)
    ijl  = ( i2 - 1 + 1 ) * ( j2 - 1 + 1 )

    allocate(emissions(i2,j2,self%km,self%nbins),  __STAT__)
    emissions = 0.0
    allocate(emissions_point, mold=delp,  __STAT__)
    emissions_point = 0.0
    allocate(emissions_surface(i2,j2,self%nbins), __STAT__)
    emissions_surface = 0.0

!   Get surface gridded emissions
!   -----------------------------
    select case (self%emission_scheme)

    case ('k14')
       allocate(ustar_, mold=U10M,    __STAT__)
       allocate(ustar_t_, mold=U10M,  __STAT__)
       allocate(ustar_ts_, mold=U10M, __STAT__)
       allocate(R_, mold=U10M,        __STAT__)
       allocate(H_w_, mold=U10M,      __STAT__)
       allocate(f_erod_, mold=U10M,   __STAT__)
       allocate(z_, mold=U10M,        __STAT__)

       z_ = 10.0 ! wind is at 10m

       call DustEmissionK14( self%km, tsoil1, wcsf, rhos,        &
                             du_z0, z_, u10n, v10n, ustar,    &
                             frland, asnow,               &
                             du_src,                       &
                             du_sand, du_silt, du_clay,             &
                             du_texture, du_veg, du_gvf,     &
                             self%f_swc, self%f_scl, self%uts_gamma, &
                             MAPL_UNDEF, MAPL_GRAV, MAPL_KARMAN,   &
                             self%clayFlag, self%Ch_DU/1.e-9,  &
                             emissions_surface,            &
                             ustar_,                       &
                             ustar_t_,                     &
                             ustar_ts_,                    &
                             R_, H_w_, f_erod_,            &
                             __RC__ )


       if (associated(DU_UST)) DU_UST = ustar_
       if (associated(DU_UST_T)) DU_UST_T = ustar_t_
       if (associated(DU_UST_T)) DU_UST_T = ustar_ts_
       if (associated(DU_DPC)) DU_DPC = R_
       if (associated(DU_SMC)) DU_SMC = H_w_
       if (associated(DU_EROD)) DU_EROD = f_erod_

    case ('fengsha')
       call DustEmissionFENGSHA (frlake, frsnow, lwi, slc, du_clay, du_sand, du_silt,       &
                                 du_ssm, du_rdrag, airdens(:,:,self%km), ustar, du_uthres,  &
                                 self%alpha, self%gamma, self%kvhmax, MAPL_GRAV,   &
                                 self%rhop, self%sdist, emissions_surface, __RC__)
    case ('ginoux')

       call DustEmissionGOCART2G(self%radius*1.e-6, frlake, wet1, lwi, u10m, v10m, &
                                 self%Ch_DU, du_src, MAPL_GRAV, &
                                 emissions_surface, __RC__)
    case default
       _ASSERT_RC(.false.,'missing dust emission scheme. Allowed: ginoux, fengsha, k14',ESMF_RC_NOT_IMPL)
    end select

!   Read point emissions file once per day
!   --------------------------------------
    if (self%doing_point_emissions) then
       if (self%day_save /= idd) then
          self%day_save = idd
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
       end if
    end if

!   Get indices for point emissions
!   -------------------------------
    if (self%nPts > 0) then
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
    end if

!   Update aerosol state
!   --------------------
    call UpdateAerosolState (emissions, emissions_surface, emissions_point, &
                             self%sfrac, self%nPts, self%km, self%CDT, MAPL_GRAV, &
                             self%nbins, delp, DU, __RC__)

    if (associated(DUEM)) then
       DUEM = sum(emissions, dim=3)
    end if

!   Clean up
!   --------
    deallocate(emissions, emissions_surface, emissions_point, __STAT__)
    if (self%nPts > 0) then
        deallocate(iPoint, jPoint, __STAT__)
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
    type (DU2G_GridComp), pointer     :: self

    integer                           :: n
    real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
    real                              :: fwet
    logical                           :: KIN

    real, parameter ::  cpd    = 1004.16

#include "DU2G_DeclarePointer___.h"

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
    call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'DU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

    associate (scheme => self%emission_scheme)
#include "DU2G_GetPointer___.h"
    end associate

    allocate(dqa, mold=wet1, __STAT__)
    allocate(drydepositionfrequency, mold=wet1, __STAT__)

!   Dust Settling
!   -------------
    do n = 1, self%nbins
       call Chem_Settling (self%km, self%klid, n, self%rhFlag, self%cdt, MAPL_GRAV, &
                           self%radius(n)*1.e-6, self%rhop(n), DU(:,:,:,n), t, airdens, &
                           rh2, zle, delp, DUSD, correctionMaring=self%maringFlag, __RC__)


    end do

!   Dust Deposition
!   ----------------
   do n = 1, self%nbins
      drydepositionfrequency = 0.
      call DryDeposition(self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
                         MAPL_KARMAN, cpd, MAPL_GRAV, z0h, drydepositionfrequency, status, &
                         self%radius(n)*1.e-6, self%rhop(n), u10m, v10m, frlake, wet1)
      VERIFY_(status)

      dqa = 0.
      dqa = max(0.0, DU(:,:,self%km,n)*(1.-exp(-drydepositionfrequency*self%cdt)))
      DU(:,:,self%km,n) = DU(:,:,self%km,n) - dqa

    if (associated(DUDP)) then
       DUDP(:,:,n) = dqa*delp(:,:,self%km)/MAPL_GRAV/self%cdt
    end if
   end do


!  Dust Large-scale Wet Removal
!  ----------------------------
   KIN = .TRUE.
   do n = 1, self%nbins
      fwet = 0.8
      call WetRemovalGOCART2G(self%km, self%klid, self%nbins, self%nbins, n, self%cdt, 'dust', &
                              KIN, MAPL_GRAV, fwet, DU(:,:,:,n), ple, t, airdens, &
                              pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, DUWT, __RC__)
   end do

!  Compute diagnostics
!  -------------------
!  Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
   call Aero_Compute_Diags (self%diag_Mie, self%km, self%klid, 1, self%nbins, self%rlow, &
                            self%rup, self%wavelengths_profile*1.0e-9, &
                            self%wavelengths_vertint*1.0e-9, DU, MAPL_GRAV, t, airdens, &
                            rh2, u, v, delp, ple,tropp, &
                            DUSMASS, DUCMASS, DUMASS, DUEXTTAU, DUSTEXTTAU, DUSCATAU,DUSTSCATAU, &
                            DUSMASS25, DUCMASS25, DUMASS25, DUEXTT25, DUSCAT25, &
                            DUFLUXU, DUFLUXV, DUCONC, DUEXTCOEF, DUSCACOEF, &
                            DUEXTTFM, DUSCATFM, DUANGSTR, DUAERIDX, NO3nFlag=.false., __RC__ )

    RETURN_(ESMF_SUCCESS)

  end subroutine Run2

!============================================================================
!BOP
! !IROUTINE: Run_data -- ExtData Dust Grid Component

! !INTERFACE:
  subroutine Run_data (GC, import, export, internal, RC)

    ! !ARGUMENTS:
    type (ESMF_GridComp), intent(inout) :: GC       ! Gridded component 
    type (ESMF_State),    intent(inout) :: import   ! Import state
    type (ESMF_State),    intent(inout) :: export   ! Export state
    type (ESMF_State),    intent(inout) :: internal ! Interal state
    integer, optional,    intent(  out) :: RC       ! Error code:

! !DESCRIPTION: Updates pointers in Internal state with fields from ExtData. 

! Locals
    character (len=ESMF_MAXSTR)        :: COMP_NAME
    type (wrap_)                       :: wrap
    type (DU2G_GridComp), pointer      :: self

    integer                            :: i
    character (len=ESMF_MAXSTR)        :: field_name

    real, pointer, dimension(:,:,:,:)  :: ptr4d_int
    real, pointer, dimension(:,:,:)    :: ptr3d_imp

    __Iam__('Run_data')

!EOP
!*****************************************************************************
!   Begin... 

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet (GC, NAME=COMP_NAME, __RC__)
    Iam = trim(COMP_NAME) //'::'//Iam

!   Get my private internal state
!   ------------------------------
    call ESMF_UserCompGetInternalState(GC, 'DU2G_GridComp', wrap, STATUS)
    VERIFY_(STATUS)
    self => wrap%ptr

!   Update interal data pointers with ExtData
!   -----------------------------------------
    call MAPL_GetPointer (internal, NAME='DU', ptr=ptr4d_int, __RC__)

    do i = 1, self%nbins
    write (field_name, '(A, I0.3)') 'du', i
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
    type(DU2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld

    real(kind=DP), dimension(:,:,:), allocatable     :: ext_s, ssa_s, asy_s  ! (lon:,lat:,lev:)
    real, dimension(:,:,:), allocatable              :: x
    integer                                          :: instance
    integer                                          :: n, nbins
    integer                                          :: i1, j1, i2, j2, km
    integer                                          :: band

    integer :: k

    __Iam__('DU2G::aerosol_optics')

!   Begin... 

!   Mie Table instance/index
!   ------------------------
    call ESMF_AttributeGet (state, name='mie_table_instance', value=instance, __RC__)

!   Radiation band
!   --------------
    band = 0
    call ESMF_AttributeGet (state, name='band_for_aerosol_optics', value=band, __RC__)

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

    allocate(ext_s(i1:i2, j1:j2, km), &
             ssa_s(i1:i2, j1:j2, km), &
             asy_s(i1:i2, j1:j2, km), &
                 x(i1:i2, j1:j2, km), __STAT__)

    call ESMF_StateGet (state, 'DU', field=fld, __RC__)
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

    call ESMF_AttributeGet (state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = ext_s(:,:,:)
    end if

    call ESMF_AttributeGet (state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
        var = ssa_s(:,:,:)
    end if

    call ESMF_AttributeGet (state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
    if (fld_name /= '') then
        call MAPL_GetPointer (state, var, trim(fld_name), __RC__)
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
    real(kind=DP), intent(  out) :: bext_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=DP), intent(  out) :: bssa_s (size(ext_s,1),size(ext_s,2),size(ext_s,3))
    real(kind=DP), intent(  out) :: basym_s(size(ext_s,1),size(ext_s,2),size(ext_s,3))
    integer,                       intent(  out) :: rc

    ! local
    integer                           :: l
    real                              :: bext (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! extinction
    real                              :: bssa (size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! SSA
    real                              :: gasym(size(ext_s,1),size(ext_s,2),size(ext_s,3))  ! asymmetry parameter

    __Iam__('DU2G::aerosol_optics::mie_')

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
    type(DU2G_GridComp), pointer                     :: self

    character (len=ESMF_MAXSTR)                      :: fld_name
    type(ESMF_Field)                                 :: fld

    real, dimension(:,:,:), allocatable              :: tau_s, tau, x ! (lon:,lat:,lev:)
    integer                                          :: instance
    integer                                          :: n, nbins, k
    integer                                          :: i1, j1, i2, j2, km, i, j
    real                                             :: wavelength

    __Iam__('DU2G::monochromatic_aerosol_optics')

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

    call ESMF_StateGet (state, 'DU', field=fld, __RC__)
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


end module DU2G_GridCompMod

