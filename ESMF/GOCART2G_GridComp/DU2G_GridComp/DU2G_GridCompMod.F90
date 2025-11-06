#include "MAPL_Generic.h"

!BOP
!MODULE: DU2G_GridCompMod - GOCART refactoring of the DU gridded component

!INTERFACE:
module DU2G_GridCompMod

   !USES:
   use ESMF
   use pflogger, only: logger_t => logger
   use mapl_ErrorHandling, only: MAPL_Verify, MAPL_Assert, MAPL_Return
   ! use MAPL
   use MAPL, only: MAPL_get_num_threads, MAPL_MetaComp
   use MAPL_MaplGrid, only: MAPL2_GridGet => MAPL_GridGet
   use mapl3g_generic, only: MAPL_GridCompGet, MAPL_GridCompGetResource, MAPL_GridCompGetInternalState
   use mapl3g_generic, only: MAPL_GridCompSetEntryPoint
   use mapl3g_generic, only: MAPL_GridCompAddSpec
   use mapl3g_generic, only: MAPL_STATEITEM_STATE, MAPL_STATEITEM_FIELDBUNDLE
   use mapl3g_generic, only: MAPL_ClockGet
   use mapl3g_VerticalStaggerLoc, only: VERTICAL_STAGGER_NONE, VERTICAL_STAGGER_CENTER, VERTICAL_STAGGER_EDGE
   use mapl3g_RestartModes, only: MAPL_RESTART_SKIP
   use mapl3g_UngriddedDim, only: UngriddedDim
   use mapl3g_State_API, only: MAPL_StateGetPointer
   use GOCART2G_MieMod
   use Chem_AeroGeneric
   use iso_c_binding, only: c_loc, c_f_pointer, c_ptr

   use GOCART2G_Process       ! GOCART2G process library
   use GA_EnvironmentMod
   use MAPL_StringTemplate, only: StrTemplate
   !$ use omp_lib

   implicit none
   private

   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData = 2

   !PUBLIC MEMBER FUNCTIONS:
   public  SetServices

   !DESCRIPTION: This module implements GOCART's Dust (DU) Gridded Component.

   !REVISION HISTORY:
   ! 4January2024   Collow - Updated call for ChemSettling
   ! 16Oct2019  Sherman, da Silva, Darmenov, Clune -  First attempt at refactoring

   !EOP

   integer, parameter :: NHRES = 6

   !Dust state
   type :: ThreadWorkspace
      integer :: day_save = -1
      integer :: nPts = -1
      integer, allocatable, dimension(:) :: pstart, pend
      real, allocatable, dimension(:) :: pLat, pLon, pBase, pTop, pEmis
   end type ThreadWorkspace

   type, extends(GA_Environment) :: DU2G_GridComp
      real, allocatable :: rlow(:)      ! particle effective radius lower bound [um]
      real, allocatable :: rup(:)       ! particle effective radius upper bound [um]
      real, allocatable :: sfrac(:)     ! fraction of total source
      real, allocatable :: sdist(:)     ! FENGSHA aerosol fractional size distribution [1]
      real, allocatable :: Ch_DU_res(:) ! resolutions used for Ch_DU
      real :: alpha                     ! FENGSHA scaling factor
      real :: gamma                     ! FENGSHA tuning exponent
      real :: kvhmax                    ! FENGSHA max. vertical/horizontal mass flux ratio [1]
      real :: f_sdl                     ! FENGSHA drylimit tuning factor
      real :: Ch_DU                     ! dust emission tuning coefficient [kg s2 m-5].
      real :: f_swc                     ! soil mosture scaling factor
      real :: f_scl                     ! clay content scaling factor
      real :: uts_gamma                 ! threshold friction velocity parameter 'gamma'
      logical :: maringFlag=.false.     ! maring settling velocity correction
      integer :: drag_opt               ! FENGSHA drag option 1 - input only, 2 - Darmenova, 3 - Leung
      integer :: distribution_opt       ! FENGSHA distribution option 1 - Kok, 2 - Kok 2021, 3 - Meng 2022
      integer :: day_save = -1
      integer :: clayFlag               ! clay and silt term in K14
      character(len=:), allocatable :: emission_scheme ! emission scheme selector
      ! Workspae for point emissions
      logical :: doing_point_emissions = .false.
      character(len=:), allocatable :: point_emissions_srcfilen ! filename for pointwise emissions
      type(ThreadWorkspace), allocatable :: workspaces(:)
   end type DU2G_GridComp

   type wrap_
      type (DU2G_GridComp), pointer :: PTR !=> null()
   end type wrap_

contains

   !BOP
   !IROUTINE: SetServices
   !INTERFACE:
   subroutine SetServices(gc, rc)

      !ARGUMENTS:
      type(ESMF_GridComp), intent(inout) :: gc
      integer, intent(out) :: rc  ! return code

      !DESCRIPTION: This version uses MAPL_GenericSetServices, which sets
      !     the Initialize and Finalize services to generic versions. It also
      !     allocates our instance of a generic state and puts it in the
      !     gridded component (gc). Here we only set the two-stage run method
      !     and declare the data services.

      !REVISION HISTORY:
      !   16oct2019   E.Sherman, A.Da Silva, A.Darmenov, T.Clune  First attempt at refactoring
      !EOP

      character(len=ESMF_MAXSTR) :: comp_name
      type(wrap_) :: wrap
      type(DU2G_GridComp), pointer :: self

      character(len=ESMF_MAXSTR) :: field_name
      character(len=:), allocatable :: emission_scheme
      real :: DEFVAL
      logical :: data_driven = .true., file_exists
      integer :: i, num_threads
      character(len=255) :: msg
      class(logger_t), pointer :: logger
      type(UngriddedDim) :: ungrd_nbins
      type(UngriddedDim) :: ungrd_wavelengths_profile, ungrd_wavelengths_vertint
      integer :: status

      call ESMF_GridCompGet(gc, name=comp_name, _RC)
      call MAPL_GridCompGet(gc, logger=logger, _RC)
      call logger%info("SetServices:: starting...")

      ! Wrap internal state for storing in gc
      allocate(self, _STAT)
      wrap%ptr => self

      num_threads = MAPL_get_num_threads()
      allocate(self%workspaces(0:num_threads-1), __STAT__)

      ! ! Load resource file
      ! cfg = ESMF_ConfigCreate (_RC)
      ! inquire(file='DU2G_instance_'//trim(comp_name)//'.rc', exist=file_exists)
      ! if (file_exists) then
      !    call ESMF_ConfigLoadFile (cfg, 'DU2G_instance_'//trim(comp_name)//'.rc', _RC)
      ! else
      !    if (mapl_am_i_root()) print*,'DU2G_instance_'//trim(comp_name)//'.rc does not exist! Loading DU2G_GridComp_DU.rc instead'
      !    call ESMF_ConfigLoadFile (cfg, 'DU2G_instance_DU.rc', _RC)
      ! end if

      ! process generic config items
      call self%GA_Environment%load_from_config(gc, _RC)

      ! Defined UngriddedDim items
      ungrd_nbins = UngriddedDim(self%nbins, name="nbins", units="1")
      ungrd_wavelengths_profile = UngriddedDim( &
           size(self%wavelengths_profile), &
           name="wavelengths_profile", &
           units="nm")
      ungrd_wavelengths_vertint = UngriddedDim( &
           size(self%wavelengths_vertint), &
           name="wavelengths_vertint", &
           units="nm")

      ! allocate(self%sfrac(self%nbins), self%rlow(self%nbins), self%rup(self%nbins), __STAT__)
      ! process DU-specific items
      call MAPL_GridCompGetResource(gc, "maringFlag", self%maringFlag, _RC)
      call MAPL_GridCompGetResource(gc, "source_fraction", self%sfrac, _RC)
      call MAPL_GridCompGetResource(gc, "Ch_DU", self%Ch_DU_res, _RC)
      _ASSERT(size(self%Ch_DU_res)==NHRES, "incorrect size of Ch_DU")

      call MAPL_GridCompGetResource(gc, "radius_lower", self%rlow, _RC)
      call MAPL_GridCompGetResource(gc, "radius_upper", self%rup, _RC)

      ! Choose Emission Scheme
      call MAPL_GridCompGetResource(gc, "emission_scheme", emission_scheme, default="ginoux", _RC)
      self%emission_scheme = ESMF_UtilStringLowerCase(trim(emission_scheme), _RC)

      ! Test if our scheme is allowed, if so, print it out
      _ASSERT(any(self%emission_scheme == [character(len=7) :: 'ginoux','k14','fengsha']), "Error. Unallowed emission scheme: "//trim(self%emission_scheme)//". Allowed: ginoux, k14, fengsha")

      ! Point Sources
      call MAPL_GridCompGetResource(gc, "point_emissions_srcfilen", self%point_emissions_srcfilen, default='/dev/null', _RC)
      if ( (index(self%point_emissions_srcfilen,'/dev/null')>0) ) then
         self%doing_point_emissions = .false. ! disable it if no file specified
      else
         self%doing_point_emissions = .true.  ! we are good to go
      end if

      ! read scheme-specific parameters
      select case (self%emission_scheme)
      case ('fengsha')
         call MAPL_GridCompGetResource(gc, "alpha", self%alpha, _RC)
         call MAPL_GridCompGetResource(gc, "gamma", self%gamma, _RC)
         call MAPL_GridCompGetResource(gc, "soil_moisture_factor", self%f_swc, _RC)
         call MAPL_GridCompGetResource(gc, "soil_drylimit_factor", self%f_sdl, _RC)
         call MAPL_GridCompGetResource(gc, "vertical_to_horizontal_flux_ratio_limit", self%kvhmax, _RC)
         call MAPL_GridCompGetResource(gc, "drag_partition_option", self%drag_opt, _RC)
      case ('k14')
         call MAPL_GridCompGetResource(gc, "clayFlag", self%clayFlag, _RC)
         call MAPL_GridCompGetResource(gc, "soil_moisture_factor", self%f_swc, _RC)
         call MAPL_GridCompGetResource(gc, "soil_clay_factor", self%f_scl, _RC)
         call MAPL_GridCompGetResource(gc, "uts_gamma", self%uts_gamma, _RC)
      case ('ginoux')
         ! nothing to do
      case default
         _ASSERT_RC(.false., "Unallowed emission scheme: "//trim(self%emission_scheme)//". Allowed: ginoux, k14, fengsha", ESMF_RC_NOT_IMPL)
      end select

      ! Is DU data driven?
      call determine_data_driven(comp_name, data_driven, _RC)

      ! Set entry points
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, _RC)
      call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, phase_name="Run1", _RC)
      if (data_driven .neqv. .true.) then
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run2, phase_name="Run2", _RC)
         call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run0, phase_name="Run0", _RC)
      end if

      DEFVAL = 0.0

      ! IMPORT STATE
      if (data_driven) then
         _FAIL("data driver section has not been activated yet")

         ! call MAPL_AddInternalSpec(gc, &
         !      short_name='DU', &
         !      long_name='Dust Mixing Ratio all bins', &
         !      units='kg kg-1', &
         !      dims=MAPL_DimsHorzVert, &
         !      vlocation=MAPL_VlocationCenter, &
         !      restart=MAPL_RestartOptional, &
         !      ungridded_dims=[self%nbins], &
         !      ! friendlyto='DYNAMICS:TURBULENCE:MOIST', &
         !      add2export=.true., _RC)

         ! ! Pressure at layer edges
         ! call MAPL_AddImportSpec(gc, &
         !      SHORT_NAME='PLE', &
         !      LONG_NAME='air_pressure', &
         !      UNITS='Pa', &
         !      DIMS=MAPL_DimsHorzVert, &
         !      VLOCATION=MAPL_VLocationEdge, &
         !      RESTART=MAPL_RestartSkip, _RC)

         ! ! RH: is between 0 and 1
         ! call MAPL_AddImportSpec(gc, &
         !      SHORT_NAME='RH2', &
         !      LONG_NAME='Rel_Hum_after_moist', &
         !      UNITS='1', &
         !      DIMS=MAPL_DimsHorzVert, &
         !      VLOCATION=MAPL_VLocationCenter, &
         !      RESTART=MAPL_RestartSkip, _RC)

         ! do i = 1, self%nbins
         !    write (field_name, '(A, I0.3)') '', i
         !    call MAPL_AddImportSpec(gc, &
         !         short_name='climdu'//trim(field_name), &
         !         long_name='Dust Mixing Ratio (bin '//trim(field_name)//')', &
         !         units='kg kg-1', &
         !         restart=MAPL_RestartSkip,  &
         !         dims=MAPL_DimsHorzVert, &
         !         vlocation=MAPL_VLocationCenter, _RC)

         !    ! dry deposition
         !    call MAPL_AddImportSpec(gc, &
         !         short_name='climDUDP'//trim(field_name), &
         !         long_name ='Dust Mixing Ratio (bin '//trim(field_name)//')', &
         !         units='kg kg-1', &
         !         dims=MAPL_DimsHorzOnly, &
         !         vlocation=MAPL_VLocationCenter, &
         !         restart=MAPL_RestartSkip, _RC)

         !    !        !wet deposition
         !    call MAPL_AddImportSpec(gc, &
         !         short_name='climDUWT'//trim(field_name), &
         !         long_name='Dust Mixing Ratio (bin '//trim(field_name)//')', &
         !         units='kg kg-1', &
         !         dims=MAPL_DimsHorzOnly, &
         !         vlocation=MAPL_VLocationCenter, &
         !         restart=MAPL_RestartSkip, _RC)

         !    !        !gravitational settling
         !    call MAPL_AddImportSpec(gc, &
         !         short_name='climDUSD'//trim(field_name), &
         !         long_name='Dust Mixing Ratio (bin '//trim(field_name)//')', &
         !         units='kg kg-1', &
         !         dims=MAPL_DimsHorzOnly, &
         !         vlocation=MAPL_VLocationCenter, &
         !         restart=MAPL_RestartSkip, _RC)

         !    !        !convective scavenging
         !    call MAPL_AddImportSpec(gc, &
         !         short_name='climDUSV'//trim(field_name), &
         !         long_name='Dust Mixing Ratio (bin '//trim(field_name)//')', &
         !         units='kg kg-1', &
         !         dims=MAPL_DimsHorzOnly, &
         !         vlocation=MAPL_VLocationCenter, &
         !         restart=MAPL_RestartSkip, _RC)
         ! end do
      end if ! (data_driven)

      ! Computational Instance
      if (.not. data_driven) then
#include "DU2G_Export___.h"
         associate (scheme => self%emission_scheme)
#include "DU2G_Import___.h"
         end associate
#include "DU2G_Internal___.h"
      end if

      ! This state holds fields needed by radiation
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name=trim(comp_name)//"_AERO", &
           standard_name="aerosols_from_"//trim(comp_name), &
           units="kg kg-1", &
           dims="xyz", &
           vstagger=VERTICAL_STAGGER_CENTER, &
           itemtype=MAPL_STATEITEM_STATE, _RC)

      ! This bundle is needed by surface for snow albedo modification
      ! by aerosol settling and deposition
      ! ~~~DEVELOPERS NOTE~~~ Change to StateItem when possible
      call MAPL_GridCompAddSpec(gc, &
           state_intent=ESMF_STATEINTENT_EXPORT, &
           short_name=trim(comp_name)//"_AERO_DP", &
           standard_name="aerosol_deposition_from_"//trim(comp_name), &
           units="kg m-2 s-1", &
           dims="xy", &
           vstagger=VERTICAL_STAGGER_CENTER, &
           itemtype=MAPL_STATEITEM_FIELDBUNDLE, _RC)

      ! Store internal state in gc
      call ESMF_UserCompSetInternalState(gc, "DU2G_GridComp", wrap, _RC)

      call logger%info("SetServices:: ...complete")
      _RETURN(_SUCCESS)

   end subroutine SetServices

   !BOP
   !IROUTINE: Initialize
   !INTERFACE:
   subroutine Initialize(gc, import, export, clock, rc)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: This initializes DU's Grid Component. It primaryily fills
      !               GOCART's AERO states with its dust fields.

      !REVISION HISTORY:
      ! 16oct2019  E.Sherman, A.da Silva, T.Clune, A.Darmenov - First attempt at refactoring
      !EOP

      type(ESMF_Grid) :: grid
      type(ESMF_Geom) :: geom
      type(ESMF_State) :: internal, aero, provider_state
      type(ESMF_FieldBundle) :: bundle_dp
      type(ESMF_Field) :: field, fld
      type(ESMF_Info) :: field_info, aero_info
      type(wrap_) :: wrap
      type(DU2G_GridComp), pointer :: self
      integer, allocatable :: mieTable_pointer(:), channels_(:)
      real :: CDT ! chemistry timestep (secs)
      real :: HDT ! model timestep (secs)
      integer :: ibin, dims(3), km, instance, nmom_, status
      logical :: data_driven, file_exists
      character(len=ESMF_MAXSTR) :: bin_index, prefix
      character(len=:), allocatable :: comp_name, file_
      class(logger_t), pointer :: logger

      call MAPL_GridCompGet (gc, name=comp_name, logger=logger, _RC)
      call logger%info("Initialize:: starting...")

      ! Get my internal private state
      call ESMF_UserCompGetInternalState(gc, "DU2G_GridComp", wrap, _RC)
      self => wrap%ptr

      ! Global dimensions are needed here for choosing tuning parameters
      call MAPL_GridCompGet(gc, grid=grid, num_levels=km, _RC)
      call MAPL2_GridGet(grid, globalCellCountPerDim=dims, _RC)
      self%km = km

      ! Dust emission tuning coefficient [kg s2 m-5]. NOT bin specific.
      ! TO DO: find a more robust way to implement resolution dependent tuning
      self%Ch_DU = Chem_UtilResVal(dims(1), dims(2), self%Ch_DU_res(:), _RC)
      self%Ch_DU = self%Ch_DU * 1.0e-9

      ! Dust emission size distribution for FENGSHA
      if (self%emission_scheme == "fengsha") then
         allocate(self%sdist(self%nbins), __STAT__)
         call DustAerosolDistributionKok(self%radius, self%rup, self%rlow, self%sdist)
      end if

      ! Get DTs
      call MAPL_ClockGet(clock, dt=HDT, _RC)
      call MAPL_GridCompGetResource(gc, "GOCART2G_DT", CDT, default=real(HDT), _RC)
      self%CDT = CDT

      ! Get parameters from generic state.
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is DU data driven?
      call determine_data_driven(comp_name, data_driven, _RC)

      ! If this is a data component, the data is provided in the import
      ! state via ExtData instead of the actual GOCART children
      provider_state = export
      prefix = ""
      if (data_driven) then
         provider_state = import
         prefix = "clim"
      end if

      ! Add attribute information for DU export. Used in NI hetergenous chemistry.
      call ESMF_StateGet(export, "DU", field, _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key="radius", values=self%radius, _RC)
      call ESMF_InfoSet(field_info, key="fnum", values=self%fnum, _RC)

      ! Add attribute information to internal state variables
      ! Fill AERO States with dust fields
      call ESMF_StateGet (export, trim(comp_name)//"_AERO"    , aero    , _RC)
      call ESMF_StateGet (export, trim(comp_name)//"_AERO_DP" , bundle_dp, _RC)

      call ESMF_StateGet (internal, "DU", field, _RC)
      ! fld = MAPL_FieldCreate (field, "DU", _RC)
      ! call MAPL_StateAdd (aero, fld, _RC) ! pchakrab: TODO - this is equivalent to fld = field
      call ESMF_StateAdd(aero, [field], _RC)
      ! call ESMF_AttributeSet(field, NAME="ScavengingFractionPerKm", value=self%fscav(1), _RC)
      call ESMF_InfoGetFromHost(field, field_info, _RC)
      call ESMF_InfoSet(field_info, key="ScavengingFractionPerKm", value=self%fscav(1), _RC)      

      if (data_driven) then
         instance = instanceData
         do ibin = 1, self%nbins
            write (bin_index, "(A, I0.3)") "", ibin
            ! Dry deposition
            call append_to_bundle("DUDP"//trim(bin_index), provider_state, prefix, bundle_dp, _RC)
            ! Wet deposition (Convective scavenging)
            call append_to_bundle("DUSV"//trim(bin_index), provider_state, prefix, bundle_dp, _RC)
            ! Wet deposition
            call append_to_bundle("DUWT"//trim(bin_index), provider_state, prefix, bundle_dp, _RC)
            ! Gravitational Settling
            call append_to_bundle("DUSD"//trim(bin_index), provider_state, prefix, bundle_dp, _RC)
         end do
      else
         instance = instanceComputational
         ! Dry deposition
         call append_to_bundle("DUDP", provider_state, prefix, bundle_dp, _RC)
         ! Wet deposition (Convective scavenging)
         call append_to_bundle("DUSV", provider_state, prefix, bundle_dp, _RC)
         ! Wet deposition
         call append_to_bundle("DUWT", provider_state, prefix, bundle_dp, _RC)
         ! Gravitational Settling
         call append_to_bundle("DUSD", provider_state, prefix, bundle_dp, _RC)
      end if
      self%instance = instance

      ! Create Radiation Mie Table
      call MAPL_GridCompGetResource(gc, "aerosol_radBands_optics_file", file_, _RC )
      self%rad_Mie = GOCART2G_Mie(trim(file_), _RC)

      ! Create Diagnostics Mie Table
      ! Get file names for the optical tables
      call MAPL_GridCompGetResource(gc, "aerosol_monochromatic_optics_file", file_, _RC )
      call MAPL_GridCompGetResource(gc, "n_moments", nmom_, default=0,  _RC)
      call MAPL_GridCompGetResource(gc, "aerosol_monochromatic_optics_wavelength_in_nm_from_LUT", channels_, _RC)
      self%diag_Mie = GOCART2G_Mie(trim(file_), channels_*1.e-9, nmom=nmom_, _RC)

      ! Mie Table instance/index
      ! call ESMF_AttributeSet (aero, name="mie_table_instance", value=instance, _RC)
      call ESMF_InfoGetFromHost(aero, aero_info, _RC)
      call ESMF_InfoSet(aero_info, key="mie_table_instance", value=instance, _RC)      

      ! Add variables to DU instance's aero state. This is used in aerosol optics calculations
      call MAPL_GridCompGet(gc, geom=geom, _RC)
      call add_aero(aero, label="air_pressure_for_aerosol_optics", label2="PLE", geom=geom, km=self%km, _RC)
      call add_aero(aero, label="relative_humidity_for_aerosol_optics", label2="RH", geom=geom, km=self%km, _RC)
      ! call ESMF_StateGet (import, "PLE", field, _RC)
      ! call MAPL2_StateAdd (aero, field, _RC)
      ! call ESMF_StateGet (import, "RH2", field, _RC)
      ! call MAPL2_StateAdd (aero, field, _RC)
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
      call add_aero( &
           aero, &
           label="monochromatic_extinction_in_air_due_to_ambient_aerosol", label2="monochromatic_EXT", &
           geom=geom, typekind=ESMF_TYPEKIND_R4,_RC)
      call add_aero(aero, label="sum_of_internalState_aerosol", label2="aerosolSum", geom=geom, km=self%km, _RC)

      call ESMF_InfoGetFromHost(aero, aero_info, _RC)
      call ESMF_InfoSet(aero_info, key="band_for_aerosol_optics", value=0, _RC)
      call ESMF_InfoSet(aero_info, key="wavelength_for_aerosol_optics", value=0., _RC)
      mieTable_pointer = transfer(c_loc(self), [1])
      call ESMF_InfoSet(aero_info, key="mieTable_pointer", values=mieTable_pointer, _RC)
      call ESMF_InfoSet(aero_info, key="internal_variable_name", value="DU", _RC)
      ! Add callback methods
      call ESMF_MethodAdd(aero, label="aerosol_optics", userRoutine=aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label="monochromatic_aerosol_optics", userRoutine=monochromatic_aerosol_optics, _RC)
      call ESMF_MethodAdd(aero, label="get_mixR", userRoutine=get_mixR, _RC)

      call logger%info("Initialize:: ...complete")
      _RETURN(_SUCCESS)

   end subroutine Initialize

   !BOP
   !IROUTINE: Run0
   !INTERFACE:
   subroutine Run0(gc, import, export, clock, RC)
      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION:  Clears klid to 0.0 for Dust
      !EOP

      type(ESMF_State) :: internal
      type(wrap_) :: wrap
      type(DU2G_GridComp), pointer :: self
      real, pointer, dimension(:,:,:) :: ple
      real, allocatable, dimension(:,:,:) :: ple0
      real, pointer, dimension(:,:,:,:) :: ptr4d_int
      integer :: i1, i2, j1, j2, km, status
      class(logger_t), pointer :: logger

      call MAPL_GridCompGet(gc, logger=logger, _RC)
      call logger%info("Run0:: starting...")

      ! Get parameters from generic state.
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, "DU2G_GridComp", wrap, _RC)
      self => wrap%ptr

      ! Set klid and Set internal values to 0 above klid
      km = self%km
      call MAPL_StateGetPointer(import, ple, "PLE", _RC)
      i1 = lbound(ple, 1); i2 = ubound(ple, 1); j1 = lbound(ple, 2); j2 = ubound(ple, 2)
      allocate(ple0(i1:i2, j1:j2, 0:km), source=ple(i1:i2, j1:j2, 1:km+1))
      call findKlid(self%klid, self%plid, ple0, _RC)
      call MAPL_StateGetPointer(internal, ptr4d_int, "DU", _RC)
      call setZeroKlid4d(self%km, self%klid, ptr4d_int)

      call logger%info("Run0:: ...complete")
      _RETURN(_SUCCESS)
   end subroutine Run0

   !BOP
   !IROUTINE: Run
   !INTERFACE:
   subroutine Run(gc, import, export, clock, rc)
      !ARGUMENTS:
      type (ESMF_GridComp) :: gc
      type (ESMF_State) :: import
      type (ESMF_State) :: export
      type (ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: Run method for the Dust Grid Component. Determines whether to
      !                 run data or computational run method.
      !EOP

      character(len=:), allocatable :: comp_name
      type(ESMF_State) :: internal
      logical :: data_driven
      integer :: status
      class(logger_t), pointer :: logger

      call MAPL_GridCompGet(gc, name=comp_name, logger=logger, _RC)
      call logger%info("Run0:: starting...")

      ! Get parameters from generic state.
      call MAPL_GridCompGetInternalState(gc, internal, _RC)

      ! Is DU data driven?
      call determine_data_driven(comp_name, data_driven, _RC)

      ! Update INTERNAL state variables with ExtData
      if (data_driven) then
         call Run_data(gc, import, export, internal, _RC)
      else
         call Run1(gc, import, export, clock, _RC)
      end if

      call logger%info("Run0:: ...complete")
      _RETURN(_SUCCESS)
   end subroutine Run

   !BOP
   !IROUTINE: Run1
   !INTERFACE:
   subroutine Run1 (gc, import, export, clock, RC)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION:  Computes emissions/sources for Dust
      !EOP

      character (len=ESMF_MAXSTR)       :: comp_name
      type (MAPL_MetaComp), pointer     :: mapl
      type (ESMF_State)                 :: internal
      type (ESMF_Grid)                  :: grid
      type (wrap_)                      :: wrap
      type (DU2G_GridComp), pointer     :: self
      type(ESMF_Time)                   :: time

      integer  :: nymd, nhms, iyr, imm, idd, ihr, imn, isc
      integer  :: import_shape(2), i2, j2
      ! real, dimension(:,:), pointer     :: emissions_surface
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
      type(ThreadWorkspace), pointer :: workspace
      integer :: thread, jstart, jend

#include "DU2G_DeclarePointer___.h"

      __Iam__('Run1')

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet (gc, NAME=comp_name, __RC__)
      Iam = trim(comp_name) //'::'// Iam

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC (gc, mapl, __RC__)

      call MAPL_Get(mapl, grid=grid, __RC__)

      ! Get parameters from generic state.
      call MAPL_Get (mapl, INTERNAL_ESMF_STATE=internal, __RC__)

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, 'DU2G_GridComp', wrap, STATUS)
      VERIFY_(STATUS)
      self => wrap%ptr

      ! Extract nymd(yyyymmdd) from clock
      call ESMF_ClockGet (clock, currTime=time, __RC__)
      call ESMF_TimeGet (time ,YY=iyr, MM=imm, DD=idd, H=ihr, M=imn, S=isc, __RC__)
      call MAPL_PackTime (nymd, iyr, imm , idd)
      call MAPL_PackTime (nhms, ihr, imn, isc)

      associate (scheme => self%emission_scheme)
#include "DU2G_GetPointer___.h"
      end associate

      ! Set du_src to 0 where undefined
      if (associated(du_src)) then
         where (1.01*du_src > MAPL_UNDEF) du_src = 0.
      endif
      ! Get dimensions
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

      ! Get surface gridded emissions
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

         call DustEmissionK14( &
              self%km, tsoil1, wcsf, rhos, &
              du_z0, z_, u10n, v10n, ustar, &
              frland, asnow, &
              du_src, &
              du_sand, du_silt, du_clay, &
              du_texture, du_veg, du_gvf, &
              self%f_swc, self%f_scl, self%uts_gamma, &
              MAPL_UNDEF, MAPL_GRAV, MAPL_KARMAN, &
              self%clayFlag, self%Ch_DU/1.e-9, &
              emissions_surface, &
              ustar_, &
              ustar_t_, &
              ustar_ts_, &
              R_, H_w_, f_erod_, __RC__ )

         if (associated(DU_UST)) DU_UST = ustar_
         if (associated(DU_UST_T)) DU_UST_T = ustar_t_
         if (associated(DU_UST_T)) DU_UST_T = ustar_ts_
         if (associated(DU_DPC)) DU_DPC = R_
         if (associated(DU_SMC)) DU_SMC = H_w_
         if (associated(DU_EROD)) DU_EROD = f_erod_

      case ('fengsha')
         call DustEmissionFENGSHA( &
              frlake, frsnow, lwi, slc, du_clay, du_sand, du_silt, &
              du_ssm, du_rdrag, airdens(:,:,self%km), ustar, du_gvf, du_lai, du_uthres, &
              self%alpha, self%gamma, self%kvhmax, MAPL_GRAV, &
              self%rhop, self%sdist, self%f_sdl, self%f_swc, self%drag_opt, emissions_surface,  __RC__)

      case ('ginoux')

         call DustEmissionGOCART2G( &
              self%radius*1.e-6, frlake, wet1, lwi, u10m, v10m, &
              self%Ch_DU, du_src, MAPL_GRAV, &
              emissions_surface, __RC__)
      case default
         _ASSERT_RC(.false.,'missing dust emission scheme. Allowed: ginoux, fengsha, k14',ESMF_RC_NOT_IMPL)
      end select

      ! Read point emissions file once per day
      thread = MAPL_get_current_thread()
      workspace => self%workspaces(thread)
      if (self%doing_point_emissions) then
         if (workspace%day_save /= idd) then
            workspace%day_save = idd
            call StrTemplate( &
                 fname, self%point_emissions_srcfilen, xid='unknown', &
                 nymd=nymd, nhms=120000 )
            inquire( file=fname, exist=fileExists)
            if (fileExists) then
               call ReadPointEmissions( &
                    nymd, fname, workspace%nPts, workspace%pLat, workspace%pLon, &
                    workspace%pBase, workspace%pTop, workspace%pEmis, workspace%pStart, &
                    workspace%pEnd, label='source', __RC__)
            else if (.not. fileExists) then
               !$omp critical (DU2G_1)
               if(mapl_am_i_root()) print*,'GOCART2G ',trim(comp_name),': ',trim(fname),' not found; proceeding.'
               !$omp end critical (DU2G_1)
               workspace%nPts = -1 ! set this back to -1 so the "if (workspace%nPts > 0)" conditional is not exercised.
            end if
         end if
      end if

      ! Get indices for point emissions
      if (workspace%nPts > 0) then
         allocate(iPoint(workspace%nPts), jPoint(workspace%nPts),  __STAT__)
         call MAPL_GetHorzIJIndex( &
              workspace%nPts, iPoint, jPoint, &
              grid=grid, &
              lon=workspace%pLon/real(MAPL_RADIANS_TO_DEGREES), &
              lat=workspace%pLat/real(MAPL_RADIANS_TO_DEGREES), &
              rc=status)
         if (status /= 0) then
            !$omp critical (DU2G_2)
            if (mapl_am_i_root()) print*, trim(Iam), ' - cannot get indices for point emissions'
            !$omp end critical (DU2G_2)
            VERIFY_(status)
         end if

         call updatePointwiseEmissions( &
              self%km, workspace%pBase, workspace%pTop, workspace%pEmis, workspace%nPts, &
              workspace%pStart, workspace%pEnd, zle, &
              area, iPoint, jPoint, nhms, emissions_point, __RC__)
      end if

      ! Update aerosol state
      call UpdateAerosolState( &
           emissions, emissions_surface, emissions_point, &
           self%sfrac, workspace%nPts, self%km, self%CDT, MAPL_GRAV, &
           self%nbins, delp, DU, __RC__)

      if (associated(DUEM)) then
         DUEM = sum(emissions, dim=3)
      end if

      ! Clean up
      deallocate(emissions, emissions_surface, emissions_point, __STAT__)
      if (workspace%nPts > 0) then
         deallocate(iPoint, jPoint, __STAT__)
      end if

      RETURN_(ESMF_SUCCESS)

   end subroutine Run1

   !BOP
   !IROUTINE: Run2
   !INTERFACE:
   subroutine Run2 (gc, import, export, clock, RC)

      !ARGUMENTS:
      type(ESMF_GridComp) :: gc
      type(ESMF_State) :: import
      type(ESMF_State) :: export
      type(ESMF_Clock) :: clock
      integer, intent(out) :: rc

      !DESCRIPTION: Run2 method for the Dust Grid Component.
      !EOP

      character (len=ESMF_MAXSTR)       :: comp_name
      type (MAPL_MetaComp), pointer     :: MAPL
      type (ESMF_State)                 :: internal
      type (wrap_)                      :: wrap
      type (DU2G_GridComp), pointer     :: self

      integer                           :: n
      real, allocatable, dimension(:,:) :: drydepositionfrequency, dqa
      real                              :: fwet
      logical                           :: KIN

      integer                           :: i1, j1, i2, j2, km
      real, dimension(3)                :: rainout_eff
      real, parameter ::  cpd    = 1004.16
      real, target, allocatable, dimension(:,:,:)   :: RH20,RH80
      real, pointer, dimension(:,:)     :: flux_ptr
      integer                           :: settling_opt
#include "DU2G_DeclarePointer___.h"

      __Iam__('Run2')

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet (gc, NAME=comp_name, __RC__)
      Iam = trim(comp_name) // '::' // Iam

      ! Get my internal MAPL_Generic state
      call MAPL_GetObjectFromGC (gc, MAPL, __RC__)

      ! Get parameters from generic state.
      call MAPL_Get (MAPL, INTERNAL_ESMF_STATE=internal, __RC__)

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, 'DU2G_GridComp', wrap, STATUS)
      VERIFY_(STATUS)
      self => wrap%ptr

      associate (scheme => self%emission_scheme)
#include "DU2G_GetPointer___.h"
      end associate

      allocate(dqa, mold=wet1, __STAT__)
      allocate(drydepositionfrequency, mold=wet1, __STAT__)

      ! Set klid and Set internal DU values to 0 above klid
      call findKlid (self%klid, self%plid, ple, __RC__)
      call setZeroKlid4d (self%km, self%klid, DU)

      ! Dust Settling
      select case (self%settling_scheme)
      case ('gocart')
         settling_opt = 1
      case ('ufs')
         settling_opt = 2
      case default
         _ASSERT_RC(.false.,'Unsupported settling scheme: '//trim(self%settling_scheme),ESMF_RC_NOT_IMPL)
      end select

      do n = 1, self%nbins
         nullify(flux_ptr)
         if (associated(DUSD)) flux_ptr => DUSD(:,:,n)
         call Chem_SettlingSimple( &
              self%km, self%klid, self%diag_Mie, n, self%cdt, MAPL_GRAV, &
              DU(:,:,:,n), t, airdens, &
              rh2, zle, delp, flux_ptr, correctionMaring=self%maringFlag, &
              settling_scheme=settling_opt, __RC__)
      end do

      ! Dust Deposition
      do n = 1, self%nbins
         drydepositionfrequency = 0.
         call DryDeposition( &
              self%km, t, airdens, zle, lwi, ustar, zpbl, sh,&
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


      ! Dust Large-scale Wet Removal
      KIN = .TRUE.
      select case (self%wet_removal_scheme)
      case ('gocart')
         do n = 1, self%nbins
            fwet = 1.0
            call WetRemovalGOCART2G( &
                 self%km, self%klid, self%nbins, self%nbins, n, self%cdt, 'dust', &
                 KIN, MAPL_GRAV, fwet, DU(:,:,:,n), ple, t, airdens, &
                 pfl_lsan, pfi_lsan, cn_prcp, ncn_prcp, DUWT, __RC__)
         end do
      case ('ufs')
         rainout_eff = 0.0
         do n = 1, self%nbins
            rainout_eff(1)   = self%fwet_ice(n)  ! remove with ice
            rainout_eff(2)   = self%fwet_snow(n) ! remove with snow
            rainout_eff(3)   = self%fwet_rain(n) ! remove with rain
            call WetRemovalUFS( &
                 self%km, self%klid, n, self%cdt, 'dust', KIN, MAPL_GRAV, &
                 self%radius(n), rainout_eff, self%washout_tuning, self%wet_radius_thr, &
                 DU(:,:,:,n), ple, t, airdens, pfl_lsan, pfi_lsan, DUWT, __RC__)
         end do
      case default
         _ASSERT_RC(.false.,'Unsupported wet removal scheme: '//trim(self%wet_removal_scheme),ESMF_RC_NOT_IMPL)
      end select

      ! Compute diagnostics
      !  Certain variables are multiplied by 1.0e-9 to convert from nanometers to meters
      call Aero_Compute_Diags( &
           self%diag_Mie, self%km, self%klid, 1, self%nbins, self%rlow, &
           self%rup, self%wavelengths_profile*1.0e-9, &
           self%wavelengths_vertint*1.0e-9, DU, MAPL_GRAV, t, airdens, &
           rh2, u, v, delp, ple,tropp, &
           DUSMASS, DUCMASS, DUMASS, DUEXTTAU, DUSTEXTTAU, DUSCATAU,DUSTSCATAU, &
           DUSMASS25, DUCMASS25, DUMASS25, DUEXTT25, DUSCAT25, &
           DUFLUXU, DUFLUXV, DUCONC, DUEXTCOEF, DUSCACOEF, &
           DUBCKCOEF,DUEXTTFM, DUSCATFM, DUANGSTR, DUAERIDX, NO3nFlag=.false., __RC__ )

      i1 = lbound(RH2, 1); i2 = ubound(RH2, 1)
      j1 = lbound(RH2, 2); j2 = ubound(RH2, 2)
      km = ubound(RH2, 3)

      allocate(RH20(i1:i2,j1:j2,km), __STAT__)
      allocate(RH80(i1:i2,j1:j2,km), __STAT__)

      RH20(:,:,:) = 0.20
      call Aero_Compute_Diags( &
           mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=self%nbins, rlow=self%rlow, &
           rup=self%rup, wavelengths_profile=self%wavelengths_profile*1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint*1.0e-9, aerosol=DU, &
           grav=MAPL_GRAV, tmpu=t, rhoa=airdens, &
           rh=rh20, u=u, v=v, delp=delp, ple=ple,tropp=tropp, &
           extcoef = DUEXTCOEFRH20, scacoef = DUSCACOEFRH20, NO3nFlag=.False., __RC__)

      RH80(:,:,:) = 0.80

      call Aero_Compute_Diags( &
           mie=self%diag_Mie, km=self%km, klid=self%klid, nbegin=1, &
           nbins=self%nbins, rlow=self%rlow, &
           rup=self%rup, wavelengths_profile=self%wavelengths_profile*1.0e-9, &
           wavelengths_vertint=self%wavelengths_vertint*1.0e-9, aerosol=DU, &
           grav=MAPL_GRAV, tmpu=t, rhoa=airdens, &
           rh=rh80, u=u, v=v, delp=delp, ple=ple,tropp=tropp, &
           extcoef = DUEXTCOEFRH80, scacoef = DUSCACOEFRH80, NO3nFlag=.False., __RC__)

      deallocate(RH20,RH80)
      RETURN_(ESMF_SUCCESS)

   end subroutine Run2

   !BOP
   !IROUTINE: Run_data -- ExtData Dust Grid Component

   !INTERFACE:
   subroutine Run_data (gc, import, export, internal, RC)

      !ARGUMENTS:
      type (ESMF_GridComp), intent(inout) :: gc       ! Gridded component
      type (ESMF_State),    intent(inout) :: import   ! Import state
      type (ESMF_State),    intent(inout) :: export   ! Export state
      type (ESMF_State),    intent(inout) :: internal ! Interal state
      integer, optional,    intent(  out) :: RC       ! Error code:

      !DESCRIPTION: Updates pointers in Internal state with fields from ExtData.

      character (len=ESMF_MAXSTR)        :: comp_name
      type (wrap_)                       :: wrap
      type (DU2G_GridComp), pointer      :: self

      integer                            :: i
      character (len=ESMF_MAXSTR)        :: field_name

      real, pointer, dimension(:,:,:,:)  :: ptr4d_int
      real, pointer, dimension(:,:,:)    :: ptr3d_imp

      __Iam__('Run_data')

      !EOP

      ! Get my name and set-up traceback handle
      call ESMF_GridCompGet (gc, NAME=comp_name, __RC__)
      Iam = trim(comp_name) //'::'//Iam

      ! Get my private internal state
      call ESMF_UserCompGetInternalState(gc, 'DU2G_GridComp', wrap, STATUS)
      VERIFY_(STATUS)
      self => wrap%ptr

      ! Update interal data pointers with ExtData
      call MAPL_GetPointer (internal, NAME='DU', ptr=ptr4d_int, __RC__)

      do i = 1, self%nbins
         write (field_name, '(A, I0.3)') 'du', i
         call MAPL_GetPointer (import,  NAME='clim'//trim(field_name), ptr=ptr3d_imp, __RC__)

         ptr4d_int(:,:,:,i) = ptr3d_imp
      end do

      RETURN_(ESMF_SUCCESS)

   end subroutine Run_data

   subroutine aerosol_optics(state, rc)

      implicit none

      !ARGUMENTS:
      type(ESMF_State) :: state
      integer, intent(out) :: rc

      !Local
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

      ! Mie Table instance/index
      call ESMF_AttributeGet (state, name='mie_table_instance', value=instance, __RC__)

      ! Radiation band
      band = 0
      call ESMF_AttributeGet (state, name='band_for_aerosol_optics', value=band, __RC__)

      ! Pressure at layer edges
      call ESMF_AttributeGet (state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
      call MAPL_GetPointer (state, ple, trim(fld_name), __RC__)

      ! call MAPL_GetPointer (state, ple, 'PLE', __RC__)

      i1 = lbound(ple, 1); i2 = ubound(ple, 1)
      j1 = lbound(ple, 2); j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      ! Relative humidity
      call ESMF_AttributeGet (state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
      call MAPL_GetPointer (state, rh, trim(fld_name), __RC__)

      ! call MAPL_GetPointer (state, rh, 'RH2', __RC__)

      allocate( &
           ext_s(i1:i2, j1:j2, km), &
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

   subroutine monochromatic_aerosol_optics(state, rc)

      implicit none

      !ARGUMENTS:
      type (ESMF_State)                                :: state
      integer,            intent(out)                  :: rc

      !Local
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

      ! Mie Table instance/index
      call ESMF_AttributeGet (state, name='mie_table_instance', value=instance, __RC__)

      ! Radiation band
      wavelength = 0.
      call ESMF_AttributeGet (state, name='wavelength_for_aerosol_optics', value=wavelength, __RC__)

      ! Pressure at layer edges
      call ESMF_AttributeGet (state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
      call MAPL_GetPointer (state, ple, trim(fld_name), __RC__)

      ! call MAPL_GetPointer (state, ple, 'PLE', __RC__)
      i1 = lbound(ple, 1); i2 = ubound(ple, 1)
      j1 = lbound(ple, 2); j2 = ubound(ple, 2)
      km = ubound(ple, 3)

      ! Relative humidity
      call ESMF_AttributeGet (state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
      call MAPL_GetPointer (state, rh, trim(fld_name), __RC__)

      ! call MAPL_GetPointer (state, rh, 'RH2', __RC__)

      allocate( &
           tau_s(i1:i2, j1:j2, km), &
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

