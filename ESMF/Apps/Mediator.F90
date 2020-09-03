#include "MAPL_Generic.h"
#include "NUOPC_ErrLog.h"

module VarSpecMod
  use MAPL_Mod
  use ESMF
  use NUOPC
  implicit none

  integer, parameter :: NUM_GEOS_LEVELS = 72
  private
  public :: VarSpec, get_var_spec_from_transferred_field

  type VarSpec
     character(len=:), allocatable :: name
     type(ESMF_TypeKind_Flag) :: typekind
     integer, allocatable :: grid_to_field_map(:), ungridded_lbound(:), ungridded_ubound(:)
     integer :: dims, vlocation
   contains
     procedure :: complete_field_in_state
     procedure :: create_field
  end type VarSpec
  
contains

   function get_var_spec_from_transferred_field(field) result(spec)
    type(VarSpec) :: spec
    type(ESMF_Field), intent(in) :: field

    integer :: count
    integer :: typekind_int, rc
    character(len=ESMF_MAXSTR) :: field_name

    call ESMF_FieldGet(field, name = field_name, rc = rc)
    VERIFY_NUOPC_(rc)
    spec%name = field_name
    
    call ESMF_AttributeGet(field, name = "TypeKind", convention = "NUOPC", &
         purpose = "Instance", value = typekind_int, rc = rc)
    VERIFY_NUOPC_(rc)

    spec%typekind = typekind_int

    call NUOPC_GetAttribute(field, name = "GridToFieldMap", itemCount = count, rc = rc)
    VERIFY_NUOPC_(rc)

    if (count > 0) then
       allocate(spec%grid_to_field_map(count))

       call ESMF_AttributeGet(field, name = "GridToFieldMap", &
            convention = "NUOPC", purpose = "Instance", &
            valueList = spec%grid_to_field_map, rc = rc)
       VERIFY_NUOPC_(rc)
    end if

    _ASSERT(allocated(spec%grid_to_field_map), "grid_to_field_map on transferred field must be present")
    

     call NUOPC_GetAttribute(field, name = "UngriddedLBound", itemCount = count, rc = rc)
    VERIFY_NUOPC_(rc)

    if (count > 0) then
       allocate(spec%ungridded_lbound(count))

       call ESMF_AttributeGet(field, name = "UngriddedLBound", &
            convention = "NUOPC", purpose = "Instance", &
            valueList = spec%ungridded_lbound, rc = rc)
       VERIFY_NUOPC_(rc)
    end if
    
    
    call NUOPC_GetAttribute(field, name = "UngriddedUBound", itemCount = count, rc = rc)
    VERIFY_NUOPC_(rc)

    if (count > 0) then
       allocate(spec%ungridded_ubound(count))

       call ESMF_AttributeGet(field, name = "UngriddedUBound", &
            convention = "NUOPC", purpose = "Instance", &
            valueList = spec%ungridded_ubound, rc = rc)
       VERIFY_NUOPC_(rc)
    end if
    
  end function get_var_spec_from_transferred_field
  
  
  subroutine complete_field_in_state(spec, state, rc)
    class(VarSpec), intent(in) :: spec
    type(ESMF_State), intent(inout) :: state
    integer, intent(out) :: rc

    type(ESMF_Field) :: field

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, field = field, itemName = spec%name, rc = rc)
    VERIFY_NUOPC_(rc)

    call complete_field(spec, field, rc)
    VERIFY_NUOPC_(rc)

  end subroutine complete_field_in_state

  
  subroutine complete_field(spec, field, rc)
    class(VarSpec), intent(in) :: spec
    type(ESMF_Field), intent(inout) :: field
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    if (allocated(spec%ungridded_lbound) .and. allocated(spec%ungridded_ubound)) then
       call ESMF_FieldEmptyComplete(field, typekind = spec%typekind, &
            gridToFieldMap = spec%grid_to_field_map, &
            ungriddedLBound = spec%ungridded_lbound, &
            ungriddedUBound = spec%ungridded_ubound, rc = rc)
       VERIFY_NUOPC_(rc)
    else
       call ESMF_FieldEmptyComplete(field, spec%typekind, &
            gridToFieldMap = spec%grid_to_field_map, rc = rc)
       VERIFY_NUOPC_(rc)
    end if

    call set_MAPL_field_attributes(spec, field, rc = rc)
    VERIFY_NUOPC_(rc)
    
  end subroutine complete_field


  function create_field(spec, grid, rc) result(field)
    class(VarSpec), intent(in) :: spec
    type(ESMF_Grid), intent(in) :: grid
    integer, intent(out) :: rc

    type(ESMF_Field) :: field

    rc = ESMF_SUCCESS

    if (allocated(spec%ungridded_lbound) .and. allocated(spec%ungridded_ubound)) then
       field = ESMF_FieldCreate(grid, name = spec%name, &
            typeKind = spec%typekind, &
            gridToFieldMap = spec%grid_to_field_map, &
            ungriddedLBound = spec%ungridded_lbound, &
            ungriddedUBound = spec%ungridded_ubound, rc = rc)
       VERIFY_NUOPC_(rc)
    else
       field = ESMF_FieldCreate(grid, spec%typekind, &
            gridToFieldMap = spec%grid_to_field_map, name = spec%name, rc = rc)
       VERIFY_NUOPC_(rc)
    end if

    call set_MAPL_field_attributes(spec, field, rc)
    VERIFY_NUOPC_(rc)

  end function create_field


  subroutine set_MAPL_field_attributes(spec, field, rc)
    type(VarSpec), intent(in) :: spec
    type(ESMF_Field), intent(inout) :: field
    integer, intent(out) :: rc

    integer :: num_levels
    rc = ESMF_SUCCESS

    if (allocated(spec%ungridded_lbound) .and. allocated(spec%ungridded_ubound)) then
       num_levels = spec%ungridded_ubound(1) - spec%ungridded_lbound(1) + 1
       call ESMF_AttributeSet(field, name = "DIMS", value = MAPL_DimsHorzVert, rc = rc)
       VERIFY_NUOPC_(rc)
       if (num_levels == NUM_GEOS_LEVELS) then
          call ESMF_AttributeSet(field, name = "VLOCATION", value = MAPL_VLocationCenter, rc = rc)
          VERIFY_NUOPC_(rc)
       else if(num_levels == NUM_GEOS_LEVELS+1) then
          call ESMF_AttributeSet(field, name = "VLOCATION", value = MAPL_VLocationEdge, rc = rc)
          VERIFY_NUOPC_(rc)
       end if
    else
       call ESMF_AttributeSet(field, name = "DIMS", value = MAPL_DimsHorzOnly, rc = rc)
       VERIFY_NUOPC_(rc)
       call ESMF_AttributeSet(field, name = "VLOCATION", value = MAPL_VLocationNone, rc = rc)
       VERIFY_NUOPC_(rc)
    end if

    call ESMF_AttributeSet(field, name = "LONG_NAME", value = spec%name, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_AttributeSet(field, name = "UNITS", value = "1", rc = rc)
    VERIFY_NUOPC_(rc)

  end subroutine set_MAPL_field_attributes


end module VarSpecMod


module agcm_ctm_mediator
  use ESMF
  use MAPL_Mod
  use NUOPC
  use NUOPC_Mediator, &
       mediator_routine_SS    => SetServices, &
       mediator_label_Advance => label_Advance, &
       mediator_label_DataInitialize  => label_DataInitialize, &
       mediator_label_CheckImport => label_CheckImport, &
       mediator_label_TimestampExport  => label_TimestampExport, &
       mediator_routine_run => routine_run, &
       mediator_label_set_run_clock => label_SetRunClock

  use Mediator_InterpolatorMod
  use VarSpecMod

  use, intrinsic :: iso_fortran_env, only: real32, real64, int64
  
  implicit none

  private

  public SetServices


  type mediator_internal_state
     type(ESMF_State) :: prev_state, regridded_import_state
     type(ESMF_Config) :: config
     logical :: use_interpolation, create_restart, use_regridding ! interpolate, or directly copy states from import to export
     type(mediator_interpolator) :: interpolator
     type(ESMF_Clock) :: slow_clock, fast_clock, regrid_clock
     type(ESMF_Grid) :: agcm_grid, ctm_grid
     class(AbstractRegridder), pointer :: regridder => null()
     integer :: regrid_method = REGRID_METHOD_BILINEAR
  end type mediator_internal_state

  type mediator_internal_state_wrapper
     type(mediator_internal_state), pointer :: ptr
  end type mediator_internal_state_wrapper

  integer :: mpi_rank, de_count

#include "mpif.h"

contains

  subroutine SetServices(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(mediator, mediator_routine_SS, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
         userRoutine = initialize_p0, phase = 0, rc = rc)
    VERIFY_NUOPC_(rc)

    call MPI_Comm_rank(MPI_COMM_WORLD, mpi_rank, rc)


    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=["IPDv05p1"], userRoutine = advertise_fields, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=["IPDv05p5"], userRoutine = modify_decomposition, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_INITIALIZE, &
         phaseLabelList=["IPDv05p6"], userRoutine = realize_fields, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_DataInitialize, &
         specRoutine = initialize_data, rc=rc)
    VERIFY_NUOPC_(rc)



    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_RUN, phaseLabelList = ["interpolate"], &
         userRoutine = mediator_routine_run, rc = rc)
    VERIFY_NUOPC_(rc)

    ! attach specializing method(s)
    call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_advance, &
         specPhaseLabel = "interpolate", specRoutine = interpolate, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_set_run_clock, &
         specPhaseLabel = "interpolate", specRoutine = set_run_clock_fast, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_MethodRemove(mediator, label=mediator_label_CheckImport, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_CheckImport, &
         specRoutine=CheckImport, specPhaseLabel = "interpolate", rc=rc)
    VERIFY_NUOPC_(rc)



    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_RUN, &
         phaseLabelList = ["update_interpolator"], userRoutine = mediator_routine_run, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_advance, &
         specPhaseLabel = "update_interpolator", specRoutine = update_interpolator, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_TimestampExport, &
         specPhaseLabel = "update_interpolator", specRoutine = timestamp_update_interpolator, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_set_run_clock, &
         specPhaseLabel = "update_interpolator", specRoutine = set_run_clock_slow, rc = rc)
    VERIFY_NUOPC_(rc)


    
    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_RUN, phaseLabelList = ["regrid_import_state"], &
         userRoutine = mediator_routine_run, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_advance, &
         specPhaseLabel = "regrid_import_state", specRoutine = regrid_import_state, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_TimestampExport, &
         specPhaseLabel = "regrid_import_state", specRoutine = timestamp_regrid, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_set_run_clock, &
         specPhaseLabel = "regrid_import_state", specRoutine = set_run_clock_regrid, rc = rc)
    VERIFY_NUOPC_(rc)



    call NUOPC_CompSetEntryPoint(mediator, ESMF_METHOD_RUN, phaseLabelList = ["write_restart"], &
         userRoutine = mediator_routine_run, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel=mediator_label_advance, &
         specPhaseLabel = "write_restart", specRoutine = write_restart, rc=rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_TimestampExport, &
         specPhaseLabel = "write_restart", specRoutine = timestamp_regrid, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompSpecialize(mediator, specLabel = mediator_label_set_run_clock, &
         specPhaseLabel = "write_restart", specRoutine = set_run_clock_fast, rc = rc)
    VERIFY_NUOPC_(rc)

    
    call NUOPC_CompSpecialize(mediator, specLabel = label_Finalize, &
         specRoutine = mediator_finalize, rc = rc)
    VERIFY_NUOPC_(rc)


  end subroutine SetServices


  subroutine initialize_p0(mediator, import_state, export_state, clock, rc)
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: import_state, export_state
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(mediator_internal_state_wrapper) :: wrapper

    type(ESMF_TimeInterval) :: fast_interval, interp_interval
    type(ESMF_Time) :: current_time
    character(len=ESMF_MAXSTR) :: grid_type_attr

    type(ESMF_Config) :: cap_config
    integer :: agcm_dt
    integer(int64) :: id

    _UNUSED_DUMMY(import_state)
    _UNUSED_DUMMY(export_state)

    call NUOPC_CompFilterPhaseMap(mediator, ESMF_METHOD_INITIALIZE, &
         acceptStringList=["IPDv05p"], rc=rc)

    allocate(wrapper%ptr)

    wrapper%ptr%prev_state = ESMF_StateCreate(name = "interpolation_fields", rc = rc)
    VERIFY_NUOPC_(rc)

    wrapper%ptr%slow_clock = ESMF_ClockCreate(clock, rc = rc)
    VERIFY_NUOPC_(rc)
    wrapper%ptr%regrid_clock = ESMF_ClockCreate(clock, rc = rc)
    VERIFY_NUOPC_(rc)
    wrapper%ptr%fast_clock = ESMF_ClockCreate(clock, rc = rc)
    VERIFY_NUOPC_(rc)

    cap_config = ESMF_ConfigCreate(rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigLoadFile(cap_config, "CTM_CAP.rc", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigGetAttribute(cap_config, agcm_dt, label = "HEARTBEAT_DT:", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_TimeIntervalSet(fast_interval, s = agcm_dt, rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ClockSet(wrapper%ptr%fast_clock, timeStep = fast_interval, rc = rc)
    VERIFY_NUOPC_(rc)


    call ESMF_ClockGet(clock, currTime = current_time, timeStep = interp_interval, rc = rc)
    VERIFY_NUOPC_(rc)

    wrapper%ptr%interpolator = mediator_interpolator(start_time = current_time, &
         interp_interval = interp_interval)

    wrapper%ptr%config = ESMF_ConfigCreate(rc=rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigLoadFile(wrapper%ptr%config, "NUOPC_run_config.txt", rc=rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(wrapper%ptr%config, wrapper%ptr%use_interpolation, &
         label = "use_time_interpolation:", default = .false., rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(wrapper%ptr%config, wrapper%ptr%create_restart, &
         label = "create_restart:", default = .false., rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(wrapper%ptr%config, wrapper%ptr%use_regridding, &
         label = "use_regridding:", default = .false., rc = rc)
    VERIFY_NUOPC_(rc)


    wrapper%ptr%agcm_grid = make_grid_from_config("GEOSCTM.rc", "AGCM", rc)
    VERIFY_NUOPC_(rc)
    wrapper%ptr%ctm_grid = make_grid_from_config("GEOSCTM.rc", "GEOSctm", rc)
    VERIFY_NUOPC_(rc)

    if (wrapper%ptr%use_regridding) then
       wrapper%ptr%regridded_import_state = ESMF_StateCreate(name = "regridded_import_state", rc = rc)
       VERIFY_NUOPC_(rc)
    end if

    call ESMF_UserCompSetInternalState(mediator, "internal_state", wrapper, rc)
    VERIFY_NUOPC_(rc)

  end subroutine initialize_p0

  !-----------------------------------------------------------------------------

  subroutine advertise_fields(mediator, import_state, export_state, clock, rc)
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: import_state, export_state
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    !type(StringVarSpecIterator) :: iter
    type(mediator_internal_state), pointer :: med_state

    _UNUSED_DUMMY(mediator)
    _UNUSED_DUMMY(clock)

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_SetAttribute(import_state, "FieldTransferPolicy", "transferAll", rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_SetAttribute(export_state, "FieldTransferPolicy", "transferAll", rc = rc)
    VERIFY_NUOPC_(rc)

    ! iter = var_specs%begin()
    ! do while(iter /= var_specs%end())
    !    call NUOPC_Advertise(import_state, standardName = iter%key(), &
    !         transferOfferGeomObject = "cannot provide", rc = rc)
    !    VERIFY_NUOPC_(rc)

    !    call NUOPC_Advertise(export_state, standardName = iter%key(), &
    !         transferOfferGeomObject = "cannot provide", rc = rc)
    !    VERIFY_NUOPC_(rc)
    !    call iter%next()
    ! end do

  end subroutine advertise_fields     



  subroutine modify_decomposition(mediator, import_state, export_state, clock, rc)  
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: import_state, export_state
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(mediator_internal_state), pointer :: med_state

    integer :: i
    type(ESMF_Field) :: import_field, export_field

    character(len=ESMF_MAXSTR), allocatable :: field_names(:)

    _UNUSED_DUMMY(clock)
    _UNUSED_DUMMY(export_state)

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    if (.not. med_state%use_regridding) then
       return
    end if

    field_names = get_state_names(import_state, rc)
    VERIFY_NUOPC_(rc)

    do i = 1, size(field_names)
       call ESMF_StateGet(import_state, field_names(i), import_field, rc = rc)
       VERIFY_NUOPC_(rc)
       call ESMF_FieldEmptySet(import_field, grid = med_state%agcm_grid, rc = rc)
       VERIFY_NUOPC_(rc)

       call ESMF_StateGet(export_state, field_names(i), export_field, rc = rc)
       VERIFY_NUOPC_(rc)
       call ESMF_FieldEmptySet(export_field, grid = med_state%ctm_grid, rc = rc)
       VERIFY_NUOPC_(rc)
    end do

  end subroutine modify_decomposition

  
  function make_grid_from_config(config_name, comp_name, rc) result(grid)
    use MAPL_GridManagerMod, only: grid_manager
    type(ESMF_Grid) :: grid
    character(len=*), intent(in) :: config_name, comp_name
    integer, intent(out) :: rc

    character(len=*), parameter :: CF_COMPONENT_SEPARATOR = '.'
    integer :: ny, nn
    type(ESMF_Config) :: config
    character(len=ESMF_MAXSTR) :: grid_name
    character(len=2) :: dateline

    config = ESMF_ConfigCreate(rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigLoadFile(config, config_name, rc = rc)
    VERIFY_NUOPC_(rc)

    call MAPL_ConfigPrepend(config, COMP_NAME, CF_COMPONENT_SEPARATOR, 'NX:', rc = rc)
    VERIFY_NUOPC_(rc)
    call MAPL_ConfigPrepend(config, COMP_NAME, CF_COMPONENT_SEPARATOR, 'NY:', rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(config, grid_name, label = COMP_NAME//CF_COMPONENT_SEPARATOR//'GRIDNAME:', rc = rc)
    VERIFY_NUOPC_(rc)
    nn = len_trim(grid_name)
    dateline = grid_name(nn-1:nn)
    if (dateline == 'CF') then
       call ESMF_ConfigGetAttribute(config, ny, label = COMP_NAME//CF_COMPONENT_SEPARATOR//'NY:', rc = rc)
       VERIFY_NUOPC_(rc)
       call MAPL_ConfigSetAttribute(config, value = ny/6, label = COMP_NAME//CF_COMPONENT_SEPARATOR//'NY:', rc = rc)
       VERIFY_NUOPC_(rc)
    end if

    grid = grid_manager%make_grid(config, prefix = COMP_NAME//CF_COMPONENT_SEPARATOR, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigDestroy(config, rc = rc)
    VERIFY_NUOPC_(rc)

  contains

    subroutine MAPL_ConfigPrepend(cf,prefix,separator,label,rc)
      type(ESMF_Config), intent(inout) :: cf
      character(len=*) , intent(in   ) :: prefix
      character(len=*) , intent(in   ) :: separator
      character(len=*) , intent(in   ) :: label
      integer, optional , intent(out  ) :: rc

      integer  :: status
      character(len=ESMF_MAXSTR) :: Iam = "MAPL_ConfigPrepend"
      integer  :: val

      call ESMF_ConfigGetAttribute( cf, val, label=trim(prefix)//trim(separator)//trim(label), rc = status )
      if (status /= ESMF_SUCCESS) then
         call ESMF_ConfigGetAttribute(CF,val,label=trim(label),rc=status)
         _VERIFY(status)
         call MAPL_ConfigSetAttribute(CF, val, label=trim(prefix)//trim(separator)//trim(label),rc=status)
         _VERIFY(status)
      end if

      _RETURN(ESMF_SUCCESS)

    end subroutine MAPL_ConfigPrepend

  end function make_grid_from_config
  

  subroutine realize_fields(mediator, import_state, export_state, clock, rc)
    type(ESMF_GridComp)  :: mediator
    type(ESMF_State)     :: import_state, export_state
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    type(mediator_internal_state), pointer :: med_state
    character(len=ESMF_MAXSTR), allocatable :: export_names(:)
    integer :: i
    type(VarSpec) :: export_spec, import_spec
    type(ESMF_Grid) :: import_grid, export_grid
    type(ESMF_Field) :: interpolation_field, regridded_field, export_field, import_field

    _UNUSED_DUMMY(clock)

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    export_names = get_state_names(export_state, rc)
    VERIFY_NUOPC_(rc)

    export_grid = get_grid_from_state(export_state, rc)
    VERIFY_NUOPC_(rc)
    import_grid = get_grid_from_state(import_state, rc)
    VERIFY_NUOPC_(rc)

    do i = 1, size(export_names)

       call ESMF_StateGet(export_state, export_names(i), export_field, rc = rc)
       VERIFY_NUOPC_(rc)
       call ESMF_StateGet(export_state, export_names(i), import_field, rc = rc)
       VERIFY_NUOPC_(rc)
       
       export_spec = get_var_spec_from_transferred_field(export_field)
       import_spec = get_var_spec_from_transferred_field(import_field)

       call import_spec%complete_field_in_state(import_state, rc)
       VERIFY_NUOPC_(rc)
       call export_spec%complete_field_in_state(export_state, rc)
       VERIFY_NUOPC_(rc)

       if (med_state%use_interpolation) then
          interpolation_field = export_spec%create_field(med_state%ctm_grid, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_StateAdd(med_state%prev_state, [interpolation_field], rc = rc)
          VERIFY_NUOPC_(rc)
       end if

       if (med_state%use_regridding) then
          regridded_field = export_spec%create_field(med_state%ctm_grid, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_StateAdd(med_state%regridded_import_state, [regridded_field], rc = rc)
          VERIFY_NUOPC_(rc)
       end if

    end do

    if (med_state%use_regridding) then       
       med_state%regridder => new_regridder_manager%make_regridder(&
            grid_in = med_state%agcm_grid, &
            grid_out = med_state%ctm_grid, &
            regrid_method = med_state%regrid_method, rc = rc)
       VERIFY_NUOPC_(rc)
    end if

  end subroutine realize_fields


  subroutine set_run_clock_slow(mediator, rc)
    type(ESMF_GridComp) :: mediator
    integer, intent(out) :: rc

    type(mediator_internal_state), pointer :: med_state
    type(ESMF_Clock) :: driver_clock

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_MediatorGet(mediator, driverClock=driver_clock, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompSet(mediator, clock = med_state%slow_clock, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompCheckSetClock(mediator, driver_clock, rc = rc)
    VERIFY_NUOPC_(rc)

  end subroutine set_run_clock_slow


  subroutine set_run_clock_regrid(mediator, rc)
    type(ESMF_GridComp) :: mediator
    integer, intent(out) :: rc

    type(mediator_internal_state), pointer :: med_state
    type(ESMF_Clock) :: driver_clock

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_MediatorGet(mediator, driverClock=driver_clock, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompSet(mediator, clock = med_state%regrid_clock, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompCheckSetClock(mediator, driver_clock, rc = rc)
    VERIFY_NUOPC_(rc)

  end subroutine set_run_clock_regrid


  subroutine set_run_clock_fast(mediator, rc)
    type(ESMF_GridComp) :: mediator
    integer, intent(out) :: rc

    type(mediator_internal_state), pointer :: med_state
    type(ESMF_Clock) :: driver_clock

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_MediatorGet(mediator, driverClock=driver_clock, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompSet(mediator, clock = med_state%fast_clock, rc = rc)
    VERIFY_NUOPC_(rc)

    call NUOPC_CompCheckSetClock(mediator, driver_clock, rc = rc)
    VERIFY_NUOPC_(rc)

  end subroutine set_run_clock_fast


  subroutine initialize_data(mediator, rc)
    type(ESMF_GridComp) :: mediator
    integer, intent(out) :: rc

    type(ESMF_State) :: import_state, export_state
    type(ESMF_Clock) :: clock

    type(ESMF_VM) :: vm
    type(ArrDescr) :: descriptor
    integer :: comm, im_world, jm_world, lm_world, nx, ny, num_readers, num_writers, is, ie, js, je
    type(ESMF_Grid) :: grid

    type(mediator_internal_state), pointer :: med_state

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompGet(mediator, clock = clock, importState = import_state, &
         exportState = export_state, vm = vm, rc = rc)
    VERIFY_NUOPC_(rc)


    call ESMF_VMGet(vm, mpicommunicator = comm, rc = rc)
    VERIFY_NUOPC_(rc)

    grid = get_grid_from_state(export_state, rc)
    VERIFY_NUOPC_(rc)

    call MAPL_grid_interior(grid, is, ie, js, je)

    num_writers = 1
    num_readers = 1

    call get_decomp(nx, ny, im_world, jm_world, lm_world)

    call ArrDescrInit(descriptor, comm, im_world, jm_world, lm_world, nx, ny, &
         num_readers, num_writers, is, ie, js, je, rc = rc)
    VERIFY_NUOPC_(rc)

    descriptor%grid = grid
    descriptor%tile = .false.

    if (med_state%use_interpolation) then
       call MAPL_VarReadNCPar(filename = "mediator_import_rst", &
            state = med_state%prev_state, arrdes = descriptor, bootstrapable = .false., rc = rc)
       VERIFY_NUOPC_(rc)
    end if

    call NUOPC_CompAttributeSet(mediator, &
         name="InitializeDataComplete", value="true", rc=rc)
    VERIFY_NUOPC_(rc)

  end subroutine initialize_data

 

  subroutine write_restart(mediator, rc)
    type(ESMF_GridComp) :: mediator
    integer, intent(out) :: rc

    type(ESMF_State) :: import_state, export_state
    type(ESMF_Clock) :: clock, write_clock

    type(ESMF_VM) :: vm
    type(ArrDescr) :: descriptor
    integer :: comm, im_world, jm_world, lm_world, nx, ny, num_readers, num_writers, is, ie, js, je
    type(ESMF_Grid) :: grid

    type(mediator_internal_state), pointer :: med_state

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompGet(mediator, clock = clock, importState = import_state, &
         exportState = export_state, vm = vm, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_VMGet(vm, mpicommunicator = comm, rc = rc)
    VERIFY_NUOPC_(rc)

    grid = get_grid_from_state(export_state, rc)
    VERIFY_NUOPC_(rc)

    call MAPL_grid_interior(grid, is, ie, js, je)

    num_writers = 1
    num_readers = 1

    call get_decomp(nx, ny, im_world, jm_world, lm_world)

    call ArrDescrInit(descriptor, comm, im_world, jm_world, lm_world, nx, ny, &
         num_readers, num_writers, is, ie, js, je, rc = rc)
    VERIFY_NUOPC_(rc)

    descriptor%grid = grid
    descriptor%tile = .false.

    write_clock = ESMF_ClockCreate(clock)
    call ESMF_ClockAdvance(write_clock)

    if (MAPL_AM_I_ROOT()) then
       print *, "Writting mediator resetart"
    end if
    
    if (.not. med_state%use_regridding) then
       call MAPL_VarWriteNCPar(filename = "mediator_import_rst", &
            state = import_state, arrdes = descriptor, clock = write_clock, rc = rc)
       VERIFY_NUOPC_(rc)
    else
       call MAPL_VarWriteNCPar(filename = "mediator_import_rst", &
            state = med_state%regridded_import_state, arrdes = descriptor, clock = write_clock, rc = rc)
       VERIFY_NUOPC_(rc)
    end if


  end subroutine write_restart


  subroutine get_decomp(nx, ny, im, jm, lm)
    integer, intent(out) :: nx, ny, im, jm, lm
    integer :: rc

    type(ESMF_Config) :: config

    config = ESMF_ConfigCreate(rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigLoadFile(config, "GEOSCTM.rc", rc=rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(config, nx, label = "NX:", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigGetAttribute(config, ny, label = "NY:", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigGetAttribute(config, im, label = "GEOSCTM_IM:", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigGetAttribute(config, jm, label = "GEOSCTM_JM:", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigGetAttribute(config, lm, label = "GEOSCTM_LM:", rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigDestroy(config, rc = rc)
    VERIFY_NUOPC_(rc)

  end subroutine get_decomp


  subroutine CheckImport(model, rc)
    type(ESMF_GridComp)   :: model
    integer, intent(out)  :: rc

    _UNUSED_DUMMY(model)

    rc = ESMF_SUCCESS

  end subroutine CheckImport


  subroutine timestamp_update_interpolator(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    _UNUSED_DUMMY(mediator)
    rc = ESMF_SUCCESS

  end subroutine timestamp_update_interpolator


  subroutine timestamp_regrid(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    !type(ESMF_Clock) :: mediator_clock
    !type(ESMF_State) :: import_state, export_state
    !type(ESMF_TimeInterval) :: time_step

    _UNUSED_DUMMY(mediator)
    rc = ESMF_SUCCESS
    
    VERIFY_NUOPC_(rc)
    
  end subroutine timestamp_regrid


  subroutine update_interpolator(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    type(mediator_internal_state), pointer :: med_state

    type(ESMF_State) :: import_state, export_state

    rc = ESMF_SUCCESS

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    _ASSERT(med_state%use_interpolation, "Can't use setup_interpolation if have use_interpolation is not enabled")

    call ESMF_GridCompGet(mediator, importState = import_state, &
         exportState = export_state, rc = rc)
    VERIFY_NUOPC_(rc)

    call med_state%interpolator%update()

    if (.not. med_state%use_regridding) then
       call copy_state(import_state, med_state%prev_state)
    else
       !call regrid_state(med_state%regridder, src_state = import_state, dst_state = med_state%prev_state, rc = rc)
       call copy_state(med_state%regridded_import_state, med_state%prev_state)
       VERIFY_NUOPC_(rc)
    end if

  end subroutine update_interpolator


  subroutine check_fields_and_get_rank_typekind(field1, field2, rank, typekind, rc)
    type(ESMF_Field), intent(in) :: field1, field2
    integer, intent(out) :: rank
    type(ESMF_TypeKind_Flag), intent(out) :: typekind
    integer, intent(out) :: rc

    type(ESMF_TypeKind_Flag) :: typekind1, typekind2
    integer :: rank1, rank2

    rc = ESMF_SUCCESS

    call ESMF_FieldGet(field1, typeKind = typekind1, rank = rank1, rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_FieldGet(field2, typeKind = typekind2, rank = rank2, rc = rc)
    VERIFY_NUOPC_(rc)

    _ASSERT(rank1 == rank2, "Field Ranks must match!")
    _ASSERT(typekind1 == typekind2, "Fields must have same typekind")

    rank = rank1
    typekind = typekind1    
  end subroutine check_fields_and_get_rank_typekind


  function get_state_names(state, rc) result(names)
    type(ESMF_State), intent(in) :: state
    character(len=ESMF_MAXSTR), allocatable :: names(:)
    integer, intent(out) :: rc
    integer :: num_fields

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemCount = num_fields, rc = rc)
    VERIFY_NUOPC_(rc)

    allocate(names(num_fields))

    call ESMF_StateGet(state, itemNameList = names, rc = rc)
    VERIFY_NUOPC_(rc)

  end function get_state_names


  subroutine interpolate(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    type(ESMF_Clock)              :: mediator_clock
    type(ESMF_State)              :: import_state, export_state

    type(mediator_internal_state), pointer :: med_state
    character(len=ESMF_MAXSTR), allocatable :: field_names(:)
    integer :: i

    type(ESMF_Field) :: import_field, export_field, prev_field
    type(ESMF_Time) :: current_time

    rc = ESMF_SUCCESS

    call NUOPC_MediatorGet(mediator, mediatorClock = mediator_clock, &
         importState = import_state, exportState = export_state, rc = rc)

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ClockGet(mediator_clock, currTime = current_time, rc = rc)
    VERIFY_NUOPC_(rc)

    field_names = get_state_names(import_state, rc)
    VERIFY_NUOPC_(rc)

    
    if (med_state%use_interpolation) then
       do i = 1, size(field_names)

          ! choose whether to use the import state directly or the re-gridded import state
          if (.not. med_state%use_regridding) then
             call ESMF_StateGet(import_state, field_names(i), import_field, rc = rc)
             VERIFY_NUOPC_(rc)
          else
             call ESMF_StateGet(med_state%regridded_import_state, field_names(i), import_field, rc = rc)
             VERIFY_NUOPC_(rc)
          end if       

          call ESMF_StateGet(med_state%prev_state, field_names(i), prev_field, rc = rc)
          VERIFY_NUOPC_(rc)

          call ESMF_StateGet(export_state, field_names(i), export_field, rc = rc)
          VERIFY_NUOPC_(rc)

          call med_state%interpolator%interpolate(current_time = current_time, &
               in_field0 = prev_field, in_field1 = import_field, out_field = export_field)
       end do
    else
       ! copy data from agcm to ctm without interpolation
       ! re-grid if they are on different grids
       if (.not. med_state%use_regridding) then
          call copy_state(import_state, export_state)
       else
          call regrid_state(med_state%regridder, import_state, export_state, rc)
          VERIFY_NUOPC_(rc)
       end if
    end if


  end subroutine interpolate


  subroutine regrid_import_state(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    type(ESMF_State)              :: import_state, export_state
    type(mediator_internal_state), pointer :: med_state

    call NUOPC_MediatorGet(mediator, &
         importState = import_state, exportState = export_state, rc = rc)

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    _ASSERT(med_state%use_regridding, "Re-gridding must be enabled using mediator label regrid_import_state")

    call regrid_state(med_state%regridder, import_state, med_state%regridded_import_state, rc)
    VERIFY_NUOPC_(rc)
  end subroutine regrid_import_state
  

  subroutine mediator_finalize(mediator, rc)
    type(ESMF_GridComp)  :: mediator
    integer, intent(out) :: rc

    type(mediator_internal_state), pointer :: med_state

    med_state => get_internal_state(mediator, rc)
    VERIFY_NUOPC_(rc)

    if (med_state%create_restart) then
       call write_restart(mediator, rc)
       VERIFY_NUOPC_(rc)
    end if

  end subroutine mediator_finalize


  subroutine regrid_state(regridder, src_state, dst_state, rc)
    class(AbstractRegridder), intent(inout) :: regridder
    type(ESMF_State), intent(in) :: src_state
    type(ESMF_State), intent(inout) :: dst_state
    integer, intent(out) :: rc

    integer :: i
    character(len=ESMF_MAXSTR), allocatable :: field_names(:)
    type(ESMF_Field) :: src_field, dst_field

    type(ESMF_Field) :: x_src_field, y_src_field, x_dst_field, y_dst_field

    character(len=ESMF_MAXSTR), allocatable :: vector_components(:)
    vector_components = ["U   ", "V   ", "U10M", "V10M", "U10N", "V10N", "MFX ", "MFY ", "CX  ", "CY  "]

    rc = ESMF_SUCCESS

    field_names = get_state_names(src_state, rc)
    VERIFY_NUOPC_(rc)

    do i = 1, size(field_names)      

       if (.not. is_vector_component(field_names(i), vector_components)) then

          call ESMF_StateGet(src_state, field_names(i), src_field, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_StateGet(dst_state, field_names(i), dst_field, rc = rc)
          VERIFY_NUOPC_(rc)

          call regrid_field_scalar(regridder, src_field, dst_field, rc)
          VERIFY_NUOPC_(rc)
       end if

    end do

    do i = 1, size(vector_components), 2
       call ESMF_StateGet(src_state, trim(vector_components(i)), x_src_field, rc = rc)
       VERIFY_NUOPC_(rc)
       call ESMF_StateGet(src_state, trim(vector_components(i+1)), y_src_field, rc = rc)
       VERIFY_NUOPC_(rc)

       call ESMF_StateGet(dst_state, trim(vector_components(i)), x_dst_field, rc = rc)
       VERIFY_NUOPC_(rc)
       call ESMF_StateGet(dst_state, trim(vector_components(i+1)), y_dst_field, rc = rc)
       VERIFY_NUOPC_(rc)

       call regrid_field_vector(regridder, x_src_field, y_src_field, x_dst_field, y_dst_field, rc = rc)
       VERIFY_NUOPC_(rc)
    end do

  end subroutine regrid_state

  
  pure logical function is_vector_component(name, vector_components) 
    character(len=*), intent(in) :: name
    character(len=ESMF_MAXSTR), intent(in) :: vector_components(:)
    integer :: i

    is_vector_component = .false.

    do i = 1, size(vector_components)
       if (trim(name) == trim(vector_components(i))) then
          is_vector_component = .true.
          exit
       end if       
     end do
    
  end function is_vector_component
  

  subroutine regrid_field_scalar(regridder, src_field, dst_field, rc)
    class(AbstractRegridder), intent(inout) :: regridder
    type(ESMF_Field), intent(in) :: src_field
    type(ESMF_Field), intent(inout) :: dst_field
    integer, intent(out) :: rc

    type(ESMF_TypeKind_Flag) :: typekind

    real(real32), pointer :: src_2d_r4(:,:), dst_2d_r4(:,:)
    real(real32), pointer :: src_3d_r4(:,:,:), dst_3d_r4(:,:,:)

    real(real64), pointer :: src_2d_r8(:,:), dst_2d_r8(:,:)
    real(real64), pointer :: src_3d_r8(:,:,:), dst_3d_r8(:,:,:)

    integer :: rank

    call check_fields_and_get_rank_typekind(src_field, dst_field, rank, typekind, rc)
    VERIFY_NUOPC_(rc)

    select case(rank)
    case(2)
       if (typekind == ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(src_field, farrayPtr = src_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = dst_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(src_2d_r4, dst_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
       else if(typekind == ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(src_field, farrayPtr = src_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = dst_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(src_2d_r8, dst_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
       else
          _ASSERT(.false., "Typekind must be real32 or real64")
       end if
    case(3)
       if (typekind == ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(src_field, farrayPtr = src_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = dst_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(src_3d_r4, dst_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
       else if(typekind == ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(src_field, farrayPtr = src_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = dst_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(src_3d_r8, dst_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
       else
          _ASSERT(.false., "Typekind must be real32 or real64")
       end if
    case default
       _ASSERT(.false., "Field rank must be 2 or 3")
    end select

  end subroutine regrid_field_scalar



  subroutine regrid_field_vector(regridder, x_src_field, y_src_field, x_dst_field, y_dst_field, rc)
    class(AbstractRegridder), intent(inout) :: regridder
    type(ESMF_Field), intent(in) :: x_src_field, y_src_field
    type(ESMF_Field), intent(inout) :: x_dst_field, y_dst_field
    integer, intent(out) :: rc

    type(ESMF_TypeKind_Flag) :: typekind

    real(real32), pointer :: x_src_2d_r4(:,:), y_src_2d_r4(:,:), x_dst_2d_r4(:,:), y_dst_2d_r4(:,:)
    real(real32), pointer :: x_src_3d_r4(:,:,:), y_src_3d_r4(:,:,:), x_dst_3d_r4(:,:,:), y_dst_3d_r4(:,:,:)

    real(real64), pointer :: x_src_2d_r8(:,:), y_src_2d_r8(:,:), x_dst_2d_r8(:,:), y_dst_2d_r8(:,:)
    real(real64), pointer :: x_src_3d_r8(:,:,:), y_src_3d_r8(:,:,:), x_dst_3d_r8(:,:,:), y_dst_3d_r8(:,:,:)

    integer :: rank

    call check_fields_and_get_rank_typekind(x_src_field, x_dst_field, rank, typekind, rc)
    VERIFY_NUOPC_(rc)

    select case(rank)
    case(2)
       if (typekind == ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(x_src_field, farrayPtr = x_src_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_src_field, farrayPtr = y_src_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)

          call ESMF_FieldGet(x_dst_field, farrayPtr = x_dst_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_dst_field, farrayPtr = y_dst_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(x_src_2d_r4, y_src_2d_r4, x_dst_2d_r4, y_dst_2d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
       else if(typekind == ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(x_src_field, farrayPtr = x_src_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_src_field, farrayPtr = y_src_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)

          call ESMF_FieldGet(x_dst_field, farrayPtr = x_dst_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_dst_field, farrayPtr = y_dst_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(x_src_2d_r8, y_src_2d_r8, x_dst_2d_r8, y_dst_2d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
       else
          _ASSERT(.false., "Typekind must be real32 or real64")
       end if
    case(3)
       if (typekind == ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(x_src_field, farrayPtr = x_src_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_src_field, farrayPtr = y_src_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)

          call ESMF_FieldGet(x_dst_field, farrayPtr = x_dst_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_dst_field, farrayPtr = y_dst_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(x_src_3d_r4, y_src_3d_r4, x_dst_3d_r4, y_dst_3d_r4, rc = rc)
          VERIFY_NUOPC_(rc)
       else if(typekind == ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(x_src_field, farrayPtr = x_src_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_src_field, farrayPtr = y_src_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)

          call ESMF_FieldGet(x_dst_field, farrayPtr = x_dst_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(y_dst_field, farrayPtr = y_dst_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)

          call regridder%regrid(x_src_3d_r8, y_src_3d_r8, x_dst_3d_r8, y_dst_3d_r8, rc = rc)
          VERIFY_NUOPC_(rc)
       else
          _ASSERT(.false., "Typekind must be real32 or real64")
       end if
    case default
       _ASSERT(.false., "Field rank must be 2 or 3")
    end select

  end subroutine regrid_field_vector


  ! copy field data from state1 to state2
  subroutine copy_state(src_state, dst_state)
    type(ESMF_State) :: src_state, dst_state

    integer :: i, rc
    character(len=ESMF_MAXSTR), allocatable :: field_names(:)
    type(ESMF_Field) :: src_field, dst_field

    field_names = get_state_names(src_state, rc)
    VERIFY_NUOPC_(rc)

    do i = 1, size(field_names)
       call ESMF_StateGet(src_state, field_names(i), src_field, rc = rc)
       VERIFY_NUOPC_(rc)

       call ESMF_StateGet(dst_state, field_names(i), dst_field, rc = rc)
       VERIFY_NUOPC_(rc)

       call copy_field(src_field, dst_field, rc)
       VERIFY_NUOPC_(rc)
    end do

  end subroutine copy_state


  subroutine copy_field(src_field, dst_field, rc)
    type(ESMF_Field) :: src_field, dst_field
    integer, intent(out) :: rc
    integer :: field_rank_in, field_rank_out
    type(ESMF_TypeKind_Flag) :: typekind_in, typekind_out

    real, pointer :: ptr2d_r4_in(:,:), ptr2d_r4_out(:,:)
    real, pointer :: ptr3d_r4_in(:,:,:), ptr3d_r4_out(:,:,:)

    real(real64), pointer :: ptr2d_r8_in(:,:), ptr2d_r8_out(:,:)
    real(real64), pointer :: ptr3d_r8_in(:,:,:), ptr3d_r8_out(:,:,:)

    rc = ESMF_SUCCESS
    ! add check to make sure that fields match
    call ESMF_FieldGet(src_field, typeKind = typekind_in, rank = field_rank_in, rc = rc)
    call ESMF_FieldGet(dst_field, typeKind = typekind_out, rank = field_rank_out, rc = rc)

    _ASSERT(field_rank_in == field_rank_out, "Field ranks must match")

    select case(field_rank_in)
    case(2)
       if (typekind_in == ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(src_field, farrayPtr = ptr2d_r4_in, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = ptr2d_r4_out, rc = rc)
          VERIFY_NUOPC_(rc)
          ptr2d_r4_out = ptr2d_r4_in
       else if (typekind_in == ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(src_field, farrayPtr = ptr2d_r8_in, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = ptr2d_r8_out, rc = rc)
          VERIFY_NUOPC_(rc)
          ptr2d_r8_out = ptr2d_r8_in
       else
          print *, "unknown 2d typekind in mediator"
          stop
       end if
    case(3)
       if (typekind_in == ESMF_TYPEKIND_R4) then
          call ESMF_FieldGet(src_field, farrayPtr = ptr3d_r4_in, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = ptr3d_r4_out, rc = rc)
          VERIFY_NUOPC_(rc)
          ptr3d_r4_out = ptr3d_r4_in
       else if (typekind_in == ESMF_TYPEKIND_R8) then
          call ESMF_FieldGet(src_field, farrayPtr = ptr3d_r8_in, rc = rc)
          VERIFY_NUOPC_(rc)
          call ESMF_FieldGet(dst_field, farrayPtr = ptr3d_r8_out, rc = rc)
          VERIFY_NUOPC_(rc)
          ptr3d_r8_out = ptr3d_r8_in
       else
          print *, "unknown type kind in mediator"
          stop
       end if
    end select
  end subroutine copy_field


  function get_internal_state(mediator, rc) result(state)
    type(ESMF_GridComp), intent(in) :: mediator
    integer, intent(out) :: rc
    type(mediator_internal_state), pointer :: state
    type(mediator_internal_state_wrapper) :: wrapper

    call ESMF_UserCompGetInternalState(mediator, "internal_state", wrapper, rc)
    VERIFY_NUOPC_(rc)

    state => wrapper%ptr
  end function get_internal_state


  function get_grid_from_state(state, rc) result(grid)
    type(ESMF_State), intent(in) :: state
    integer, intent(out) :: rc
    type(ESMF_Grid) :: grid

    integer :: item_count
    character(len=ESMF_MAXSTR), allocatable :: item_names(:)
    type(ESMF_Field) :: field

    rc = ESMF_SUCCESS

    call ESMF_StateGet(state, itemCount = item_count, rc = rc)
    VERIFY_NUOPC_(rc)

    _ASSERT(item_count /= 0, "State has no fields to get state from")

    allocate(item_names(item_count))

    call ESMF_StateGet(state, itemnamelist = item_names, rc = rc)
    VERIFY_NUOPC_(rc)

    ! take the  grid from the first item in the state
    call ESMF_StateGet(state, item_names(1), field, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_FieldGet(field, grid = grid, rc = rc)
    VERIFY_NUOPC_(rc)

  end function get_grid_from_state

end module agcm_ctm_mediator

