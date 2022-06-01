#include "MAPL_ErrLog.h"

module Aerosol_Cap

  use ESMF
  use NUOPC
  use NUOPC_Model, only : &
    NUOPC_ModelGet, &
    model_routine_SS            => SetServices,          &
    model_routine_Run           => routine_Run,          &
    model_label_Advance         => label_Advance,        &
    model_label_DataInitialize  => label_DataInitialize, &
    model_label_Finalize        => label_Finalize

  use Aerosol_Comp_mod
  use Aerosol_Diag_mod
  use Aerosol_Internal_mod
  use Aerosol_Logger_mod
  use Aerosol_Shared_mod
  use Aerosol_Tracer_mod

  use Aerosol_GridComp_mod, only : &
    aerosol_routine_SS          => SetServices

  use MAPL

  use mpi

  implicit none

  ! -- model name
  character(len=*), parameter :: modelName = "UFS Aerosols"

  ! -- import fields
  character(len=*), dimension(*), parameter :: &
    importFieldNames = [ &
      "inst_pres_interface                  ", &
      "inst_pres_levels                     ", &
      "inst_geop_interface                  ", &
      "inst_geop_levels                     ", &
      "inst_temp_levels                     ", &
      "inst_zonal_wind_levels               ", &
      "inst_merid_wind_levels               ", &
      "inst_cloud_frac_levels               ", &
      "inst_ice_nonconv_tendency_levels     ", &
      "inst_liq_nonconv_tendency_levels     ", &
      "inst_tracer_mass_frac                ", &
      "inst_pbl_height                      ", &
      "surface_cell_area                    ", &
      "inst_convective_rainfall_amount      ", &
      "inst_rainfall_amount                 ", &
      "inst_zonal_wind_height10m            ", &
      "inst_merid_wind_height10m            ", &
      "inst_friction_velocity               ", &
      "inst_land_sea_mask                   ", &
      "inst_temp_height_surface             ", &
      "inst_up_sensi_heat_flx               ", &
      "inst_surface_roughness               ", &
      "inst_surface_soil_wetness            ", &
      "inst_soil_moisture_content           ", &
      "ice_fraction_in_atm                  ", &
      "lake_fraction                        ", &
      "ocean_fraction                       ", &
      "surface_snow_area_fraction           "  &
    ]
  ! -- export fields
  character(len=*), dimension(*), parameter :: &
    exportFieldNames = [ &
      "inst_tracer_mass_frac                ", &
      "inst_tracer_up_surface_flx           ", &
      "inst_tracer_down_surface_flx         "  &
    ]


  private

  public :: SetServices

contains

  subroutine SetServices(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! begin
    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! switch to IPDv03
    call ESMF_GridCompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      userRoutine=ModelInitializeP0, phase=0, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=ModelInitializeP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! - set component's initialization status
    call NUOPC_CompSpecialize(model, &
      specLabel=model_label_DataInitialize, specRoutine=ModelDataInitialize, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! - finalize method
    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
     specRoutine=ModelFinalize, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine SetServices

  subroutine ModelInitializeP0(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState
    type(ESMF_State)     :: exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! local variables
    integer                    :: verbosity
    character(len=ESMF_MAXSTR) :: name

    ! local parameters
    character(len=*), parameter :: rName = "ModelInitializeP0"

    ! begin
    rc = ESMF_SUCCESS

    ! startup
    call AerosolLog(modelName//': Initializing ...', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

   ! get component's info
    call NUOPC_CompGet(model, name=name, verbosity=verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(model, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelInitializeP0

  subroutine ModelInitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    ! begin
    rc = ESMF_SUCCESS

    ! -- advertise imported fields
    call NUOPC_Advertise(importState, importFieldNames, &
      TransferOfferGeomObject="cannot provide", &
      SharePolicyField="share", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! -- advertise exported fields
    call NUOPC_Advertise(exportState, exportFieldNames, &
      TransferOfferGeomObject="cannot provide", &
      SharePolicyField="share", &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelInitializeP1

  subroutine ModelDataInitialize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer          :: stat, urc
    integer          :: localDeCount
    integer          :: modelComm
    integer          :: status
    integer          :: lb(2), ub(2)
    integer          :: item, nlev, rank
    type(ESMF_Clock) :: clock
    type(ESMF_Grid)  :: grid
    type(ESMF_State) :: importState, exportState
    type(ESMF_VM)    :: vm
    type(ESMF_Info)  :: tracerInfo
    type(ESMF_Field), pointer :: fieldList(:)
    type(MAPL_Cap),   pointer :: cap
    type(MAPL_CapOptions)     :: maplCapOptions
    type(Aerosol_InternalState_T) :: is
    type(Aerosol_Tracer_T), pointer :: trp

    ! begin
    rc = ESMF_SUCCESS

    ! retrieve import/export states and clock from NUOPC component and
    ! pass them to the GOCART component
    call NUOPC_ModelGet(model, modelClock=clock, &
      importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! retrieve member list from import state, if any
    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, &
      nestedFlag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! retrieve model's metadata from imported tracer field
    call AerosolModelGet(model, grid=grid, numLevels=nlev, tracerInfo=tracerInfo, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg="Unable to retrieve model grid and metadata ", &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! retrieve model's MPI communicator
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, mpiCommunicator=modelComm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! allocate memory for the internal state and store it into component
    allocate(is % wrap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Unable to allocate internal state", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(is % wrap % maplCap)
    nullify(is % wrap % tracers)

    ! set MAPL options
    maplCapOptions = MAPL_CapOptions(_RC)

    maplCapOptions % use_comm_world = .false.
    maplCapOptions % comm           = modelComm
    maplCapOptions % logging_config = ""

    call MPI_Comm_size(modelComm, maplCapOptions % npes_model, urc)
    if (urc /= MPI_SUCCESS) then
      call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! startup MAPL
    allocate(is % wrap % maplCap, stat=stat, source = MAPL_Cap("MAPL Driver", aerosol_routine_SS, cap_options=maplCapOptions, rc=rc))
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Unable to allocate MAPL cap", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    cap => is % wrap % maplCap

    ! initialize MAPL I/O on the same PETs as the model component
    call cap % initialize_io_clients_servers(cap % get_comm_world(), _RC)

    ! create aerosol grid component
    call cap % initialize_cap_gc(_RC)

    call cap % cap_gc % set_services(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! provide model grid to MAPL
    call cap % cap_gc % set_grid(grid, lm=nlev, _RC)

    ! provide model clock to MAPL
    call cap % cap_gc % set_clock(clock, _RC)

    ! initialize aerosol grid component
    call cap % cap_gc % initialize(rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! create maps linking imported aerosol tracers to MAPL fields
    ! initialize tracer datatype containing tracer names and mapping information
    allocate(is % wrap % tracers, stat=stat, source = AerosolTracer(cap % get_cap_rc_file(), tracerInfo, rc=rc))
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Unable to allocate tracers data structure", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    trp => is % wrap % tracers

    ! print tracer maps
    call AerosolTracerPrint(trp, 'GOCART2G Tracer Map')

    ! set component's internal state
    call ESMF_GridCompSetInternalState(model, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(model, &
      name="InitializeDataComplete", value="true", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! log end of model's initialize
    call AerosolLog(modelName//': Model initialized', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelDataInitialize

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer                       :: diagnostic
    character(len=ESMF_MAXSTR)    :: name
    character(len=ESMF_MAXSTR)    :: msgString, timeString
    type(ESMF_State)              :: importState
    type(ESMF_State)              :: exportState
    type(ESMF_Clock)              :: clock
    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep
    type(MAPL_Cap),       pointer :: cap
    type(Aerosol_InternalState_T) :: is
    type(Aerosol_InternalData_T), pointer :: this

    ! local parameters
    character(len=*), parameter :: rName = "ModelAdvance"

    ! begin
    rc = ESMF_SUCCESS

    ! get component's information
    call NUOPC_CompGet(model, name=name, diagnostic=diagnostic, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! log time step
    call AerosolLogStep(clock, msg=modelName, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(model, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    nullify(cap)
    if (associated(is % wrap)) cap => is % wrap % maplCap

    ! advance MAPL component
    if (associated(cap)) then
      ! -- import
      call AerosolStateUpdate(model, cap, "import", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! -- run cap
      call cap % cap_gc % run(rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! -- export tracers
      call AerosolStateUpdate(model, cap, "export", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! -- export diagnostics
      call AerosolDiagUpdate(model, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    ! print field diagnostics
    if (btest(diagnostic,17)) then
      call AerosolFieldDiagnostics(model, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      if (associated(cap)) then
        call MAPLFieldDiagnostics(model, cap % cap_gc % import_state, "import", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call MAPLFieldDiagnostics(model, cap % cap_gc % export_state, "export", rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end if
    end if

  end subroutine ModelAdvance

  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer :: status
    type(MAPL_Cap),       pointer :: cap
    type(Aerosol_InternalState_T) :: is
    type(Aerosol_InternalData_T), pointer :: this

    ! begin
    rc = ESMF_SUCCESS

    ! finalize
    call AerosolLog(modelName//': Finalizing ...', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! get component's internal state
    call ESMF_GridCompGetInternalState(model, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    nullify(this, cap)
    if (associated(is % wrap)) then
      this => is % wrap
      cap  => this % maplCap
      if (associated(cap)) then
        ! finalize MAPL component
        call cap % cap_gc % finalize(_RC)

        ! finalize I/O
        call cap % finalize_io_clients_servers(_RC)

        ! finalize MAPL framework
        call MAPL_Finalize(comm=cap % get_comm_world(), _RC)

        ! deallocate MAPL Cap
        deallocate(this % maplCap, stat=status)
        if (ESMF_LogFoundDeallocError(statusToCheck=status, &
          msg="Failed to free MAPL Cap.", &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        nullify(this % maplCap, cap)
      end if
      ! finally, deallocate internal state
      deallocate(is % wrap, stat=status)
      if (ESMF_LogFoundDeallocError(statusToCheck=status, &
        msg="Failed to free internal state memory.", &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      nullify(is % wrap, this)
    end if

    ! log completion of model's finalize
    call AerosolLog(modelName//': Model finalized', rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelFinalize

end module Aerosol_Cap
