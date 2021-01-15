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

  use Aerosol_GridComp_mod, only : &
    aerosol_routine_SS          => SetServices

  use MAPL

  use mpi

  implicit none

  ! -- import fields
  integer, parameter :: importFieldCount = 25
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "inst_pres_interface                  ", &
      "inst_pres_levels                     ", &
      "inst_geop_interface                  ", &
      "inst_geop_levels                     ", &
      "inst_temp_levels                     ", &
      "inst_zonal_wind_levels               ", &
      "inst_merid_wind_levels               ", &
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
      "inst_sensi_heat_flx                  ", &
      "inst_surface_roughness               ", &
      "inst_surface_soil_wetness            ", &
      "ice_fraction                         ", &
      "lake_fraction                        ", &
      "ocean_fraction                       "  &
    /)
  ! -- export fields
  integer, parameter :: exportFieldCount = 3
  character(len=*), dimension(exportFieldCount), parameter :: &
    exportFieldNames = (/ &
      "inst_tracer_mass_frac                ", &
      "inst_tracer_up_surface_flx           ", &
      "inst_tracer_down_surface_flx         "  &
    /)

  type Aerosol_InternalState_T
    type(MAPL_Cap), pointer :: maplCap
  end type Aerosol_InternalState_T

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
    if (importFieldCount > 0) then
      call NUOPC_Advertise(importState, importFieldNames, &
        TransferOfferGeomObject="cannot provide", &
        SharePolicyField="share", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

    ! -- advertise exported fields
    if (exportFieldCount > 0) then
      call NUOPC_Advertise(exportState, exportFieldNames, &
        TransferOfferGeomObject="cannot provide", &
        SharePolicyField="share", &
        rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

  end subroutine ModelInitializeP1

  subroutine ModelDataInitialize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer          :: stat, urc
    integer          :: localDeCount
    integer          :: localPet, petCount, activePetCount
    integer          :: maplComm, modelComm
    integer          :: status
    integer          :: lb(2), ub(2)
    integer          :: item, nlev, rank
    integer, dimension(:), allocatable :: petList, recvBuffer
    type(ESMF_Clock) :: clock
    type(ESMF_Field), pointer :: fieldList(:)
    type(ESMF_Grid)  :: grid
    type(ESMF_State) :: importState, exportState
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_VM)    :: vm
    type(Aerosol_InternalState_T) :: is
    type(MAPL_Cap), pointer :: this
    type(MAPL_CapOptions) :: maplCapOptions

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

    ! select PETs carrying data payloads from imported fields
    nlev = 1
    if (associated(fieldList)) then
      do item = 1, size(fieldList)
        call ESMF_FieldGet(fieldList(item), rank=rank, localDeCount=localDeCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
        if (localDeCount /= 1) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="localDeCount must be 1", &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)
        end if
        if (rank == 4) then
          call ESMF_FieldGet(fieldList(item), grid=grid, &
            ungriddedLBound=lb, ungriddedUBound=ub, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
          nlev = ub(1) - lb(1) + 1
        end if
      end do
      deallocate(fieldList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Unable to deallocate internal memory", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
      nullify(fieldList)
    end if

    ! select PETs carrying data payloads from imported fields
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=modelComm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    nullify(is % maplCap)

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
    allocate(is % maplCap, stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Unable to allocate internal memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    is % maplCap = MAPL_Cap("MAPL Driver", aerosol_routine_SS, cap_options=maplCapOptions, _RC)

    this => is % maplCap

    call this % nuopc_fill_mapl_comm(_RC)

    call this % initialize_cap_gc(this % get_mapl_comm())

    call this % cap_gc % set_services(_RC)

    ! provide model grid to MAPL
    call this % cap_gc % set_grid(grid, lm=nlev, _RC)

    ! provide model clock to MAPL
    call this % cap_gc % set_clock(clock, _RC)

    ! initialize aerosol grid component
    call this % cap_gc % initialize(_RC)

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

  end subroutine ModelDataInitialize

  subroutine ModelAdvance(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer                       :: status
    integer                       :: diagnostic
    integer                       :: item
    character(len=ESMF_MAXSTR)    :: name
    type(ESMF_Field)              :: ifield, ofield
    type(ESMF_State)              :: importState
    type(ESMF_State)              :: exportState
    type(ESMF_Clock)              :: clock
    type(ESMF_Time)               :: currTime
    type(ESMF_TimeInterval)       :: timeStep
    type(Aerosol_InternalState_T) :: is

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

    ! get component's internal state
    nullify(is % maplCap)
    call ESMF_GridCompGetInternalState(model, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! print tinestep details
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing Aerosol from: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_TimePrint(currTime + timeStep, &
      preString="---------------------> to: ", rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__)) &
      return  ! bail out

    ! advance MAPL component
    if (associated(is % maplCap)) then
      ! -- import
      call AerosolStateUpdate(model, is % maplCap, "import", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out

      ! -- run cap
      call is % maplCap % cap_gc % run(_RC)

      ! -- export tarcers
      call AerosolStateUpdate(model, is % maplCap, "export", rc=rc)
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
      call MAPLFieldDiagnostics(model, is % maplCap % cap_gc % import_state, "import", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
      call MAPLFieldDiagnostics(model, is % maplCap % cap_gc % export_state, "export", rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

  end subroutine ModelAdvance

  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer :: status
    type(Aerosol_InternalState_T) :: is

    ! begin
    rc = ESMF_SUCCESS

    ! get component's internal state
    nullify(is % maplCap)
    call ESMF_GridCompGetInternalState(model, is, rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    if (associated(is % maplCap)) then

      ! finalize MAPL component
      call is % maplCap % cap_gc % finalize(_RC)

      ! finalize I/O
      call is % maplCap % finalize_io_clients_servers(_RC)

      ! finalize MAPL framework
      call MAPL_Finalize(comm=is % maplCap % get_comm_world(), _RC)

      ! finally, deallocate internal state
      deallocate(is % maplCap, stat=status)
      if (ESMF_LogFoundDeallocError(statusToCheck=status, &
        msg="Failed to free internal state memory.", &
        line=__LINE__,  &
        file=__FILE__)) &
        return  ! bail out
      nullify(is % maplCap)
    end if

  end subroutine ModelFinalize

end module Aerosol_Cap
