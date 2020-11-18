#include "MAPL_ErrLog.h"
#include "NUOPC_ErrLog.h"

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
  integer, parameter :: importFieldCount = 22
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "inst_pres_interface                  ", &
      "inst_pres_levels                     ", &
      "inst_geop_interface                  ", &
      "inst_geop_levels                     ", &
      "inst_temp_levels                     ", &
      "inst_zonal_wind_levels               ", &
      "inst_merid_wind_levels               ", &
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
      "ice_fraction                         ", &
      "lake_fraction                        ", &
      "ocean_fraction                       "  &
    /)
  ! -- export fields
  integer, parameter :: exportFieldCount = 1
  character(len=*), dimension(exportFieldCount), parameter :: &
    exportFieldNames = (/ &
      "inst_tracer_mass_frac                "  &
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

    integer :: status

    ! begin
    rc = ESMF_SUCCESS

    ! the NUOPC model component will register the generic methods
    call NUOPC_CompDerive(model, model_routine_SS, rc=status)
    VERIFY_NUOPC_(status)

    ! switch to IPDv03
    call ESMF_GridCompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      userRoutine=ModelInitializeP0, phase=0, rc=status)
    VERIFY_NUOPC_(status)

    ! set entry point for methods that require specific implementation
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
      phaseLabelList=(/"IPDv03p1"/), userRoutine=ModelInitializeP1, rc=status)
    VERIFY_NUOPC_(status)

    ! - set component's initialization status
    call NUOPC_CompSpecialize(model, &
      specLabel=model_label_DataInitialize, specRoutine=ModelDataInitialize, &
      rc=status)
    VERIFY_NUOPC_(status)

    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specRoutine=ModelAdvance, rc=status)
    VERIFY_NUOPC_(status)

    ! - finalize method
    call NUOPC_CompSpecialize(model, specLabel=model_label_Finalize, &
     specRoutine=ModelFinalize, rc=status)
    VERIFY_NUOPC_(status)

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

    integer :: status

    ! begin
    rc = ESMF_SUCCESS

   ! get component's info
    call NUOPC_CompGet(model, name=name, verbosity=verbosity, rc=status)
    VERIFY_NUOPC_(status)

    ! intro
    call NUOPC_LogIntro(name, rName, verbosity, rc=status)
    VERIFY_NUOPC_(status)

    ! switch to IPDv01 by filtering all other phaseMap entries
    call NUOPC_CompFilterPhaseMap(model, ESMF_METHOD_INITIALIZE, &
      acceptStringList=(/"IPDv03p"/), rc=status)
    VERIFY_NUOPC_(status)

    ! extro
    call NUOPC_LogExtro(name, rName, verbosity, rc=status)
    VERIFY_NUOPC_(status)

  end subroutine ModelInitializeP0

  subroutine ModelInitializeP1(model, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: model
    type(ESMF_State)     :: importState, exportState
    type(ESMF_Clock)     :: clock
    integer, intent(out) :: rc

    integer :: status

    ! begin
    rc = ESMF_SUCCESS

    ! -- advertise imported fields
    if (importFieldCount > 0) then
      call NUOPC_Advertise(importState, importFieldNames, &
        TransferOfferGeomObject="cannot provide", &
        SharePolicyField="share", &
        rc=status)
      VERIFY_NUOPC_(status)
    end if

    ! -- advertise exported fields
    if (exportFieldCount > 0) then
      call NUOPC_Advertise(exportState, exportFieldNames, &
        TransferOfferGeomObject="cannot provide", &
        SharePolicyField="share", &
        rc=status)
      VERIFY_NUOPC_(status)
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
      importState=importState, exportState=exportState, rc=status)
    VERIFY_NUOPC_(status)

    ! retrieve member list from import state, if any
    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, &
      nestedFlag=.true., rc=status)
    VERIFY_NUOPC_(status)

    ! select PETs carrying data payloads from imported fields
    nlev = 1
    if (associated(fieldList)) then
      do item = 1, size(fieldList)
        call ESMF_FieldGet(fieldList(item), rank=rank, localDeCount=localDeCount, rc=status)
        VERIFY_NUOPC_(status)
        if (localDeCount /= 1) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="localDeCount must be 1", &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)
        end if
        if (rank == 4) then
          call ESMF_FieldGet(fieldList(item), grid=grid, &
            ungriddedLBound=lb, ungriddedUBound=ub, rc=status)
          VERIFY_NUOPC_(status)
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
    call ESMF_GridCompGet(model, vm=vm, rc=status)
    VERIFY_NUOPC_(status)

    call ESMF_VMGet(vm, localPet=localPet, mpiCommunicator=modelComm, rc=status)
    VERIFY_NUOPC_(status)

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

    ! add public API to MAPL_CapGridCompMod to set grid and clock?
    call this % cap_gc % set_grid(grid, lm=nlev, _RC)

    call this % cap_gc % initialize(_RC)

    ! set component's internal state
    call ESMF_GridCompSetInternalState(model, is, status)
    VERIFY_NUOPC_(status)

    ! -> set InitializeDataComplete Component Attribute to "true", indicating
    ! to the driver that this Component has fully initialized its data
    call NUOPC_CompAttributeSet(model, &
      name="InitializeDataComplete", value="true", rc=status)
    VERIFY_NUOPC_(status)

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
    call NUOPC_CompGet(model, name=name, diagnostic=diagnostic, rc=status)
    VERIFY_NUOPC_(status)

    ! query the Component for its clock, importState and exportState
    call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
      exportState=exportState, rc=status)
    VERIFY_NUOPC_(status)

    ! get component's internal state
    nullify(is % maplCap)
    call ESMF_GridCompGetInternalState(model, is, status)
    VERIFY_NUOPC_(status)

    ! print tinestep details
    call ESMF_ClockPrint(clock, options="currTime", &
      preString="------>Advancing Aerosol from: ", rc=status)
    VERIFY_NUOPC_(status)

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=status)
    VERIFY_NUOPC_(status)

    call ESMF_TimePrint(currTime + timeStep, &
      preString="---------------------> to: ", rc=status)
    VERIFY_NUOPC_(status)

    ! advance MAPL component
    if (associated(is % maplCap)) then
      ! -- import
      call AerosolStateUpdate(model, is % maplCap, "import", rc=status)
      VERIFY_NUOPC_(status)

      ! -- run cap
      call is % maplCap % cap_gc % run(_RC)

      ! -- export tarcers
      call AerosolStateUpdate(model, is % maplCap, "export", rc=status)
      VERIFY_NUOPC_(status)
    end if

    ! print field diagnostics
    if (btest(diagnostic,17)) then
      call AerosolFieldDiagnostics(model, rc=status)
      VERIFY_NUOPC_(status)

      call MAPLFieldDiagnostics(model, is % maplCap % cap_gc % import_state, "import", rc=status)
      VERIFY_NUOPC_(status)

      call MAPLFieldDiagnostics(model, is % maplCap % cap_gc % export_state, "export", rc=status)
      VERIFY_NUOPC_(status)
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
    call ESMF_GridCompGetInternalState(model, is, status)
    VERIFY_NUOPC_(status)

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
