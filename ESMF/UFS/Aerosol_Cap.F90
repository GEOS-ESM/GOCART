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
  integer, parameter :: importFieldCount = 26
  character(len=*), dimension(importFieldCount), parameter :: &
    importFieldNames = (/ &
      "inst_pres_interface                  ", &
      "inst_pres_levels                     ", &
      "inst_geop_interface                  ", &
      "inst_geop_levels                     ", &
      "inst_temp_levels                     ", &
      "inst_zonal_wind_levels               ", &
      "inst_merid_wind_levels               ", &
      "inst_omega_levels                    ", &
      "inst_tracer_mass_frac                ", &
      "soil_type                            ", &
      "inst_pbl_height                      ", &
      "surface_cell_area                    ", &
      "inst_convective_rainfall_amount      ", &
      "inst_exchange_coefficient_heat_levels", &
      "inst_spec_humid_conv_tendency_levels ", &
      "inst_friction_velocity               ", &
      "inst_rainfall_amount                 ", &
      "inst_soil_moisture_content           ", &
      "inst_down_sw_flx                     ", &
      "inst_land_sea_mask                   ", &
      "inst_temp_height_surface             ", &
      "inst_up_sensi_heat_flx               ", &
      "inst_lwe_snow_thickness              ", &
      "vegetation_type                      ", &
      "inst_vegetation_area_frac            ", &
      "inst_surface_roughness               "  &
    /)
  ! -- export fields
  integer, parameter :: exportFieldCount = 5
  character(len=*), dimension(exportFieldCount), parameter :: &
    exportFieldNames = (/ &
      "inst_tracer_mass_frac                ", &
      "inst_tracer_up_surface_flx           ", &
      "inst_tracer_down_surface_flx         ", &
      "inst_tracer_clmn_mass_dens           ", &
      "inst_tracer_anth_biom_flx            "  &
    /)

  ! -- MAPL field map
  integer, parameter :: importMapSize = 3
  character(len=*), dimension(importMapSize, 2), parameter :: &
    importFieldMap = reshape((/ &
      "inst_temp_levels                     ", "T                                    ", &
      "inst_zonal_wind_levels               ", "U                                    ", &
      "inst_merid_wind_levels               ", "V                                    "  &
      /), (/importMapSize, 2/))

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


    integer :: item
    integer                   :: counts(ESMF_MAXDIM)
    real(ESMF_KIND_R8), pointer :: r8d2(:,:)

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
    localDeCount = 1
    if (associated(fieldList)) then
      if (size(fieldList) > 0) then
        call ESMF_FieldGet(fieldList(1), grid=grid, localDeCount=localDeCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__)) &
          return  ! bail out
      end if
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
    print *,'maplCap: n_members, npes_member, npes = ',maplCapOptions % npes_model, &
      maplCapOptions % n_members ; flush 6

    call this % initialize_cap_gc(this % get_mapl_comm())

    call this % cap_gc % set_services(_RC)

    ! add public API to MAPL_CapGridCompMod to set grid and clock?
    call this % cap_gc % set(grid=grid, _RC)

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
      ! set imports
!     do item = 1, importMapSize
      do item = 1, 0
        call is % maplCap % cap_gc % get_field_from_state(trim(importFieldMap(item,2)), &
          'Cap_Imports', ESMF_STATEINTENT_IMPORT, ofield, _RC)
        call ESMF_StateGet(importState, trim(importFieldMap(item,1)), ifield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
        call ESMF_FieldCopy(ofield, ifield, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__)) &
          return  ! bail out
      end do
      call is % maplCap % cap_gc % run(_RC)
      ! get exports
    end if

    ! print field diagnostics
    if (btest(diagnostic,17)) then
      call AerosolFieldDiagnostics(model, rc=rc)
      if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__)) &
        return  ! bail out
    end if

   ! temporarily reroute tracers (no update)
   call AerosolTracerReroute(model, rc=rc)
   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
     line=__LINE__, &
     file=__FILE__)) &
     return  ! bail out

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
