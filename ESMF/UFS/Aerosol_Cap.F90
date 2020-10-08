module Aerosol_Cap

  use ESMF
  use NUOPC
  use NUOPC_Model, only : &
    NUOPC_ModelGet, &
    model_routine_SS            => SetServices,          &
    model_routine_Run           => routine_Run,          &
    model_label_Advance         => label_Advance,        &
    model_label_CheckImport     => label_CheckImport,    &
    model_label_DataInitialize  => label_DataInitialize, &
    model_label_Finalize        => label_Finalize,       &
    model_label_TimestampExport => label_TimestampExport

  use GOCART2G_GridCompMod, only : &
    gocart_routine_SS           => SetServices

  implicit none


  type(ESMF_GridComp) :: gocartComp

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
!   call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
!     phaseLabelList=(/"IPDv03p1"/), userRoutine=ModelInitializeP1, rc=rc)
!   if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!     line=__LINE__, &
!     file=__FILE__)) &
!     return  ! bail out

    ! - set component's initialization status
    call NUOPC_CompSpecialize(model, &
      specLabel=model_label_DataInitialize, specRoutine=ModelDataInitialize, &
      rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    
    ! attach specializing method(s)
    ! - run methods
    call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_RUN, &
      phaseLabelList=(/ "phase1", "phase2" /), &
      userRoutine=model_routine_Run, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    ! - advance phase 1
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specPhaseLabel="phase1", specRoutine=ModelAdvanceP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    ! - advance phase 2
    call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
      specPhaseLabel="phase2", specRoutine=ModelAdvanceP2, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! support multi-phase run sequences by properly advancing the clock
    ! only during last run phase
    call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
      specPhaseLabel="phase1", specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_TimestampExport, &
      specPhaseLabel="phase1", specRoutine=TimestampExportP1, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

    call NUOPC_CompSpecialize(model, specLabel=model_label_CheckImport, &
      specPhaseLabel="phase2", specRoutine=NUOPC_NoOp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

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

  ! -- no advertise/realize yet

  subroutine ModelDataInitialize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer          :: stat, urc
    integer          :: localDeCount
    integer          :: localPet, petCount, activePetCount
    integer, dimension(:), allocatable :: petList, recvBuffer
    type(ESMF_Clock) :: clock
    type(ESMF_Field), pointer :: fieldList(:)
    type(ESMF_State) :: importState, exportState
    type(ESMF_VM)    :: vm
    

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
    call ESMF_GridCompGet(model, vm=vm, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    call ESMF_VMGet(vm, localPet=localPet, petCount=petCount, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out


    localDeCount = 1
    if (associated(fieldList)) then
      if (size(fieldList) > 0) then
        call ESMF_FieldGet(fieldList(1), localDeCount=localDeCount, rc=rc)
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

    allocate(recvBuffer(petCount), petList(petCount), stat=stat)
    if (ESMF_LogFoundAllocError(statusToCheck=stat, &
      msg="Unable to allocate internal memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    petList = 0
    petList(localPet+1) = localDeCount
    recvBuffer = 0

    call ESMF_VMAllReduce(vm, petList, recvBuffer, petCount, ESMF_REDUCE_SUM, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    activePetCount = 0
    petList  = 0

    do localPet = 0, petCount - 1
      if (recvBuffer(localPet + 1) > 0) then
        activePetCount = activePetCount + 1
        petList(activePetCount) = localPet
      end if
    end do

    ! create GOCART component on same PET list
    gocartComp = ESMF_GridCompCreate(name="GOCART", petList=petList(1:activePetCount), rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

    ! free up memory
    deallocate(recvBuffer, petList, stat=stat)
    if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
      msg="Unable to deallocate internal memory", &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! specialize GOCART component
    call ESMF_GridCompSetServices(gocartComp, gocart_routine_SS, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! initialize GOCART component
    call ESMF_GridCompInitialize(gocartComp, importState=importState, &
      exportState=exportState, clock=clock, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
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

  subroutine ModelAdvanceP1(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer          :: urc
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: clock

    ! local parameters
    character(len=*), parameter :: rName = "ModelAdvanceP1"

    ! begin
    rc = ESMF_SUCCESS

    ! retrieve model states
    call NUOPC_ModelGet(model, modelClock=clock, &
      importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
   
    ! advance GOCART component
    call ESMF_GridCompRun(gocartComp, importState=importState, &
      exportState=exportState, clock=clock, phase=1, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    
  end subroutine ModelAdvanceP1

  subroutine ModelAdvanceP2(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer          :: urc
    type(ESMF_State) :: importState
    type(ESMF_State) :: exportState
    type(ESMF_Clock) :: clock

    ! local parameters
    character(len=*), parameter :: rName = "ModelAdvanceP2"

    ! begin
    rc = ESMF_SUCCESS

    ! retrieve model states
    call NUOPC_ModelGet(model, modelClock=clock, &
      importState=importState, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
   
    ! advance GOCART component
    call ESMF_GridCompRun(gocartComp, importState=importState, &
      exportState=exportState, clock=clock, phase=2, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out
    
  end subroutine ModelAdvanceP2

  subroutine ModelFinalize(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    integer :: urc

    ! begin
    rc = ESMF_SUCCESS

    ! finalize GOCART component
    call ESMF_GridCompFinalize(gocartComp, userRc=urc, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out
    if (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) &
      return  ! bail out

    ! free up memory
    call ESMF_GridCompDestroy(gocartComp, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__)) &
      return  ! bail out

  end subroutine ModelFinalize

  !-----------------------------------------------------------------------------

  subroutine TimestampExportP1(model, rc)
    type(ESMF_GridComp)  :: model
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Clock)     :: driverClock, modelClock
    type(ESMF_State)     :: exportState

    ! begin
    rc = ESMF_SUCCESS

    ! get driver and model clock
    call NUOPC_ModelGet(model, driverClock=driverClock, &
      modelClock=modelClock, exportState=exportState, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

    ! reset model clock to initial time
    call NUOPC_CheckSetClock(modelClock, driverClock, &
      forceCurrTime=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

    ! update timestamp on export Fields
    call NUOPC_SetTimestamp(exportState, modelClock, rc=rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, file=__FILE__)) return  ! bail out

  end subroutine TimestampExportP1

end module Aerosol_Cap
