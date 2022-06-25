module Aerosol_Shared_Mod

  use ESMF
  use NUOPC
  use NUOPC_Model, only: NUOPC_ModelGet

  implicit none

  ! -- Shared physical constants
  ! gravity (m/s2)
  real(ESMF_KIND_R8), parameter :: con_g = 9.80665e+0_ESMF_KIND_R8
  ! inverse gravity (m/s2)
  real(ESMF_KIND_R8), parameter :: onebg = 1._ESMF_KIND_R8 / con_g

  ! -- Shared methods
  interface AerosolGetPtr
    module procedure AerosolGetPtr2D
    module procedure AerosolGetPtr3D
    module procedure AerosolGetPtr4D
  end interface

  private

  public :: con_g, onebg

  public :: AerosolGetPtr
  public :: AerosolFieldDiagnostics
  public :: AerosolModelGet
  public :: MAPLFieldDiagnostics

contains

  ! -- Field data methods

  subroutine AerosolGetPtr2D(state, fieldName, farrayPtr, localDe, rc)
    type(ESMF_State)                :: state
    character(len=*),   intent(in)  :: fieldName
    real(ESMF_KIND_R8), pointer     :: farrayPtr(:,:)
    integer, optional,  intent(in)  :: localDe
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer          :: localrc
    type(ESMF_Field) :: field

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    nullify(farrayPtr)

    call ESMF_StateGet(state, fieldName, field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=farrayPtr, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolGetPtr2D

  subroutine AerosolGetPtr3D(state, fieldName, farrayPtr, localDe, rc)
    type(ESMF_State)                :: state
    character(len=*),   intent(in)  :: fieldName
    real(ESMF_KIND_R8), pointer     :: farrayPtr(:,:,:)
    integer, optional,  intent(in)  :: localDe
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer          :: localrc
    type(ESMF_Field) :: field

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    nullify(farrayPtr)

    call ESMF_StateGet(state, fieldName, field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=farrayPtr, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolGetPtr3D

  subroutine AerosolGetPtr4D(state, fieldName, farrayPtr, localDe, rc)
    type(ESMF_State)                :: state
    character(len=*),   intent(in)  :: fieldName
    real(ESMF_KIND_R8), pointer     :: farrayPtr(:,:,:,:)
    integer, optional,  intent(in)  :: localDe
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer          :: localrc
    type(ESMF_Field) :: field

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    nullify(farrayPtr)

    call ESMF_StateGet(state, fieldName, field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=farrayPtr, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolGetPtr4D

  ! -- Diagnostic methods

  subroutine AerosolFieldDiagnostics(model, rc)
    type(ESMF_GridComp)            :: model
    integer, optional, intent(out) :: rc

    ! -- local variables
    logical :: isConnected
    integer :: localrc
    integer :: item
    integer :: localDe, localDeCount, rank
    integer :: fieldCount, maxLength
    character(len=ESMF_MAXSTR)          :: msgString
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR), pointer :: connectedList(:)
    character(len=ESMF_MAXSTR), pointer :: fieldNames(:)
    real(ESMF_KIND_R8)                              :: minValue, maxValue
    real(ESMF_KIND_R8), dimension(:),       pointer :: localMin, localMax, globalMin, globalMax
    real(ESMF_KIND_R8), dimension(:),       pointer :: fp1d
    real(ESMF_KIND_R8), dimension(:,:),     pointer :: fp2d
    real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: fp3d
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: fp4d
    type(ESMF_Field), pointer :: fieldList(:)
    type(ESMF_State)  :: importState
    type(ESMF_VM)     :: vm

    ! -- local parameters
    character(len=*), parameter :: rName = "diagnostics"

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- get component information
    call NUOPC_CompGet(model, name=name, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- query the Component for its VM and importState
    call ESMF_GridCompGet(model, vm=vm, importState=importState, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(fieldNames, connectedList, fieldList)
    call NUOPC_GetStateMemberLists(importState, StandardNameList=fieldNames, &
      ConnectedList=connectedList, fieldList=fieldList, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    isConnected = .false.
    if (associated(fieldNames)) isConnected = any(connectedList == "true")

    nullify(localMin, localMax, globalMin, globalMax)

    ! -- check values of imported fields, if requested
    if (isConnected) then

      fieldCount = size(fieldList)

      allocate(localMin(fieldCount), localMax(fieldCount), &
        globalMin(fieldCount), globalMax(fieldCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

      localMin = huge(0._ESMF_KIND_R8)
      localMax = -localMin
      globalMin = 0._ESMF_KIND_R8
      globalMax = 0._ESMF_KIND_R8

      ! -- find longest field name for formatting purpose
      maxLength = 0
      do item = 1, fieldCount
        maxLength = max(maxLength, len_trim(fieldNames(item)))
      end do

      do item = 1, fieldCount
        if (connectedList(item) == "true") then
          ! --- get field data
          call ESMF_FieldGet(fieldList(item), rank=rank, &
            localDeCount=localDeCount, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out

          do localDe = 0, localDeCount - 1
            minValue = localMin(item)
            maxValue = localMax(item)
            select case(rank)
              case(1)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp1d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp1d)
                maxValue = maxval(fp1d)
              case(2)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp2d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp2d)
                maxValue = maxval(fp2d)
              case(3)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp3d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp3d)
                maxValue = maxval(fp3d)
              case(4)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp4d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp4d)
                maxValue = maxval(fp4d)
              case default
                call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="Field rank not implemented.", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return ! bail out
            end select
            localMin(item) = min(minValue, localMin(item))
            localMax(item) = max(maxValue, localMax(item))
          end do
        end if
      end do

      ! -- compute global min and max values
      call ESMF_VMAllReduce(vm, localMin, globalMin, fieldCount, ESMF_REDUCE_MIN, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call ESMF_VMAllReduce(vm, localMax, globalMax, fieldCount, ESMF_REDUCE_MAX, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      ! -- log results
      do item = 1, fieldCount
        if (connectedList(item) == "true") then
          write(msgString,'(a,": ",a,": ",a,"[",i0,"]: local/global min/max =",4g20.8)') &
            trim(name), rName, fieldNames(item)(1:maxLength), item, localMin(item), &
            localMax(item), globalMin(item), globalMax(item)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
        end if
      end do

      deallocate(localMin, localMax, globalMin, globalMax, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

      if (associated(fieldNames)) then
        deallocate(fieldNames, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if
      if (associated(connectedList)) then
        deallocate(connectedList, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if
      if (associated(fieldList)) then
        deallocate(fieldList, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if
    end if

  end subroutine AerosolFieldDiagnostics

  subroutine MAPLFieldDiagnostics(model, state, label, rc)
    type(ESMF_GridComp)            :: model
    type(ESMF_State)               :: state
    character(len=*),  intent(in)  :: label
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, item
    integer :: localDe, localDeCount, rank
    integer :: fieldCount, maxLength
    character(len=ESMF_MAXSTR)          :: msgString
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR), pointer :: shapeList(:)
    character(len=ESMF_MAXSTR), pointer :: fieldNames(:)
    real(ESMF_KIND_R4)                              :: minValue, maxValue
    real(ESMF_KIND_R4), dimension(:),       pointer :: localMin, localMax, globalMin, globalMax
    real(ESMF_KIND_R4), dimension(:),       pointer :: fp1d
    real(ESMF_KIND_R4), dimension(:,:),     pointer :: fp2d
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer :: fp3d
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: fp4d
    type(ESMF_Field) :: field
    type(ESMF_VM)    :: vm
    type(ESMF_FieldStatus_Flag) :: fieldStatus
    type(ESMF_StateItem_Flag), pointer :: itemTypeList(:)

    ! -- local parameters
    character(len=*), parameter :: rName = "MAPL: diagnostic"

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- get component information
    call NUOPC_CompGet(model, name=name, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- query the Component for its VM and importState
    call ESMF_GridCompGet(model, vm=vm, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(state, itemCount=fieldCount, nestedFlag=.true., rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    write(msgString,'(a,": ",a,": ",a,": found ",i0," items")') &
      trim(name), rName, trim(label), fieldCount
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    allocate(fieldNames(fieldCount), itemTypeList(fieldCount), shapeList(fieldCount))
    fieldNames = ""
    shapeList = ""

    call ESMF_StateGet(state, itemNameList=fieldNames, &
      itemTypeList=itemTypeList, nestedFlag=.true., rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(localMin, localMax, globalMin, globalMax)

    ! -- check values of imported fields, if requested
    if (fieldCount > 0) then

      allocate(localMin(fieldCount), localMax(fieldCount), &
        globalMin(fieldCount), globalMax(fieldCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

      localMin = huge(0._ESMF_KIND_R4)
      localMax = -localMin
      globalMin = 0._ESMF_KIND_R4
      globalMax = 0._ESMF_KIND_R4

      ! -- find longest field name for formatting purpose
      maxLength = 0
      do item = 1, fieldCount
        maxLength = max(maxLength, len_trim(fieldNames(item)))
      end do

      do item = 1, fieldCount
        if (itemTypeList(item) == ESMF_STATEITEM_FIELD) then
          ! --- get field data
          call ESMF_StateGet(state, trim(fieldNames(item)), field, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          localDeCount = 0
          rank = 0
          call ESMF_FieldGet(field, status=fieldStatus, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          if (fieldStatus == ESMF_FIELDSTATUS_COMPLETE) then
            call ESMF_FieldGet(field, rank=rank, &
              localDeCount=localDeCount, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
          else
            shapeList(item) = "incomplete"
          end if

          do localDe = 0, localDeCount - 1
            minValue = localMin(item)
            maxValue = localMax(item)
            select case(rank)
              case(1)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp1d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp1d)
                maxValue = maxval(fp1d)
                write(shapeList(item), '(i0,":",i0)') lbound(fp1d,1), ubound(fp1d,1)
              case(2)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp2d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp2d)
                maxValue = maxval(fp2d)
                write(shapeList(item), '(i0,":",i0,",",i0,":",i0)') (lbound(fp2d,i), ubound(fp2d,i), i=1,2)
              case(3)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp3d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp3d)
                maxValue = maxval(fp3d)
                write(shapeList(item), '(i0,":",i0,2(",",i0,":",i0))') (lbound(fp3d,i), ubound(fp3d,i), i=1,3)
              case(4)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp4d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp4d)
                maxValue = maxval(fp4d)
                write(shapeList(item), '(i0,":",i0,3(",",i0,":",i0))') (lbound(fp4d,i), ubound(fp4d,i), i=1,4)
              case default
                call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="Field rank not implemented.", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return ! bail out
            end select
            localMin(item) = min(minValue, localMin(item))
            localMax(item) = max(maxValue, localMax(item))
          end do
        end if
      end do

      ! -- compute global min and max values
      call ESMF_VMAllReduce(vm, localMin, globalMin, fieldCount, ESMF_REDUCE_MIN, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call ESMF_VMAllReduce(vm, localMax, globalMax, fieldCount, ESMF_REDUCE_MAX, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      ! -- log results
      do item = 1, fieldCount
        if (itemTypeList(item) == ESMF_STATEITEM_FIELD) then
          write(msgString,'(a,": ",a,": ",a,"[",i0,"]: [",a,"] local/global min/max =",4g20.8)') &
            trim(name), rName, fieldNames(item)(1:maxLength), item, trim(shapeList(item)), &
            localMin(item), localMax(item), globalMin(item), globalMax(item)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
        end if
      end do

      deallocate(localMin, localMax, globalMin, globalMax, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

    if (associated(fieldNames)) then
      deallocate(fieldNames, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if
    if (associated(itemTypeList)) then
      deallocate(itemTypeList, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if
    if (associated(shapeList)) then
      deallocate(shapeList, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

  end subroutine MAPLFieldDiagnostics

  ! -- Metadata methods

  subroutine AerosolModelGet(model, grid, numLevels, tracerInfo, rc)

    ! -- arguments
    type(ESMF_GridComp)                     :: model
    type(ESMF_Grid),  optional, intent(out) :: grid
    integer,          optional, intent(out) :: numLevels
    type(ESMF_Info),  optional, intent(out) :: tracerInfo
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer          :: localrc, stat
    integer          :: item, rank, localDeCount
    integer, dimension(2) :: lb, ub
    type(ESMF_Array) :: array
    type(ESMF_State) :: importState
    type(ESMF_Field), pointer :: fieldList(:)
    type(ESMF_Info)  :: info

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- initialize output variables
    if (present(numLevels))  numLevels  = 0

    ! retrieve import state from gridded component
    call NUOPC_ModelGet(model, importState=importState, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! retrieve member list from import state, if any
    nullify(fieldList)
    call NUOPC_GetStateMemberLists(importState, fieldList=fieldList, &
      nestedFlag=.true., rc=rc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! retrieve number of vertical levels from imported fields
    if (associated(fieldList)) then
      do item = 1, size(fieldList)
        call ESMF_FieldGet(fieldList(item), rank=rank, localDeCount=localDeCount, rc=rc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
        ! -- validate field data decomposition
        if (localDeCount /= 1) then
          call ESMF_LogSetError(ESMF_RC_INTNRL_BAD, msg="localDeCount must be 1", &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)
        end if
        if (rank == 4) then
          call ESMF_FieldGet(fieldList(item), array=array, grid=grid, &
            ungriddedLBound=lb, ungriddedUBound=ub, rc=rc)
          if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__)) &
            return  ! bail out
          ! -- populate remaining output arguments
          if (present(numLevels)) numLevels = ub(1) - lb(1) + 1
          if (present(tracerInfo)) then
            call ESMF_InfoGetFromHost(array, tracerInfo, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__)) &
              return  ! bail out
          end if
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

  end subroutine AerosolModelGet

end module Aerosol_Shared_Mod
