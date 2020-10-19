module Aerosol_Comp_Mod

  use ESMF
  use NUOPC

  implicit none

  private

  public :: AerosolFieldDiagnostics
  public :: AerosolTracerReroute

contains

  subroutine AerosolTracerReroute(model, rc)
    type(ESMF_GridComp)            :: model
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: localDe, localDeCount
    type(ESMF_Field) :: iField, eField
    type(ESMF_State) :: importState, exportState
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: fpImp, fpExp

    logical, save :: first = .true.

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- query the Component for its VM and importState
    call ESMF_GridCompGet(model, importState=importState, &
      exportState=exportState, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- retrieve
    call ESMF_StateGet(importState, "inst_tracer_mass_frac", iField, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(exportState, "inst_tracer_mass_frac", eField, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (first) then
      call ESMF_FieldFill(eField, dataFillScheme="one", const1=1.e-05_ESMF_KIND_R8, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
      first = .false.
      return
    end if

    call ESMF_FieldGet(eField, localDeCount=localDeCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    do localDe = 0, localDeCount - 1
      nullify(fpExp, fpImp)
      call ESMF_FieldGet(eField, localDe=localDe, farrayPtr=fpImp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return  ! bail out
      call ESMF_FieldGet(iField, localDe=localDe, farrayPtr=fpExp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return  ! bail out
      fpExp = fpImp
    end do

  end subroutine AerosolTracerReroute

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

end module Aerosol_Comp_Mod
