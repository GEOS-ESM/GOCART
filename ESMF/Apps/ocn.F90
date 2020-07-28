!==============================================================================
! Earth System Modeling Framework
! Copyright 2002-2020, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!==============================================================================
#include "NUOPC_ErrLog.h"

module OCN

!-----------------------------------------------------------------------------
! OCN Component.
!-----------------------------------------------------------------------------

    use ESMF
    use NUOPC
    use NUOPC_Model, &
            model_routine_SS      => SetServices, &
            model_label_SetClock  => label_SetClock, &
            model_label_Advance   => label_Advance

    implicit none
    private

    public SetServices

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

    subroutine SetServices(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        type(ESMF_VM) :: vm
        integer       :: localPet

        rc = ESMF_SUCCESS

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"OCN: ", localPet, " start SetServices"

        ! the NUOPC model component will register the generic methods
        print*, "OCN: ", localPet, " register Generic SetServices"
        call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
        VERIFY_NUOPC_(rc)

        ! set entry point for methods that require specific implementation
        print*, "OCN: ", localPet, " attach InitializeP1"
        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
                phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "OCN: ", localPet, " attach InitializeP2"
        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
                phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
        VERIFY_NUOPC_(rc)

        ! attach specializing method(s)
        print*, "OCN: ", localPet, " set SetClock"
        call NUOPC_CompSpecialize(model, specLabel=model_label_SetClock, &
                specRoutine=SetClock, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "OCN: ", localPet, " set ModelAdvance"
        call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
                specRoutine=ModelAdvance, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"OCN: ", localPet, " finish SetServices"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine InitializeP1(model, importState, exportState, clock, rc)
        type(ESMF_GridComp)  :: model
        type(ESMF_State)     :: importState, exportState
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        type(ESMF_VM) :: vm
        integer       :: localPet

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        rc = ESMF_SUCCESS

        print*,"OCN: ", localPet, " start InitializeP1"

        ! Disabling the following macro, e.g. renaming to WITHIMPORTFIELDS_disable,
        ! will result in a model component that does not advertise any importable
        ! Fields. Use this if you want to drive the model independently.
#define WITHIMPORTFIELDS
#ifdef WITHIMPORTFIELDS
        ! importable field: air_pressure_at_sea_level
        print*,"OCN: ", localPet, " advertise IMPORTS"
        call NUOPC_Advertise(importState, &
                StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
        VERIFY_NUOPC_(rc)

        ! importable field: surface_net_downward_shortwave_flux
        call NUOPC_Advertise(importState, &
                StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=rc)
        VERIFY_NUOPC_(rc)
#endif

        ! exportable field: sea_surface_temperature
        print*,"OCN: ", localPet, " advertise EXPORTS"
        call NUOPC_Advertise(exportState, &
                StandardName="sea_surface_temperature", name="sst", rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"OCN: ", localPet, " finish InitializeP1"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine InitializeP2(model, importState, exportState, clock, rc)
        type(ESMF_GridComp)  :: model
        type(ESMF_State)     :: importState, exportState
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_TimeInterval) :: stabilityTimeStep
        type(ESMF_Field)        :: field
        type(ESMF_Grid)         :: gridIn
        type(ESMF_Grid)         :: gridOut

        type(ESMF_VM) :: vm
        integer       :: localPet

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        rc = ESMF_SUCCESS

        print*,"OCN: ", localPet, " start InitializeP2"

        ! create a Grid object for Fields
        print*,"OCN: ", localPet, " create grids"
        gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/100, 10/), &
                minCornerCoord=(/10._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
                maxCornerCoord=(/100._ESMF_KIND_R8, 200._ESMF_KIND_R8/), &
                coordSys=ESMF_COORDSYS_CART, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
                rc=rc)
        VERIFY_NUOPC_(rc)
        gridOut = gridIn ! for now out same as in

#ifdef WITHIMPORTFIELDS
        ! importable field: air_pressure_at_sea_level
        print*,"OCN: ", localPet, " create import field: pmsl"
        field = ESMF_FieldCreate(name="pmsl", grid=gridIn, &
                typekind=ESMF_TYPEKIND_R8, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"OCN: ", localPet, " realize pmsl"
        call NUOPC_Realize(importState, field=field, rc=rc)
        VERIFY_NUOPC_(rc)

        ! importable field: surface_net_downward_shortwave_flux
        print*,"OCN: ", localPet, " create import field: rsns"
        field = ESMF_FieldCreate(name="rsns", grid=gridIn, &
                typekind=ESMF_TYPEKIND_R8, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"OCN realize rsns"
        call NUOPC_Realize(importState, field=field, rc=rc)
        VERIFY_NUOPC_(rc)
#endif

        ! exportable field: sea_surface_temperature
        print*,"OCN: ", localPet, " create export field: sst"
        field = ESMF_FieldCreate(name="sst", grid=gridOut, &
                typekind=ESMF_TYPEKIND_R8, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"OCN: ", localPet, " realize sst"
        call NUOPC_Realize(exportState, field=field, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"OCN: ", localPet, " finish InitializeP2"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine SetClock(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Clock)              :: clock
        type(ESMF_TimeInterval)       :: stabilityTimeStep

        type(ESMF_VM) :: vm
        integer       :: localPet

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        rc = ESMF_SUCCESS

        print*,"OCN: ", localPet, " start SetClock"

        ! query the Component for its clock, importState and exportState
        print*,"OCN: ", localPet, " query GridComp"
        call NUOPC_ModelGet(model, modelClock=clock, rc=rc)
        VERIFY_NUOPC_(rc)

        ! initialize internal clock
        ! here: parent Clock and stability timeStep determine actual model timeStep
        print*,"OCN: ", localPet, " set time interval"
        call ESMF_TimeIntervalSet(stabilityTimeStep, m=5, rc=rc) ! 5 minute steps
        VERIFY_NUOPC_(rc)
        print*,"OCN: ", localPet, " set clock"
        call NUOPC_CompSetClock(model, clock, stabilityTimeStep, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"OCN: ", localPet, " finish SetClock"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine ModelAdvance(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Clock)            :: clock
        type(ESMF_State)            :: importState, exportState
        type(ESMF_Time)             :: currTime
        type(ESMF_TimeInterval)     :: timeStep
        character(len=160)          :: msgString

        type(ESMF_Field)      :: field
        real(kind=8), pointer :: ptr_pmsl(:,:), ptr_rsns(:,:)

        type(ESMF_VM) :: vm
        integer       :: localPet

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"OCN: ", localPet, " start ModelAdvance"

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
        print*,"OCN: ", localPet, " enter ESMF tracing Region"
        call ESMF_TraceRegionEnter("OCN:ModelAdvance")
#endif

        rc = ESMF_SUCCESS

        ! query the Component for its clock, importState and exportState
        print*,"OCN: ", localPet, " query GridComp"
        call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
                exportState=exportState, rc=rc)
        VERIFY_NUOPC_(rc)

        ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

        ! Because of the way that the internal Clock was set in SetClock(),
        ! its timeStep is likely smaller than the parent timeStep. As a consequence
        ! the time interval covered by a single parent timeStep will result in 
        ! multiple calls to the ModelAdvance() routine. Every time the currTime
        ! will come in by one internal timeStep advanced. This goes until the
        ! stopTime of the internal Clock has been reached.

        print*,"OCN: ", localPet, " print current time"
        call ESMF_ClockPrint(clock, options="currTime", &
                preString="------>Advancing OCN from: ", unit=msgString, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Get: ", localPet, " clock times"
        call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"OCN: ", localPet, " get the pmsl field"
        call ESMF_StateGet(importState, field=field, &
                itemName="pmsl", rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"OCN: ", localPet, " get the pmsl value"
        call ESMF_FieldGet(field, localDE=0, farrayPtr=ptr_pmsl, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"OCN: ", localPet, " current value of pmsl:", minval(ptr_pmsl), maxval(ptr_pmsl)

        print*,"OCN: ", localPet, " get the rsns field"
        call ESMF_StateGet(importState, field=field, &
                itemName="rsns", rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"OCN: ", localPet, " get the rsns value"
        call ESMF_FieldGet(field, localDE=0, farrayPtr=ptr_rsns, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"OCN: ", localPet, " current value of rsns:", minval(ptr_rsns), maxval(ptr_rsns)

        print*,"OCN: ", localPet, " print stop (next) time"
        call ESMF_TimePrint(currTime + timeStep, &
                preString="---------------------> to: ", unit=msgString, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
        VERIFY_NUOPC_(rc)

#ifdef NUOPC_TRACE
        print*,"OCN: ", localPet, " exit ESMF tracing Region"
        call ESMF_TraceRegionExit("OCN:ModelAdvance")
#endif
        print*,"OCN: ", localPet, " finish ModelAdvance"
    end subroutine
end module
