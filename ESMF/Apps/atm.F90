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

module ATM

!-----------------------------------------------------------------------------
! ATM Component.
!-----------------------------------------------------------------------------

    use ESMF
    use NUOPC
    use NUOPC_Model, &
            model_routine_SS    => SetServices, &
            model_label_Advance => label_Advance

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

        print*,"ATM: ", localPet, " start SetServices"

        ! the NUOPC model component will register the generic methods
        print*, "ATM: ", localPet, " register Generic SetServices"
        call NUOPC_CompDerive(model, model_routine_SS, rc=rc)
        VERIFY_NUOPC_(rc)

        ! set entry point for methods that require specific implementation
        print*, "ATM: ", localPet, " attach InitializeP1"
        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
                phaseLabelList=(/"IPDv00p1"/), userRoutine=InitializeP1, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "ATM: ", localPet, " attach InitializeP2"
        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE, &
                phaseLabelList=(/"IPDv00p2"/), userRoutine=InitializeP2, rc=rc)
        VERIFY_NUOPC_(rc)

        ! attach specializing method(s)
        print*, "ATM: ", localPet, " set ModelAdvance"
        call NUOPC_CompSpecialize(model, specLabel=model_label_Advance, &
                specRoutine=ModelAdvance, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"ATM: ", localPet, " finish SetServices"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine InitializeP1(model, importState, exportState, clock, rc)
        type(ESMF_GridComp)  :: model
        type(ESMF_State)     :: importState, exportState
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        type(ESMF_VM) :: vm
        integer       :: localPet

        rc = ESMF_SUCCESS

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"ATM: ", localPet, " start InitializeP1"

        ! Disabling the following macro, e.g. renaming to WITHIMPORTFIELDS_disable,
        ! will result in a model component that does not advertise any importable
        ! Fields. Use this if you want to drive the model independently.
#define WITHIMPORTFIELDS
#ifdef WITHIMPORTFIELDS
        ! importable field: sea_surface_temperature
        print*,"ATM: ", localPet, " advertise IMPORTS"
        call NUOPC_Advertise(importState, &
                StandardName="sea_surface_temperature", name="sst", rc=rc)
        VERIFY_NUOPC_(rc)
#endif

        print*,"ATM: ", localPet, " advertise EXPORTS"
        ! exportable field: air_pressure_at_sea_level
        call NUOPC_Advertise(exportState, &
                StandardName="air_pressure_at_sea_level", name="pmsl", rc=rc)
        VERIFY_NUOPC_(rc)

        ! exportable field: surface_net_downward_shortwave_flux
        call NUOPC_Advertise(exportState, &
                StandardName="surface_net_downward_shortwave_flux", name="rsns", rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"ATM: ", localPet, " finish InitializeP1"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine InitializeP2(model, importState, exportState, clock, rc)
        type(ESMF_GridComp)  :: model
        type(ESMF_State)     :: importState, exportState
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Field)        :: field
        type(ESMF_Grid)         :: gridIn
        type(ESMF_Grid)         :: gridOut

        type(ESMF_VM) :: vm
        integer       :: localPet

        rc = ESMF_SUCCESS

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"ATM: ", localPet, " start InitializeP2"

        print*,"ATM: ", localPet, " create grids"
        ! create a Grid object for Fields
        gridIn = ESMF_GridCreateNoPeriDimUfrm(maxIndex=(/10, 100/), &
                minCornerCoord=(/10._ESMF_KIND_R8, 20._ESMF_KIND_R8/), &
                maxCornerCoord=(/100._ESMF_KIND_R8, 200._ESMF_KIND_R8/), &
                coordSys=ESMF_COORDSYS_CART, staggerLocList=(/ESMF_STAGGERLOC_CENTER/), &
                rc=rc)
        VERIFY_NUOPC_(rc)
        gridOut = gridIn ! for now out same as in

#ifdef WITHIMPORTFIELDS
        print*,"ATM: ", localPet, " create import field: sst"
        ! importable field: sea_surface_temperature
        field = ESMF_FieldCreate(name="sst", grid=gridIn, &
                typekind=ESMF_TYPEKIND_R8, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " realize psml"
        call NUOPC_Realize(importState, field=field, rc=rc)
        VERIFY_NUOPC_(rc)
#endif

        ! exportable field: air_pressure_at_sea_level
#ifdef CREATE_AND_REALIZE
        ! This branch shows the standard procedure of creating a complete field
        ! with Grid and memory allocation, and then calling Realize() for it.
        print*,"ATM: ", localPet, " create export field: pmsl, with allocation"
        field = ESMF_FieldCreate(name="pmsl", grid=gridOut, &
                typekind=ESMF_TYPEKIND_R8, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " realize psml"
        call NUOPC_Realize(exportState, field=field, rc=rc)
        VERIFY_NUOPC_(rc)
#else
        ! This branch shows the alternative way of "realizing" an advertised field.
        ! It accesses the empty field that was created during advertise, and
        ! finishes it, setting a Grid on it, and then calling FieldEmptyComplete().
        ! No formal Realize() is then needed.
        print*,"ATM: ", localPet, " create export field: pmsl, from advertised field"
        call ESMF_StateGet(exportState, field=field, itemName="pmsl", rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " set psml grid"
        call ESMF_FieldEmptySet(field, grid=gridOut, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " complete field pmsl"
        call ESMF_FieldEmptyComplete(field, typekind=ESMF_TYPEKIND_R8, rc=rc)
        VERIFY_NUOPC_(rc)
#define WITH_FORMAL_REALIZE
#ifdef WITH_FORMAL_REALIZE
        ! There is not need to formally call Realize() when completing the
        ! adverised field directly. However, calling Realize() also works.
        print*,"ATM: ", localPet, " realize psml formally"
        call NUOPC_Realize(exportState, field=field, rc=rc)
        VERIFY_NUOPC_(rc)
#endif
#endif
        ! exportable field: surface_net_downward_shortwave_flux
        print*,"ATM: ", localPet, " create export field: rsns"
        field = ESMF_FieldCreate(name="rsns", grid=gridOut, &
                typekind=ESMF_TYPEKIND_R8, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " realize rsns"
        call NUOPC_Realize(exportState, field=field, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"ATM: ", localPet, " finish InitializeP2"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine ModelAdvance(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Clock)            :: clock
        type(ESMF_State)            :: importState, exportState
        character(len=160)          :: msgString

        type(ESMF_Field)      :: field
        real(kind=8), pointer :: ptr_pmsl(:,:), ptr_rsns(:,:)

        type(ESMF_VM) :: vm
        integer       :: localPet

        call ESMF_GridCompGet(model, vm=vm, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, localPet=localPet, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"ATM: ", localPet, " start ModelAdvance"

#define NUOPC_TRACE__OFF
#ifdef NUOPC_TRACE
        print*,"ATM: ", localPet, " enter ESMF tracing Region"
        call ESMF_TraceRegionEnter("ATM:ModelAdvance")
#endif

        rc = ESMF_SUCCESS

        ! query the Component for its clock, importState and exportState
        print*,"ATM: ", localPet, " query GridComp"
        call NUOPC_ModelGet(model, modelClock=clock, importState=importState, &
                exportState=exportState, rc=rc)
        VERIFY_NUOPC_(rc)

        ! HERE THE MODEL ADVANCES: currTime -> currTime + timeStep

        ! Because of the way that the internal Clock was set by default,
        ! its timeStep is equal to the parent timeStep. As a consequence the
        ! currTime + timeStep is equal to the stopTime of the internal Clock
        ! for this call of the ModelAdvance() routine.

        print*,"ATM: ", localPet, " print current time"
        call ESMF_ClockPrint(clock, options="currTime", &
                preString="------>Advancing ATM from: ", unit=msgString, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"ATM: ", localPet, " get the pmsl field"
        call ESMF_StateGet(exportState, field=field, &
                itemName="pmsl", rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " get the pmsl value"
        call ESMF_FieldGet(field, localDE=0, farrayPtr=ptr_pmsl, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " current value of pmsl:", minval(ptr_pmsl), maxval(ptr_pmsl)

        ptr_pmsl = ptr_pmsl + 1.d0
        print*,"ATM: ", localPet, " new value of pmsl:", minval(ptr_pmsl), maxval(ptr_pmsl)

        print*,"ATM: ", localPet, " get the rsns field"
        call ESMF_StateGet(exportState, field=field, &
                itemName="rsns", rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " get the rsns value"
        call ESMF_FieldGet(field, localDE=0, farrayPtr=ptr_rsns, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"ATM: ", localPet, " current value of rsns:", minval(ptr_rsns), maxval(ptr_rsns)

        ptr_rsns = ptr_rsns + 10.d0
        print*,"ATM: ", localPet, " new value of rsns:", minval(ptr_rsns), maxval(ptr_rsns)

        print*,"ATM: ", localPet, " print stop (next) time"
        call ESMF_ClockPrint(clock, options="stopTime", &
                preString="---------------------> to: ", unit=msgString, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_LogWrite(msgString, ESMF_LOGMSG_INFO, rc=rc)
        VERIFY_NUOPC_(rc)

#ifdef NUOPC_TRACE
        print*,"ATM: ", localPet, " exit ESMF tracing Region"
        call ESMF_TraceRegionExit("ATM:ModelAdvance")
#endif
        print*,"ATM: ", localPet, " finish ModelAdvance"
    end subroutine
end module
