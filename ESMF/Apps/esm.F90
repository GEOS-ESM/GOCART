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

module ESM

    !-----------------------------------------------------------------------------
    ! Code that specializes generic ESM Component code.
    !-----------------------------------------------------------------------------

    use ESMF
    use NUOPC
    use NUOPC_Driver, &
            driver_routine_SS             => SetServices, &
            driver_label_SetModelServices => label_SetModelServices

    use ATM, only: atmSS => SetServices
    use OCN, only: ocnSS => SetServices

    use NUOPC_Connector, only: cplSS => SetServices

    use NUOPC_MAPLcapMod, &
            NUOPC_CAP_SS => SetServices

    implicit none

    private

    public SetServices

!-----------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------

    subroutine SetServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS

        print*, "ESM start SetServices"

        print*, "ESM register Generic SetServices"
        ! NUOPC_Driver registers the generic methods
        call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
        VERIFY_NUOPC_(rc)

        ! attach specializing method(s)
        print*, "ESM attach SetModelServices"
        call NUOPC_CompSpecialize(driver, &
                specLabel=driver_label_SetModelServices, &
                specRoutine=SetModelServices, rc=rc)
        VERIFY_NUOPC_(rc)

        ! set driver verbosity
        print*, "ESM set attribute"
        call NUOPC_CompAttributeSet(driver, name="Verbosity", value="low", rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "ESM finish SetServices"
    end subroutine

!-----------------------------------------------------------------------------

    subroutine SetModelServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        ! local variables
        type(ESMF_Grid)               :: grid
        type(ESMF_Field)              :: field
        type(ESMF_Time)               :: startTime
        type(ESMF_Time)               :: stopTime
        type(ESMF_TimeInterval)       :: timeStep
        type(ESMF_Clock)              :: internalClock
        type(ESMF_GridComp)           :: child, test
        type(ESMF_CplComp)            :: connector

        rc = ESMF_SUCCESS

        print*, "ESM start SetModelServices"

        ! SetServices for ATM
        print*, "ESM add ATM comp"
        call NUOPC_DriverAddComp(driver, "ATM", atmSS, comp=child, rc=rc)
        VERIFY_NUOPC_(rc)
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="low", rc=rc)
        VERIFY_NUOPC_(rc)

        call NUOPC_DriverAddComp(driver, "test", NUOPC_CAP_SS, comp=child, rc=rc)
        VERIFY_NUOPC_(rc)

        ! SetServices for OCN
        print*, "ESM add OCN comp"
        call NUOPC_DriverAddComp(driver, "OCN", ocnSS, comp=child, rc=rc)
        VERIFY_NUOPC_(rc)
        call NUOPC_CompAttributeSet(child, name="Verbosity", value="low", rc=rc)
        VERIFY_NUOPC_(rc)

        ! Disabling the following macro, e.g. renaming to WITHCONNECTORS_disable,
        ! will result in a driver that does not call connectors between the model
        ! components. This mode can be used if all model components are driven 
        ! as independent models. However, even for independent models the
        ! connectors can be set here, but will turn into no-ops.
#define WITHCONNECTORS
#ifdef WITHCONNECTORS
        ! SetServices for atm2ocn
        print*, "ESM Connect ATM -> OCN"
        call NUOPC_DriverAddComp(driver, srcCompLabel="ATM", dstCompLabel="OCN", &
                compSetServicesRoutine=cplSS, comp=connector, rc=rc)
        VERIFY_NUOPC_(rc)
        call NUOPC_CompAttributeSet(connector, name="Verbosity", value="low", &
                rc=rc)
        VERIFY_NUOPC_(rc)

        ! SetServices for ocn2atm
        print*, "ESM Connect OCN -> ATM"
        call NUOPC_DriverAddComp(driver, srcCompLabel="OCN", dstCompLabel="ATM", &
                compSetServicesRoutine=cplSS, comp=connector, rc=rc)
        VERIFY_NUOPC_(rc)
        call NUOPC_CompAttributeSet(connector, name="Verbosity", value="low", &
                rc=rc)
        VERIFY_NUOPC_(rc)
#endif

        ! set the driver clock
        print*, "ESM Setup Time Interval"
        call ESMF_TimeIntervalSet(timeStep, m=15, rc=rc) ! 15 minute steps
        VERIFY_NUOPC_(rc)

        print*, "ESM Setup Start Time"
        call ESMF_TimeSet(startTime, yy=2010, mm=6, dd=1, h=0, m=0, &
                calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "ESM Setup Stop Time"
        call ESMF_TimeSet(stopTime, yy=2010, mm=6, dd=1, h=1, m=0, &
                calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "ESM Create clock"
        internalClock = ESMF_ClockCreate(name="Application Clock", &
                timeStep=timeStep, startTime=startTime, stopTime=stopTime, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "ESM Set Clock"
        call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "ESM finish SetModelServices"
    end subroutine
end module
