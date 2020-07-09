#include "NUOPC_ErrLog.h"
#include "MAPL_Generic.h"

module MaplNuopc_driver
    use ESMF
    use NUOPC

    use NUOPC_Driver, &
            driverSS => SetServices, &
            modelSS  => label_SetModelServices, &
            runSS    => label_SetRunSequence

    use MAPL
    use MAPL_NUOPCWrapperMod, only: wrapperSS => SetServices, init_wrapper

    use MAPL_Provide_GridComp, only: provideSS => SetServices
    use MAPL_Recieve_GridComp, only: recieveSS => SetServices

    implicit none
    private

    public SetServices

contains
    subroutine SetServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        type(ESMF_Config) :: config

        rc = ESMF_SUCCESS

        ! NUOPC_Driver registers the generic methods
        call NUOPC_CompDerive(driver, driverSS, rc=rc)
        VERIFY_ESMF_(rc)

        ! attach specializing method(s)
        call NUOPC_CompSpecialize(driver, specLabel=modelSS, &
                specRoutine=SetModelServices, rc=rc)
        VERIFY_ESMF_(rc)

        call NUOPC_CompSpecialize(driver, specLabel=runSS, &
                specRoutine=SetRunSequence, rc=rc)
        VERIFY_ESMF_(rc)

        config = ESMF_ConfigCreate(rc=rc)
        VERIFY_ESMF_(rc)
        call ESMF_ConfigLoadFile(config, "NUOPC_run_config.txt", rc=rc)
        VERIFY_ESMF_(rc)
        call ESMF_GridCompSet(driver, config=config, rc=rc)
        VERIFY_ESMF_(rc)
    end subroutine SetServices

    subroutine SetModelServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        type(ESMF_VM)       :: vm
        type(ESMF_Config)   :: config
        type(ESMF_GridComp) :: provide, recieve

        logical              :: seq
        integer              :: i, n_pes, n_provide_pes, n_recieve_pes
        integer, allocatable :: provide_petlist(:), recieve_petlist(:)

        rc = ESMF_SUCCESS

        call set_clock(driver)

        call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
        VERIFY_ESMF_(rc)
        call ESMF_VMGet(vm, petCount=n_pes, rc=rc)
        VERIFY_ESMF_(rc)

        call ESMF_ConfigGetAttribute(config, seq, label="sequential:", rc=rc)
        VERIFY_ESMF_(rc)
        call ESMF_ConfigGetAttribute(config, n_provide_pes, label="provide_pets:", rc=rc)
        VERIFY_ESMF_(rc)
        call ESMF_ConfigGetAttribute(config, n_recieve_pes, label="recieve_pets:", rc=rc)
        VERIFY_ESMF_(rc)

        allocate(provide_petlist(n_provide_pes), recieve_petlist(n_recieve_pes))

        if (seq) then
            _ASSERT((n_provide_pes == n_recieve_pes), "provide_pets must be equal to recieve_pets in sequential")
            _ASSERT((n_provide_pes == n_pes), "provide_pets must be equal to number of pets in sequential")

            provide_petlist = [(i, i = 0, n_provide_pes - 1)]
            recieve_petlist = [(i, i = 0, n_recieve_pes - 1)]
        else
            _ASSERT(((n_provide_pes + n_recieve_pes) == n_pes), "provide_pets + recieve_pets must be equal to number of pets")

            provide_petlist = [(i, i = 0, n_provide_pes - 1)]
            recieve_petlist = [(i, i = n_provide_pes, n_pes - 1)]
        end if


        call NUOPC_DriverAddComp(driver, "provide", wrapperSS, comp=provide, &
                petlist=provide_petlist, rc=rc)
        VERIFY_ESMF_(rc)
        call init_wrapper(wrapper_gc=provide, name="provide", &
                cap_rc_file="PROVIDE_CAP.rc", root_set_services=provideSS, rc=rc)
        VERIFY_ESMF_(rc)

        call NUOPC_DriverAddComp(driver, "recieve", wrapperSS, comp=recieve, &
                petlist=recieve_petlist, rc=rc)
        VERIFY_ESMF_(rc)
        call init_wrapper(wrapper_gc=recieve, name="recieve", &
                cap_rc_file="RECIEVE_CAP.rc", root_set_services=recieveSS, rc=rc)
        VERIFY_ESMF_(rc)
    end subroutine SetModelServices

    subroutine set_clock(driver)
        type(ESMF_GridComp), intent(inout) :: driver

        type(ESMF_Time)         :: startTime
        type(ESMF_Time)         :: stopTime
        type(ESMF_TimeInterval) :: timeStep
        type(ESMF_Clock)        :: internalClock
        type(ESMF_Config)       :: config

        integer :: start_date_and_time(2), end_date_and_time(2), dt, file_unit, yy, mm, dd, h, m, s, rc

        ! Read the start time
        open(newunit = file_unit, file = "cap_restart", form = 'formatted', &
                status = 'old', action = 'read')
        read(file_unit, *) start_date_and_time
        close(file_unit)

        ! Set the start time
        call UnpackDateTime(start_date_and_time, yy, mm, dd, h, m, s)
        call ESMF_TimeSet(startTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
                calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
        VERIFY_ESMF_(rc)

        ! Read the end time
        call ESMF_GridCompGet(driver, config = config, rc = rc)
        VERIFY_ESMF_(rc)
        call ESMF_ConfigGetAttribute(config, end_date_and_time(1), label = "end_date:", rc = rc)
        VERIFY_ESMF_(rc)
        call ESMF_ConfigGetAttribute(config, end_date_and_time(2), label = "end_time:", rc = rc)
        VERIFY_ESMF_(rc)

        ! Set the end time
        call UnpackDateTime(end_date_and_time, yy, mm, dd, h, m, s)
        call ESMF_TimeSet(stopTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
                calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
        VERIFY_ESMF_(rc)

        ! Read the interpolation time interval
        call ESMF_ConfigGetAttribute(config, dt, label = "interpolation_dt:", rc = rc)
        VERIFY_ESMF_(rc)

        ! Create the driver clock
        call ESMF_TimeIntervalSet(timeStep, s=dt, rc=rc)
        VERIFY_ESMF_(rc)

        internalClock = ESMF_ClockCreate(name="Driver Clock", timeStep = timeStep, &
                startTime=startTime, stopTime=stopTime, rc=rc)
        VERIFY_ESMF_(rc)

        ! set the driver clock
        call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
        VERIFY_ESMF_(rc)
    contains
        subroutine UnpackDateTime(DATETIME, YY, MM, DD, H, M, S)
            integer, intent(IN)  :: DATETIME(:)
            integer, intent(OUT) :: YY, MM, DD, H, M, S

            YY =     datetime(1)/10000
            MM = mod(datetime(1),10000)/100
            DD = mod(datetime(1),100)
            H  =     datetime(2)/10000
            M  = mod(datetime(2),10000)/100
            S  = mod(datetime(2),100)
            return
        end subroutine UnpackDateTime
    end subroutine set_clock

    subroutine SetRunSequence(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS
    end subroutine SetRunSequence
end module MaplNuopc_driver