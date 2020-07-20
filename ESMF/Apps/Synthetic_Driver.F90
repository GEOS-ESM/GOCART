#include "MAPL_Generic.h"
#include "NUOPC_ErrLog.h"

module synthetic_driver
    use ESMF
    use NUOPC
    use NUOPC_Driver, &
            driverSS => SetServices, &
            modelSS  => label_SetModelServices, &
            runSS    => label_SetRunSequence

    use MAPL
    use MAPL_NUOPCWrapperMod, only: wrapperSS => SetServices, init_wrapper
    use NUOPC_Connector,      only: cplSS     => SetServices

    use Provider_GridCompMod, only: providerSS => SetServices
    use Reciever_GridCompMod, only: recieverSS => SetServices

    use UFS_Testing_Cap, only: ufsSS => SetServices

    implicit none
    private

    public SetServices
contains
    subroutine SetServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc
        type(ESMF_Config)    :: config

        print*, "Driver start Set Services"

        ! NUOPC_Driver registers the generic methods
        print*,"Driver add Generic Set Services"
        call NUOPC_CompDerive(driver, driverSS, rc=rc)
        VERIFY_NUOPC_(rc)

        ! attach specializing method(s)
        print*,"Driver add Set Model Services"
        call NUOPC_CompSpecialize(driver, specLabel=modelSS, &
                specRoutine=SetModelServices, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver add Set Run Sequence"
        call NUOPC_CompSpecialize(driver, specLabel=runSS, &
                specRoutine=SetRunSequence, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver create config"
        config = ESMF_ConfigCreate(rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"Driver read config"
        call ESMF_ConfigLoadFile(config, "NUOPC_run_config.txt", rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"Driver add config"
        call ESMF_GridCompSet(driver, config=config, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "Driver finish Set Services"

        _RETURN(_SUCCESS)
    end subroutine SetServices

    subroutine SetModelServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        type(ESMF_GridComp) :: provider, reciever, mediator
        type(ESMF_CplComp)  :: connector
        type(ESMF_VM)       :: vm
        type(ESMF_Config)   :: config

        logical              :: seq
        integer              :: i, npes, n_provider_pes, n_reciever_pes
        integer, allocatable :: provider_petlist(:), reciever_petlist(:)

        print*, "Driver start Set Model Services"

        print*,"Driver set clock"
        call set_clock(driver)

        print*,"Driver read config"
        call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"Driver read npes"
        call ESMF_VMGet(vm, petCount=npes, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver read from config"
        call ESMF_ConfigGetAttribute(config, seq, label="sequential:", rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_ConfigGetAttribute(config, n_provider_pes, &
                label="provider_pets:", rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_ConfigGetAttribute(config, n_reciever_pes, &
                label="reciever_pets:", rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver create pet lists"
        allocate(provider_petlist(n_provider_pes), reciever_petlist(n_reciever_pes))

        if (seq) then
            _ASSERT((n_provider_pes == n_reciever_pes), "provider_pets must be equal to reciever_pets in sequential")
            _ASSERT((n_provider_pes == npes), "provider_pets must be equal to number of pets in sequential")

            provider_petlist = [(i, i = 0, n_provider_pes - 1)]
            reciever_petlist = [(i, i = 0, n_reciever_pes - 1)]
        else
            _ASSERT(((n_provider_pes + n_reciever_pes) == npes), "provider_pets + reciever_pets must be equal to number of pets")

            provider_petlist = [(i, i = 0, n_provider_pes - 1)]
            reciever_petlist = [(i, i = n_provider_pes, npes - 1)]
        end if

        print*,"Driver add provider"
        call NUOPC_DriverAddComp(driver, "provider", wrapperSS, comp=provider, &
                petlist=provider_petlist, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"Driver wrap provider MAPL"
        call init_wrapper(wrapper_gc=provider, name="provider", &
                cap_rc_file="AGCM_CAP.rc", root_set_services=providerSS, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver add reciever"
        call NUOPC_DriverAddComp(driver, "reciever", wrapperSS, comp=reciever, &
                petlist=reciever_petlist, rc=rc)
        VERIFY_NUOPC_(rc)
        print*,"Driver wrap reciever MAPL"
        call init_wrapper(wrapper_gc=reciever, name="reciever", &
                cap_rc_file="CTM_CAP.rc", root_set_services=recieverSS, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver connect compoinents"
        call NUOPC_DriverAddComp(driver, srcCompLabel="provider", dstCompLabel="reciever", &
                compSetServicesRoutine=cplSS, comp=connector, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "Driver finish Set Model Services"

        _RETURN(_SUCCESS)
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
        VERIFY_NUOPC_(rc)

        call ESMF_GridCompGet(driver, config = config, rc = rc)
        VERIFY_NUOPC_(rc)

        ! Read the end time
        call ESMF_ConfigGetAttribute(config, end_date_and_time(1), &
                label="end_date:", rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_ConfigGetAttribute(config, end_date_and_time(2), &
                label="end_time:", rc=rc)
        VERIFY_NUOPC_(rc)

        ! Set the end time
        call UnpackDateTime(end_date_and_time, yy, mm, dd, h, m, s)
        call ESMF_TimeSet(stopTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
                calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
        VERIFY_NUOPC_(rc)

        ! Read the interpolation time interval
        call ESMF_ConfigGetAttribute(config, dt, label="interpolation_dt:", rc = rc)
        VERIFY_NUOPC_(rc)

        ! Set the time interval
        call ESMF_TimeIntervalSet(timeStep, s=dt, rc=rc)
        VERIFY_NUOPC_(rc)
        internalClock = ESMF_ClockCreate(name="Driver Clock", timeStep=timeStep, &
                startTime=startTime, stopTime=stopTime, rc=rc)
        VERIFY_NUOPC_(rc)

        ! Set the driver clock
        call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
        VERIFY_NUOPC_(rc)
    contains
        subroutine UnpackDateTime(DATETIME, YY, MM, DD, H, M, S)
            integer, intent(IN   ) :: DATETIME(:)
            integer, intent(  OUT) :: YY, MM, DD, H, M, S

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

        ! local variables
        type(ESMF_Time)         :: startTime, stopTime
        type(ESMF_TimeInterval) :: timeStep
        type(ESMF_Config)       :: config
        type(NUOPC_FreeFormat)  :: run_sequence_ff

        print*, "Driver start Set Run Sequence"

        print*,"Driver read config"
        call ESMF_GridCompGet(driver, config=config, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver read run sequence"
        run_sequence_ff = NUOPC_FreeFormatCreate(config, label="run_sequence::", rc=rc)
        VERIFY_NUOPC_(rc)

        ! ingest FreeFormat run sequence
        print*,"Driver ingest run sequence"
        call NUOPC_DriverIngestRunSequence(driver, run_sequence_ff, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"Driver destroy run sequence"
        call NUOPC_FreeFormatDestroy(run_sequence_ff, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "Driver finish Set Run Sequence"

        _RETURN(_SUCCESS)
    end subroutine SetRunSequence
end module synthetic_driver
