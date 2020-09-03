#include "MAPL_Generic.h"
#include "NUOPC_ErrLog.h"

module aerosol_driver
    use ESMF
    use NUOPC

    use NUOPC_Driver, &
            driverSS => SetServices, &
            modelSS  => label_SetModelServices, &
            runSS    => label_SetRunSequence

    use MAPL
    use MAPL_NUOPCWrapperMod, only: wrapperSS => SetServices, init_wrapper
    use NUOPC_Connector,      only: cplSS     => SetServices

    use Aerosol_GridComp_mod, only: aerosolSS => SetServices

    implicit none
    private

    public SetServices
contains
    subroutine SetServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        type(ESMF_Config) :: config

        ! Register NUOPC generic methods
        call NUOPC_CompDerive(driver, driverSS, rc=rc)
        VERIFY_NUOPC_(rc)

        ! Attach specializing methods
        call NUOPC_CompSpecialize(driver, specLabel=modelSS, &
                specRoutine=SetModelServices, rc=rc)
        VERIFY_NUOPC_(rc)

!        call NUOPC_CompSpecialize(driver, specLabel=runSS, &
!                specRoutine=SetRunSequence, rc=rc)
!        VERIFY_NUOPC_(rc)

        ! Set the NUOPC configuration file
        config = ESMF_ConfigCreate(rc=rc)
        VERIFY_NUOPC_(rc)

        call ESMF_ConfigLoadFile(config, "NUOPC_run_config.txt", rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_GridCompSet(driver, config=config, rc=rc)
        VERIFY_NUOPC_(rc)

        _RETURN(_SUCCESS)
    end subroutine SetServices

    subroutine SetModelServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        type(ESMF_VM)       :: vm
        type(ESMF_Config)   :: config
        type(ESMF_GridComp) :: comp

        integer              :: i, n_pes
        integer, allocatable :: petlist(:)

        call set_clock(driver)

        ! Read GridComp configuration
        call ESMF_GridCompGet(driver, vm=vm, config=config, rc=rc)
        VERIFY_NUOPC_(rc)
        call ESMF_VMGet(vm, petCount=n_pes, rc=rc)
        VERIFY_NUOPC_(rc)

        allocate(petlist(n_pes))
        petlist = [(i, i = 0, n_pes - 1)]

        ! Create the MAPL grid_comp
        call NUOPC_DriverAddComp(driver, "aerosol", wrapperSS, comp=comp, &
                petlist=petlist, rc=rc)
        VERIFY_NUOPC_(rc)
        call init_wrapper(wrapper_gc=comp, name="aerosol", &
                cap_rc_file="AEROSOL_CAP.rc", root_set_services=aerosolSS, rc=rc)
        VERIFY_NUOPC_(rc)

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

        ! Read the end time
        call ESMF_GridCompGet(driver, config = config, rc = rc)
        VERIFY_NUOPC_(rc)
        call ESMF_ConfigGetAttribute(config, end_date_and_time(1), label = "end_date:", rc = rc)
        VERIFY_NUOPC_(rc)
        call ESMF_ConfigGetAttribute(config, end_date_and_time(2), label = "end_time:", rc = rc)
        VERIFY_NUOPC_(rc)

        ! Set the end time
        call UnpackDateTime(end_date_and_time, yy, mm, dd, h, m, s)
        call ESMF_TimeSet(stopTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
                calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
        VERIFY_NUOPC_(rc)

        ! Read the interpolation time interval
        call ESMF_ConfigGetAttribute(config, dt, label = "interpolation_dt:", rc = rc)
        VERIFY_NUOPC_(rc)

        ! Create the driver clock
        call ESMF_TimeIntervalSet(timeStep, s=dt, rc=rc)
        VERIFY_NUOPC_(rc)

        internalClock = ESMF_ClockCreate(name="Driver Clock", timeStep = timeStep, &
                startTime=startTime, stopTime=stopTime, rc=rc)
        VERIFY_NUOPC_(rc)

        ! set the driver clock
        call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
        VERIFY_NUOPC_(rc)

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
!
!    subroutine SetRunSequence(driver, rc)
!        type(ESMF_GridComp)  :: driver
!        integer, intent(out) :: rc
!
!        rc = ESMF_SUCCESS
!
!    end subroutine SetRunSequence
end module aerosol_driver