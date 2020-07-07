#include "NUOPC_ErrLog.h"

module aerosol_driver
    use ESMF
    use NUOPC

    use NUOPC_Driver, &
        driverSS => SetServices, &
        modelSS  => label_SetModelServices, &
        runSS    => label_SetRunSequence

    implicit none
    private

    public SetServices
contains
    subroutine SetServices(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        type(ESMF_Config) :: config

        rc = ESMF_SUCCESS

        ! Register NUOPC generic methods
        call NUOPC_CompDerive(driver, driverSS, rc=rc)
        VERIFY_ESMF_(rc)

        ! Attach specializing methods
        call NUOPC_CompSpecialize(driver, specLabel=modelSS, &
                specRoutine=SetModelServices, rc=rc)
        VERIFY_ESMF_(rc)

        call NUOPC_CompSpecialize(driver, specLabel=runSS, &
                specRoutine=SetRunSequence, rc=rc)
        VERIFY_ESMF_(rc)

        ! Set the NUOPC configuration file
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

        rc = ESMF_SUCCESS

    end subroutine SetModelServices

    subroutine SetRunSequence(driver, rc)
        type(ESMF_GridComp)  :: driver
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS

    end subroutine SetRunSequence
end module aerosol_driver