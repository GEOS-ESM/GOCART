#include "NUOPC_ErrLog.h"

program AerosolMain
    use ESMF
    use NUOPC
    use, intrinsic :: iso_fortran_env, only: INT64

    use aerosol_driver, only: driverSS => SetServices
    use MAPL_Profiler,  only: BaseProfiler, TimeProfiler, get_global_time_profiler

    implicit none
    include "mpif.h"

    integer             :: rank, file_unit, rc, urc
    integer(kind=INT64) :: t0, t1, count_rate
    real                :: elapsed_time

    type(ESMF_GridComp) :: driver_gridcomp

    call mpi_init(rc)
    VERIFY_ESMF_(rc)

    call mpi_comm_rank(MPI_COMM_WORLD, rank, rc)
    VERIFY_ESMF_(rc)

    if (rank == 0) then
        call system_clock(t0)
    end if

    ! Initialize ESMF
    call ESMF_Initialize(LogKindFlag=ESMF_LOGKIND_MULTI, rc=rc)
    VERIFY_ESMF_(rc)

    call ESMF_LogWrite("Aerosol Test App STARTING", ESMF_LOGMSG_INFO, rc=rc)
    VERIFY_ESMF_(rc)

    ! Create the model
    driver_gridcomp = ESMF_GridCompCreate(name="driver", rc=rc)
    VERIFY_ESMF_(rc)

    call ESMF_GridCompSetServices(driver_gridcomp, driverSS, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    ! Initialize the model
    call ESMF_GridCompInitialize(driver_gridcomp, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    ! Run the Model
    call ESMF_GridCompRun(driver_gridcomp, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    ! Get the run time
    if (rank == 0) then
        call system_clock(t1, count_rate)
        open(newunit=file_unit, file="elapsed_time.txt", &
                status="replace", action="write")
        elapsed_time = (t1 - t0) / real(count_rate)
        write(file_unit, '(f0.0)') elapsed_time
        close(file_unit)
    end if

    ! Finalize the model
    call ESMF_GridCompFinalize(driver_gridcomp, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    ! Destroy the model
    call ESMF_GridCompDestroy(driver_gridcomp, rc=rc)
    VERIFY_ESMF_(rc)

    ! Write MPI errors
    call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, rc)
    VERIFY_ESMF_(rc)

    call ESMF_LogWrite("Aerosol Test App FINISHED", ESMF_LOGMSG_INFO, rc=rc)
    VERIFY_ESMF_(rc)

    ! Finalize ESMF
    call ESMF_Finalize()
end program AerosolMain
