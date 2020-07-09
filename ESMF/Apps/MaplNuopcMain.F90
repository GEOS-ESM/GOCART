#include "NUOPC_ErrLog.h"

program MaplNuopcMain
    use ESMF
    use NUOPC
    use, intrinsic :: iso_fortran_env, only: INT64

    use MaplNuopc_driver, only: driverSS => SetServices
    use MAPL_Profiler,  only: BaseProfiler, TimeProfiler, get_global_time_profiler

    implicit none
    include "mpif.h"

    integer                      :: rank, file_unit, rc, urc
    integer(kind=INT64)          :: t0, t1, count_rate
    real                         :: elapsed_time
    class(BaseProfiler), pointer :: t_p

    type(ESMF_GridComp) :: prototype_gridcomp

    call mpi_init(rc)
    VERIFY_ESMF_(rc)

    call mpi_comm_rank(MPI_COMM_WORLD, rank, rc)
    VERIFY_ESMF_(rc)

    if (rank == 0) then
        call system_clock(t0)
    end if

    t_p => get_global_time_profiler()
    t_p = TimeProfiler('All', comm_world=MPI_COMM_WORLD)
    call t_p%start()

    ! Initialize ESMF
    call ESMF_Initialize(LogKindFlag=ESMF_LOGKIND_MULTI, rc=rc)
    VERIFY_ESMF_(rc)

    call ESMF_LogWrite("MAPL-NUOPC Test STARTING", ESMF_LOGMSG_INFO, rc=rc)
    VERIFY_ESMF_(rc)

    ! Create the prototype
    prototype_gridcomp = ESMF_GridCompCreate(name="prototype", rc=rc)
    VERIFY_ESMF_(rc)

    ! Run set services for prototype
    call ESMF_GridCompSetServices(prototype_gridcomp, driverSS, userRc=urc, rc=rc)
    VERIFY_NUOPC_(rc, urc)

    ! Initialize the prototype
    call ESMF_GridCompInitialize(prototype_gridcomp, userRc=urc, rc=rc)
    VERIFY_NUOPC_(rc, urc)

    ! Run the the prototype
    call ESMF_GridCompRun(prototype_gridcomp, userRc=urc, rc=rc)
    VERIFY_NUOPC_(rc, urc)

    ! Get the run time
    if (rank == 0) then
        call system_clock(t1, count_rate)
        open(newunit=file_unit, file="elapsed_time.txt", &
                status="replace", action="write")
        elapsed_time = (t1 - t0) / real(count_rate)
        write(file_unit, '(f0.0)') elapsed_time
        close(file_unit)
    end if

    ! Finalize the prototype
    call ESMF_GridCompFinalize(prototype_gridcomp, userRc=urc, rc=rc)
    VERIFY_NUOPC_(rc, urc)

    ! Destroy the prototype
    call ESMF_GridCompDestroy(prototype_gridcomp, rc=rc)
    VERIFY_ESMF_(rc)

    ! Write MPI errors
    call MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_RETURN, rc)
    VERIFY_ESMF_(rc)

    call ESMF_LogWrite("MAPL-NUOPC Test FINISHED", ESMF_LOGMSG_INFO, rc=rc)
    VERIFY_ESMF_(rc)

    ! Finalize ESMF
    call ESMF_Finalize()

    call t_p%stop()
end program MaplNuopcMain