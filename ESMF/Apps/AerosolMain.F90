#include "NUOPC_ErrLog.h"

program AerosolMain
    use ESMF
    use NUOPC

    implicit none
    include "mpif.h"

    integer :: rc


    ! Initialize ESMF
    call ESMF_Initialize(LogKindFlag=ESMF_LOGKIND_MULTI, rc=rc)
!    VERIFY_ESMF_(rc)

    call ESMF_LogWrite("Aerosol Test App Starting", ESMF_LOGMSG_INFO, rc=rc)
!    VERIFY_ESMF_(rc)
end program AerosolMain