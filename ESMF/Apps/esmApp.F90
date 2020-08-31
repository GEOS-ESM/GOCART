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

program esmApp

    !-----------------------------------------------------------------------------
    ! Generic ESM application driver
    !-----------------------------------------------------------------------------

    use ESMF
    use NUOPC
    use ESM, only: esmSS => SetServices

    implicit none

    integer                 :: rc, urc
    type(ESMF_GridComp)     :: esmComp

    print*,"START Program"

    print*,"Initialize ESMF"
    ! Initialize ESMF
    call ESMF_Initialize(logkindflag=ESMF_LOGKIND_MULTI, rc=rc)
    VERIFY_ESMF_(rc)
    call ESMF_LogWrite("esmApp STARTING", ESMF_LOGMSG_INFO, rc=rc)
    VERIFY_ESMF_(rc)

    print*,"Create ESM"
    ! Create the earth system Component
    esmComp = ESMF_GridCompCreate(name="esm", rc=rc)
    VERIFY_ESMF_(rc)

    print*,"ESM call SetServices"
    ! SetServices for the earth system Component
    call ESMF_GridCompSetServices(esmComp, esmSS, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    print*,"ESM set profiling"
    ! Set Profiling Attribute
    call NUOPC_CompAttributeSet(esmComp, name="Profiling", value="0", rc=rc)
    VERIFY_ESMF_(rc)

    print*,"ESM call Initialize"
    ! Call Initialize for the earth system Component
    call ESMF_GridCompInitialize(esmComp, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    print*,"ESM call Run"
    ! Call Run  for earth the system Component
    call ESMF_GridCompRun(esmComp, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    print*,"ESM call Finalize"
    ! Call Finalize for the earth system Component
    call ESMF_GridCompFinalize(esmComp, userRc=urc, rc=rc)
    VERIFY_ALL_ESMF_(rc, urc)

    print*,"Destroy ESM"
    ! Destroy the earth system Component
    call ESMF_GridCompDestroy(esmComp, rc=rc)
    VERIFY_ESMF_(rc)

    call ESMF_LogWrite("esmApp FINISHED", ESMF_LOGMSG_INFO, rc=rc)
    VERIFY_ESMF_(rc)

    print*,"Finalize ESMF"
    ! Finalize ESMF
    call ESMF_Finalize()

    print*,"FINISH Program"
end program
