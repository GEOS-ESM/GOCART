#include "NUOPC_ErrLog.h"
#include "MAPL_Generic.h"

module UFS_Testing_Cap
    use ESMF
    use NUOPC
    use NUOPC_Model, &
        ufsSS => SetServices, &
        ufsLA => label_Advance

    use MAPL

    implicit none
    private

    public SetServices

!    integer, parameter :: ImportFieldCount = 0
!    character(len=*), dimension(ImportFieldCount), parameter :: &
!            ImportFieldNames

    integer, parameter :: ExportFieldCount = 1
    character(len=*), dimension(ExportFieldCount), parameter :: &
            ExportFieldNames = [ &
                "inst_temp_levels" &
            ]
contains
    subroutine SetServices(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS

        ! Register NUOPC generic methods
        call NUOPC_CompDerive(model, ufsSS, rc=rc)
        VERIFY_ESMF_(rc)

        ! set entry point for methods that require specific implementation
        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE,&
                phaseLabelList=["IPDv05p1"], userRoutine=AdvertiseFields, rc=rc)
        VERIFY_ESMF_(rc)

        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE,&
                phaseLabelList=["IPDv05p4"], userRoutine=RealizeFields, rc=rc)
        VERIFY_ESMF_(rc)

        ! attach specializing method
        call NUOPC_CompSpecialize(model, specLabel=ufsLA, &
                specRoutine=ModelAdvance, rc=rc)
    end subroutine SetServices

    subroutine AdvertiseFields(model, import_state, export_state, clock, rc)
        type(ESMF_GridComp)  :: model
        type(ESMF_State)     :: import_state, export_state
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        _UNUSED_DUMMY(clock)

        rc = ESMF_SUCCESS

        if (ExportFieldCount > 0) then
            call NUOPC_Advertise(export_state, ExportFieldNames, &
                    TransferOfferGeomObject="cannot provide", &
                    SharePolicyField="share", rc=rc)
            VERIFY_ESMF_(rc)
        end if
    end subroutine AdvertiseFields

    subroutine RealizeFields(model, import_state, export_state, clock, rc)
        type(ESMF_GridComp)  :: model
        type(ESMF_State)     :: import_state, export_state
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        _UNUSED_DUMMY(clock)

        rc = ESMF_SUCCESS
    end subroutine RealizeFields

    subroutine ModelAdvance(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS
    end subroutine ModelAdvance
end module UFS_Testing_Cap