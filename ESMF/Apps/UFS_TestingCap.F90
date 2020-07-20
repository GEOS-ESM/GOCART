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

    integer, parameter :: ImportFieldCount = 1
    character(len=*), dimension(ImportFieldCount), parameter :: &
            ImportFieldNames = [ &
                    "var1" &
                    ]

    integer, parameter :: ExportFieldCount = 1
    character(len=*), dimension(ExportFieldCount), parameter :: &
            ExportFieldNames = [ &
                "var2" &
            ]
contains
    subroutine SetServices(model, rc)
        type(ESMF_GridComp)  :: model
        integer, intent(out) :: rc

        rc = ESMF_SUCCESS

        print*,"UFS start SetServices"


        ! Register NUOPC generic methods
        call NUOPC_CompDerive(model, ufsSS, rc=rc)
        VERIFY_NUOPC_(rc)

        ! set entry point for methods that require specific implementation
        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE,&
                phaseLabelList=["IPDv05p1"], userRoutine=AdvertiseFields, rc=rc)
        VERIFY_NUOPC_(rc)

        call NUOPC_CompSetEntryPoint(model, ESMF_METHOD_INITIALIZE,&
                phaseLabelList=["IPDv05p4"], userRoutine=RealizeFields, rc=rc)
        VERIFY_NUOPC_(rc)

        ! attach specializing method
        call NUOPC_CompSpecialize(model, specLabel=ufsLA, &
                specRoutine=ModelAdvance, rc=rc)
        VERIFY_NUOPC_(rc)

        print*,"UFS finish SetServices"
    end subroutine SetServices

    subroutine AdvertiseFields(model, import_state, export_state, clock, rc)
        type(ESMF_GridComp)  :: model
        type(ESMF_State)     :: import_state, export_state
        type(ESMF_Clock)     :: clock
        integer, intent(out) :: rc

        _UNUSED_DUMMY(clock)

        rc = ESMF_SUCCESS
        print*,"UFS start advertise"

        if (ImportFieldCount > 0) then
            call NUOPC_Advertise(import_state, ImportFieldNames, &
                    TransferOfferGeomObject="cannot provide", &
                    SharePolicyField="share", rc=rc)
            VERIFY_NUOPC_(rc)
        end if

        if (ExportFieldCount > 0) then
            call NUOPC_Advertise(export_state, ExportFieldNames, &
                    TransferOfferGeomObject="cannot provide", &
                    SharePolicyField="share", rc=rc)
            VERIFY_NUOPC_(rc)
        end if

        print*,"UFS finish advertise"
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

        integer          :: i, num_items_import, num_items_export
        type(ESMF_State) :: import_state, export_state
        type(ESMF_Field) :: field
        real, pointer    :: ptr2d_in(:,:), ptr2d_out(:,:)

        character(len=ESMF_MAXSTR), allocatable :: item_names_import(:), item_names_export(:)

        rc = ESMF_SUCCESS

        print*,"UFS start advance"

        print*, "UFS get import and export states"
        call NUOPC_ModelGet(model, importState=import_state, &
                exportState=export_state, rc=rc)
        VERIFY_NUOPC_(rc)

        print*, "UFS get number of import items"
        call ESMF_StateGet(import_state, itemcount=num_items_import, rc=rc)
        VERIFY_NUOPC_(rc)
        print*, "UFS num import:", num_items_import

        allocate(item_names_import(num_items_import))
        print*, "UFS get import names"
        call ESMF_StateGet(import_state, itemnamelist=item_names_import, rc=rc)
        VERIFY_NUOPC_(rc)
        print*, "UFS import names are:", item_names_import

        print*, "UFS get number of export items"
        call ESMF_StateGet(export_state, itemcount=num_items_export, rc=rc)
        VERIFY_NUOPC_(rc)
        print*, "UFS num export:", num_items_export

        allocate(item_names_export(num_items_export))
        print*, "UFS get export names"
        call ESMF_StateGet(export_state, itemnamelist=item_names_export, rc=rc)
        VERIFY_NUOPC_(rc)
        print*, "UFS export names are:", item_names_export

!        print*, "UFS get imports"
!        do i=1, ImportFieldCount
!            print*, "UFS get field"
!            call ESMF_StateGet(import_state, field=field, &
!                    itemName=trim(ImportFieldNames(i)), rc=rc)
!            VERIFY_NUOPC_(rc)
!            print*, "UFS get ptr"
!            call ESMF_FieldGet(field, localDE=0, farrayPtr=ptr2d_in, rc=rc)
!            VERIFY_NUOPC_(rc)
!
!            print*, 'The value of ', trim(ImportFieldNames(i)), ' is:', minval(ptr2d_in), maxval(ptr2d_in)
!        end do
        print*,"UFS finish advance"
    end subroutine ModelAdvance
end module UFS_Testing_Cap
