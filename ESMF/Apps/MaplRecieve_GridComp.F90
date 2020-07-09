#include "MAPL_Generic.h"

module MAPL_Recieve_GridComp
    use ESMF
    use MAPL

    implicit none
    private

    public SetServices

    include "mpif.h"
contains
    subroutine SetServices(gc, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        integer,             intent(  out) :: rc

        character(len=ESMF_MAXSTR) :: comp_name

        __Iam__('SetServices')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, __RC__)
        call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, __RC__)

        call MAPL_AddImportSpec(gc, &
                short_name='var1', &
                long_name='na', &
                units='na', &
                dims=MAPL_DimsHorzOnly, &
                vlocation=MAPL_VLocationNone, __RC__)

        call MAPL_GenericSetServices(gc, __RC__)

        _RETURN(_SUCCESS)
    end subroutine SetServices

    subroutine Initialize(gc, import, export, clock, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        type(ESMF_State),    intent(inout) :: import
        type(ESMF_State),    intent(inout) :: export
        type(ESMF_Clock),    intent(inout) :: clock
        integer, optional,   intent(  out) :: rc

        character(len=ESMF_MAXSTR) :: comp_name

        __Iam__('Initialize')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GridCreate(gc, __RC__)

        call MAPL_GenericInitialize(gc, import, export, clock, __RC__)

        _RETURN(_SUCCESS)
    end subroutine Initialize

    subroutine Run(gc, import, export, clock, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        type(ESMF_State),    intent(inout) :: import
        type(ESMF_State),    intent(inout) :: export
        type(ESMF_Clock),    intent(inout) :: clock
        integer, optional,   intent(  out) :: rc

        character(len=ESMF_MAXSTR) :: comp_name

        real, pointer :: ptr2d(:,:)

        __Iam__('Run')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GetPointer(import, ptr2d, 'var1', __RC__)

        print*, 'The value of var1 is:', minval(ptr2d), maxval(ptr2d)

        _RETURN(_SUCCESS)
    end subroutine Run
end module MAPL_Recieve_GridComp