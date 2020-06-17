#include "MAPL_Generic.h"

module Aerosols_GridComp_mod
    use ESMF
    use MAPL

    implicit none
    private

    public SetServices

    type Aerosol_GridComp
        integer :: nbins
    end type Aerosol_GridComp

contains
    subroutine SetServices(gc, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        integer, optional,   intent(  out) :: rc

        type(Aerosol_GridComp), pointer :: self

        character(len=ESMF_MAXSTR) :: comp_name

        __Iam__('SetServices')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, __RC__)
        call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, __RC__)

#include "Aerosol_Internal___.h"
#include "Aerosol_Export___.h"
#include "Aerosol_Import___.h"

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

        type(MAPL_MetaComp), pointer :: MAPL

        __Iam__('Initialize')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GetObjectFromgc(gc, MAPL, __RC__)

        _RETURN(_SUCCESS)
    end subroutine Initialize

    subroutine Run(gc, import, export, clock, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        type(ESMF_State),    intent(inout) :: import
        type(ESMF_State),    intent(inout) :: export
        type(ESMF_Clock),    intent(inout) :: clock
        integer, optional,   intent(  out) :: rc

        character(len=ESMF_MAXSTR) :: comp_name

        type(MAPL_MetaComp), pointer :: MAPL

#include "Aerosol_DeclarePointer___.h"

        __Iam__('Run')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GetObjectFromgc(gc, MAPL, __RC__)

        _RETURN(_SUCCESS)
    end subroutine Run
end module Aerosols_GridComp_mod