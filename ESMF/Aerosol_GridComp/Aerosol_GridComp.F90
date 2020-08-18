#include "MAPL_Generic.h"

module Aerosol_GridComp_mod
    use ESMF
    use MAPL

    implicit none
    private

    public SetServices

    type Aerosol_GridCompState
        integer :: nbins
    end type Aerosol_GridCompState

    type :: Aerosol_GridCompState_Wrapper
        type(Aerosol_GridCompState), pointer :: ptr => null()
    end type Aerosol_GridCompState_Wrapper

    character(*), parameter :: internal_name = "Aerosol_GridCompState"

contains
    subroutine SetServices(gc, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        integer,             intent(  out) :: rc

        type(Aerosol_GridCompState_Wrapper)  :: wrap
        type(Aerosol_GridCompState), pointer :: self

        character(len=ESMF_MAXSTR) :: comp_name

        __Iam__('SetServices')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        allocate(self, __STAT__)
        wrap%ptr => self

        call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_INITIALIZE, Initialize, __RC__)
        call MAPL_GridCompSetEntryPoint(gc, ESMF_METHOD_RUN, Run, __RC__)

        call ESMF_UserCompSetInternalState(gc, internal_name, wrap, status)
        VERIFY_(status)

        ! declare imports and exports based on StateSpecs.rc tables
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

        type(Aerosol_GridCompState_Wrapper)  :: wrap
        type(Aerosol_GridCompState), pointer :: self
        type(MAPL_MetaComp),         pointer :: MAPL

        __Iam__('Initialize')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GetObjectFromgc(gc, MAPL, __RC__)

        call ESMF_UserCompGetInternalState(gc, internal_name, wrap, status)
        VERIFY_(status)
        self => wrap%ptr

        call MAPL_GridCreate(gc, __RC__)

        call MAPL_GenericInitialize(gc, import, export, clock, __RC__)

        ! Temporary measure
        call ForceAllocation(export, __RC__)

        _RETURN(_SUCCESS)
    end subroutine Initialize

    subroutine Run(gc, import, export, clock, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        type(ESMF_State),    intent(inout) :: import
        type(ESMF_State),    intent(inout) :: export
        type(ESMF_Clock),    intent(inout) :: clock
        integer, optional,   intent(  out) :: rc

        character(len=ESMF_MAXSTR) :: comp_name

        type(ESMF_State)             :: internal
        type(MAPL_MetaComp), pointer :: MAPL

        ! declare pointers to Aerosol_StateSpecs table entries
#include "Aerosol_DeclarePointer___.h"

        __Iam__('Run')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GetObjectFromGc(gc, MAPL, __RC__)
        call MAPL_Get(mapl, internal_ESMF_State=internal, __RC__)

        ! Get/set field information from MAPL for fields in Aerosol_StateSpecs table entries
#include "Aerosol_GetPointer___.h"

        !!! Do run step

        _RETURN(_SUCCESS)
    end subroutine Run

    subroutine ForceAllocation(state, rc)
        type(ESMF_State),  intent(inout) :: state
        integer, optional, intent(  out) :: rc

        integer                                 :: itemCount, i, dims
        character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
        type(ESMF_StateItem_FLAG),  allocatable :: itemTypeList(:)

        type(ESMF_Field) :: field
        real, pointer    :: ptr2d(:,:), ptr3d(:,:,:)

        __Iam__('ForceAllocation')

        call ESMF_StateGet(state, itemCount=itemCount, __RC__)
        allocate(itemNameList(itemCount), stat=status)
        VERIFY_(status)
        allocate(itemTypeList(itemCount), stat=status)
        VERIFY_(status)

        call ESMF_StateGet(state, itemNameList=itemNameList, &
                itemTypeList=itemTypeList, __RC__)

        if (itemCount /= 0) then
            do i=1, itemCount
                if (itemTypeList(i) == ESMF_STATEITEM_FIELD) then
                    call ESMF_StateGet(state, trim(itemNameList(i)), field, __RC__)
                    call ESMF_AttributeGet(field, name='DIMS', value=dims, __RC__)
                    if (dims == MAPL_DimsHorzOnly) then
                        call MAPL_GetPointer(state, ptr2d, trim(itemNameList(i)), &
                                alloc=.true., __RC__)
                    else if (dims == MAPL_DimsHorzVert) then
                        call MAPL_GetPointer(state, ptr3d, trim(itemNameList(i)), &
                                alloc=.true., __RC__)
                    else
                        _ASSERT((1 == 2), "invalid field for test defined")
                    end if
                end if
            end do
        end if

        _RETURN(_SUCCESS)
    end subroutine ForceAllocation
end module Aerosol_GridComp_mod
