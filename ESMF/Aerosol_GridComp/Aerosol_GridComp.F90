#include "MAPL_Generic.h"

module Aerosol_GridComp_mod
    use ESMF
    use MAPL

    use GOCART2G_GridCompMod, only : GOCART2G_SetServices => SetServices

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
    character(*), parameter :: num_bands    = 'NUM_BANDS:'

contains
    subroutine SetServices(gc, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        integer,             intent(  out) :: rc

        type(Aerosol_GridCompState_Wrapper)  :: wrap
        type(Aerosol_GridCompState), pointer :: self

        integer :: gocart
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

        ! add children
        gocart = MAPL_AddChild(gc, name='GOCART2G', ss=GOCART2G_SetServices, __RC__)

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

        logical :: gridIsPresent

        __Iam__('Initialize')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        call MAPL_GetObjectFromgc(gc, MAPL, __RC__)

        call get_NumBands(gc, __RC__)

        call ESMF_UserCompGetInternalState(gc, internal_name, wrap, status)
        VERIFY_(status)
        self => wrap%ptr

        call ESMF_GridCompGet(gc, gridIsPresent=gridIsPresent, __RC__)

        if (.not.gridIsPresent) then
          call MAPL_GridCreate(gc, __RC__)
        end if

        call MAPL_GenericInitialize(gc, import, export, clock, __RC__)

        ! Temporary measure
        call ForceAllocation(import, __RC__)
        call ForceAllocation(export, __RC__)

        _RETURN(_SUCCESS)

    end subroutine Initialize

    subroutine Run(gc, import, export, clock, rc)
        type(ESMF_GridComp), intent(inout) :: gc
        type(ESMF_State),    intent(inout) :: import
        type(ESMF_State),    intent(inout) :: export
        type(ESMF_Clock),    intent(inout) :: clock
        integer, optional,   intent(  out) :: rc

        integer :: urc
        integer :: item, gcCount
        integer :: phase, runPhaseCount
        character(len=ESMF_MAXSTR) :: comp_name
        character(len=ESMF_MAXSTR) :: itemName

        type(ESMF_GridComp), pointer :: gcs(:)
        type(ESMF_State),    pointer :: gim(:)
        type(ESMF_State),    pointer :: gex(:)
        type(ESMF_State)             :: internal
        type(MAPL_MetaComp), pointer :: MAPL, maplObj

        ! declare pointers to Aerosol_StateSpecs table entries
#include "Aerosol_DeclarePointer___.h"

        __Iam__('Run')

        call ESMF_GridCompGet(gc, name=comp_name, __RC__)
        Iam = trim(comp_name) //'::'// Iam

        ! get info from MAPL
        call MAPL_GetObjectFromGc(gc, MAPL, __RC__)
        nullify(gcs, gim, gex)
        call MAPL_Get(mapl, gcs=gcs, gim=gim, gex=gex, internal_ESMF_State=internal, __RC__)

        ! Get/set field information from MAPL for fields in Aerosol_StateSpecs table entries
#include "Aerosol_GetPointer___.h"

        !!! Do run step

        ! Run children
        gcCount = 0
        if (associated(gcs)) gcCount = size(gcs)

        do item = 1, gcCount
          call ESMF_GridCompGet(gcs(item), name=itemName, __RC__)
          call MAPL_GetObjectFromGC(gcs(item), maplObj, __RC__)
          call MAPL_Get(maplObj, NumRunPhases=runPhaseCount, __RC__)
          do phase = 1, runPhaseCount
            call MAPL_TimerOn(mapl, trim(itemName))
            call ESMF_GridCompRun(gcs(item), &
                    importState = gim(item), &
                    exportState = gex(item), &
                          clock = clock,     &
                          phase = phase,     &
                         userRc = urc,       &
                                  __RC__)
            _ASSERT(urc==ESMF_SUCCESS,'Failed to run child component')
            call MAPL_TimerOff(mapl, trim(itemName))
          end do
        end do

        _RETURN(_SUCCESS)

    end subroutine Run

    subroutine ForceAllocation(state, rc)
        type(ESMF_State),  intent(inout) :: state
        integer, optional, intent(  out) :: rc

        integer                                 :: itemCount, i
        character(len=ESMF_MAXSTR), allocatable :: itemNameList(:)
        type(ESMF_StateItem_FLAG),  allocatable :: itemTypeList(:)

        type(ESMF_Field)            :: field
        type(ESMF_FieldStatus_Flag) :: fieldStatus


        __Iam__('ForceAllocation')

        call ESMF_StateGet(state, itemCount=itemCount, __RC__)

        if (itemCount > 0) then
           allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=status)
           VERIFY_(status)

           call ESMF_StateGet(state, itemNameList=itemNameList, &
                itemTypeList=itemTypeList, __RC__)

           do i=1, itemCount
              if (itemTypeList(i) == ESMF_STATEITEM_FIELD) then
                 call ESMF_StateGet(state, trim(itemNameList(i)), field, __RC__)
                 call ESMF_FieldGet(field, status=fieldStatus, __RC__)
                 if (fieldStatus == ESMF_FIELDSTATUS_GRIDSET) then
                    call MAPL_AllocateCoupling(field, __RC__)
                 end if
              end if
           end do

           deallocate(itemNameList, itemTypeList, stat=status)
           VERIFY_(status)
        end if

        _RETURN(_SUCCESS)
    end subroutine ForceAllocation

    subroutine get_NumBands(gc, rc)
       type(ESMF_GridComp), intent(inout) :: gc
       integer, optional,   intent(  out) :: rc

       type(MAPL_MetaComp), pointer :: MAPL
       type(ESMF_Config)            :: cfg
       type(ESMF_Grid)              :: grid
       logical                      :: is_present_in_grid, is_present_in_config
       integer                      :: n_bands
       integer                      :: status

       call MAPL_GetObjectFromGc(gc, MAPL, __RC__)
       call MAPL_Get(MAPL, cf=cfg, __RC__)

       call ESMF_GridCompGet(gc, grid=grid, __RC__)
       call ESMF_AttributeGet(grid, name=num_bands, isPresent=is_present_in_grid, __RC__)

       if (is_present_in_grid) then
          call ESMF_AttributeGet(grid, value=n_bands, name=num_bands, __RC__)
          call MAPL_ConfigSetAttribute(cfg, n_bands, label=num_bands, __RC__)

       else
          call ESMF_ConfigFindLabel(cfg, label=num_bands, isPresent=is_present_in_config, __RC__)
          if (is_present_in_config) then
             if (MAPL_Am_I_Root()) print*, "WARNING: falling back on MAPL NUM_BANDS"
          else
             _ASSERT(0==1, "NUM_BANDS must be an attribute of grid or defined in config file of Aerosol_GridComp")
          end if
       end if


       _RETURN(_SUCCESS)
    end subroutine get_NumBands

end module Aerosol_GridComp_mod
