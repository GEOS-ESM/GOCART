#include "MAPL_Generic.h"

module LinearFields_mod
    use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

    use ESMF
    use MAPL

    implicit none
    private

    public scale_field
    public shift_field

    type, abstract :: MAPL_Field
        private
        type(ESMF_Field) :: field
    contains
        procedure(i_scale_real32), deferred :: scale_real32
        procedure(i_scale_real64), deferred :: scale_real64
        generic                             :: scale => scale_real32, scale_real64

        procedure(i_shift_real32), deferred :: shift_real32
        procedure(i_shift_real64), deferred :: shift_real64
        generic                             :: shift => shift_real32, shift_real64
    end type MAPL_Field

    type, extends(MAPL_Field) :: MAPL_Real32_2DField
    contains
        procedure :: scale_real32 => scale_real32_2D_field_real32
        procedure :: scale_real64 => scale_real32_2D_field_real64

        procedure :: shift_real32 => shift_real32_2D_field_real32
        procedure :: shift_real64 => shift_real32_2D_field_real64
    end type MAPL_Real32_2DField

    type, extends(MAPL_Field) :: MAPL_Real32_3DField
    contains
        procedure :: scale_real32 => scale_real32_3D_field_real32
        procedure :: scale_real64 => scale_real32_3D_field_real64

        procedure :: shift_real32 => shift_real32_3D_field_real32
        procedure :: shift_real64 => shift_real32_3D_field_real64
    end type MAPL_Real32_3DField

    type, extends(MAPL_Field) :: MAPL_Real64_2DField
    contains
        procedure :: scale_real32 => scale_real64_2D_field_real32
        procedure :: scale_real64 => scale_real64_2D_field_real64

        procedure :: shift_real32 => shift_real64_2D_field_real32
        procedure :: shift_real64 => shift_real64_2D_field_real64
    end type MAPL_Real64_2DField

    type, extends(MAPL_Field) :: MAPL_Real64_3DField
    contains
        procedure :: scale_real32 => scale_real64_3D_field_real32
        procedure :: scale_real64 => scale_real64_3D_field_real64

        procedure :: shift_real32 => shift_real64_3D_field_real32
        procedure :: shift_real64 => shift_real64_3D_field_real64
    end type MAPL_Real64_3DField

    abstract interface
        subroutine i_scale_real32(this, scale_factor, rc)
            use, intrinsic :: iso_fortran_env, only: REAL32
            import MAPL_Field
            class(MAPL_Field), intent(inout) :: this
            real(kind=REAL32), intent(in   ) :: scale_factor
            integer, optional, intent(  out) :: rc
        end subroutine i_scale_real32

        subroutine i_scale_real64(this, scale_factor, rc)
            use, intrinsic :: iso_fortran_env, only: REAL64
            import MAPL_Field
            class(MAPL_Field), intent(inout) :: this
            real(kind=REAL64), intent(in   ) :: scale_factor
            integer, optional, intent(  out) :: rc
        end subroutine i_scale_real64

        subroutine i_shift_real32(this, shift_factor, rc)
            use, intrinsic :: iso_fortran_env, only: REAL32
            import MAPL_Field
            class(MAPL_Field), intent(inout) :: this
            real(kind=REAL32), intent(in   ) :: shift_factor
            integer, optional, intent(  out) :: rc
        end subroutine i_shift_real32

        subroutine i_shift_real64(this, shift_factor, rc)
            use, intrinsic :: iso_fortran_env, only: REAL64
            import MAPL_Field
            class(MAPL_Field), intent(inout) :: this
            real(kind=REAL64), intent(in   ) :: shift_factor
            integer, optional, intent(  out) :: rc
        end subroutine i_shift_real64
    end interface

    interface scale_field
        module procedure scale_field_real32
        module procedure scale_field_real64
    end interface scale_field

    interface shift_field
        module procedure shift_field_real32
        module procedure shift_field_real64
    end interface shift_field
contains
    subroutine wrap_field(field, wrapper, rc)
        type(ESMF_Field),               intent(in   ) :: field
        class(MAPL_Field), allocatable, intent(  out) :: wrapper
        integer, optional,              intent(  out) :: rc

        type(ESMF_TypeKind_Flag)       :: typekind
        integer                        :: status, field_rank

        call ESMF_FieldGet(field, typekind=typekind, rank=field_rank, __RC__)

        if (typekind == ESMF_TYPEKIND_R4) then
            select case(field_rank)
            case (2)
                wrapper = MAPL_Real32_2DField(field)
            case (3)
                wrapper = MAPL_Real32_3DField(field)
            case default
                _FAIL("Unsupported ESMF_Field Rank")
            end select
        elseif (typekind == ESMF_TYPEKIND_R8) then
            select case(field_rank)
            case (2)
                wrapper = MAPL_Real64_2DField(field)
            case (3)
                wrapper = MAPL_Real64_3DField(field)
            case default
                _FAIL("Unsupported ESMF_Field Rank")
            end select
        else
            _FAIL("Unsupported ESMF_TYPEKIND")
        end if

        _RETURN(_SUCCESS)
    end subroutine wrap_field

    subroutine scale_field_real32(field, scale_factor, rc)
        type(ESMF_Field),  intent(in   ) :: field
        real(kind=REAL32), intent(in   ) :: scale_factor
        integer, optional, intent(  out) :: rc

        class(MAPL_Field), allocatable :: wrapper
        integer                        :: status

        call wrap_field(field, wrapper, __RC__)

        call wrapper%scale(scale_factor, __RC__)

        _RETURN(_SUCCESS)
    end subroutine scale_field_real32

    subroutine scale_field_real64(field, scale_factor, rc)
        type(ESMF_Field),  intent(in   ) :: field
        real(kind=REAL64), intent(in   ) :: scale_factor
        integer, optional, intent(  out) :: rc

        class(MAPL_Field), allocatable :: wrapper
        integer                        :: status

        call wrap_field(field, wrapper, __RC__)

        call wrapper%scale(scale_factor, __RC__)

        _RETURN(_SUCCESS)
    end subroutine scale_field_real64

    subroutine scale_real32_2D_field_real32(this, scale_factor, rc)
        class(MAPL_Real32_2DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)

        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real32_2D_field_real32

    subroutine scale_real32_2D_field_real64(this, scale_factor, rc)
        class(MAPL_Real32_2DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real32_2D_field_real64

    subroutine scale_real32_3D_field_real32(this, scale_factor, rc)
        class(MAPL_Real32_3DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real32_3D_field_real32

    subroutine scale_real32_3D_field_real64(this, scale_factor, rc)
        class(MAPL_Real32_3DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real32_3D_field_real64

    subroutine scale_real64_2D_field_real32(this, scale_factor, rc)
        class(MAPL_Real64_2DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real64_2D_field_real32

    subroutine scale_real64_2D_field_real64(this, scale_factor, rc)
        class(MAPL_Real64_2DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real64_2D_field_real64

    subroutine scale_real64_3D_field_real32(this, scale_factor, rc)
        class(MAPL_Real64_3DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real64_3D_field_real32

    subroutine scale_real64_3D_field_real64(this, scale_factor, rc)
        class(MAPL_Real64_3DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: scale_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array*scale_factor

        _RETURN(_SUCCESS)
    end subroutine scale_real64_3D_field_real64

    subroutine shift_field_real32(field, shift_factor, rc)
        type(ESMF_Field),  intent(in   ) :: field
        real(kind=REAL32), intent(in   ) :: shift_factor
        integer, optional, intent(  out) :: rc

        class(MAPL_Field), allocatable :: wrapper
        integer                        :: status

        call wrap_field(field, wrapper, __RC__)

        call wrapper%shift(shift_factor, __RC__)

        _RETURN(_SUCCESS)
    end subroutine shift_field_real32

    subroutine shift_field_real64(field, shift_factor, rc)
        type(ESMF_Field),  intent(in   ) :: field
        real(kind=REAL64), intent(in   ) :: shift_factor
        integer, optional, intent(  out) :: rc

        class(MAPL_Field), allocatable :: wrapper
        integer                        :: status

        call wrap_field(field, wrapper, __RC__)

        call wrapper%shift(shift_factor, __RC__)

        _RETURN(_SUCCESS)
    end subroutine shift_field_real64

    subroutine shift_real32_2D_field_real32(this, shift_factor, rc)
        class(MAPL_Real32_2DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)

        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real32_2D_field_real32

    subroutine shift_real32_2D_field_real64(this, shift_factor, rc)
        class(MAPL_Real32_2DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real32_2D_field_real64

    subroutine shift_real32_3D_field_real32(this, shift_factor, rc)
        class(MAPL_Real32_3DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real32_3D_field_real32

    subroutine shift_real32_3D_field_real64(this, shift_factor, rc)
        class(MAPL_Real32_3DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL32), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real32_3D_field_real64

    subroutine shift_real64_2D_field_real32(this, shift_factor, rc)
        class(MAPL_Real64_2DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real64_2D_field_real32

    subroutine shift_real64_2D_field_real64(this, shift_factor, rc)
        class(MAPL_Real64_2DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real64_2D_field_real64

    subroutine shift_real64_3D_field_real32(this, shift_factor, rc)
        class(MAPL_Real64_3DField), intent(inout) :: this
        real(kind=REAL32),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real64_3D_field_real32

    subroutine shift_real64_3D_field_real64(this, shift_factor, rc)
        class(MAPL_Real64_3DField), intent(inout) :: this
        real(kind=REAL64),          intent(in   ) :: shift_factor
        integer, optional,          intent(  out) :: rc

        real(kind=REAL64), pointer :: array(:,:,:)
        integer                    :: status

        call ESMF_FieldGet(this%field, localDE=0, farrayPtr=array, __RC__)
        array = array+shift_factor

        _RETURN(_SUCCESS)
    end subroutine shift_real64_3D_field_real64
end module LinearFields_mod