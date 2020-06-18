#include "MAPL_Generic.h"

module ConvertUnits_private_mod
    use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

    use ESMF
    use MAPL
    use gFTL_StringReal32Map
    use gFTL_StringReal64Map

    use LinearFields_mod

    implicit none
    private

    public ScaleMapReal32
    public ScaleMapReal64

    public read_scale_config
    public read_scale_config_real32
    public read_scale_config_real64

    public scale_field_in_real32
    public scale_field_in_real64
    public scale_field_out_real32
    public scale_field_out_real64

    interface read_scale_config
        module procedure read_scale_config_real32
        module procedure read_scale_config_real64
    end interface read_scale_config

    type, extends(StringReal32Map) :: ScaleMapReal32
    contains
        procedure :: scale_in  => scale_in_real32
        procedure :: scale_out => scale_out_real32
    end type ScaleMapReal32

    type, extends(StringReal64Map) :: ScaleMapReal64
    contains
        procedure :: scale_in  => scale_in_real64
        procedure :: scale_out => scale_out_real64
    end type ScaleMapReal64
contains
    subroutine read_scale_config_real32(config, scale_map, rc)
        type(ESMF_Config),    intent(inout) :: config
        type(ScaleMapReal32), intent(  out) :: scale_map
        integer, optional,    intent(  out) :: rc

        integer                    :: i, num_fields, columns, status
        real(kind=REAL32)          :: scale_factor
        character(len=ESMF_MaxStr) :: name

        call ESMF_ConfigGetDim(config, label="fields_to_scale:", lineCount=num_fields, &
                columnCount=columns, __RC__)
        call ESMF_ConfigFindLabel(config, label="fields_to_scale:", __RC__)
        call ESMF_ConfigNextLine(config, __RC__)

        do i=1, num_fields
            call ESMF_ConfigGetAttribute(config, value=name, __RC__)
            call ESMF_ConfigGetAttribute(config, value=scale_factor, __RC__)

            call scale_map%insert(trim(name), scale_factor)

            call ESMF_ConfigNextLine(config, __RC__)
        end do

        _RETURN(_SUCCESS)
    end subroutine read_scale_config_real32

    subroutine read_scale_config_real64(config, scale_map, rc)
        type(ESMF_Config),    intent(inout) :: config
        type(ScaleMapReal64), intent(  out) :: scale_map
        integer, optional,    intent(  out) :: rc

        integer                    :: i, num_fields, columns, status
        real(kind=REAL64)          :: scale_factor
        character(len=ESMF_MaxStr) :: name

        call ESMF_ConfigGetDim(config, label="fields_to_scale:", lineCount=num_fields, &
                columnCount=columns, __RC__)
        call ESMF_ConfigFindLabel(config, label="fields_to_scale:", __RC__)
        call ESMF_ConfigNextLine(config, __RC__)

        do i=1, num_fields
            call ESMF_ConfigGetAttribute(config, value=name, __RC__)
            call ESMF_ConfigGetAttribute(config, value=scale_factor, __RC__)

            call scale_map%insert(trim(name), scale_factor)

            call ESMF_ConfigNextLine(config, __RC__)
        end do

        _RETURN(_SUCCESS)
    end subroutine read_scale_config_real64

    subroutine scale_field_in_real32(state, scale_iter, rc)
        type(ESMF_State),              intent(in   ) :: state
        type(StringReal32MapIterator), intent(in   ) :: scale_iter
        integer, optional,             intent(  out) :: rc

        type(ESMF_Field) :: field
        integer          :: status

        call ESMF_StateGet(state, scale_iter%key(), field, __RC__)
        call scale_field(field, scale_iter%value(), __RC__)

        _RETURN(_SUCCESS)
    end subroutine scale_field_in_real32

    subroutine scale_in_real32(this, state, rc)
        class(ScaleMapReal32), intent(inout) :: this
        type(ESMF_State),      intent(in   ) :: state
        integer, optional,     intent(  out) :: rc

        type(StringReal32MapIterator) :: iter
        integer                       :: status

        iter = this%begin()
        do while (iter /= this%end())
            call scale_field_in_real32(state, iter, __RC__)
            call iter%next()
        end do

        _RETURN(_SUCCESS)
    end subroutine scale_in_real32

    subroutine scale_field_in_real64(state, scale_iter, rc)
        type(ESMF_State),              intent(in   ) :: state
        type(StringReal64MapIterator), intent(in   ) :: scale_iter
        integer, optional,             intent(  out) :: rc

        type(ESMF_Field) :: field
        integer          :: status

        call ESMF_StateGet(state, scale_iter%key(), field, __RC__)
        call scale_field(field, scale_iter%value(), __RC__)

        _RETURN(_SUCCESS)
    end subroutine scale_field_in_real64

    subroutine scale_in_real64(this, state, rc)
        class(ScaleMapReal64), intent(inout) :: this
        type(ESMF_State),      intent(in   ) :: state
        integer, optional,     intent(  out) :: rc

        type(StringReal64MapIterator) :: iter
        integer                       :: status

        iter = this%begin()
        do while (iter /= this%end())
            call scale_field_in_real64(state, iter, __RC__)
            call iter%next()
        end do

        _RETURN(_SUCCESS)
    end subroutine scale_in_real64

    subroutine scale_field_out_real32(state, scale_iter, rc)
        type(ESMF_State),              intent(in   ) :: state
        type(StringReal32MapIterator), intent(in   ) :: scale_iter
        integer, optional,             intent(  out) :: rc

        type(ESMF_Field) :: field
        integer          :: status

        call ESMF_StateGet(state, scale_iter%key(), field, __RC__)
        call scale_field(field, 1.e0/scale_iter%value(), __RC__)

        _RETURN(_SUCCESS)
    end subroutine scale_field_out_real32

    subroutine scale_out_real32(this, state, rc)
        class(ScaleMapReal32), intent(inout) :: this
        type(ESMF_State),      intent(in   ) :: state
        integer, optional,     intent(  out) :: rc

        type(StringReal32MapIterator) :: iter
        integer                       :: status

        iter = this%begin()
        do while (iter /= this%end())
            call scale_field_out_real32(state, iter, __RC__)
            call iter%next()
        end do

        _RETURN(_SUCCESS)
    end subroutine scale_out_real32

    subroutine scale_field_out_real64(state, scale_iter, rc)
        type(ESMF_State),              intent(in   ) :: state
        type(StringReal64MapIterator), intent(in   ) :: scale_iter
        integer, optional,             intent(  out) :: rc

        type(ESMF_Field) :: field
        integer          :: status

        call ESMF_StateGet(state, scale_iter%key(), field, __RC__)
        call scale_field(field, 1.d0/scale_iter%value(), __RC__)

        _RETURN(_SUCCESS)
    end subroutine scale_field_out_real64

    subroutine scale_out_real64(this, state, rc)
        class(ScaleMapReal64), intent(inout) :: this
        type(ESMF_State),      intent(in   ) :: state
        integer, optional,     intent(  out) :: rc

        type(StringReal64MapIterator) :: iter
        integer                       :: status

        iter = this%begin()
        do while (iter /= this%end())
            call scale_field_out_real64(state, iter, __RC__)
            call iter%next()
        end do

        _RETURN(_SUCCESS)
    end subroutine scale_out_real64
end module ConvertUnits_private_mod