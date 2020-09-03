#include "MAPL_Generic.h"

module NOAA_TracerMap_mod
    use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

    use ESMF
    use MAPL
    use gFTL_StringIntegerMap

    implicit none
    private

    public :: TracerMap

    type, extends(StringIntegerMap) :: TracerMap
    contains
        procedure, nopass :: read_tracer_name
        procedure, nopass :: remove_first_quote
        procedure, nopass :: remove_second_quote

        procedure :: read_field_table_line
        procedure :: read_field_table

        procedure :: create_tracer_real32
        procedure :: create_tracer_real64
        procedure :: create_tracer
    end type TracerMap
contains
    function read_tracer_name(str) result(tracer_name)
        character(*), intent(in)  :: str
        character(:), allocatable :: tracer_name

        integer :: str_index

        str_index = scan(str, ',', back=.true.)
        tracer_name = str(str_index:)
    end function read_tracer_name

    function remove_first_quote(str) result(tracer_name)
        character(*), intent(in)  :: str
        character(:), allocatable :: tracer_name

        integer :: str_index

        str_index = scan(str, '"')
        tracer_name = str(str_index + 1:)
    end function remove_first_quote

    function remove_second_quote(str) result(tracer_name)
        character(*), intent(in)  :: str
        character(:), allocatable :: tracer_name

        integer :: str_index

        str_index = scan(str, '"', back=.true.)
        tracer_name = str(:str_index - 1)
    end function remove_second_quote

    subroutine read_field_table_line(this, str, tracer_index)
        class(TracerMap), intent(inout) :: this
        character(*),     intent(in   ) :: str
        integer,          intent(inout) :: tracer_index

        character(:), allocatable :: tracer_name

        if (index(str, 'TRACER') > 0) then
            tracer_index = tracer_index + 1

            tracer_name = this%remove_second_quote(&
                  this%remove_first_quote(this%read_tracer_name(str)))
            call this%insert(tracer_name, tracer_index)
        end if
    end subroutine read_field_table_line

    subroutine read_field_table(this, filename)
        class(TracerMap), intent(inout) :: this
        character(*),     intent(in   ) :: filename

        character(len=ESMF_MAXSTR) :: str

        integer :: file_unit, iostat, tracer_index, str_index

        open(newunit=file_unit, file=filename, form='formatted', access='sequential', status='old')

        tracer_index = 0
        do
            read(file_unit, '(a)', iostat=iostat) str
            if (iostat /= 0) exit

            call this%read_field_table_line(str, tracer_index)
        end do

        close(unit=file_unit)
    end subroutine read_field_table

    subroutine create_tracer_real32(this, field, name, tracer, rc)
        class(TracerMap),  intent(inout) :: this
        type(ESMF_Field),  intent(in   ) :: field
        character(*),      intent(in   ) :: name
        type(ESMF_Field),  intent(  out) :: tracer
        integer, optional, intent(  out) :: rc

        type(ESMF_Grid)            :: grid
        real(kind=REAL32), pointer :: field_array(:,:,:,:)
        real(kind=REAL32), pointer :: tracer_array(:,:,:)
        integer                    :: tracer_size(4)
        integer                    :: idx, status

        idx = this%at(name)

        call ESMF_FieldGet(field,  grid=grid, __RC__)
        call ESMF_FieldGet(field,  localDE=0, farrayPtr=field_array, __RC__)

        ! check tracer index
        tracer_size = shape(field_array)
        _ASSERT(idx <= tracer_size(4), "invalid tracer index")

        tracer_array => field_array(:,:,:,idx)
        tracer = ESMF_FieldCreate(grid, tracer_array, name=name, __RC__)

        _RETURN(_SUCCESS)
    end subroutine create_tracer_real32

    subroutine create_tracer_real64(this, field, name, tracer, rc)
        class(TracerMap),  intent(inout) :: this
        type(ESMF_Field),  intent(in   ) :: field
        character(*),      intent(in   ) :: name
        type(ESMF_Field),  intent(  out) :: tracer
        integer, optional, intent(  out) :: rc

        type(ESMF_Grid)            :: grid
        real(kind=REAL64), pointer :: field_array(:,:,:,:)
        real(kind=REAL64), pointer :: tracer_array(:,:,:)
        integer                    :: tracer_size(4)
        integer                    :: idx, status

        idx = this%at(name)

        call ESMF_FieldGet(field,  grid=grid, __RC__)
        call ESMF_FieldGet(field,  localDE=0, farrayPtr=field_array, __RC__)

        ! check tracer index
        tracer_size = shape(field_array)
        _ASSERT(idx <= tracer_size(4), "invalid tracer index")

        tracer_array => field_array(:,:,:,idx)
        tracer = ESMF_FieldCreate(grid, tracer_array, name=name, __RC__)

        _RETURN(_SUCCESS)
    end subroutine create_tracer_real64

    subroutine create_tracer(this, field, name, tracer, rc)
        class(TracerMap),  intent(inout) :: this
        type(ESMF_Field),  intent(in   ) :: field
        character(*),      intent(in   ) :: name
        type(ESMF_Field),  intent(  out) :: tracer
        integer, optional, intent(  out) :: rc

        type(ESMF_TypeKind_Flag) :: typekind
        integer                  :: status

        call ESMF_FieldGet(field, typekind=typekind, __RC__)

        if (typekind == ESMF_TYPEKIND_R4) then
            call this%create_tracer_real32(field, name, tracer, __RC__)
        elseif (typekind == ESMF_TYPEKIND_R8) then
            call this%create_tracer_real64(field, name, tracer, __RC__)
        else
            _FAIL("Unsupported ESMF_TYPEKIND")
        end if

        _RETURN(_SUCCESS)
    end subroutine create_tracer
end module NOAA_TracerMap_mod
