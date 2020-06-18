#include "MAPL_Generic.h"

module ChemTracers_private_mod
    use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

    use ESMF
    use MAPL
    use gFTL_StringIntegerMap

    implicit none
    private

    public TracerMap

    public read_tracer_config
    public create_chem_tracer_real32
    public create_chem_tracer_real64
    public create_chem_tracer
    public create_tracer_bundle

    type, extends(StringIntegerMap) :: TracerMap
    contains
        procedure :: create_tracer_bundle
    end type TracerMap
contains
    subroutine read_tracer_config_option(config, option, rc)
        type(ESMF_Config),         intent(inout) :: config
        character(:), allocatable, intent(  out) :: option
        integer, optional,         intent(  out) :: rc

        character(len=ESMF_MaxStr) :: value
        integer                    :: status

        call ESMF_ConfigGetAttribute(config, value=value, label="NOAA_tracer_config:", __RC__)
        option = trim(value)//".inst_tracer_mass_frac:"

        _RETURN(_SUCCESS)
    end subroutine read_tracer_config_option

    subroutine read_tracer_config(config, tracer_map, rc)
        type(ESMF_Config), intent(inout) :: config
        type(TracerMap),   intent(  out) :: tracer_map
        integer, optional, intent(  out) :: rc

        character(:), allocatable  :: option
        integer                    :: i, num_tracers, columns, idx, status
        character(len=ESMF_MaxStr) :: name

        call read_tracer_config_option(config, option, __RC__)

        call ESMF_ConfigGetDim(config, label=option, lineCount=num_tracers, &
                columnCount=columns, __RC__)
        call ESMF_ConfigFindLabel(config, label=option, __RC__)
        call ESMF_ConfigNextLine(config, __RC__)

        do i=1, num_tracers
            call ESMF_ConfigGetAttribute(config, value=name, __RC__)
            call ESMF_ConfigGetAttribute(config, value=idx,  __RC__)

            call tracer_map%insert(trim(name), idx)

            call ESMF_ConfigNextLine(config, __RC__)
        end do

        _RETURN(_SUCCESS)
    end subroutine read_tracer_config

    subroutine create_chem_tracer_real32(field, name, idx, tracer, rc)
        type(ESMF_Field),  intent(in   ) :: field
        character(*),      intent(in   ) :: name
        integer,           intent(in   ) :: idx
        type(ESMF_Field),  intent(  out) :: tracer
        integer, optional, intent(  out) :: rc

        type(ESMF_Grid)            :: grid
        real(kind=REAL32), pointer :: field_array(:,:,:,:)
        real(kind=REAL32), pointer :: tracer_array(:,:,:)
        integer                    :: tracer_size(4)
        integer                    :: status

        call ESMF_FieldGet(field, grid=grid, __RC__)
        call ESMF_FieldGet(field, localDE=0, farrayPtr=field_array, __RC__)

        ! check tracer index
        tracer_size = shape(field_array)
        _ASSERT(idx <= tracer_size(4), "Invalid tracer index")

        tracer_array => field_array(:,:,:,idx)
        tracer =  ESMF_FieldCreate(grid, tracer_array, name=name, __RC__)

        _RETURN(_SUCCESS)
    end subroutine create_chem_tracer_real32

    subroutine create_chem_tracer_real64(field, name, idx, tracer, rc)
        type(ESMF_Field),  intent(in   ) :: field
        character(*),      intent(in   ) :: name
        integer,           intent(in   ) :: idx
        type(ESMF_Field),  intent(  out) :: tracer
        integer, optional, intent(  out) :: rc

        type(ESMF_Grid)            :: grid
        real(kind=REAL64), pointer :: field_array(:,:,:,:)
        real(kind=REAL64), pointer :: tracer_array(:,:,:)
        integer                    :: tracer_size(4)
        integer                    :: status

        call ESMF_FieldGet(field, grid=grid, __RC__)
        call ESMF_FieldGet(field, localDE=0, farrayPtr=field_array, __RC__)

        ! check tracer index
        tracer_size = shape(field_array)
        _ASSERT(idx <= tracer_size(4), "Invalid tracer index")

        tracer_array => field_array(:,:,:,idx)
        tracer =  ESMF_FieldCreate(grid, tracer_array, name=name, __RC__)

        _RETURN(_SUCCESS)
    end subroutine create_chem_tracer_real64

    subroutine create_chem_tracer(field, name, idx, tracer, rc)
        type(ESMF_Field),  intent(in   ) :: field
        character(*),      intent(in   ) :: name
        integer,           intent(in   ) :: idx
        type(ESMF_Field),  intent(  out) :: tracer
        integer, optional, intent(  out) :: rc

        type(ESMF_TypeKind_Flag) :: typekind
        integer                  :: status

        call ESMF_FieldGet(field, typekind=typekind, __RC__)

        if (typekind==ESMF_TYPEKIND_R4) then
            call create_chem_tracer_real32(field, name, idx, tracer, __RC__)
        elseif (typekind==ESMF_TYPEKIND_R8) then
            call create_chem_tracer_real64(field, name, idx, tracer, __RC__)
        else
            _FAIL("Unsupported ESMF_TYPEKIND")
        end if

        _RETURN(_SUCCESS)
    end subroutine create_chem_tracer

    subroutine create_tracer_bundle(this, field, bundle, rc)
        class(TracerMap),       intent(inout) :: this
        type(ESMF_Field),       intent(in   ) :: field
        type(ESMF_FieldBundle), intent(  out) :: bundle
        integer, optional,      intent(  out) :: rc

        type(ESMF_Field), allocatable  :: tracers(:)
        type(StringIntegerMapIterator) :: iter
        integer                        :: status, i

        allocate(tracers(this%size()))
        i    = 1
        iter = this%begin()
        do while (iter /= this%end())
            call create_chem_tracer(field, iter%key(), iter%value(), tracers(i), __RC__)
            i = i + 1
            call iter%next()
        end do

        bundle = ESMF_FieldBundleCreate(fieldList=tracers, name="inst_mass_tracers", __RC__)

        _RETURN(_SUCCESS)
    end subroutine create_tracer_bundle
end module ChemTracers_private_mod