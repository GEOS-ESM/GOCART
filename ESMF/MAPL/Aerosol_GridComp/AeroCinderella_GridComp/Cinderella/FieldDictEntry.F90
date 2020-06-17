module FieldDictEntry_mod
    use gFTL_StringVector

    implicit none
    private

    public FieldDictEntry
    public init_FieldDictEntry

    type FieldDictEntry
        character(:), allocatable :: standard_name
        character(:), allocatable :: units

        character(:), allocatable :: long_name
        character(:), allocatable :: description
!        type(StringVector)        :: alias
    contains
        procedure :: add_long_name
        procedure :: add_description
!        procedure :: add_alias
    end type FieldDictEntry
contains
    function init_FieldDictEntry(standard_name, units) result(dict_entry)
        character(*), intent(in) :: standard_name
        character(*), intent(in) :: units
        type(FieldDictEntry)     :: dict_entry

        dict_entry%standard_name = standard_name
        dict_entry%units         = units
        dict_entry%long_name     = standard_name

!        call dict_entry%alias%push_back(standard_name)
    end function init_FieldDictEntry

    subroutine add_long_name(this, long_name, success)
        class(FieldDictEntry), intent(inout) :: this
        character(*),          intent(in   ) :: long_name
        logical,               intent(  out) :: success

        if (long_name /= this%long_name) then
            this%long_name = long_name
!            call this%alias%push_back(long_name)

            success = .true.
        else
            success = .false.
        end if
    end subroutine add_long_name

    subroutine add_description(this, description)
        class(FieldDictEntry), intent(inout) :: this
        character(*),          intent(in   ) :: description

        if (.not. allocated(this%description)) then
            this%description = description
        end if
    end subroutine add_description

!    subroutine add_alias(this, alias, success)
!        class(FieldDictEntry), intent(inout) :: this
!        character(*),          intent(in   ) :: alias
!        logical,               intent(  out) :: success
!
!        type(StringVectorIterator) :: iter
!
!        success = .true.
!        iter    = this%alias%begin()
!        do while (iter /= this%alias%end())
!            if (iter%get() == alias) then
!                success = .false.
!                exit
!            end if
!            call iter%next()
!        end do
!
!        if (success) then
!            call this%alias%push_back(alias)
!        end if
!    end subroutine add_alias
end module FieldDictEntry_mod