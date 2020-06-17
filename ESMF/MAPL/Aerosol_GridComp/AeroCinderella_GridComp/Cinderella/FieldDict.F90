module FieldDict_mod
    use FieldDictEntry_mod
    use FieldDictMap_mod

    use yaFyaml
    use gFTL_StringVector

    implicit none
    private

    public FieldDict
    public read_FieldDict

!    type, extends(FieldDictMap) :: FieldDict
!    end type FieldDict
contains
    subroutine read_base(config, name, field_entry, names)
        type(Configuration),  intent(in   ) :: config
        character(*),         intent(in   ) :: name
        type(FieldDictEntry), intent(  out) :: field_entry
        type(StringVector),   intent(  out) :: names

        character(:), allocatable :: units

        call names%push_back(name)

        units       = config%at('units')
        field_entry = init_FieldDictEntry(name, units)
    end subroutine read_base

    subroutine read_long_name(config, field_entry, names)
        type(Configuration),  intent(in   ) :: config
        type(FieldDictEntry), intent(inout) :: field_entry
        type(StringVector),   intent(inout) :: names

        character(:), allocatable :: long_name
        logical                   :: success

        long_name = config%at('long_name')

        if (long_name /= '') then
            call field_entry%add_long_name(long_name, success)

            if (success) call names%push_back(long_name)
        end if
    end subroutine read_long_name

    subroutine read_description(config, field_entry)
        type(Configuration),  intent(in   ) :: config
        type(FieldDictEntry), intent(inout) :: field_entry

        character(:), allocatable :: description

        description = config%at('description')

        if (description /= '') then
            call field_entry%add_description(description)
        end if

    end subroutine read_description

!    subroutine read_alias(config, field_entry, names)
!        type(Configuration),  intent(in   ) :: config
!        type(FieldDictEntry), intent(inout) :: field_entry
!        type(StringVector),   intent(inout) :: names
!
!        type(Configuration)         :: alias_config
!        type(ConfigurationIterator) :: iter
!        character(:), allocatable   :: alias
!        logical                     :: success
!
!        call config%get(alias_config, 'alias')
!        iter = alias_config%begin()
!        do while(iter /= alias_config%end())
!            alias = iter%get()
!
!            if (alias /= '') then
!                call field_entry%add_alias(alias, success)
!                if (success) call names%push_back(alias)
!            end if
!
!            call iter%next()
!        end do
!    end subroutine read_alias

    subroutine add_to_FieldDict(field_dict, field_entry, names)
        type(FieldDict),      intent(inout) :: field_dict
        type(FieldDictEntry), intent(in   ) :: field_entry
        type(StringVector),   intent(in   ) :: names

        type(StringVectorIterator) :: iter

        iter = names%begin()
        do while(iter /= names%end())
            call field_dict%insert(iter%get(), field_entry)
            call iter%next()
        end do

    end subroutine add_to_FieldDict

    subroutine read_FieldDict(file_name, field_dict)
        character(*),    intent(in   ) :: file_name
        type(FieldDict), intent(  out) :: field_dict

        type(Parser)                :: p
        type(Configuration)         :: main_config
        type(Configuration)         :: field_config
        type(Configuration)         :: entries_config
        type(Configuration)         :: entry_config
        type(Configuration)         :: alias_config
        type(ConfigurationIterator) :: iter, alias_iter

        type(StringVector)        :: names
        type(FieldDictEntry)      :: field_entry
        character(:), allocatable :: units, long_name, description, alias
        logical                   :: defined

        p           = Parser('core')
        main_config = p%load(FileStream(file_name))

        call main_config%get(field_config, 'field_dictionary')
        call field_config%get(entries_config, 'entries')

        iter = entries_config%begin()
        do while(iter /=entries_config%end())
            call entries_config%get(entry_config, iter%key())

            call read_base(entries_config, iter%key(), field_entry, names)
            call read_long_name(entries_config, field_entry, names)
            call read_description(entries_config, field_entry)
!            call read_alias(entries_config, field_entry, names)
            call add_to_FieldDict(field_dict, field_entry, names)

            call iter%next()
        end do

    end subroutine read_FieldDict
end module FieldDict_mod