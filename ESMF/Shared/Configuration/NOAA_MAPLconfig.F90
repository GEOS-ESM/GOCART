#include "MAPL_Generic.h"
#include "NUOPC_ErrLog.h"

module NOAA_MAPLpolicies
   use ESMF
   use NUOPC
   use MAPL
   use yaFyaml

   implicit none
   private

   public :: PoliciesConfig

   character(*), parameter :: TransferOfferGeomObject = 'TransferOfferGeomObject'
   character(*), parameter :: SharePolicyField        = 'SharePolicyField'
   character(*), parameter :: SharePolicyGeomObject   = 'SharePolicyGeomObject'

   character(*), parameter :: default_TransferOfferGeomObject = 'will provide'
   character(*), parameter :: default_SharePolicyField        = 'not share'

   type :: PoliciesConfig
      character(:), allocatable :: TransferOfferGeomObject
      character(:), allocatable :: SharePolicyField
      character(:), allocatable :: SharePolicyGeomObject
   contains
      procedure :: fill_defaults
      procedure :: read_policies_config
      procedure :: advertise_to_state
   end type PoliciesConfig
contains
   subroutine fill_defaults(this)
      class(PoliciesConfig), intent(inout) :: this

      if (.not. allocated(this%TransferOfferGeomObject)) this%TransferOfferGeomObject = default_TransferOfferGeomObject
      if (.not. allocated(this%SharePolicyField))        this%SharePolicyField        = default_SharePolicyField
      if (.not. allocated(this%SharePolicyGeomObject))   this%SharePolicyGeomObject   = this%SharePolicyField
   end subroutine fill_defaults

   subroutine read_policies_config(this, config)
      class(PoliciesConfig), intent(inout) :: this
      type(Configuration),   intent(in   ) :: config

      type(ConfigurationIterator) :: iter
      character(:), pointer       :: key

      iter = config%begin()
      do while(iter /= config%end())
         key => iter%key()

         select case(key)
         case(TransferOfferGeomObject)
            this%TransferOfferGeomObject = iter%value()
         case(SharePolicyField)
            this%SharePolicyField        = iter%value()
         case(SharePolicyGeomObject)
            this%SharePolicyGeomObject   = iter%value()
         end select

         call iter%next()
      end do

      call this%fill_defaults()
   end subroutine read_policies_config

   subroutine advertise_to_state(this, state, name, rc)
      class(PoliciesConfig), intent(inout) :: this
      type(ESMF_State),      intent(inout) :: state
      character(*),          intent(in   ) :: name
      integer, optional,     intent(  out) :: rc

      integer :: status

      call NUOPC_Advertise(state, standardName=name, &
         TransferOfferGeomObject=this%TransferOfferGeomObject, &
         SharePolicyField=this%SharePolicyField, &
         SharePolicyGeomObject=this%SharePolicyGeomObject, &
         rc=status)
      VERIFY_NUOPC_(status)

      _RETURN(_SUCCESS)
   end subroutine advertise_to_state
end module NOAA_MAPLpolicies

! module NOAA_MAPLoptions
!    use, intrinsic :: iso_fortran_env, only: REAL32, REAL64
!    use ESMF
!    use NUOPC
!    use MAPL
!    use yaFyaml

!    implicit none
!    private

!    public :: OptionsConfig

!    type :: OptionsConfig
!       character(:), allocatable :: state

!       logical :: reverse_levels
!       logical :: reverse_mask

!       real(kind=REAL64) :: convert2MAPL
!    end type OptionsConfig
! end module NOAA_MAPLoptions

! module NOAA_MAPLfieldConfig
!    use ESMF
!    use NUOPC
!    use MAPL
!    use yaFyaml

!    implicit none
!    private

!    public FieldConfig
!    public create_field_config

!    character(*), parameter :: standard_name           = 'standard_name'
!    character(*), parameter :: units                   = 'units'
!    character(*), parameter :: TransferOfferGeomObject = 'TransferOfferGeomObject'
!    character(*), parameter :: SharePolicyField        = 'SharePolicyField'
!    character(*), parameter :: SharePolicyGeomObject   = 'SharePolicyGeomObject'

!    character(*), parameter :: default_units                   = '1'
!    character(*), parameter :: default_TransferOfferGeomObject = 'will provide'
!    character(*), parameter :: default_SharePolicyField        = 'not share'

!    character(*), parameter :: cannot_provide = 'cannot provide'
!    character(*), parameter :: share          = 'share'

!    type :: FieldConfig
!       character(:), allocatable :: short_name
!       character(:), allocatable :: standard_name
!       character(:), allocatable :: units
!       character(:), allocatable :: TransferOfferGeomObject
!       character(:), allocatable :: SharePolicyField
!       character(:), allocatable :: SharePolicyGeomObject
!       logical                   :: doNotAllocate
!    contains
!       procedure :: add_to_field_dictionary
!       procedure :: add_name_to_field_dictionary
!       procedure :: create_synonyms
!       procedure :: advertise_to_state
!       procedure :: set_doNotAllocate
!       procedure :: fill_defaults
!       procedure :: read_field_config
!    end type FieldConfig

! contains
!    subroutine add_name_to_field_dictionary(this, name, rc)
!       class(FieldConfig), intent(inout) :: this
!       character(*),       intent(in   ) :: name
!       integer, optional,  intent(  out) :: rc

!       logical :: has_entry
!       integer :: status

!       has_entry = NUOPC_FieldDictionaryHasEntry(standardName=name, rc=status)
!       VERIFY_NUOPC_(status)

!       if (.not. has_entry) then
!          call NUOPC_FieldDictionaryAddEntry(standardName=name, canonicalUnits=units, rc=status)
!          VERIFY_NUOPC_(status)
!       end if

!       _RETURN(_SUCCESS)
!    end subroutine add_name_to_field_dictionary

!    subroutine create_synonyms(this, rc)
!       class(FieldConfig), intent(inout) :: this
!       integer, optional,  intent(  out) :: rc

!       logical :: are_synonymns
!       integer :: status

!       are_synonymns = NUOPC_FieldDictionaryMatchSyno(standardName1=this%short_name, &
!          standardName2=this%standard_name, rc=status)
!       VERIFY_NUOPC_(status)

!       if (.not. are_synonymns) then
!          call NUOPC_FieldDictionarySetSyno([this%short_name, this%standard_name], rc=status)
!          VERIFY_NUOPC_(status)
!       end if

!       _RETURN(_SUCCESS)
!    end subroutine create_synonyms

!    subroutine add_to_field_dictionary(this, rc)
!       class(FieldConfig), intent(inout) :: this
!       integer, optional,  intent(  out) :: rc

!       integer :: status

!       call this%add_name_to_field_dictionary(this%short_name, __RC__)
!       call this%add_name_to_field_dictionary(this%standard_name, __RC__)
!       call this%create_synonyms(__RC__)

!       _RETURN(_SUCCESS)
!    end subroutine add_to_field_dictionary

!    subroutine advertise_to_state(this, state, rc)
!       class(FieldConfig), intent(inout) :: this
!       type(ESMF_State),   intent(inout) :: state
!       integer, optional,  intent(  out) :: rc

!       integer :: status

!       call NUOPC_Advertise(state, standardName=this%short_name, &
!          TransferOfferGeomObject=this%TransferOfferGeomObject, &
!          SharePolicyField=this%SharePolicyField, &
!          SharePolicyGeomObject=this%SharePolicyGeomObject, &
!          rc=status)
!       VERIFY_NUOPC_(status)

!       _RETURN(_SUCCESS)
!    end subroutine advertise_to_state

!    subroutine set_doNotAllocate(this)
!       class(FieldConfig), intent(inout) :: this

!       if ((this%TransferOfferGeomObject == cannot_provide) .and. (this%SharePolicyField == share)) then
!          this%doNotAllocate = .true.
!       else
!          this%doNotAllocate = .false.
!       end if
!    end subroutine set_doNotAllocate

!    subroutine fill_defaults(this)
!       class(FieldConfig), intent(inout) :: this

!       if (.not. allocated(this%standard_name))           this%standard_name           = this%short_name
!       if (.not. allocated(this%units))                   this%units                   = default_units
!       if (.not. allocated(this%TransferOfferGeomObject)) this%TransferOfferGeomObject = default_TransferOfferGeomObject
!       if (.not. allocated(this%SharePolicyField))        this%SharePolicyField        = default_SharePolicyField
!       if (.not. allocated(this%SharePolicyGeomObject))   this%SharePolicyGeomObject   = this%SharePolicyField

!       call this%set_doNotAllocate()
!    end subroutine fill_defaults

!    subroutine read_field_config(this, short_name, config)
!       class(FieldConfig),  intent(inout) :: this
!       character(*),        intent(in   ) :: short_name
!       type(Configuration), intent(inout) :: config

!       type(ConfigurationIterator) :: iter
!       character(:), pointer       :: key

!       this%short_name = short_name

!       iter = config%begin()
!       do while (iter /= config%end())
!          key => iter%key()

!          select case(key)
!          case(standard_name)
!             this%standard_name = iter%value()
!          case(units)
!             this%units = iter%value()
!          case(TransferOfferGeomObject)
!             this%TransferOfferGeomObject = iter%value()
!          case(SharePolicyField)
!             this%SharePolicyField = iter%value()
!          case(SharePolicyGeomObject)
!             this%SharePolicyGeomObject = iter%value()
!          end select

!          call iter%next()
!       end do

!       call this%fill_defaults()
!    end subroutine read_field_config

!    function create_field_config(name, config) result(field_config)
!       type(FieldConfig) :: field_config
!       character(*),        intent(in   ) :: name
!       type(Configuration), intent(inout) :: config

!       call field_config%read_field_config(name, config)
!    end function create_field_config
! end module NOAA_MAPLfieldConfig

! module NOAA_MAPLfieldConfigMap
!    use NOAA_MAPLfieldConfig

! #include "types/key_deferredLengthString.inc"
! #define _value type(FieldConfig)

! #define _map FieldConfigMap
! #define _iterator FieldConfigMapIterator
! #define _alt
! #include "templates/map.inc"
! end module NOAA_MAPLfieldConfigMap

! module NOAA_MAPLconfigMod
!    use ESMF
!    use MAPL_Mod
!    use yaFyaml

!    use NOAA_MAPLfieldConfig
!    use NOAA_MAPLfieldConfigMap

!    implicit none
!    private

!    public NOAA_MAPLconfig

!    character(*), parameter :: NOAA_imports = 'NOAA_imports'
!    character(*), parameter :: NOAA_exports = 'NOAA_exports'

!    type :: NOAA_MAPLconfig
!       type(FieldConfigMap) :: imports
!       type(FieldConfigMap) :: exports
!    contains
!       procedure, nopass :: read_from_config
!       procedure         :: read_config
!       procedure, nopass :: advertise_map
!       procedure         :: advertise
!       procedure, nopass :: add_map_to_field_dictionary
!       procedure         :: add_to_field_dictionary
!    end type NOAA_MAPLconfig

! contains
!    subroutine add_map_to_field_dictionary(field_config_map, rc)
!       type(FieldConfigMap), intent(inout) :: field_config_map
!       integer, optional,    intent(  out) :: rc

!       type(FieldConfigMapIterator) :: iter
!       type(FieldConfig)            :: field_config

!       integer :: status

!       iter = field_config_map%begin()
!       do while(iter /= field_config_map%end())
!          field_config = iter%value()

!          call field_config%add_to_field_dictionary(__RC__)

!          call iter%next()
!       end do

!       _RETURN(_SUCCESS)
!    end subroutine add_map_to_field_dictionary

!    subroutine add_to_field_dictionary(this, rc)
!       class(NOAA_MAPLconfig), intent(inout) :: this
!       integer, optional,      intent(  out) :: rc

!       integer :: status

!       call this%add_map_to_field_dictionary(this%imports, __RC__)
!       call this%add_map_to_field_dictionary(this%exports, __RC__)

!       _RETURN(_SUCCESS)
!    end subroutine add_to_field_dictionary

!    subroutine advertise_map(field_config_map, state, rc)
!       class(FieldConfigMap), intent(inout) :: field_config_map
!       type(ESMF_State),      intent(inout) :: state
!       integer, optional,     intent(  out) :: rc

!       type(FieldConfigMapIterator) :: iter
!       type(FieldConfig)            :: field_config

!       integer :: status

!       iter = field_config_map%begin()
!       do while(iter /= field_config_map%end())
!          field_config = iter%value()

!          call field_config%advertise_to_state(state, __RC__)

!          call iter%next()
!       end do

!       _RETURN(_SUCCESS)
!    end subroutine advertise_map

!    subroutine advertise(this, import_state, export_state, rc)
!       class(NOAA_MAPLconfig), intent(inout) :: this
!       type(ESMF_State),       intent(inout) :: import_state
!       type(ESMF_State),       intent(inout) :: export_state
!       integer, optional,      intent(  out) :: rc

!       integer :: status

!       call this%advertise_map(this%imports, import_state, __RC__)
!       call this%advertise_map(this%exports, export_state, __RC__)

!       _RETURN(_SUCCESS)
!    end subroutine advertise

!    function read_from_config(config) result(field_config_map)
!       type(FieldConfigMap)               :: field_config_map
!       type(Configuration), intent(inout) :: config

!       type(Configuration)         :: sub_config
!       type(ConfigurationIterator) :: iter

!       character(:), pointer :: name

!       iter = config%begin()
!       do while(iter /= config%end())
!          name       => iter%key()
!          sub_config =  iter%value()

!          call field_config_map%insert(name, create_field_config(name, sub_config))
!          call iter%next()
!       end do
!    end function read_from_config

!    subroutine read_config(this, filename)
!       class(NOAA_MAPLconfig), intent(inout) :: this
!       character(*),            intent(in   ) :: filename

!       type(Parser)                :: p
!       type(FileStream)            :: file_stream
!       type(Configuration)         :: config, sub_config
!       type(ConfigurationIterator) :: iter

!       character(:), pointer :: key

!       p           = Parser('core')
!       file_stream = FileStream(filename)
!       config      = p%load(file_stream)

!       iter = config%begin()
!       do while (iter /= config%end())
!          key => iter%key()

!          select case(key)
!          case (NOAA_imports)
!             sub_config = iter%value()
!             this%imports = this%read_from_config(sub_config)
!          case (NOAA_exports)
!             sub_config = iter%value()
!             this%exports= this%read_from_config(sub_config)
!          end select

!          call iter%next()
!       end do

!       call file_stream%close()
!    end subroutine read_config
! end module NOAA_MAPLconfigMod
