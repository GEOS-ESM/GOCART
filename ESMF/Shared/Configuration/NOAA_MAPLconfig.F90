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

module NOAA_MAPLoptions
   use, intrinsic :: iso_fortran_env, only: REAL32, REAL64
   use ESMF
   use NUOPC
   use MAPL
   use yaFyaml

   use LinearFields

   implicit none
   private

   public :: OptionsConfig

   character(*), parameter :: reverseLevels   = 'reverseLevels'
   character(*), parameter :: reverseMask     = 'reverseMask'
   character(*), parameter :: convert2MAPL    = 'convert2MAPL'

   integer,           parameter :: default_reverseLevels = 0
   logical,           parameter :: default_reverseMask   = .false.
   real(kind=REAL32), parameter :: default_convert2MAPL  = 1.0_REAL32

   type :: OptionsConfig
      integer,           allocatable :: reverseLevels
      logical,           allocatable :: reverseMask
      real(kind=REAL32), allocatable :: convert2MAPL
   contains
      procedure :: fill_defaults
      procedure :: read_options_config

      procedure :: reverse_levels
      procedure :: reverse_mask
      procedure :: convert_to_MAPL
      procedure :: convert_from_MAPL

      procedure :: to_MAPL
      procedure :: from_MAPL
   end type OptionsConfig
contains
   subroutine fill_defaults(this)
      class(OptionsConfig), intent(inout) :: this

      if (.not. allocated(this%reverseLevels)) this%reverseLevels = default_reverseLevels
      if (.not. allocated(this%reverseMask))   this%reverseMask   = default_reverseMask
      if (.not. allocated(this%convert2MAPL))  this%convert2MAPL  = default_convert2MAPL
   end subroutine fill_defaults

   subroutine read_options_config(this, config)
      class(OptionsConfig), intent(inout) :: this
      type(Configuration),  intent(in   ) :: config

      type(ConfigurationIterator) :: iter
      character(:), pointer       :: key

      integer           :: integer_tmp
      logical           :: logical_tmp
      real(kind=REAL32) :: real_tmp

      iter = config%begin()
      do while(iter /= config%end())
         key => iter%key()

         select case(key)
         case(reverseLevels)
            integer_tmp        = iter%value()
            this%reverseLevels = integer_tmp
         case(reverseMask)
            logical_tmp        = iter%value()
            this%reverseMask   = logical_tmp
         case(convert2MAPL)
            real_tmp           = iter%value()
            this%convert2MAPL  = real_tmp
         end select

         call iter%next()
      end do

      call this%fill_defaults()
   end subroutine read_options_config

   subroutine reverse_levels(this, field, rc)
      class(OptionsConfig), intent(inout) :: this
      type(ESMF_Field),     intent(in   ) :: field
      integer, optional,    intent(  out) :: rc

      integer :: status

      if (this%reverseLevels > 0) then
         call reverse_field(field, this%reverseLevels, __RC__)
      end if

      _RETURN(_SUCCESS)
   end subroutine reverse_levels

   subroutine reverse_mask(this, field, rc)
      class(OptionsConfig), intent(inout) :: this
      type(ESMF_Field),     intent(in   ) :: field
      integer, optional,    intent(  out) :: rc

      integer :: status

      if (this%reverseMask) then
         call shift_field(field, -1.0_REAL64, __RC__)
         call scale_field(field, -1.0_REAL64, __RC__)
      end if

      _RETURN(_SUCCESS)
   end subroutine reverse_mask

   subroutine convert_to_MAPL(this, field, rc)
      class(OptionsConfig), intent(inout) :: this
      type(ESMF_Field),     intent(in   ) :: field
      integer, optional,    intent(  out) :: rc

      integer :: status

      if (this%convert2MAPL /= 1.0_REAL32) then
         call scale_field(field, this%convert2MAPL, __RC__)
      end if

      _RETURN(_SUCCESS)
   end subroutine convert_to_MAPL

   subroutine convert_from_MAPL(this, field, rc)
      class(OptionsConfig), intent(inout) :: this
      type(ESMF_Field),     intent(in   ) :: field
      integer, optional,    intent(  out) :: rc

      integer :: status

      if (this%convert2MAPL /= 1.0_REAL32) then
         call scale_field(field, (1.0_REAL32 / this%convert2MAPL), __RC__)
      end if

      _RETURN(_SUCCESS)
   end subroutine convert_from_MAPL

   subroutine to_MAPL(this, field, rc)
      class(OptionsConfig), intent(inout) :: this
      type(ESMF_Field),     intent(in   ) :: field
      integer, optional,    intent(  out) :: rc

      integer :: status

      call this%reverse_levels( field, __RC__)
      call this%reverse_mask(   field, __RC__)
      call this%convert_to_MAPL(field, __RC__)

      _RETURN(_SUCCESS)
   end subroutine to_MAPL

   subroutine from_MAPL(this, field, rc)
      class(OptionsConfig), intent(inout) :: this
      type(ESMF_Field),     intent(in   ) :: field
      integer, optional,    intent(  out) :: rc

      integer :: status

      call this%reverse_levels( field, __RC__)
      call this%reverse_mask(   field, __RC__)
      call this%convert_from_MAPL(field, __RC__)

      _RETURN(_SUCCESS)
   end subroutine from_MAPL
end module NOAA_MAPLoptions

module NOAA_MAPLfield
   use, intrinsic :: iso_fortran_env, only: REAL32, REAL64
   use ESMF
   use NUOPC
   use MAPL
   use yaFyaml
   use LinearFields

   use NOAA_TracerEntry
   use NOAA_TracerMap
   use NOAA_MAPLpolicies
   use NOAA_MAPLoptions

   implicit none
   private

   public :: FieldConfig
   public :: create_FieldConfig

   character(*), parameter :: standardName = 'standardName'
   character(*), parameter :: units        = 'units'

   character(*), parameter :: MAPL_state  = 'MAPL_state'
   character(*), parameter :: NUOPC_state = 'NUOPC_state'

   character(*), parameter :: policies     = 'policies'
   character(*), parameter :: options      = 'options'
   character(*), parameter :: tracer       = 'tracer'

   character(*), parameter :: default_units = '1'

   type :: FieldConfig
      character(:), allocatable :: MAPL_name
      character(:), allocatable :: units
      character(:), allocatable :: standardName

      character(:), allocatable :: MAPL_state
      character(:), allocatable :: NUOPC_state

      type(PoliciesConfig), allocatable :: policies_config
      type(OptionsConfig),  allocatable :: options_config
      type(TracerEntry),    allocatable :: tracer_entry

      type(ESMF_Field), pointer :: MAPL_field  => null()
      type(ESMF_Field), pointer :: NUOPC_field => null()
   contains
      procedure :: fill_defaults
      procedure :: read_field_config

      procedure :: add_to_field_dictionary
      procedure :: add_name_to_field_dictionary
      procedure :: create_synonyms

      ! Convert, cast, recast...
      procedure :: copy_NUOPC_to_MAPL
      procedure :: copy_MAPL_to_NUOPC

      procedure :: advertise_as_MAPL
      procedure :: advertise_as_NUOPC

      ! Remove fields
      procedure :: initialize_fields
      procedure :: initialize_MAPL_field
      procedure :: initialize_NUOPC_field
   end type FieldConfig
contains
   function create_FieldConfig(MAPL_name, config, default_state) result(field_config)
      type(FieldConfig)                  :: field_config
      character(*),        intent(in   ) :: MAPL_name
      type(Configuration), intent(in   ) :: config
      character(*),        intent(in   ) :: default_state

      call field_config%read_field_config(MAPL_name, config, default_state)
   end function create_FieldConfig

   subroutine fill_defaults(this, default_state)
      class(FieldConfig), intent(inout) :: this
      character(*),       intent(in   ) :: default_state

      character(:), allocatable :: name

      type(PoliciesConfig) :: policies_config
      type(OptionsConfig)  :: options_config

      if (.not. allocated(this%units)) this%units = default_units

      if (.not. allocated(this%MAPL_state))  this%MAPL_state  = default_state
      if (.not. allocated(this%NUOPC_state)) this%NUOPC_state = default_state

      if (.not. allocated(this%policies_config)) then
         call policies_config%fill_defaults()
         this%policies_config = policies_config
      end if

      if (.not. allocated(this%options_config)) then
         call options_config%fill_defaults()
         this%options_config = options_config
      end if

      if (allocated(this%tracer_entry)) then
         name = this%tracer_entry%name
      else
         name = this%MAPL_name
      end if

      if (.not. allocated(this%standardName)) this%standardName = name
   end subroutine fill_defaults

   subroutine read_field_config(this, MAPL_name, config, default_state)
      class(FieldConfig),  intent(inout) :: this
      character(*),        intent(in   ) :: MAPL_name
      type(Configuration), intent(in   ) :: config
      character(*),        intent(in   ) :: default_state

      type(PoliciesConfig) :: policies_config
      type(OptionsConfig)  :: options_config
      type(TracerEntry)    :: tracer_entry

      type(Configuration)         :: sub_config
      type(ConfigurationIterator) :: iter
      character(:), pointer       :: key

      this%MAPL_name = MAPL_name

      iter = config%begin()
      do while(iter /= config%end())
         key => iter%key()

         select case(key)
         case (standardName)
            this%standardName = iter%value()
         case (units)
            this%units        = iter%value()
         case (MAPL_State)
            this%MAPL_state   = iter%value()
         case (NUOPC_State)
            this%NUOPC_state  = iter%value()
         case (policies)
            sub_config = iter%value()
            call policies_config%read_policies_config(sub_config)
            this%policies_config = policies_config
         case (options)
            sub_config = iter%value()
            call options_config%read_options_config(sub_config)
            this%options_config = options_config
         case (tracer)
            sub_config = iter%value()
            call tracer_entry%read_tracer_entry_config(sub_config)
            this%tracer_entry = tracer_entry
         end select

         call iter%next()
      end do

      call this%fill_defaults(default_state)
   end subroutine read_field_config

   subroutine add_name_to_field_dictionary(this, name, rc)
      class(FieldConfig), intent(inout) :: this
      character(*),       intent(in   ) :: name
      integer, optional,  intent(  out) :: rc

      logical :: has_entry
      integer :: status

      has_entry = NUOPC_FieldDictionaryHasEntry(standardName=name, rc=status)
      VERIFY_NUOPC_(status)

      if (.not. has_entry) then
         call NUOPC_FieldDictionaryAddEntry(standardName=name, canonicalUnits=this%units, rc=status)
         VERIFY_NUOPC_(status)
      end if

      _RETURN(_SUCCESS)
   end subroutine add_name_to_field_dictionary

   subroutine create_synonyms(this, rc)
      class(FieldConfig), intent(inout) :: this
      integer, optional,  intent(  out) :: rc

      logical :: are_synonymns
      integer :: status

      if (this%MAPL_name == this%standardName) then
         are_synonymns = .true.
      else
         are_synonymns = NUOPC_FieldDictionaryMatchSyno(standardName1=this%MAPL_name, &
            standardName2=this%standardName, rc=status)
         VERIFY_NUOPC_(status)
      end if

      if (.not. are_synonymns) then
         call NUOPC_FieldDictionarySetSyno([this%MAPL_name, this%standardName], rc=status)
         VERIFY_NUOPC_(status)
      end if

      _RETURN(_SUCCESS)
   end subroutine create_synonyms

   subroutine add_to_field_dictionary(this, rc)
      class(FieldConfig), intent(inout) :: this
      integer, optional,  intent(  out) :: rc

      integer :: status

      call this%add_name_to_field_dictionary(this%MAPL_name, __RC__)
      call this%add_name_to_field_dictionary(this%standardName, __RC__)
      call this%create_synonyms(__RC__)

      _RETURN(_SUCCESS)
   end subroutine add_to_field_dictionary

   subroutine copy_NUOPC_to_MAPL(this, rc)
      class(FieldConfig), intent(inout) :: this
      integer, optional,  intent(  out) :: rc

      integer :: status

      call copy_field(this%NUOPC_field, this%MAPL_field, __RC__)
      call this%options_config%to_MAPL(this%MAPL_field, __RC__)

      _RETURN(_SUCCESS)
   end subroutine copy_NUOPC_to_MAPL

   subroutine copy_MAPL_to_NUOPC(this, rc)
      class(FieldConfig), intent(inout) :: this
      integer, optional,  intent(  out) :: rc

      integer :: status

      call copy_field(this%MAPL_field, this%NUOPC_field, __RC__)
      call this%options_config%from_MAPL(this%NUOPC_field, __RC__)

      _RETURN(_SUCCESS)
   end subroutine copy_MAPL_to_NUOPC

   subroutine advertise_as_MAPL(this, state, rc)
      class(FieldConfig), intent(inout) :: this
      type(ESMF_State),   intent(inout) :: state
      integer, optional,  intent(  out) :: rc

      integer :: status

      call this%policies_config%advertise_to_state(state, this%MAPL_name, __RC__)

      _RETURN(_SUCCESS)
   end subroutine advertise_as_MAPL

   subroutine advertise_as_NUOPC(this, state, rc)
      class(FieldConfig), intent(inout) :: this
      type(ESMF_State),   intent(inout) :: state
      integer, optional,  intent(  out) :: rc

      integer :: status

      call this%policies_config%advertise_to_state(state, this%standardName, __RC__)

      _RETURN(_SUCCESS)
   end subroutine advertise_as_NUOPC

   subroutine initialize_MAPL_field(this, import_state, export_state, rc)
      class(FieldConfig), intent(inout) :: this
      type(ESMF_State),   intent(inout) :: import_state
      type(ESMF_State),   intent(inout) :: export_state
      integer, optional,  intent(  out) :: rc

      type(ESMF_Field), pointer :: field
      integer                   :: status

      select case(this%MAPL_state)
      case ('import')
         call ESMF_StateGet(import_state, itemName=this%MAPL_name, field=field, __RC__)
      case ('export')
         call ESMF_StateGet(export_state, itemName=this%MAPL_name, field=field, __RC__)
      case default
         _FAIL('Invalid MAPL_state given')
      end select

      this%MAPL_field => field

      _RETURN(_SUCCESS)
   end subroutine initialize_MAPL_field

   subroutine initialize_NUOPC_field(this, import_state, export_state, tracers, rc)
      class(FieldConfig), intent(inout) :: this
      type(ESMF_State),   intent(inout) :: import_state
      type(ESMF_State),   intent(inout) :: export_state
      type(TracerMap),    intent(inout) :: tracers
      integer, optional,  intent(  out) :: rc

      type(ESMF_Field), pointer :: field
      type(ESMF_Field), target  :: tracer
      integer                   :: status

      select case(this%NUOPC_state)
      case ('import')
         call ESMF_StateGet(import_state, itemName=this%standardName, field=field, __RC__)
      case ('export')
         call ESMF_StateGet(export_state, itemName=this%standardName, field=field, __RC__)
      case default
         _FAIL('Invalid NUOPC_state given')
      end select

      if (allocated(this%tracer_entry)) then
         call tracers%create_tracer(field, this%MAPL_name, this%tracer_entry, tracer, __RC__)
         this%NUOPC_field => tracer
      else
         this%NUOPC_field => field
      end if

      _RETURN(_SUCCESS)
   end subroutine initialize_NUOPC_field

   subroutine initialize_fields(this, import_state, export_state, tracers, rc)
      class(FieldConfig), intent(inout) :: this
      type(ESMF_State),   intent(inout) :: import_state
      type(ESMF_State),   intent(inout) :: export_state
      type(TracerMap),    intent(inout) :: tracers
      integer, optional,  intent(  out) :: rc

      integer :: status

      call this%initialize_MAPL_field( import_state, export_state,          __RC__)
      call this%initialize_NUOPC_field(import_state, export_state, tracers, __RC__)

      _RETURN(_SUCCESS)
   end subroutine initialize_fields
end module NOAA_MAPLfield

module NOAA_MAPLfieldMap
   use NOAA_MAPLfield

#include "types/key_deferredLengthString.inc"
#define _value type(FieldConfig)

#define _map FieldConfigMap
#define _iterator FieldConfigMapIterator
#define _alt
#include "templates/map.inc"
end module NOAA_MAPLfieldMap

module NOAA_MAPLconfigMod
   use ESMF
   use MAPL
   use yaFyaml

   use NOAA_MAPLfield
   use NOAA_MAPLfieldMap
   use NOAA_TracerMap

   implicit none
   private

   public :: NOAA_MAPLconfig
   public :: create_NOAA_MAPLconfig

   character(*), parameter :: imports = 'imports'
   character(*), parameter :: exports = 'exports'

   type :: NOAA_MAPLconfig
      type(FieldConfigMap) :: imports
      type(FieldConfigMap) :: exports

      type(TracerMap)      :: tracers
   contains
      procedure :: read_field_map_config
      procedure :: read_config

      procedure :: add_to_field_dictionary

      procedure :: advertise_as_MAPL
      procedure :: advertise_as_NUOPC

      procedure :: initialize_fields

      procedure :: copy_imports_NUOPC_to_MAPL
      procedure :: copy_imports_MAPL_to_NUOPC

      procedure :: copy_exports_NUOPC_to_MAPL
      procedure :: copy_exports_MAPL_to_NUOPC
   end type NOAA_MAPLconfig
contains
   function create_NOAA_MAPLconfig(config_filename, field_table_filename, rc) result(NOAA_MAPL_config)
      type(NOAA_MAPLconfig) :: NOAA_MAPL_config
      character(*),      intent(in   ) :: config_filename
      character(*),      intent(in   ) :: field_table_filename
      integer, optional, intent(  out) :: rc

      integer :: status

      type(Parser)                :: p
      type(FileStream)            :: file_stream
      type(Configuration)         :: config

      p           = Parser('core')
      file_stream = FileStream(config_filename)
      config      = p%load(file_stream)

      call NOAA_MAPL_config%read_config(config, __RC__)
      call NOAA_MAPL_config%tracers%read_field_table(field_table_filename)

      call file_stream%close()

      _RETURN(_SUCCESS)
   end function

   subroutine read_field_map_config(this, config, default_state, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      type(Configuration),    intent(in   ) :: config
      character(*),           intent(in   ) :: default_state
      integer, optional,      intent(  out) :: rc

      type(Configuration)         :: sub_config
      type(ConfigurationIterator) :: iter
      character(:), pointer       :: MAPL_name

      type(FieldConfig) :: field_config
      integer           :: status

      iter = config%begin()
      do while(iter /= config%end())
         MAPL_name  => iter%key()
         sub_config =  iter%value()

         field_config = create_FieldConfig(MAPL_name, sub_config, default_state)

         select case(default_state)
         case ('import')
            call this%imports%insert(MAPL_name, field_config)
         case ('export')
            call this%exports%insert(MAPL_name, field_config)
         case default
            _FAIL('Invalid default_state given')
         end select

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine read_field_map_config

   subroutine read_config(this, config, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      type(Configuration),    intent(in   ) :: config
      integer, optional,      intent(  out) :: rc

      type(Configuration)         :: sub_config
      type(ConfigurationIterator) :: iter
      character(:), pointer       :: key

      integer :: status

      iter = config%begin()
      do while(iter /= config%end())
         key => iter%key()

         select case(key)
         case (imports)
            sub_config = iter%value()
            call this%read_field_map_config(sub_config, 'import', __RC__)
         case (exports)
            sub_config = iter%value()
            call this%read_field_map_config(sub_config, 'export', __RC__)
         end select

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine read_config

   subroutine add_to_field_dictionary(this, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%imports%begin()
      do while(iter /= this%imports%end())
         field_config = iter%value()
         call field_config%add_to_field_dictionary(__RC__)

         call iter%next()
      end do

      iter = this%exports%begin()
      do while(iter /= this%exports%end())
         field_config = iter%value()
         call field_config%add_to_field_dictionary(__RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine add_to_field_dictionary

   subroutine advertise_as_MAPL(this, import_state, export_state, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      type(ESMF_state),       intent(inout) :: import_state
      type(ESMF_state),       intent(inout) :: export_state
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%imports%begin()
      do while(iter /= this%imports%end())
         field_config = iter%value()
         call field_config%advertise_as_MAPL(import_state, __RC__)

         call iter%next()
      end do

      iter = this%exports%begin()
      do while(iter /= this%exports%end())
         field_config = iter%value()
         call field_config%advertise_as_MAPL(export_state, __RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine advertise_as_MAPL

   subroutine advertise_as_NUOPC(this, import_state, export_state, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      type(ESMF_state),       intent(inout) :: import_state
      type(ESMF_state),       intent(inout) :: export_state
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%imports%begin()
      do while(iter /= this%imports%end())
         field_config = iter%value()
         call field_config%advertise_as_MAPL(import_state, __RC__)

         call iter%next()
      end do

      iter = this%exports%begin()
      do while(iter /= this%exports%end())
         field_config = iter%value()
         call field_config%advertise_as_MAPL(export_state, __RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine advertise_as_NUOPC

   subroutine initialize_fields(this, import_state, export_state, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      type(ESMF_state),       intent(inout) :: import_state
      type(ESMF_state),       intent(inout) :: export_state
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%imports%begin()
      do while(iter /= this%imports%end())
         field_config = iter%value()
         call field_config%initialize_fields(import_state, export_state, this%tracers, __RC__)

         call iter%next()
      end do

      iter = this%exports%begin()
      do while(iter /= this%exports%end())
         field_config = iter%value()
         call field_config%initialize_fields(import_state, export_state, this%tracers, __RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine initialize_fields

   subroutine copy_imports_NUOPC_to_MAPL(this, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%imports%begin()
      do while(iter /= this%imports%end())
         field_config = iter%value()
         call field_config%copy_NUOPC_to_MAPL(__RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine copy_imports_NUOPC_to_MAPL

   subroutine copy_exports_NUOPC_to_MAPL(this, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%exports%begin()
      do while(iter /= this%exports%end())
         field_config = iter%value()
         call field_config%copy_NUOPC_to_MAPL(__RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine copy_exports_NUOPC_to_MAPL

   subroutine copy_imports_MAPL_to_NUOPC(this, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%imports%begin()
      do while(iter /= this%imports%end())
         field_config = iter%value()
         call field_config%copy_MAPL_to_NUOPC(__RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine copy_imports_MAPL_to_NUOPC

   subroutine copy_exports_MAPL_to_NUOPC(this, rc)
      class(NOAA_MAPLconfig), intent(inout) :: this
      integer, optional,      intent(  out) :: rc

      type(FieldConfig)            :: field_config
      type(FieldConfigMapIterator) :: iter
      integer                      :: status

      iter = this%exports%begin()
      do while(iter /= this%exports%end())
         field_config = iter%value()
         call field_config%copy_MAPL_to_NUOPC(__RC__)

         call iter%next()
      end do

      _RETURN(_SUCCESS)
   end subroutine copy_exports_MAPL_to_NUOPC
end module NOAA_MAPLconfigMod
