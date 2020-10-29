#include "MAPL_Generic.h"

module NOAA_TracerEntryMod
   use ESMF
   use MAPL

   use yaFyaml
   use gFTL_StringVector

   implicit none
   private

   public :: TracerEntry

   character(*), parameter :: tracers_name   = 'name'
   character(*), parameter :: tracer_entries = 'entries'

   type :: TracerEntry
      character(:), allocatable :: name
      type(StringVector)        :: entries
   contains
      procedure :: read_tracer_entries
      procedure :: read_tracer_entry_config
   end type TracerEntry

contains
   subroutine read_tracer_entries(this, config)
      class(TracerEntry),  intent(inout) :: this
      type(Configuration), intent(inout) :: config

      type(ConfigurationIterator) :: iter
      character(:), allocatable   :: entry_name

      iter = config%begin()
      do while(iter /= config%end())
         entry_name = iter%get()
         call this%entries%push_back(entry_name)

         print*, "here"
         call iter%next()
      end do
   end subroutine read_tracer_entries

   subroutine read_tracer_entry_config(this, config)
      class(TracerEntry),  intent(inout) :: this
      type(Configuration), intent(inout) :: config

      type(Configuration)         :: sub_config
      type(ConfigurationIterator) :: iter
      character(:), pointer       :: key

      iter = config%begin()
      do while(iter /= config%end())
         key => iter%key()

         select case(key)
         case (tracers_name)
            this%name = iter%value()
         case (tracer_entries)
            sub_config = iter%value()

            call this%read_tracer_entries(sub_config)
         end select

         call iter%next()
      end do
   end subroutine read_tracer_entry_config
end module NOAA_TracerEntryMod

module NOAA_TracerMap
   use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

   use ESMF
   use MAPL
   use gFTL_StringIntegerMap
   use gFTL_StringVector

   use NOAA_TracerEntryMod

   implicit none
   private

   public :: TracerMap

   type :: TracerMap
      type(StringIntegerMap) :: tracer_map
   contains
      procedure, nopass :: read_tracer_name
      procedure, nopass :: remove_first_quote
      procedure, nopass :: remove_second_quote

      procedure :: read_field_table_line
      procedure :: read_field_table

      procedure :: get_tracer_array_real32
      procedure :: create_tracer_real32_3D
      procedure :: create_tracer_real32_4D
      procedure :: create_tracer_real32

      procedure :: get_tracer_array_real64
      procedure :: create_tracer_real64_3D
      procedure :: create_tracer_real64_4D
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
         call this%tracer_map%insert(tracer_name, tracer_index)
      end if
   end subroutine read_field_table_line

   subroutine read_field_table(this, filename)
      class(TracerMap), intent(inout) :: this
      character(*),     intent(in   ) :: filename

      character(len=ESMF_MAXSTR) :: str

      integer :: file_unit, iostat, tracer_index, str_index

      open(newunit=file_unit, file=filename, form='formatted',&
         access='sequential', status='old')

      tracer_index = 0
      do
         read(file_unit, '(a)', iostat=iostat) str
         if (iostat /= 0) exit

         call this%read_field_table_line(str, tracer_index)
      end do

      close(unit=file_unit)
   end subroutine read_field_table

   subroutine get_tracer_array_real32(this, field_array, entry_name, tracer_array, rc)
      class(TracerMap),           intent(inout) :: this
      real(kind=REAL32), pointer, intent(in   ) :: field_array(:,:,:,:)
      character(*),               intent(in   ) :: entry_name
      real(kind=REAL32), pointer, intent(  out) :: tracer_array(:,:,:)
      integer, optional,          intent(  out) :: rc

      integer :: tracer_size(4)
      integer :: idx, status

      idx =this%tracer_map%at(entry_name)

      ! check tracer index
      tracer_size = shape(field_array)
      _ASSERT(idx <= tracer_size(4), "invalid tracer index")

      tracer_array => field_array(:,:,:,idx)

      _RETURN(_SUCCESS)
   end subroutine get_tracer_array_real32

   subroutine create_tracer_real32_3D(this, field, name, entry_name, tracer, rc)
      class(TracerMap),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      character(*),      intent(in   ) :: entry_name
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      type(ESMF_Grid)            :: grid
      real(kind=REAL32), pointer :: field_array(:,:,:,:)
      real(kind=REAL32), pointer :: tracer_array(:,:,:)
      integer                    :: tracer_size(4)
      integer                    :: idx, status

      call ESMF_FieldGet(field,  grid=grid, __RC__)
      call ESMF_FieldGet(field,  localDE=0, farrayPtr=field_array, __RC__)

      call this%get_tracer_array_real32(field_array, entry_name, tracer_array, __RC__)

      tracer = ESMF_FieldCreate(grid, tracer_array, name=name, __RC__)

      _RETURN(_SUCCESS)
   end subroutine create_tracer_real32_3D

   subroutine create_tracer_real32_4D(this, field, name, entries, tracer, rc)
      class(TracerMap),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      type(TracerEntry), intent(in   ) :: entries
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      type(StringVectorIterator) :: iter

      type(ESMF_Grid)            :: grid
      real(kind=REAL32), pointer :: field_array(:,:,:,:)
      real(kind=REAL32), pointer :: tracer_array(:,:,:)
      real(kind=REAL32), pointer :: tracer_4D_array(:,:,:,:)
      integer                    :: idx, status

      call ESMF_FieldGet(field,  grid=grid, __RC__)
      call ESMF_FieldGet(field,  localDE=0, farrayPtr=field_array, __RC__)

      idx  = 1
      iter = entries%entries%begin()

      do while(iter /= entries%entries%end())
         tracer_array => tracer_4D_array(:,:,:,idx)
         call this%get_tracer_array_real32(field_array, iter%get(), tracer_array, __RC__)

         idx = idx + 1
         call iter%next()
      end do

      tracer = ESMF_FieldCreate(grid, tracer_4D_array, name=name, __RC__)

      _RETURN(_SUCCESS)
   end subroutine create_tracer_real32_4D

   subroutine create_tracer_real32(this, field, name, entries, tracer, rc)
      class(TracerMap),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      type(TracerEntry), intent(in   ) :: entries
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      character(:), allocatable :: entry_name
      integer                   :: status

      if (entries%entries%size() == 1) then
         entry_name = entries%entries%front()
         call this%create_tracer_real32_3D(field, name, entry_name, tracer, __RC__)
      else if (entries%entries%size() > 1) then
         call this%create_tracer_real32_4D(field, name, entries, tracer, __RC__)
      else
         _ASSERT(1 == 2, "No tracer entries found")
      end if

      _RETURN(_SUCCESS)
   end subroutine create_tracer_real32

   subroutine get_tracer_array_real64(this, field_array, entry_name, tracer_array, rc)
      class(TracerMap),           intent(inout) :: this
      real(kind=REAL64), pointer, intent(in   ) :: field_array(:,:,:,:)
      character(*),               intent(in   ) :: entry_name
      real(kind=REAL64), pointer, intent(  out) :: tracer_array(:,:,:)
      integer, optional,          intent(  out) :: rc

      integer :: tracer_size(4)
      integer :: idx, status

      idx =this%tracer_map%at(entry_name)

      ! check tracer index
      tracer_size = shape(field_array)
      _ASSERT(idx <= tracer_size(4), "invalid tracer index")

      tracer_array => field_array(:,:,:,idx)

      _RETURN(_SUCCESS)
   end subroutine get_tracer_array_real64

   subroutine create_tracer_real64_3D(this, field, name, entry_name, tracer, rc)
      class(TracerMap),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      character(*),      intent(in   ) :: entry_name
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      type(ESMF_Grid)            :: grid
      real(kind=REAL64), pointer :: field_array(:,:,:,:)
      real(kind=REAL64), pointer :: tracer_array(:,:,:)
      integer                    :: tracer_size(4)
      integer                    :: idx, status

      call ESMF_FieldGet(field,  grid=grid, __RC__)
      call ESMF_FieldGet(field,  localDE=0, farrayPtr=field_array, __RC__)

      call this%get_tracer_array_real64(field_array, entry_name, tracer_array, __RC__)

      tracer = ESMF_FieldCreate(grid, tracer_array, name=name, __RC__)

      _RETURN(_SUCCESS)
   end subroutine create_tracer_real64_3D

   subroutine create_tracer_real64_4D(this, field, name, entries, tracer, rc)
      class(TracerMap),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      type(TracerEntry), intent(in   ) :: entries
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      type(StringVectorIterator) :: iter

      type(ESMF_Grid)            :: grid
      real(kind=REAL64), pointer :: field_array(:,:,:,:)
      real(kind=REAL64), pointer :: tracer_array(:,:,:)
      real(kind=REAL64), pointer :: tracer_4D_array(:,:,:,:)
      integer                    :: idx, status

      call ESMF_FieldGet(field,  grid=grid, __RC__)
      call ESMF_FieldGet(field,  localDE=0, farrayPtr=field_array, __RC__)

      idx  = 1
      iter = entries%entries%begin()

      do while(iter /= entries%entries%end())
         tracer_array => tracer_4D_array(:,:,:,idx)
         call this%get_tracer_array_real64(field_array, iter%get(), tracer_array, __RC__)

         idx = idx + 1
         call iter%next()
      end do

      tracer = ESMF_FieldCreate(grid, tracer_4D_array, name=name, __RC__)

      _RETURN(_SUCCESS)
   end subroutine create_tracer_real64_4D

   subroutine create_tracer_real64(this, field, name, entries, tracer, rc)
      class(TracerMap),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      type(TracerEntry), intent(in   ) :: entries
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      character(:), allocatable :: entry_name
      integer                   :: status

      if (entries%entries%size() == 1) then
         entry_name = entries%entries%front()
         call this%create_tracer_real64_3D(field, name, entry_name, tracer, __RC__)
      else if (entries%entries%size() > 1) then
         call this%create_tracer_real64_4D(field, name, entries, tracer, __RC__)
      else
         _ASSERT(1 == 2, "No tracer entries found")
      end if

      _RETURN(_SUCCESS)
   end subroutine create_tracer_real64

   subroutine create_tracer(this, field, name, entries, tracer, rc)
      class(TracerMap),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      type(TracerEntry), intent(in   ) :: entries
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      type(ESMF_TypeKind_Flag) :: typekind
      integer                  :: status

      call ESMF_FieldGet(field, typekind=typekind, __RC__)

      if (typekind == ESMF_TYPEKIND_R4) then
         call this%create_tracer_real32(field, name, entries, tracer, __RC__)
      elseif (typekind == ESMF_TYPEKIND_R8) then
         call this%create_tracer_real64(field, name, entries, tracer, __RC__)
      else
         _FAIL("Unsupported ESMF_TYPEKIND")
      end if

      _RETURN(_SUCCESS)
   end subroutine create_tracer
end module NOAA_TracerMap

module NOAA_TracersMap
   use NOAA_TracerMap

#include "types/key_deferredLengthString.inc"
#define _value type(TracerMap)

#define _map TracersMap
#define _iterator TracersMapIterator
#define _alt
#include "templates/map.inc"
end module NOAA_TracersMap

module NOAA_TracersMod
   use ESMF
   use MAPL
   use yaFyaml

   use NOAA_TracerEntryMod
   use NOAA_TracerMap
   use NOAA_TracersMap

   implicit none
   private

   character(*),  parameter :: field_table = 'field_table'

   type :: NOAA_Tracers
      type(TracersMap) :: tracer_map

   contains
      procedure :: read_tracer_config
      procedure :: create_tracer
   end type NOAA_Tracers

contains
   subroutine read_tracer_config(this, name, config)
      class(NOAA_Tracers), intent(inout) :: this
      character(*),        intent(in   ) :: name
      type(Configuration), intent(inout) :: config

      type(ConfigurationIterator) :: iter
      character(:), pointer       :: key
      type(TracerMap)             :: tracer_map
      character(:), allocatable   :: filename

      iter = config%begin()
      do while(iter /= config%end())
         key => iter%key()

         select case(key)
         case (field_table)
            filename = iter%value()
            call tracer_map%read_field_table(filename)
         end select

         call this%tracer_map%insert(name, tracer_map)

      call iter%next()
      end do
   end subroutine read_tracer_config

   subroutine create_tracer(this, field, name, entries, tracer, rc)
      class(NOAA_Tracers),  intent(inout) :: this
      type(ESMF_Field),  intent(in   ) :: field
      character(*),      intent(in   ) :: name
      type(TracerEntry), intent(in   ) :: entries
      type(ESMF_Field),  intent(  out) :: tracer
      integer, optional, intent(  out) :: rc

      type(TracerMap)           :: tracer_map
      character(:), allocatable :: tracers_name
      integer                   :: status

      tracers_name = entries%name
      tracer_map   = this%tracer_map%at(tracers_name)

      call tracer_map%create_tracer(field, name, entries, tracer, __RC__)

      _RETURN(_SUCCESS)
   end subroutine create_tracer
end module NOAA_TracersMod
