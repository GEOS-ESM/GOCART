#include "MAPL_Generic.h"

module NOAA_GOCARTtracers_mod
   use, intrinsic :: iso_fortran_env, only: REAL32, REAL64

   use ESMF
   use MAPL
   use yaFyaml
   use gFTL_StringVector
   use NOAA_TracerMap_mod

   implicit none
   private

   public :: GOCARTtracers

   character(*), parameter :: rc_label     = 'GOCART_tracer_config:'
   character(*), parameter :: bundle_label = 'inst_mass_tracers'

   type, extends(StringVector) :: GOCARTtracers
      type(TracerMap) :: tracer_map
   contains
      procedure, nopass :: read_filename_from_config
      procedure :: parse_yaml
      procedure :: initialize

      procedure :: create_tracer_bundle
   end type GOCARTtracers
contains
   subroutine read_filename_from_config(config, filename, rc)
      type(ESMF_Config),         intent(inout) :: config
      character(:), allocatable, intent(  out) :: filename
      integer, optional,         intent(  out) :: rc

      character(len=ESMF_MaxStr) :: value
      integer                    :: status

      call ESMF_ConfigGetAttribute(config, value=value, label=rc_label, __RC__)
      filename = trim(value)

      _RETURN(_SUCCESS)
   end subroutine read_filename_from_config

   subroutine parse_yaml(this, filename)
      class(GOCARTtracers), intent(inout) :: this
      character(*),         intent(in   ) :: filename

      type(Parser)                :: p
      type(FileStream)            :: file_stream
      type(Configuration)         :: config, sub_config, sub_sub_config
      type(ConfigurationIterator) :: iter, sub_iter

      character(:), pointer     :: key
      character(:), allocatable :: field_table, field_name

      p      = Parser('core')
      file_stream = FileStream(filename)
      config = p%load(file_stream)

      iter = config%begin()
      do while (iter /= config%end())
         key => iter%key()

         select case(key)
         case('NOAA_field_table')
            field_table = iter%value()
            call this%tracer_map%read_field_table(field_table)
         case('GOCART_Tracers')
            sub_config = iter%value()
            sub_iter = sub_config%begin()
            do while(sub_iter /= sub_config%end())
               field_name = sub_iter%get()
               call this%push_back(field_name)

               call sub_iter%next()
            end do
         end select

         call iter%next()
      end do

      call file_stream%close()
   end subroutine parse_yaml

   subroutine initialize(this, config, rc)
      class(GOCARTtracers), intent(inout) :: this
      type(ESMF_Config),    intent(inout) :: config
      integer, optional,    intent(  out) :: rc

      character(:), allocatable :: filename
      integer :: status

      call this%read_filename_from_config(config, filename, __RC__)
      call this%parse_yaml(filename)
   end subroutine initialize

   subroutine create_tracer_bundle(this, field, bundle, rc)
      class(GOCARTtracers),   intent(inout) :: this
      type(ESMF_Field),       intent(in   ) :: field
      type(ESMF_FieldBundle), intent(  out) :: bundle
      integer, optional,      intent(  out) :: rc

      type(ESMF_Field), allocatable :: tracers(:)
      type(StringVectorIterator)    :: iter
      integer                       :: status, i

      allocate(tracers(this%size()))

      i    = 1
      iter = this%begin()
      do while(iter /= this%end())
         call this%tracer_map%create_tracer(field, iter%get(), tracers(i), __RC__)

         i = i + 1
         call iter%next()
      end do

      bundle = ESMF_FieldBundleCreate(fieldList=tracers, name=bundle_label, __RC__)

      _RETURN(_SUCCESS)
   end subroutine create_tracer_bundle
end module NOAA_GOCARTtracers_mod