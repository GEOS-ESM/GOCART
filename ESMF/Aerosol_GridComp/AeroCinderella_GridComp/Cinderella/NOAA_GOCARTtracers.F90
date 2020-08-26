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

   type, extends(StringVector) :: GOCARTtracers
      type(TracerMap) :: tracer_map
   contains
      procedure :: parse_yaml
   end type GOCARTtracers
contains
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

      bundle = ESMF_FieldBundleCreate(fieldList=tracers, name="inst_mass_tracers", __RC__)

      _RETURN(_SUCCESS)
   end subroutine create_tracer_bundle

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
               sub_sub_config = sub_iter%value()
!               call sub_sub_config%get(field_name)
!               field_name = sub_iter%value()
!               call this%push_back(field_name)

               call sub_iter%next()
            end do
         end select

         call iter%next()
      end do

      call file_stream%close()
   end subroutine parse_yaml
end module NOAA_GOCARTtracers_mod