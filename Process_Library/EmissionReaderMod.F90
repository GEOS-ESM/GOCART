#define __SUCCESS__ 0
#define __VERIFY__(x) if(x/=0) then; if(present(rc)) rc=x; return; endif
#define __RC__ rc=status); __VERIFY__(status
#define __STAT__ stat=status); __VERIFY__(status
#define __IOSTAT__ iostat=status); __VERIFY__(status
#define __RETURN__(x) if (present(rc)) rc=x; return
#define __ASSERT__(expr) if(.not. (expr)) then; if (present(rc)) rc=-1; return; endif

module EmissionReaderMod
   use, intrinsic :: iso_fortran_env, only: IOSTAT_END

   implicit none
   private

   public EmissionReader

   type :: EmissionReader
      private
      integer, allocatable :: unit
   contains
      procedure :: open
      procedure :: close
      procedure :: rewind => rewind_reader
      procedure :: is_end_marker
      procedure :: read_table
      procedure :: next_line
      procedure :: count_words
      procedure :: scan_to_label
      procedure :: get_dims
   end type EmissionReader
contains
   subroutine open(this, filename, rc)
      class(EmissionReader), intent(inout) :: this
      character(*), intent(in) :: filename
      integer, optional, intent(out) :: rc

      integer :: status

      __ASSERT__(.not. allocated(this%unit))
      allocate(this%unit)

      open(newunit=this%unit, file=filename,  &
           form='formatted', access = 'sequential', status='old', &
           action='read', __IOSTAT__)

      __RETURN__(__SUCCESS__)
   end subroutine open


   subroutine close(this, rc)
      class(EmissionReader), intent(inout) :: this
      integer, optional, intent(out) :: rc

      integer :: status

      __ASSERT__(allocated(this%unit))
      close(this%unit, __IOSTAT__)
      deallocate(this%unit)

   end subroutine close


   subroutine rewind_reader(this, rc)
      class(EmissionReader), intent(in) :: this
      integer, optional, intent(out) :: rc

      integer :: status

      __ASSERT__(allocated(this%unit))
      rewind(this%unit, __IOSTAT__)

      __RETURN__(__SUCCESS__)
   end subroutine rewind_reader

   function get_dims(this, label, rc) result(dims)
      integer :: dims(2)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: label
      integer, optional, intent(out) :: rc

      integer :: status
      logical :: eof
      character(:), allocatable :: line
      integer :: n_words

      call this%rewind(__RC__)
      call this%scan_to_label(label, __RC__)
!      print*,__FILE__,__LINE__, ' found label'

      dims = 0
      do
         line = this%next_line(eof=eof, __RC__)
         __ASSERT__(.not. eof)
         if (this%is_end_marker(line)) exit

         dims(2) = dims(2) + 1

         n_words = this%count_words(line)
         dims(1) = max(dims(1), n_words)
      end do

      __RETURN__(__SUCCESS__)
   end function get_dims

   integer function count_words(this, line) result(n_words)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: line

      integer :: idx, i0

      n_words = 0
      i0 = 0
      do
         ! scan to start of next word
         idx = verify(line(i0+1:), ' ')

         n_words = n_words + 1
         i0 = i0 + idx

         ! scan to end of current word
         idx = index(line(i0+1:), ' ')
         i0 = i0 + idx
         if (idx == 0) exit

      end do

      return
   end function count_words

   logical function is_end_marker(this, line)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: line

      is_end_marker = (line == '::')

   end function is_end_marker

   function read_table(this, label, rc) result(table)
      class(EmissionReader), intent(in) :: this
      real, allocatable :: table(:,:)
      character(*), intent(in) :: label
      integer, optional, intent(out) :: rc

      integer :: i, j
      integer :: dims(2)
      integer :: status
      logical :: eof
      character(:), allocatable :: line

      dims = this%get_dims(label, __RC__)
      call this%scan_to_label(label, __RC__)

      associate (n_words => dims(1), n_lines => dims(2))
        allocate(table(n_words, n_lines), __STAT__)

        do j = 1, n_lines
           line = this%next_line(eof=eof)
           __ASSERT__(.not. eof)

           read(line,*, iostat=status) (table(i,j),i=1,n_words)
           __VERIFY__(status)
        end do

      end associate

   end function read_table

   function next_line(this, eof, rc) result(line)
      character(:), allocatable :: line
      class(EmissionReader), intent(in) :: this
      logical, intent(out) :: eof
      integer, optional, intent(out) :: rc

      integer, parameter :: MAX_LINE_LEN=1024
      character(len=MAX_LINE_LEN) :: buffer
      integer :: idx
      integer :: status

      eof = .false.
      do

         read(this%unit,'(a)', iostat=status) buffer
         if (status == IOSTAT_END) then
            eof = .true.
            __RETURN__(__SUCCESS__)
         end if
         __VERIFY__(status)

         idx = index(buffer, '#')
         if (idx == 0) idx = len(buffer)

         line = trim(buffer(:idx-1))
         if (line /= '')  exit

      end do

      __RETURN__(__SUCCESS__)
   end function next_line

   subroutine scan_to_label(this, label, rc)
      class(EmissionReader), intent(in) :: this
      character(*), intent(in) :: label
      integer, optional, intent(out) :: rc

      integer :: status
      logical :: eof
      character(:), allocatable :: line

      call this%rewind(__RC__)
      do
         line = this%next_line(eof=eof, __RC__)
         if (line == label // '::') exit
      end do

      __RETURN__(__SUCCESS__)
   end subroutine scan_to_label

end module EmissionReaderMod
