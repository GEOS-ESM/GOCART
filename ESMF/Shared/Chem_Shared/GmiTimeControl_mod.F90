!--------------------------------------------------------------------------------
!
! Mar 30, 2017: Moved this file from GmiSupportingModules/ to Chem_Shared/ for TR
!
!--------------------------------------------------------------------------------
!BOP
!
! !MODULE: GmiTimeControl_mod
!
! !INTERFACE:
!
  module GmiTimeControl_mod
!
! !USES:
!
  implicit none
!
! !PUBLIC MEMBER FUNCTIONS:
!
  private
  public  ::  Set_leapYearFlag
! public  ::  GmiAdvanceClock
  public  ::  GetSecondsFromJanuary1
  public  ::  GetDaysFromJanuary1
  public  ::  ConvertTimeToSeconds
  public  ::  ConvertSecondstoTime
  public  ::  ConvertDateToSeconds
  public  ::  GmiSplitDateTime
!
  public  :: Get_begGmiDate
  public  :: Get_begGmiTime
  public  :: Get_endGmiDate
  public  :: Get_endGmiTime
  public  :: Get_curGmiDate
  public  :: Get_curGmiTime
! public  :: Get_gmiTimeStep
  public  :: Get_gmiSeconds
  public  :: Get_numTimeSteps
  public  :: Set_begGmiDate
  public  :: Set_begGmiTime
  public  :: Set_endGmiDate
  public  :: Set_endGmiTime
  public  :: Set_curGmiDate
  public  :: Set_curGmiTime
! public  :: Set_gmiTimeStep
  public  :: Set_gmiSeconds
  public  :: Set_numTimeSteps
!
! !PUBLIC DATA MEMBERS:
!
  public  ::  SIZTUSML
  public  ::  leapYearFlag
  public  :: t_GmiClock

  type t_GmiClock
     integer  :: begGmiDate    ! beginning date of the experiment year/month/day (YYYYMMDD)
     integer  :: begGmiTime    ! beginning time of the experiment hour/min/sec   (HHMMSS)
     integer  :: endGmiDate    ! end       date of the experiment year/month/day (YYYYMMDD)
     integer  :: endGmiTime    ! end       time of the experiment hour/min/sec   (HHMMSS)
!    real*8   :: gmiTimeStep   ! GMI model time step (s)
     integer  :: curGmiDate    ! current   date of the experiment year/month/day (YYYYMMDD)
                               ! This is updated at each iteration
     integer  :: curGmiTime    ! current   time of the experiment hour/min/sec   (HHMMSS)
                               ! This is updated at each iteration
     real*8   :: gmiSeconds    ! current integration time in seconds
     integer  :: numTimeSteps  ! current number of time steps
  end type t_GmiClock

  integer, parameter :: SIZTUSML = 16
  integer :: leapYearFlag ! Used to determine how leap years are derived
!
! !DESCRIPTION:
!
! !AUTHOR:
!   John Tannahill (jrt@llnl.gov) and Jules Kouatchou (kouatchou@gsfc.nasa.gov)
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!
      CONTAINS
!
!-------------------------------------------------------------------------
      subroutine Get_begGmiDate (self, ymd)
          integer         , intent(out) :: ymd
          type(t_GmiClock), intent(in ) :: self
          ymd = self%begGmiDate
          return
      end subroutine Get_begGmiDate
!-------------------------------------------------------------------------
      subroutine Set_begGmiDate (self, ymd)
          integer         , intent(in   ) :: ymd
          type(t_GmiClock), intent(inout) :: self
           self%begGmiDate = ymd
          return
      end subroutine Set_begGmiDate
!-------------------------------------------------------------------------
      subroutine Get_begGmiTime (self, hms)
          integer         , intent(out) :: hms
          type(t_GmiClock), intent(in ) :: self
          hms = self%begGmiTime
          return
      end subroutine Get_begGmiTime
!-------------------------------------------------------------------------
      subroutine Set_begGmiTime (self, hms)
          integer         , intent(in   ) :: hms
          type(t_GmiClock), intent(inout) :: self
           self%begGmiTime = hms
          return
      end subroutine Set_begGmiTime
!-------------------------------------------------------------------------
      subroutine Get_endGmiDate (self, ymd)
          integer         , intent(out) :: ymd
          type(t_GmiClock), intent(in ) :: self
          ymd = self%endGmiDate
          return
      end subroutine Get_endGmiDate
!-------------------------------------------------------------------------
      subroutine Set_endGmiDate (self, ymd)
          integer         , intent(in   ) :: ymd
          type(t_GmiClock), intent(inout) :: self
           self%endGmiDate = ymd
          return
      end subroutine Set_endGmiDate
!-------------------------------------------------------------------------
      subroutine Get_endGmiTime (self, hms)
          integer         , intent(out) :: hms
          type(t_GmiClock), intent(in ) :: self
          hms = self%endGmiTime
          return
      end subroutine Get_endGmiTime
!-------------------------------------------------------------------------
      subroutine Set_endGmiTime (self, hms)
          integer         , intent(in   ) :: hms
          type(t_GmiClock), intent(inout) :: self
           self%endGmiTime = hms
          return
      end subroutine Set_endGmiTime
!-------------------------------------------------------------------------
      subroutine Get_curGmiDate (self, ymd)
          integer         , intent(out) :: ymd
          type(t_GmiClock), intent(in ) :: self
          ymd = self%curGmiDate
          return
      end subroutine Get_curGmiDate
!-------------------------------------------------------------------------
      subroutine Set_curGmiDate (self, ymd)
          integer         , intent(in   ) :: ymd
          type(t_GmiClock), intent(inout) :: self
           self%curGmiDate = ymd
          return
      end subroutine Set_curGmiDate
!-------------------------------------------------------------------------
      subroutine Get_curGmiTime (self, hms)
          integer         , intent(out) :: hms
          type(t_GmiClock), intent(in ) :: self
          hms = self%curGmiTime
          return
      end subroutine Get_curGmiTime
!-------------------------------------------------------------------------
      subroutine Set_curGmiTime (self, hms)
          integer         , intent(in   ) :: hms
          type(t_GmiClock), intent(inout) :: self
           self%curGmiTime = hms
          return
      end subroutine Set_curGmiTime
!-------------------------------------------------------------------------
      subroutine Get_numTimeSteps (self, numTS)
          integer         , intent(out) :: numTS
          type(t_GmiClock), intent(in ) :: self
          numTS = self%numTimeSteps
          return
      end subroutine Get_numTimeSteps
!-------------------------------------------------------------------------
      subroutine Set_numTimeSteps (self, numTS)
          integer         , intent(in   ) :: numTS
          type(t_GmiClock), intent(inout) :: self
           self%numTimeSteps = numTS
          return
      end subroutine Set_numTimeSteps
!-------------------------------------------------------------------------
!     subroutine Get_gmiTimeStep (self, tdt)
!         real*8          , intent(out) :: tdt
!         type(t_GmiClock), intent(in ) :: self
!         tdt = self%gmiTimeStep
!         return
!     end subroutine Get_gmiTimeStep
!-------------------------------------------------------------------------
!     subroutine Set_gmiTimeStep (self, tdt)
!         real*8          , intent(in   ) :: tdt
!         type(t_GmiClock), intent(inout) :: self
!          self%gmiTimeStep = tdt
!         return
!     end subroutine Set_gmiTimeStep
!-------------------------------------------------------------------------
      subroutine Get_gmiSeconds (self, gmi_secs)
          real*8          , intent(out) :: gmi_secs
          type(t_GmiClock), intent(in ) :: self
          gmi_secs = self%gmiSeconds
          return
      end subroutine Get_gmiSeconds
!-------------------------------------------------------------------------
      subroutine Set_gmiSeconds (self, gmi_secs)
          real*8          , intent(in   ) :: gmi_secs
          type(t_GmiClock), intent(inout) :: self
           self%gmiSeconds = gmi_secs
          return
      end subroutine Set_gmiSeconds
!-------------------------------------------------------------------------
      subroutine Set_leapYearFlag (leap_year_flag)
          integer, intent(in) :: leap_year_flag
          leapYearFlag = leap_year_flag
          return
      end subroutine Set_leapYearFlag
!-------------------------------------------------------------------------
!! !BOP
!! !
!! ! !IROUTINE: GmiAdvanceClock
!! !
!! ! !INTERFACE:
!! !
!!   subroutine GmiAdvanceClock (self)
!! !
!!   implicit none
!! 
!! # include "gmi_time_constants.h"
!! 
!!    type(t_GmiClock), intent(inout) :: self
!! !
!! ! !LOCAL VARIABLES:
!!   integer :: nsec
!!   real*8  :: tdt             ! model time step (s)
!!   integer :: num_time_steps  ! number of time steps
!!   integer :: nymd            ! year/month/day (YYYYMMDD)
!!   integer :: nhms            ! hour/min/sec   (HHMMSS)
!!   real*8  :: gmi_sec         ! total Gmimod seconds (s)
!! !
!! ! !DESCRIPTION:
!! !  This routine steps the Gmimod clock one time step.
!! !
!! ! !AUTHOR:
!! !
!! ! !REVISION HISTORY:
!! !  Initial code.
!! !EOP
!! !-------------------------------------------------------------------------
!! !BOC
!! 
!!   call Get_curGmiDate  (self, nymd          )
!!   call Get_curGmiTime  (self, nhms          )
!!   call Get_numTimeSteps(self, num_time_steps)
!!   call Get_gmiSeconds  (self, gmi_sec       )
!!   call Get_gmiTimeStep (self, tdt           )   WARNING: tdt differs in Run1 vs Run2
!!  
!!   num_time_steps = num_time_steps + 1
!!   gmi_sec        = num_time_steps * tdt
!! 
!!   nsec = ConvertTimeToSeconds (nhms) + Nint (tdt)
!! 
!!   if (nsec >= SECPDY) then
!!      nsec = nsec - SECPDY
!!      nymd = IncrementDate (nymd, 1)
!!   end if
!!  
!!   if (nsec < 00000) then
!!      nsec = SECPDY + nsec
!!      nymd = IncrementDate (nymd, -1)
!!   end if
!!  
!!   nhms = ConvertSecondstoTime (nsec)
!!  
!!   call Set_curGmiDate  (self, nymd          )
!!   call Set_curGmiTime  (self, nhms          )
!!   call Set_numTimeSteps(self, num_time_steps)
!!   call Set_gmiSeconds  (self, gmi_sec       )
!! 
!!   return
!!  
!!   end subroutine GmiAdvanceClock
!! !
!! !EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetSecondsFromJanuary1
!
! !INTERFACE:
!
  subroutine GetSecondsFromJanuary1 (nsec_jan1, nymd, nhms)
!
  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in)   :: nymd      ! year/month/day (YYYYMMDD)
  integer, intent(in)   :: nhms      ! hour/min/sec   (HHMMSS)
!
! !OUTPUT PARAMETERS:
  integer, intent(out)  :: nsec_jan1 ! seconds from Jan. 1st (s)
!
! !DESCRIPTION:
!  This routine returns the offset in seconds from Jan. 1st, 
!  given any nymd/nhms.
!
! !LOCAL VARIABLES:
  integer :: nsec_ymd
  integer :: nsec_hms
!
! !AUTHOR:
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
  nsec_ymd = ConvertDateToSeconds (nymd)

  nsec_hms = ConvertTimeToSeconds (nhms) 

  nsec_jan1 = nsec_ymd + nsec_hms

  return

  end subroutine GetSecondsFromJanuary1
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GetDaysFromJanuary1
!
  subroutine GetDaysFromJanuary1 (nday_jan1, nymd)
!
  implicit none

# include "gmi_time_constants.h"
!
! !INPUT PARAMETERS:
  integer, intent(in)  :: nymd           ! year/month/day (YYYYMMDD)
!
! !OUTPUT PARAMETERS:
  integer, intent(out) :: nday_jan1      ! days from Jan. 1st (days)
!
! !DESCRIPTION:
!   This routine returns the offset in days from Jan. 1st, given any
!   nymd/nhms.
!
! !LOCAL VARIABLES:
  integer :: nsec_ymd
!
! !AUTHOR:
! !HISTORY:
!
!EOP 
!-----------------------------------------------------------------------------
!BOC
      nsec_ymd  = ConvertDateToSeconds (nymd)

      nday_jan1 = (nsec_ymd / SECPDY) + 1

      return

      end subroutine GetDaysFromJanuary1
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IncrementDate


!-----------------------------------------------------------------------------
!
! ROUTINE
!   IncrementDate
!
! DESCRIPTION
!   This routine increments/decrements nymd by one day.
!
! ARGUMENTS
!   nymd  : year/month/day (YYYYMMDD)
!   dyinc : the day increment/decrement (1 or -1)
!
!-----------------------------------------------------------------------------

  function IncrementDate (nymd, dyinc)

  use GmiPrintError_mod   , only : GmiPrintError
 
  implicit none
 
# include "gmi_time_constants.h"

! ----------------------
! Argument declarations.
! ----------------------

  integer :: nymd
  integer :: dyinc


! ----------------------
! Function declarations.
! ----------------------

  integer :: IncrementDate
 
! ----------------------
! Variable declarations.
! ----------------------

  character (len=75) :: err_msg

  logical :: is_ldy
  logical :: is_lyr

  integer :: ndy, nmon, nyr

! ----------------
! Begin execution.
! ----------------
 
  if ((dyinc /= -1) .and. (dyinc /= 1)) then
    err_msg = 'dyinc range error in IncrementDate.'
    call GmiPrintError (err_msg, .true., 1, dyinc, 0, 0, 0.0d0, 0.0d0)
  end if

  nyr  = nymd / 10000
  nmon = Mod (nymd, 10000) / 100
  ndy  = Mod (nymd, 100) + dyinc
 
  if (ndy == 0) then

!   -------------------------
!   Return to previous month.
!   -------------------------

    nmon = nmon - 1
 
    if (nmon == 0) then
      nmon = MONTHS_PER_YEAR
      nyr = nyr - 1
!      if (nyr == 0) then
!        nyr = 99
!      else
!        nyr = nyr - 1
!      end if
    end if
 
    is_lyr = IsLeapYear (nyr)

    if (is_lyr) then
      ndy = DAYS_PER_MONTH_LY(nmon)
    else
      ndy = DAYS_PER_MONTH(nmon)
    end if

  else

    is_lyr = IsLeapYear (nyr)

    if (is_lyr .and. (nmon == 2) .and. (ndy == DAYS_PER_MONTH_LY(nmon))) then
      is_ldy = .true.
    else
      is_ldy = .false.
    end if

    if (.not. is_ldy) then

      if (ndy > DAYS_PER_MONTH(nmon)) then

!       ----------------------
!       Proceed to next month.
!       ----------------------
 
        ndy  = 1
        nmon = nmon + 1
 
        if (nmon > MONTHS_PER_YEAR) then
          nmon = 1
          nyr = nyr + 1
!          if (nyr == 99) then
!            nyr = 0
!          else
!            nyr = nyr + 1
!          end if
        end if
 
      end if

    end if

  end if

  IncrementDate = (nyr * 10000) + (nmon * 100) + ndy
 
  return
 
  end function IncrementDate
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: IsLeapYear
 

!-----------------------------------------------------------------------------
!
! ROUTINE
!   IsLeapYear
!
! DESCRIPTION
!   This routine always returns false if leapYearFlag < 0,
!                always returns true  if leapYearFlag > 0,
!                and calculates whether or not it is really a leap year if
!                  leapYearFlag = 0.
!
! ARGUMENTS
!   the_year : the year to check (last 2 digits)
!
!-----------------------------------------------------------------------------

      function IsLeapYear (the_year)

      implicit none

!#     include "gmi_time_utils.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: the_year
 
!     ----------------------
!     Function declarations.
!     ----------------------

      logical :: IsLeapYear

!     ----------------------
!     Variable declarations.
!     ----------------------

!c    integer, save :: the_year00 = 1900
      integer, save :: the_year00 = 2000


!     ----------------
!     Begin execution.
!     ----------------

      IsLeapYear = .false.


      if (leapYearFlag > 0) then

        IsLeapYear = .true.

      else if (leapYearFlag == 0) then

        if (Mod (the_year, 4) == 0) then

          if ((the_year /= 0) .or. (Mod (the_year00, 400) == 0)) then

            IsLeapYear = .true.

          end if

        end if

      end if

      return

      end function IsLeapYear
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertTimeToSeconds


!-----------------------------------------------------------------------------
!
! ROUTINE
!   ConvertTimeToSeconds
!
! DESCRIPTION
!   This routine converts nhms to total seconds.
!
! ARGUMENTS
!   nhms : hour/min/sec (HHMMSS)
!
!-----------------------------------------------------------------------------

      function ConvertTimeToSeconds (nhms)

     use GmiPrintError_mod   , only : GmiPrintError
 
      implicit none

# include "gmi_time_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nhms

 
!     ----------------------
!     Function declarations.
!     ----------------------

      integer :: ConvertTimeToSeconds


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: nsec

!     ----------------
!     Begin execution.
!     ----------------

      nsec = (nhms / 10000) * SECPHR +              &
             (Mod (nhms, 10000) / 100) * SECPMN +   &
             Mod (nhms, 100)

      if (nsec >= SECPDY) then

        err_msg = 'nsec too big in ConvertTimeToSeconds.'
        call GmiPrintError (err_msg, .true., 1, nsec, 0, 0, 0.0d0, 0.0d0)

      end if

      ConvertTimeToSeconds = nsec

      return

      end function ConvertTimeToSeconds
 

!-----------------------------------------------------------------------------
!
! ROUTINE
!   ConvertSecondstoTime
!
! DESCRIPTION
!   This routine converts total seconds to hour/min/sec.
!
! ARGUMENTS
!   nsec : total seconds
!
!-----------------------------------------------------------------------------

  function ConvertSecondstoTime (nsec)

    use GmiPrintError_mod   , only : GmiPrintError
 
      implicit none

# include "gmi_time_constants.h"

!     ----------------------
!     Argument declarations.
!     ----------------------

      integer :: nsec

 
!     ----------------------
!     Function declarations.
!     ----------------------

      integer :: ConvertSecondstoTime


!     ----------------------
!     Variable declarations.
!     ----------------------

      character (len=75) :: err_msg

      integer :: nhms


!     ----------------
!     Begin execution.
!     ----------------

      if (nsec >= SECPDY) then

        err_msg = 'nsec too big in ConvertSecondstoTime.'
        call GmiPrintError (err_msg, .true., 1, nsec, 0, 0, 0.0d0, 0.0d0)

      end if

      nhms = (nsec / SECPHR) * 10000 +                 &
             (Mod (nsec, SECPHR) / SECPMN) * 100 +     &
             Mod (nsec, SECPMN)


      ConvertSecondstoTime = nhms

      return

      end function ConvertSecondstoTime
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: ConvertDateToSeconds
!
! !INTERFACE:
!
  function ConvertDateToSeconds (nymd)
!
! !USES:
  use GmiPrintError_mod   , only : GmiPrintError
!
  implicit none

# include "gmi_time_constants.h"
!
! !INPUT PARAMETERS:
  integer, intent(in)   :: nymd   ! year/month/day (YYYYMMDD)
!
! !RETURNED VALUE:
  integer :: ConvertDateToSeconds
!
! !DESCRIPTION:
!  This routine converts nymd to total seconds.
!
! !LOCAL VARIABLES:
  character (len=75) :: err_msg
  logical            :: is_lyr
  integer            :: idays, isecs
  integer            :: iyrsec
  integer            :: ndy, nmon, nyr
!
! !AUTHOR:
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC

  nyr  = nymd / 10000
  nmon = Mod (nymd, 10000) / 100
  ndy  = Mod (nymd, 100)
 
  is_lyr = IsLeapYear (nyr)

  if (is_lyr) then
     idays  = (START_DAY_OF_MONTH_LY(nmon) - 1) + (ndy - 1)
     iyrsec = SECPYR_LY
  else
     idays  = (START_DAY_OF_MONTH(nmon)    - 1) + (ndy - 1)
     iyrsec = SECPYR
  end if
 
  isecs = idays * SECPDY

  if (isecs >= iyrsec) then
     err_msg = 'isecs too big in ConvertDateToSeconds.'
     call GmiPrintError (err_msg, .true., 1, isecs, 0, 0, 0.0d0, 0.0d0)
  end if

  ConvertDateToSeconds = isecs

  return

  end function ConvertDateToSeconds
!EOC
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: GmiSplitDateTime 
!
! !INTERFACE:
!
  subroutine GmiSplitDateTime (nlmr, left2, middle2, right2)
!
  implicit none
!
! !INPUT PARAMETERS:
  integer, intent(in)   :: nlmr    ! nymd or nhms ("lmr" = "left" "middle" "right")
                                   ! (YYYYMMDD or HHMMSS)
!
! !OUTPUT PARAMETERS:
  integer, intent(out)  :: left2   ! left   field (year  or hours)
  integer, intent(out)  :: middle2 ! middle field (month or minutes)
  integer, intent(out)  :: right2  ! right  field (day   or seconds)
!
! !DESCRIPTION:
!  This routine extracts the three fields from either nymd 
!  (year, month, day) or nhms (hours, minutes, seconds).
!
! !AUTHOR:
!
! !REVISION HISTORY:
!  Initial code.
!
!EOP
!-------------------------------------------------------------------------
!BOC
 
  left2   = nlmr / 10000
  middle2 = Mod (nlmr, 10000) / 100
  right2  = Mod (nlmr, 100)

  return

  end subroutine GmiSplitDateTime
!EOC
!-------------------------------------------------------------------------

  end module GmiTimeControl_mod
 
