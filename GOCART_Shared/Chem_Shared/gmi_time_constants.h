
!=============================================================================
!
! $Id$
!
! CODE DEVELOPER
!   John Tannahill, LLNL
!   jrt@llnl.gov
!
! FILE
!   gmi_time_constants.h
!
! DESCRIPTION
!   This include file contains some time constants.
!
!   Note that "LY" = "leap year".
!
!   Mar 30, 2017: Moved this file from GmiInclude/ to Chem_Shared/ for TR
!
!=============================================================================


!     ------------------
!     Integer constants.
!     ------------------

      integer, parameter ::  &
     &  MINS_PER_HOUR    =  60,  &
     &  HRS_PER_DAY      =  24,  &
     &  DAYS_PER_WEEK    =   7,  &
     &  DAYS_PER_YEAR    = 365,  &
     &  DAYS_PER_YEAR_LY = 366

      integer, parameter ::  &
     &  MONTHS_PER_YEAR  =  12

      integer, parameter ::  &
     &  SECPMN    = 60,                          & ! seconds per minute
     &  SECPHR    = SECPMN * MINS_PER_HOUR,      & ! seconds per hour
     &  SECPDY    = SECPHR * HRS_PER_DAY,        & ! seconds per day
     &  SECPWK    = SECPDY * DAYS_PER_WEEK,      & ! seconds per week
     &  SECPYR    = SECPDY * DAYS_PER_YEAR,      & ! seconds per year
     &  SECPYR_LY = SECPDY * DAYS_PER_YEAR_LY  ! seconds per leap year


      integer, parameter :: DAYS_PER_MONTH   (MONTHS_PER_YEAR) =  &
     &  (/ 31, 28, 31, 30, 31, 30,  &
     &     31, 31, 30, 31, 30, 31 /)

      integer, parameter :: DAYS_PER_MONTH_LY(MONTHS_PER_YEAR) =  &
     &  (/ 31, 29, 31, 30, 31, 30,  &
     &     31, 31, 30, 31, 30, 31 /)


      integer, parameter :: START_DAY_OF_MONTH   (MONTHS_PER_YEAR) =  &
     &  (/   1,  32,  60,  91, 121, 152,  &
     &     182, 213, 244, 274, 305, 335 /)

      integer, parameter :: START_DAY_OF_MONTH_LY(MONTHS_PER_YEAR) =  &
     &  (/   1,  32,  61,  92, 122, 153,  &
     &     183, 214, 245, 275, 306, 336 /)


!     --------------------
!     Character constants.
!     --------------------

      character (len=3), parameter :: MONTH_NAME(MONTHS_PER_YEAR) =  &
     &  (/ 'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',  &
     &     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec' /)

