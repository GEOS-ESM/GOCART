#include "MAPL_Generic.h"

module ReplenishAlarm
   use ESMF
   use MAPL


   implicit none
   private

   public :: createReplenishAlarm

   contains

   function createReplenishAlarm(gc,clock,freq,rc) result(alarm)
      type(ESMF_GridComp), intent(inout) :: gc
      type(ESMF_Clock), intent(in) :: clock
      integer, intent(in) :: freq
      integer, optional, intent(out) :: rc

      type(ESMF_Alarm) :: alarm

      integer :: status
      logical :: run_at_interval_start
      integer :: year, month, day
      integer :: hh, mm, ss, yyyymmdd, hhmmss
      integer :: reference_date, reference_time
      character(len=ESMF_MAXSTR) :: comp_name
      type(ESMF_Time) :: RingTime, currTime
      type(ESMF_TimeInterval) :: timeint, tstep

      ! this section mimics MAPL2 way to create the run alarm
      ! the goal is to have a consistent way of setting the proper
      ! offset, so that the alarm would run when the parent calls the children

      call MAPL_GridCompGetResource(gc, "RUN_AT_INTERVAL_START", run_at_interval_start, &
           default=.false., _RC)
      call ESMF_GridCompGet(gc, name=comp_name, _RC)

      call ESMF_ClockGet(clock, currTime=currTime, timestep=tstep, _RC)
      call ESMF_TimeIntervalSet(timeint, h=freq/10000, m=mod(freq,10000)/100, s=mod(freq,100), _RC)

      ! get current time from clock and create a referenace time with optonal override
      call ESMF_TimeGet(currTime, YY=YEAR, MM=MONTH, DD=DAY, &
           H=HH, M=MM, S=SS, _RC)

      yyyymmdd = year*10000 + month*100 + day
      hhmmss   = HH*10000 + MM*100 + SS

      !  Get Alarm reference date and time from resource, it defaults to midnight of the current day
      call MAPL_GridCompGetResource(gc, "REFERENCE_DATE", reference_date, &
           default=yyyymmdd, _RC)

      call MAPL_GridCompGetResource(gc, "REFERENCE_TIME", reference_time, &
           default=0, _RC)

      YEAR = reference_date/10000
      MONTH = mod(reference_date,10000)/100
      DAY = mod(reference_date,100)

      HH = reference_time/10000
      MM = mod(reference_time,10000)/100
      SS = mod(reference_time,100)

      call ESMF_TimeSet( ringTime, YY=YEAR, MM=MONTH, DD=DAY, &
           H=HH, M=MM, S=SS, _RC)

      if (ringTime > currTime) then
         ringTime = ringTime - (INT((ringTime - currTime)/TIMEINT)+1)*TIMEINT
      end if

      ! we back off current time with clock's dt since
      ! we advance the clock AFTER run method
      if (.not.run_at_interval_start) ringTime = ringTime-TSTEP

      ! make sure that ringTime is not in the past
      do while (ringTime < currTime)
         ringTime = ringTime + timeint
      end do

      alarm = ESMF_AlarmCreate(Clock = clock, &
           name = trim(comp_name) // "_ReplenishAlarm" , &
           RingInterval = timeint  ,  &
           RingTime     = ringTime,  &
           sticky       = .false.  ,  &
           _RC)
      if(ringTime == currTime) then
         call ESMF_AlarmRingerOn(alarm, _RC)
      end if

      _RETURN(ESMF_SUCCESS)
    end function createReplenishAlarm

end module ReplenishAlarm
