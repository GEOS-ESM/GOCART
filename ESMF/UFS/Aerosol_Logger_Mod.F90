module Aerosol_Logger_Mod

  use ESMF
  use MAPL, only: isAerosolLoggerOn => MAPL_am_I_root

  implicit none

  private

  public :: AerosolLog
  public :: AerosolLogStep
  public :: isAerosolLoggerOn

contains

  subroutine AerosolLog(msg, logUnit, rc)
    ! -- arguments
    character(len=*),  intent(in)  :: msg
    integer, optional, intent(in)  :: logUnit
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer                    :: localrc
    integer                    :: unit
    
    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    unit = 6
    if (present(logUnit)) unit = logUnit

    if (isAerosolLoggerOn()) then
      write(unit, '(a)') trim(msg)
    end if

    call ESMF_LogWrite(trim(msg), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolLog

  subroutine AerosolLogStep(clock, msg, logUnit, rc)
    ! -- arguments
    type(ESMF_Clock),           intent(in)  :: clock
    character(len=*), optional, intent(in)  :: msg
    integer,          optional, intent(in)  :: logUnit
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer                    :: localrc
    integer                    :: unit
    character(len=ESMF_MAXSTR) :: timeString
    character(:), allocatable  :: msgString
    type(ESMF_Time)            :: currTime
    type(ESMF_TimeInterval)    :: timeStep
    
    ! print timestep details

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    unit = 6
    if (present(logUnit)) unit = logUnit

    call ESMF_ClockGet(clock, currTime=currTime, timeStep=timeStep, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return  ! bail out

    call ESMF_TimeGet(currTime, timeStringISOFrac=timeString, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return  ! bail out

    msgString = timeString
    call ESMF_TimeGet(currTime + timeStep, timeStringISOFrac=timeString, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return  ! bail out

    msgString = 'Advancing from '//trim(msgString) // ' to ' // timeString

    if (present(msg)) msgString = trim(msg) // ': ' // trim(msgString)

    if (isAerosolLoggerOn()) then
      write(unit, '(a)') msgString
    end if

    call ESMF_LogWrite(msgString, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolLogStep

end module Aerosol_Logger_Mod
