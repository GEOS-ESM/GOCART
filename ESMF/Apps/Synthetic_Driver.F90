#include "MAPL_Generic.h"
#include "NUOPC_ErrLog.h"

module synthetic_driver

  !-----------------------------------------------------------------------------
  ! Code that specializes generic ESM Component code.
  !-----------------------------------------------------------------------------
  
  use ESMF
  use NUOPC
  use NUOPC_Driver, &
       driver_routine_SS             => SetServices, &
       driver_label_SetModelServices => label_SetModelServices, &
       driver_label_SetRunSequence   => label_SetRunSequence

  use MAPL
  use MAPL_NUOPCWrapperMod, only: wrapper_ss => SetServices, init_wrapper
  use NUOPC_Connector, only: cplSS => SetServices

  use agcm_ctm_mediator, only: mediator_set_services => SetServices
  use Provider_GridCompMod, only: provider_set_services => SetServices
  use Reciever_GridCompMod, only: reciever_set_services => SetServices

  use gFTL_StringVector
  
  implicit none

  private

  public SetServices

contains

  subroutine SetServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc
    type(ESMF_Config) :: config

    type(StringVector) :: agcm_exports
    type(StringVectorIterator) :: iter
    rc = ESMF_SUCCESS

    print*, "Driver start Set Services"

    ! NUOPC_Driver registers the generic methods
    print*,"Driver add Generic Set Services"
    call NUOPC_CompDerive(driver, driver_routine_SS, rc=rc)
    VERIFY_NUOPC_(rc)

    ! attach specializing method(s)
    print*,"Driver add Set Model Services"
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetModelServices, &
         specRoutine=SetModelServices, rc=rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver add Set Run Sequence"
    call NUOPC_CompSpecialize(driver, specLabel=driver_label_SetRunSequence, &
         specRoutine=SetRunSequence, rc=rc)
    VERIFY_NUOPC_(rc)


    print*,"Driver create config"
    config = ESMF_ConfigCreate(rc=rc)
    VERIFY_NUOPC_(rc)
    print*,"Driver read config"
    call ESMF_ConfigLoadFile(config, "NUOPC_run_config.txt", rc=rc)
    VERIFY_NUOPC_(rc)
    print*,"Driver add config"
    call ESMF_GridCompSet(driver, config=config, rc=rc)
    VERIFY_NUOPC_(rc)

    print*, "Driver finish Set Services"
  end subroutine SetServices

  !-----------------------------------------------------------------------------

  subroutine SetModelServices(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    type(ESMF_GridComp)           :: agcm, ctm, mediator
    type(ESMF_CplComp)            :: connector, connector2
    type(ESMF_VM) :: vm
    type(ESMF_Config) :: config

    logical :: seq
    integer :: i, npes, n_agcm_pes, n_ctm_pes
    integer, allocatable :: agcm_petlist(:), ctm_petlist(:)

    print*, "Driver start Set Model Services"
    rc = ESMF_SUCCESS

    print*,"Driver set clock"
    call set_clock(driver)

    print*,"Driver read config"
    call ESMF_GridCompGet(driver, vm = vm, config = config, rc = rc)
    VERIFY_NUOPC_(rc)
    print*,"Driver read npes"
    call ESMF_VMGet(vm, petCount  = npes, rc = rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver read from config"
    call ESMF_ConfigGetAttribute(config, seq, label = "sequential:", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigGetAttribute(config, n_agcm_pes, label = "agcm_pets:", rc = rc)
    VERIFY_NUOPC_(rc)
    call ESMF_ConfigGetAttribute(config, n_ctm_pes, label = "ctm_pets:", rc = rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver create pet lists"
    allocate(agcm_petlist(n_agcm_pes), ctm_petlist(n_ctm_pes))

    if (seq) then
        _ASSERT((n_agcm_pes == n_ctm_pes), "agcm_pets must be equal to ctm_pets in sequential")
        _ASSERT((n_agcm_pes == npes), "agcm_pets must be equal to number of pets in sequential")

        agcm_petlist = [(i, i = 0, n_agcm_pes - 1)]
        ctm_petlist = [(i, i = 0, n_ctm_pes - 1)]
    else
        _ASSERT(((n_agcm_pes + n_ctm_pes) == npes), "agcm_pets + ctm_pets must be equal to number of pets")

        agcm_petlist = [(i, i = 0, n_agcm_pes - 1)]
        ctm_petlist = [(i, i = n_agcm_pes, npes - 1)]
    end if

    print*,"Driver add provider"
    call NUOPC_DriverAddComp(driver, "agcm", wrapper_ss, comp = agcm, &
         petlist = agcm_petlist, rc = rc)
    VERIFY_NUOPC_(rc)
    print*,"Driver wrap provider MAPL"
    call init_wrapper(wrapper_gc = agcm, name = "agcm", &
         cap_rc_file = "AGCM_CAP.rc", root_set_services =provider_set_services, rc = rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver add reciever"
    call NUOPC_DriverAddComp(driver, "ctm", wrapper_ss, comp = ctm, &
         petlist = ctm_petlist, rc = rc)
    VERIFY_NUOPC_(rc)
    print*,"Driver wrap reciever MAPL"
    call init_wrapper(wrapper_gc = ctm, name = "ctm", &
         cap_rc_file = "CTM_CAP.rc", root_set_services = reciever_set_services, rc = rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver add mediator"
    call NUOPC_DriverAddComp(driver, "mediator", mediator_set_services, comp = mediator, &
         petlist = ctm_petlist, rc = rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver connect provider to mediator"
    call NUOPC_DriverAddComp(driver, srcCompLabel = "agcm", dstCompLabel = "mediator", &
         compSetServicesRoutine = cplSS, comp = connector, rc = rc)
    VERIFY_NUOPC_(rc)
    print*,"Driver connect mediator to reciever"
    call NUOPC_DriverAddComp(driver, srcCompLabel = "mediator", dstCompLabel = "ctm", &
         compSetServicesRoutine = cplSS, comp = connector, rc = rc)
    VERIFY_NUOPC_(rc)

    print*, "Driver finish Set Model Services"
  end subroutine SetModelServices

  subroutine set_clock(driver)
    type(ESMF_GridComp), intent(inout) :: driver

    type(ESMF_Time)               :: startTime
    type(ESMF_Time)               :: stopTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Clock)              :: internalClock
    type(ESMF_Config) :: config
    integer :: start_date_and_time(2), end_date_and_time(2), dt, file_unit, yy, mm, dd, h, m, s, rc


    open(newunit = file_unit, file = "cap_restart", form = 'formatted', &
         status = 'old', action = 'read')
    read(file_unit, *) start_date_and_time
    close(file_unit)

    call UnpackDateTime(start_date_and_time, yy, mm, dd, h, m, s)

    call ESMF_TimeSet(startTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
         calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompGet(driver, config = config, rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(config, end_date_and_time(1), label = "end_date:", rc = rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(config, end_date_and_time(2), label = "end_time:", rc = rc)
    VERIFY_NUOPC_(rc)

    call UnpackDateTime(end_date_and_time, yy, mm, dd, h, m, s)

    call ESMF_TimeSet(stopTime, yy=yy, mm=mm, dd=dd, h=h, m=m, s=s, &
         calkindflag=ESMF_CALKIND_GREGORIAN, rc=rc)
    VERIFY_NUOPC_(rc)

    call ESMF_ConfigGetAttribute(config, dt, label = "interpolation_dt:", rc = rc)
    VERIFY_NUOPC_(rc)

    ! set the driver clock
    call ESMF_TimeIntervalSet(timeStep, s=dt, rc=rc)
    VERIFY_NUOPC_(rc)

    internalClock = ESMF_ClockCreate(name="Driver Clock", timeStep = timeStep, &
         startTime=startTime, stopTime=stopTime, rc=rc)
    VERIFY_NUOPC_(rc)

    call ESMF_GridCompSet(driver, clock=internalClock, rc=rc)
    VERIFY_NUOPC_(rc)

  contains

    subroutine UnpackDateTime(DATETIME, YY, MM, DD, H, M, S)
      integer, intent(IN   ) :: DATETIME(:)
      integer, intent(  OUT) :: YY, MM, DD, H, M, S

      YY =     datetime(1)/10000
      MM = mod(datetime(1),10000)/100
      DD = mod(datetime(1),100)
      H  =     datetime(2)/10000
      M  = mod(datetime(2),10000)/100
      S  = mod(datetime(2),100)
      return
    end subroutine UnpackDateTime

  end subroutine set_clock
  
  subroutine SetRunSequence(driver, rc)
    type(ESMF_GridComp)  :: driver
    integer, intent(out) :: rc

    ! local variables
    type(ESMF_Time)               :: startTime
    type(ESMF_Time)               :: stopTime
    type(ESMF_TimeInterval)       :: timeStep
    type(ESMF_Config) :: config
    type(NUOPC_FreeFormat) :: run_sequence_ff

    rc = ESMF_SUCCESS

    print*, "Driver start Set Run Sequence"

    print*,"Driver read config"
    call ESMF_GridCompGet(driver, config=config, rc=rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver read run sequence"
    run_sequence_ff = NUOPC_FreeFormatCreate(config, label="run_sequence::", rc=rc)
    VERIFY_NUOPC_(rc)

    ! ingest FreeFormat run sequence
    print*,"Driver ingest run sequence"
    call NUOPC_DriverIngestRunSequence(driver, run_sequence_ff, rc=rc)
    VERIFY_NUOPC_(rc)

    print*,"Driver destroy run sequence"
    call NUOPC_FreeFormatDestroy(run_sequence_ff, rc=rc)
    VERIFY_NUOPC_(rc)

    print*, "Driver finish Set Run Sequence"
  end subroutine SetRunSequence

end module synthetic_driver
