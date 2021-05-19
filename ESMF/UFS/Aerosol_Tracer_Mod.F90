module Aerosol_Tracer_Mod

  use ESMF
  use gFTL_StringIntegerMap
  use gFTL_StringStringMap
  use Aerosol_Logger_mod, only: isAerosolLoggerOn

  implicit none

  type Aerosol_Tracer_T
    type(StringIntegerMap) :: indexMap
    character(len=ESMF_MAXSTR), allocatable :: names(:)
    character(len=ESMF_MAXSTR), allocatable :: units(:)
  end type

  interface AerosolTracer
    module procedure AerosolTracerCreate
  end interface
  
  private

  ! -- types
  public :: Aerosol_Tracer_T

  ! -- methods
  public :: AerosolTracer
  public :: AerosolTracerPrint
  public :: AerosolTracerGetUnitsConv

contains

  ! ----- Type constructor -----

  function AerosolTracerCreate(fileName, info, rc) result(tracers)
    type(Aerosol_Tracer_T) :: tracers
    ! -- interface variables
    character(len=*),           intent(in)  :: fileName
    type(ESMF_Info),            intent(in)  :: info
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: item, stat
    type(StringIntegerMap) :: infoMap
    type(StringStringMap)  :: configMap

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! - import tracer names
    call TracerInfoGet(info, 'tracerNames', tracers % names, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (.not.allocated(tracers % names)) then
      call ESMF_LogWrite("Unable to retrieve imported tracer list", &
        ESMF_LOGMSG_WARNING, &
        line=__LINE__, &
        file=__FILE__, &
        rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
      return
    end if

    ! - import tracer units if available
    call TracerInfoGet(info, 'tracerUnits', tracers % units, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (.not.allocated(tracers % units)) then
      allocate(tracers % units(infoMap % size()), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Unable to allocate internal workspace", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
      tracers % units = 'n/a'
    end if

    ! -- create imported tracer map
    do item = 1, size(tracers % names)
      call infoMap % insert(trim(tracers % names(item)), item)
    end do

    ! -- create map between internal fields and imported tracers
    ! -- according to the MAPL Cap configuration file
    configMap = TracerConfigMapCreate(fileName, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (configMap%empty()) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_CANNOT_GET, &
        msg="Field not mapped to imported tracers in "//trim(fileName), &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- map internal fields to tracer index in exported tracer array
    tracers % indexMap = TracerIndexMapCreate(infoMap, configMap, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- free up memory
    call infoMap%clear()
    call configMap%clear()

  end function AerosolTracerCreate

  ! ----- Type print method -----

  subroutine AerosolTracerPrint(tracers, header, unit)
    ! -- arguments
    type(Aerosol_Tracer_T), intent(in) :: tracers
    character(len=*),       intent(in) :: header
    integer, optional,      intent(in) :: unit

    ! -- local variables
    integer               :: logUnit
    integer               :: nspc
    integer,      pointer :: idx
    character(:), pointer :: fld
    character(:), allocatable :: name, units
    type(StringIntegerMapIterator) :: iter

    ! -- local parameters
    character(len=*), parameter :: unknown = 'unknown'

    ! -- begin
    if (isAerosolLoggerOn()) then

      logUnit = 6
      if (present(unit)) logUnit = unit

      name  = unknown
      units = unknown

      nspc = max(0, (37 - len_trim(header))/2)
      write(logUnit, '(38("-")/a/38("-"))') repeat(" ",nspc)//trim(header)
      iter = tracers%indexMap%begin()
      do while (iter /= tracers%indexMap%end())
        fld => iter%key()
        idx => iter%value()
        if (allocated(tracers%names)) name  = trim(tracers % names(idx))
        if (allocated(tracers%units)) units = trim(tracers % units(idx))
        write(logUnit, '(a15," -> ",a10,2x,"[",a,"]")') fld, name, units
        call iter%next()
      end do
      write(logUnit, '(38("-"))')

    end if

  end subroutine AerosolTracerPrint

  ! ----- Tracer metadata methods -----

  subroutine TracerInfoGet(info, key, values, rc)
    ! -- interface variables
    type(ESMF_Info),               intent(in)  :: info
    character(len=*),              intent(in)  :: key
    character(len=*), allocatable, intent(out) :: values(:)
    integer,          optional,    intent(out) :: rc

    ! -- local variables
    integer :: localrc
    logical :: isKeyFound

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- check if metadata key is present in ESMF_Info object
    isKeyFound = ESMF_InfoIsPresent(info, trim(key), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    if (isKeyFound) then
      isKeyFound = ESMF_InfoIsSet(info, trim(key), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

    if (isKeyFound) then
      call ESMF_InfoGetAlloc(info, trim(key), values, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

  end subroutine TracerInfoGet

  ! ----- Maps create -----

  function TracerInfoMapCreate(info, key, rc) result(map)
    type(StringIntegerMap) :: map
    ! -- interface variables
    type(ESMF_Info),            intent(in)  :: info
    character(len=*), optional, intent(in)  :: key
    integer,          optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: item
    logical :: isKeyFound
    character(len=ESMF_MAXSTR), allocatable :: values(:)

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- check if metadata key is present in ESMF_Info object
    isKeyFound = ESMF_InfoIsPresent(info, trim(key), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__, &
      file=__FILE__, &
      rcToReturn=rc)) return

    if (isKeyFound) then
      isKeyFound = ESMF_InfoIsSet(info, trim(key), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

    if (isKeyFound) then
      call ESMF_InfoGetAlloc(info, trim(key), values, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
      do item = 1, size(values)
        call map%insert(trim(values(item)), item)
      end do
    end if

  end function TracerInfoMapCreate

  function TracerConfigMapCreate(fileName, rc) result(map)
    type(StringStringMap) :: map
    ! -- argument variables
    character(len=*),  intent(in)  :: fileName
    integer, optional, intent(out) :: rc 

    ! -- local variables
    integer :: localrc
    integer :: idx, item, table, rowCount, columnCount
    logical :: eolFlag
    character(len=ESMF_MAXSTR) :: key, value
    type(ESMF_Config) :: config

    ! -- local parameters
    character(len=*), dimension(*), parameter :: fieldTables = &
      [ 'CAP_EXPORTS    ', 'CAP_DIAGNOSTICS' ]

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    config = ESMF_ConfigCreate(rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_ConfigLoadFile(config, trim(fileName), rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    do table = 1, size(fieldTables)
      call ESMF_ConfigGetDim(config, rowCount, columnCount, &
        label=trim(fieldTables(table)), rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      if (columnCount > 1) then
        call ESMF_ConfigFindLabel(config, trim(fieldTables(table))//':', rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
        do item = 1, rowCount
          call ESMF_ConfigNextLine(config, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          call ESMF_ConfigGetAttribute(config, value, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          idx = index(value, ',')
          if (idx > 0) then
            key = value(1:idx-1)
          else
            key = value
          end if
          call ESMF_ConfigGetAttribute(config, value, eolFlag=eolFlag, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          if (.not.eolFlag) call map % insert(key, value)
        end do
      end if
    end do

    call ESMF_ConfigDestroy(config, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out
    
  end function TracerConfigMapCreate

  function TracerIndexMapCreate(infoMap, fieldMap, rc) result(map)
    type(StringIntegerMap) :: map
    ! -- interface variables
    type(StringIntegerMap), intent(in)  :: infoMap
    type(StringStringMap),  intent(in)  :: fieldMap
    integer, optional,      intent(out) :: rc
 
    ! -- local variables
    integer                       :: lky
    integer, pointer              :: pos
    character(len=ESMF_MAXSTR)    :: key
    type(StringStringMapIterator) :: iter
 
    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS
 
    iter = fieldMap%begin()
    do while (iter /= fieldMap%end())
      key =  fieldMap%at(trim(iter%key()))
      lky = len_trim(key)
      if (key(lky:lky) == "*") key(lky:lky) = "1"
      pos => infoMap%at(trim(key))
      if (associated(pos)) call map%insert(trim(iter%key()), pos)
      call iter%next()
    end do

  end function TracerIndexMapCreate

  ! ----- Units conversion -----

  real(ESMF_KIND_R8) function AerosolTracerGetUnitsConv(fromUnits, toUnits)

    ! Return conversion factors from 'fromUnits' to 'toUnits'.

    ! -- arguments
    character(len=*), target, intent(in) :: fromUnits
    character(len=*), target, intent(in) :: toUnits

    ! -- local variables
    integer :: n
    integer :: ij(2)
    character(len=:), pointer :: pUnits

    ! -- local parameters
    real(ESMF_KIND_R8), dimension(3,3), parameter :: convTable = &
      reshape( [ &
      ! ---------------------------------------------------------------------------
      !        kg/kg        |       ug/kg        |       ppm              | From/To
      ! ---------------------------------------------------------------------------
        1.e+00_ESMF_KIND_R8, 1.e+09_ESMF_KIND_R8, 1.e+06_ESMF_KIND_R8, &  ! kg/kg
        1.e-09_ESMF_KIND_R8, 1.e+00_ESMF_KIND_R8, 1.e-03_ESMF_KIND_R8, &  ! ug/kg
        1.e-06_ESMF_KIND_R8, 1.e+03_ESMF_KIND_R8, 1.e+00_ESMF_KIND_R8  &  ! ppm
      ! ---------------------------------------------------------------------------
      ], [3,3], order=[2,1])
    
    ! -- begin
    AerosolTracerGetUnitsConv = 1._ESMF_KIND_R8

    ij = 0
    pUnits => fromUnits
    do n = 1, 2
      select case (trim(pUnits))
        case ("kg kg-1", "kg/kg")
          ij(n) = 1
        case ("ug kg-1", "ug/kg")
          ij(n) = 2
        case ("ppm")
          ij(n) = 3
        case default
          ij(1) = -1
          exit
      end select
      pUnits => toUnits
    end do

    if (ij(1) > 0) AerosolTracerGetUnitsConv = convTable(ij(1),ij(2))

  end function AerosolTracerGetUnitsConv

  ! ----- Map print methods -----

  subroutine TracerStringStringMapPrint(map, header, unit)

    ! -- arguments
    type(StringStringMap), intent(in) :: map
    character(len=*),      intent(in) :: header
    integer, optional,     intent(in) :: unit

    ! -- local variables
    integer :: ounit
    type(StringStringMapIterator) :: iter

    ! -- begin
    if (isAerosolLoggerOn()) then

      ounit = 6
      if (present(unit)) ounit = unit

      write(ounit, '(a/32("-"))') trim(header)
      iter = map%begin()
      do while (iter /= map%end())
        write(ounit, '(a19," -> ",a)') trim(iter%key()), trim(map%at(trim(iter%key())))
        call iter%next()
      end do
      write(ounit, '(32("-"))')

    end if
    
  end subroutine TracerStringStringMapPrint

  subroutine TracerStringStringMapPairPrint(map1, map2, header, unit)

    ! -- arguments
    type(StringStringMap), intent(in) :: map1
    type(StringStringMap), intent(in) :: map2
    character(len=*),      intent(in) :: header
    integer, optional,     intent(in) :: unit

    ! -- local variables
    integer :: ounit
    character(len=:), pointer :: pval1, pval2
    type(StringStringMapIterator) :: iter

    ! -- begin
    if (isAerosolLoggerOn()) then

      ounit = 6
      if (present(unit)) ounit = unit

      write(ounit, '(a/32("-"))') trim(header)
      iter = map1%begin()
      do while (iter /= map1%end())
        nullify(pval1,pval2)
        pval1 => map1%at(trim(iter%key()))
        pval2 => map2%at(trim(pval1))
        if (associated(pval2)) then
          write(ounit, '(a19," -> ",a10,2x,"[",a,"]")') trim(iter%key()), pval1, pval2
        else
          write(ounit, '(a19," -> ",a10,2x,"[n/a]")') trim(iter%key()), pval1
        end if
        call iter%next()
      end do
      write(ounit, '(32("-"))')

    end if

  end subroutine TracerStringStringMapPairPrint

  subroutine TracerStringIntegerMapPrint(map, header, unit)

    ! -- arguments
    type(StringIntegerMap), intent(in) :: map
    character(len=*),       intent(in) :: header
    integer, optional,      intent(in) :: unit

    ! -- local variables
    integer :: ounit
    type(StringIntegerMapIterator) :: iter

    ! -- begin
    if (isAerosolLoggerOn()) then

      ounit = 6
      if (present(unit)) ounit = unit

      write(ounit, '(a/32("-"))') trim(header)
      iter = map%begin()
      do while (iter /= map%end())
        write(ounit, '(a19," -> ",i0)') trim(iter%key()), map%at(trim(iter%key()))
        call iter%next()
      end do
      write(ounit, '(32("-"))')

    end if
    
  end subroutine TracerStringIntegerMapPrint

end module Aerosol_Tracer_Mod
