module Aerosol_Tracer_Mod

  use ESMF
  use gFTL_StringIntegerMap
  use gFTL_StringStringMap
  use MAPL, only: am_I_root => MAPL_am_I_root

  implicit none

  interface AerosolTracerMap
    module procedure AerosolTracerMapCreate
  end interface
  
  interface AerosolTracerMapPrint
    module procedure TracerStringStringMapPrint
    module procedure TracerStringIntegerMapPrint
  end interface

  private
  public :: AerosolTracerMap
  public :: AerosolTracerMapPrint
  public :: AerosolTracerGetUnitConv

contains

  function AerosolTracerInfoMapCreate(tracerInfo, separator, rc) result(map)
   type(StringIntegerMap) :: map
   ! -- interface variables
   character(len=*),           intent(in)  :: tracerInfo
   character(len=1), optional, intent(in)  :: separator
   integer,          optional, intent(out) :: rc

   ! -- local variables
   integer :: ib, ie, ix, item
   character(len=1) :: sep

   ! -- local parameter
   character(len=*), parameter :: default_separator = ':'

   ! -- begin
   if (present(rc)) rc = ESMF_SUCCESS

   if (present(separator)) then
     sep = separator
   else
     sep = default_separator
   end if

   ! -- create tracer->index map from input tracer info
   ib = 1
   ie = 0
   ix = index(tracerInfo(ib:), sep)
   item = 0
   do while (ix > 0)
     ie = ib + ix - 2
     if (ie > ib) then
       item = item + 1
       call map%insert(tracerInfo(ib:ie), item)
     end if
     ib = ie + 2
     ix = index(tracerInfo(ib:), sep)
   end do

  end function AerosolTracerInfoMapCreate


  function AerosolTracerFieldMapCreate(fileName, rc) result(map)
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
        call ESMF_ConfigFindLabel(config, trim(fieldTables(table)), rc=localrc)
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
    
  end function AerosolTracerFieldMapCreate


  function AerosolTracerIndexMapCreate(infoMap, fieldMap, rc) result(map)
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

  end function AerosolTracerIndexMapCreate


  function AerosolTracerMapCreate(fileName, tracerInfo, rc) result(map)
    type(StringIntegerMap) :: map
    ! -- interface variables
    character(len=*),  intent(in)  :: fileName
    character(len=*),  intent(in)  :: tracerInfo
    integer, optional, intent(out) :: rc 
 
    ! -- local variables
    integer :: localrc
    type(StringIntegerMap) :: infoMap
    type(StringStringMap)  :: fieldMap

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- create imported tracer map
    infoMap = AerosolTracerInfoMapCreate(tracerInfo, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (infoMap%empty()) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_CANNOT_GET, &
        msg="Imported tracer list is unavailable", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- create map between internal fields and imported tracers from the
    ! -- MAPL Cap configuration file
    fieldMap = AerosolTracerFieldMapCreate(fileName, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (fieldMap%empty()) then
      call ESMF_LogSetError(rcToCheck=ESMF_RC_CANNOT_GET, &
        msg="Field not mapped to imported tracers in "//trim(fileName), &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)
      return  ! bail out
    end if

    ! -- map internal fields to tracer index in exported tracer array
    map = AerosolTracerIndexMapCreate(infoMap, fieldMap, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call AerosolTracerMapPrint(fieldMap, "GOCART2G Field Name    Export To")

    ! -- free up memory
    call infoMap%clear()
    call fieldMap%clear()

  end function AerosolTracerMapCreate

  real(ESMF_KIND_R8) function AerosolTracerGetUnitConv(name, option)

    ! Return conversion factors to/from kg kg-1 according to
    ! the value of the option argument:
    ! option = 0: (ppm, ug kg-1) to kg kg-1
    ! option = 1: kg kg-1 to (ppm, ug kg-1)

    ! -- arguments
    character(len=*), intent(in) :: name
    integer,          intent(in) :: option

    ! -- local parameters
    real(ESMF_KIND_R8), dimension(0:1), parameter :: &
      ppm = [ 1.e-06_ESMF_KIND_R8, 1.e+06_ESMF_KIND_R8 ], &
      ug  = [ 1.e-09_ESMF_KIND_R8, 1.e+09_ESMF_KIND_R8 ]
    
    ! -- begin
    AerosolTracerGetUnitConv = 1._ESMF_KIND_R8

    if (option < 0 .or. option > 1) return

    select case (trim(name))
      case ("DMS", "MSA", "SO2")
        AerosolTracerGetUnitConv = ppm(option)
      case default
        AerosolTracerGetUnitConv = ug(option)
    end select
    
  end function AerosolTracerGetUnitConv

  subroutine TracerStringStringMapPrint(map, header, unit)

    ! -- arguments
    type(StringStringMap), intent(in) :: map
    character(len=*),      intent(in) :: header
    integer, optional,     intent(in) :: unit

    ! -- local variables
    integer :: ounit
    type(StringStringMapIterator) :: iter

    ! -- begin
    if (am_I_Root()) then

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

  subroutine TracerStringIntegerMapPrint(map, header, unit)

    ! -- arguments
    type(StringIntegerMap), intent(in) :: map
    character(len=*),       intent(in) :: header
    integer, optional,      intent(in) :: unit

    ! -- local variables
    integer :: ounit
    type(StringIntegerMapIterator) :: iter

    ! -- begin
    if (am_I_Root()) then

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
