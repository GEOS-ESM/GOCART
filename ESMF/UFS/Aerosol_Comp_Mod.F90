module Aerosol_Comp_Mod

  use ESMF
  use NUOPC
  use MAPL

  use Aerosol_Internal_mod
  use Aerosol_Shared_mod
  use Aerosol_Tracer_mod, only: AerosolTracerGetUnitsConv, &
                                Aerosol_Tracer_T

  implicit none

  character(len=*), dimension(*), parameter :: &
    fieldPairList = [ &
      "FROCEAN                         ", "ocean_fraction                  ", &
      "FRACI                           ", "ice_fraction_in_atm             ", &
      "FRLAKE                          ", "lake_fraction                   ", &
      "FRSNOW                          ", "surface_snow_area_fraction      ", &
      "LWI                             ", "inst_land_sea_mask              ", &
      "WET1                            ", "inst_surface_soil_wetness       ", &
      "U10M                            ", "inst_zonal_wind_height10m       ", &
      "V10M                            ", "inst_merid_wind_height10m       ", &
      "USTAR                           ", "inst_friction_velocity          ", &
      "TS                              ", "inst_temp_height_surface        ", &
      "AREA                            ", "surface_cell_area               ", &
      "ZPBL                            ", "inst_pbl_height                 ", &
      "SH                              ", "inst_up_sensi_heat_flx          ", &
      "T                               ", "inst_temp_levels                ", &
      "PLE                             ", "inst_pres_interface             ", &
      "U                               ", "inst_zonal_wind_levels          ", &
      "V                               ", "inst_merid_wind_levels          ", &
      "PFI_LSAN                        ", "inst_ice_nonconv_tendency_levels", &
      "PFL_LSAN                        ", "inst_liq_nonconv_tendency_levels", &
      "FCLD                            ", "inst_cloud_frac_levels          "  &
    ]

  integer, parameter :: fieldMapSize = size(fieldPairList)/2

  character(len=*), dimension(fieldMapSize,2), parameter :: &
    fieldMap = reshape(fieldPairList, [fieldMapSize,2], order=[2,1])

  private

  public :: AerosolStateUpdate

  public :: AerosolFieldDiagnostics
  public :: MAPLFieldDiagnostics


contains

  subroutine AerosolImportsUpdate(model, cap, rc)
    type(ESMF_GridComp)            :: model
    type(MAPL_Cap)                 :: cap
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc, stat
    integer :: imap, item, itemCount
    integer :: rank
    integer :: i,j, k, kk, n, ni, nj, nk, nlev, offset, v
    integer, dimension(1) :: plb, pub, rlb, rub
    real(ESMF_KIND_R4) :: blkevap, blkesat
    real(ESMF_KIND_R8) :: dt, tmp
    real(ESMF_KIND_R4), dimension(:,:),     pointer :: fp2dr
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer :: fp3dr
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: fp4dr
    real(ESMF_KIND_R8), dimension(:,:),     pointer :: fp2dp
    real(ESMF_KIND_R8), dimension(:,:),     pointer :: rain, rainc, zorl
    real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: fp3dp, slc
    real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: phii, phil, prsi, prsl, temp
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: q
    type(ESMF_Clock) :: clock
    type(ESMF_Field) :: pfield, rfield
    type(ESMF_State) :: state
    type(ESMF_TimeInterval) :: timeStep
    type(ESMF_StateItem_flag),  pointer :: itemTypeList(:)
    character(len=ESMF_MAXSTR), pointer :: itemNameList(:)

    ! -- local parameters

    real(ESMF_KIND_R8), parameter :: rd = 2.8705E+2_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: rv = 4.6150e+2_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: fv = rv / rd - 1._ESMF_KIND_R8

    real(ESMF_KIND_R8), parameter :: mwwat = 18.0153_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: mwair = 28.9628_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: epswater = mwwat / mwair
    real(ESMF_KIND_R8), parameter :: stdtemp = 273.15_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: al = 610.94_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: bl = 17.625_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: cl = 243.04_ESMF_KIND_R8

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_GridCompGet(model, importState=state, clock=clock, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_ClockGet(clock, timeStep=timeStep, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_TimeIntervalGet(timeStep, s_r8=dt, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(cap % cap_gc % import_state, itemCount=itemCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (itemCount > 0) then

      allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Unable to allocate internal workspace", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call ESMF_StateGet(cap % cap_gc % import_state, itemNameList=itemNameList, &
        itemTypeList=itemTypeList, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      ! -- retrieve pointers to field data from import state
      call AerosolGetPtr(state, "inst_pres_interface", prsi, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_pres_levels",    prsl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_geop_interface",    phii, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_geop_levels",       phil, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_temp_levels",       temp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_rainfall_amount", rain, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_convective_rainfall_amount", rainc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_surface_roughness", zorl, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_soil_moisture_content", slc, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call AerosolGetPtr(state, "inst_tracer_mass_frac",  q, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      ni = size(q, 1)
      nj = size(q, 2)
      nk = size(q, 3)

      do item = 1, itemCount
        if (itemTypeList(item) == ESMF_STATEITEM_FIELD) then

          call ESMF_StateGet(cap % cap_gc % import_state, trim(itemNameList(item)), &
            rfield, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out

          call ESMF_FieldGet(rfield, rank=rank, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out

          nullify(fp2dr, fp3dr, fp4dr)
          select case (rank)
            case (2)
              call ESMF_FieldGet(rfield, farrayPtr=fp2dr, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) return  ! bail out
            case (3)
              call ESMF_FieldGet(rfield, ungriddedLBound=rlb, ungriddedUBound=rub, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) return  ! bail out
              call ESMF_FieldGet(rfield, farrayPtr=fp3dr, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) return  ! bail out
            case (4)
              call ESMF_FieldGet(rfield, farrayPtr=fp4dr, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) return  ! bail out
          end select

          imap = 0
          do while (imap < fieldMapSize)
            imap = imap + 1
            if (trim(fieldMap(imap,1)) == trim(itemNameList(item))) then

              call ESMF_StateGet(state, trim(fieldMap(imap,2)), pfield, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) return  ! bail out

              select case (rank)
                case (2)
                  nullify(fp2dp)
                  call ESMF_FieldGet(pfield, farrayPtr=fp2dp, rc=localrc)
                  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                    line=__LINE__,  &
                    file=__FILE__,  &
                    rcToReturn=rc)) return  ! bail out
                  fp2dr = fp2dp
                case (3)
                  nullify(fp3dp)
                  call ESMF_FieldGet(pfield, ungriddedLBound=plb, ungriddedUBound=pub, rc=localrc)
                  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                    line=__LINE__,  &
                    file=__FILE__,  &
                    rcToReturn=rc)) return  ! bail out
                  call ESMF_FieldGet(pfield, farrayPtr=fp3dp, rc=localrc)
                  if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                    line=__LINE__,  &
                    file=__FILE__,  &
                    rcToReturn=rc)) return  ! bail out

                  ! -- map provider field levels to receiver field levels in reverse order
                  ! -- NOTE: if provider field has fewer vertical levels than the receiver field,
                  ! -- the remaining receiver field levels are filled by replicating values from
                  ! -- the closest available level in the provider field.
                  kk = plb(1)
                  do k = rub(1), rlb(1), -1
                    do j = 1, nj
                      do i = 1, ni
                        fp3dr(i,j,k) = fp3dp(i,j,kk)
                      end do
                    end do
                    kk = min(pub(1), kk + 1)
                  end do
                case default
                ! -- nothing to do
              end select
              imap = fieldMapSize + 1
            end if
          end do

          if (imap > fieldMapSize) cycle

          select case (trim(itemNameList(item)))
            case ("AIRDENS")
              do k = 1, nk
                kk = nk - k + 1
                do j = 1, nj
                  do i = 1, ni
                    fp3dr(i,j,kk) = prsl(i,j,k) / (rd * temp(i,j,k) * ( 1._ESMF_KIND_R8 + fv * q(i,j,k,1) ) )
                  end do
                end do
              end do
            case ("DELP")
              do k = 1, nk
                kk = nk - k + 1
                do j = 1, nj
                  do i = 1, ni
                    fp3dr(i,j,kk) = prsi(i,j,k) - prsi(i,j,k+1)
                  end do
                end do
              end do
            case ("DZ")
              fp2dr = onebg * phil(:,:,1)
            case ("CN_PRCP")
              fp2dr = rainc / dt
            case ("NCN_PRCP")
              fp2dr = max(0._ESMF_KIND_R8, rain - rainc) / dt
            case ("RH2")
              do k = 1, nk
                kk = nk - k + 1
                do j = 1, nj
                  do i = 1, ni
                    tmp = temp(i,j,k) - stdtemp
                    blkesat = al * exp( bl * ( tmp / ( tmp + cl ) ) )
                    blkevap = prsl(i,j,k) * q(i,j,k,1) / ( epswater + q(i,j,k,1) )
                    fp3dr(i,j,kk) = max( 0.005, min( 0.99, blkevap / blkesat ) )
                  end do
                end do
              end do
            case ("SLC")
              fp2dr = slc(:,:,1)
            case ("Z0H")
              fp2dr = 0.01_ESMF_KIND_R4 * zorl
            case ("ZLE")
              do k = 1, nk + 1
                kk = nk - k + 1
                do j = 1, nj
                  do i = 1, ni
                    fp3dr(i,j,kk) = onebg * phii(i,j,k)
                  end do
                end do
              end do
          end select
        end if

      end do
      deallocate(itemNameList, itemTypeList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Unable to deallocate internal workspace", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

    end if

  end subroutine AerosolImportsUpdate

  subroutine AerosolStateUpdate(model, cap, action, rc)
    type(ESMF_GridComp)            :: model
    type(MAPL_Cap)                 :: cap
    character(len=*),  intent(in)  :: action
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- update imported fields if required
    if (trim(action) == "import") then
      call AerosolImportsUpdate(model, cap, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
    end if

    ! -- update tracers
    call AerosolTracerUpdate(model, action, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolStateUpdate

  subroutine AerosolTracerUpdate(model, action, rc)
    type(ESMF_GridComp)            :: model
    character(len=*),  intent(in)  :: action
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc, stat
    integer :: item, itemCount
    integer :: rank
    integer :: i, j, k, kk, n, ni, nj, nk, v
    integer :: intentOption
    integer, pointer :: tracerId
    character(len=:), pointer :: pUnits
    character(len=ESMF_MAXSTR) :: fieldUnits, tracerUnits
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer :: fp3d
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: fp4d
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: q
    real(ESMF_KIND_R8) :: scalefac
    type(ESMF_Field) :: field
    type(ESMF_FieldStatus_flag)         :: fieldStatus
    type(ESMF_State) :: state
    type(ESMF_StateItem_flag),  pointer :: itemTypeList(:)
    character(len=ESMF_MAXSTR), pointer :: itemNameList(:)
    type(MAPL_Cap),             pointer :: cap
    type(Aerosol_Tracer_T),     pointer :: trp
    type(Aerosol_InternalState_T)       :: is

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- retrieve model import/export states
    select case (trim(action))
      case ("import")
        intentOption = 0
        call ESMF_GridCompGet(model, importState=state, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
      case ("export")
        intentOption = 1
        call ESMF_GridCompGet(model, exportState=state, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out
      case default
        ! -- action undefined, bailout
        return
    end select

    ! -- retrieve model internal state
    call ESMF_GridCompGetInternalState(model, is, localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(cap, trp)
    if (.not.associated(is % wrap)) return
    cap  => is % wrap % maplCap
    trp  => is % wrap % tracers
    if (.not.associated(cap)) return
    if (.not.associated(trp)) return

    nullify(q)
    call AerosolGetPtr(state, "inst_tracer_mass_frac", q, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(cap % cap_gc % export_state, itemCount=itemCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (itemCount > 0) then

      allocate(itemNameList(itemCount), itemTypeList(itemCount), stat=stat)
      if (ESMF_LogFoundAllocError(statusToCheck=stat, &
        msg="Unable to allocate internal workspace", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call ESMF_StateGet(cap % cap_gc % export_state, itemNameList=itemNameList, &
        itemTypeList=itemTypeList, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      ni = size(q, 1)
      nj = size(q, 2)
      nk = size(q, 3)

      do item = 1, itemCount
        if (itemTypeList(item) == ESMF_STATEITEM_FIELD) then

          call ESMF_StateGet(cap % cap_gc % export_state, trim(itemNameList(item)), &
            field, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          rank = 0
          call ESMF_FieldGet(field, status=fieldStatus, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out

          if (fieldStatus == ESMF_FIELDSTATUS_COMPLETE) then
            call ESMF_FieldGet(field, rank=rank, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
          end if

          nullify(fp3d, fp4d)
          select case (rank)
            case (3)
              call ESMF_FieldGet(field, farrayPtr=fp3d, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) return  ! bail out
            case (4)
              call ESMF_FieldGet(field, farrayPtr=fp4d, rc=localrc)
              if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                line=__LINE__,  &
                file=__FILE__,  &
                rcToReturn=rc)) return  ! bail out
            case default
              cycle
          end select

          ! -- get tracer location in coupled field
          tracerId => trp % indexMap % at(trim(itemNameList(item)))

          if (associated(tracerId)) then

            ! -- get units for automatic conversion
            ! -- (a) tracer units
            tracerUnits = trp % units(tracerId)
            ! -- (b) MAPL field units
            call ESMF_AttributeGet(field, name='UNITS', value=fieldUnits, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out

            select case (intentOption)
              case (0)
                ! -- import
                scalefac = AerosolTracerGetUnitsConv(tracerUnits, fieldUnits)
                if (associated(fp3d)) then
                  v = tracerId
                  do k = 1, nk
                    kk = nk - k + 1
                    do j = 1, nj
                      do i = 1, ni
                        fp3d(i,j,kk) = scalefac * max(0._ESMF_KIND_R8, q(i,j,k,v))
                      end do
                    end do
                  end do
                else if (associated(fp4d)) then
                  v = tracerId
                  do n = 1, size(fp4d, 4)
                    do k = 1, nk
                      kk = nk - k + 1
                      do j = 1, nj
                        do i = 1, ni
                          fp4d(i,j,kk,n) = scalefac * max(0._ESMF_KIND_R8, q(i,j,k,v))
                        end do
                      end do
                    end do
                    v = v + 1
                  end do
                end if
              case (1)
                ! -- export
                scalefac = AerosolTracerGetUnitsConv(fieldUnits, tracerUnits)
                if (associated(fp3d)) then
                  v = tracerId
                  do k = 1, nk
                    kk = nk - k + 1
                    do j = 1, nj
                      do i = 1, ni
                        q(i,j,k,v) = scalefac * max(0._ESMF_KIND_R4, fp3d(i,j,kk))
                      end do
                    end do
                  end do
                else if (associated(fp4d)) then
                  v = tracerId
                  do n = 1, size(fp4d, 4)
                    do k = 1, nk
                      kk = nk - k + 1
                      do j = 1, nj
                        do i = 1, ni
                          q(i,j,k,v) = scalefac * max(0._ESMF_KIND_R4, fp4d(i,j,kk,n))
                        end do
                      end do
                    end do
                    v = v + 1
                  end do
                end if
            end select
          end if

        end if
      end do

      deallocate(itemNameList, itemTypeList, stat=stat)
      if (ESMF_LogFoundDeallocError(statusToCheck=stat, &
        msg="Unable to deallocate internal workspace", &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

    end if

  end subroutine AerosolTracerUpdate

end module Aerosol_Comp_Mod
