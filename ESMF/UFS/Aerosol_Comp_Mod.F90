module Aerosol_Comp_Mod

  use ESMF
  use NUOPC
  use MAPL

  use Aerosol_Shared_mod

  implicit none

  character(len=*), dimension(*), parameter :: &
    fieldPairList = [ &
      "FROCEAN                         ", "ocean_fraction                  ", &
      "FRACI                           ", "ice_fraction                    ", &
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
      "SH                              ", "inst_sensi_heat_flx             ", &
      "T                               ", "inst_temp_levels                ", &
      "PLE                             ", "inst_pres_interface             ", &
      "U                               ", "inst_zonal_wind_levels          ", &
      "V                               ", "inst_merid_wind_levels          ", &
      "PFI_LSAN                        ", "inst_ice_nonconv_tendency_levels", &
      "PFL_LSAN                        ", "inst_liq_nonconv_tendency_levels", &
      "FCLD                            ", "inst_cloud_frac_levels          "  &
    ]

  character(len=*), dimension(*,2), parameter :: &
    fieldMap = reshape(fieldPairList, [size(fieldPairList)/2,2], order=[2,1])

  ! tracer map
  ! - prognostic section (advected)
  integer, parameter :: p_o3     = 7
  integer, parameter :: p_so2    = p_o3 + 1
  integer, parameter :: p_sulf   = p_o3 + 2
  integer, parameter :: p_dms    = p_o3 + 3
  integer, parameter :: p_msa    = p_o3 + 4
  integer, parameter :: p_bc_1   = p_o3 + 5
  integer, parameter :: p_bc_2   = p_o3 + 6
  integer, parameter :: p_oc_1   = p_o3 + 7
  integer, parameter :: p_oc_2   = p_o3 + 8
  integer, parameter :: p_dust_1 = p_o3 + 9
  integer, parameter :: p_dust_2 = p_o3 + 10
  integer, parameter :: p_dust_3 = p_o3 + 11
  integer, parameter :: p_dust_4 = p_o3 + 12
  integer, parameter :: p_seas_1 = p_o3 + 14
  integer, parameter :: p_seas_2 = p_o3 + 15
  integer, parameter :: p_seas_3 = p_o3 + 16
  integer, parameter :: p_seas_4 = p_o3 + 17
  integer, parameter :: p_nh3    = p_o3 + 19
  integer, parameter :: p_nh4a   = p_o3 + 20
  integer, parameter :: p_no3an1 = p_o3 + 21
  integer, parameter :: p_no3an2 = p_o3 + 22
  integer, parameter :: p_no3an3 = p_o3 + 23
  ! - diagnostic section (not advected)
  integer, parameter :: p_pm25   = p_o3 + 24
  integer, parameter :: p_pm10   = p_o3 + 25

  integer, parameter :: p_seas_d = p_seas_1 - p_so2 - 3

  private

  public :: AerosolDiagUpdate
  public :: AerosolStateUpdate

  public :: AerosolFieldDiagnostics
  public :: MAPLFieldDiagnostics


contains

  subroutine AerosolStateUpdate(model, cap, action, rc)
    type(ESMF_GridComp)            :: model
    type(MAPL_Cap)                 :: cap
    character(len=*),  intent(in)  :: action
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc, stat
    integer :: imap, item, itemCount
    integer :: rank
    integer :: fieldMapSize
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

    select case (trim(action))
      case ("import")

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
          fieldMapSize = size(fieldMap, 1)

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

        ! -- import tracer mixing ratios
        call AerosolTracerUpdate(state, cap, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out

      case ("export")

        call ESMF_GridCompGet(model, exportState=state, clock=clock, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out

        ! -- import tracer mixing ratios
        call AerosolTracerUpdate(state, cap, rc=localrc)
        if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__,  &
          file=__FILE__,  &
          rcToReturn=rc)) return  ! bail out

    end select

  end subroutine AerosolStateUpdate

  subroutine AerosolTracerUpdate(state, cap, rc)
    type(ESMF_State)               :: state
    type(MAPL_Cap)                 :: cap
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc, stat
    integer :: item, itemCount
    integer :: rank
    integer :: i, idiagn, j, k, kk, n, ni, nj, nk, offset, v
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer :: fp3d
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: fp4d
    real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: qu
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: q, qd
    real(ESMF_KIND_R8) :: iscale, escale
    type(ESMF_Field) :: field
    type(ESMF_FieldStatus_flag)         :: fieldStatus
    type(ESMF_StateIntent_flag)         :: stateintent
    type(ESMF_StateItem_flag),  pointer :: itemTypeList(:)
    character(len=ESMF_MAXSTR), pointer :: itemNameList(:)

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    call ESMF_StateGet(state, stateintent=stateintent, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(q)
    call AerosolGetPtr(state, "inst_tracer_mass_frac", q, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (stateintent == ESMF_STATEINTENT_EXPORT) then
      nullify(qu)
      call AerosolGetPtr(state, "inst_tracer_up_surface_flx", qu, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      nullify(qd)
      call AerosolGetPtr(state, "inst_tracer_down_surface_flx", qd, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
    end if

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

          idiagn = 0
          offset = 0
          iscale = 1._ESMF_KIND_R8
          escale = 1._ESMF_KIND_R8
          select case (trim(itemNameList(item)))
            case ("DMS")
              offset = p_dms - 1
              iscale = ppm2one
              escale = one2ppm
            case ("MSA")
              offset = p_msa - 1
              iscale = ppm2one
              escale = one2ppm
            case ("SO2")
              offset = p_so2 - 1
              iscale = ppm2one
              escale = one2ppm
            case ("SO4")
              offset = p_sulf - 1
              iscale = ug2kg
              escale = kg2ug
            case ("DU")
              offset = p_dust_1 - 1
              iscale = ug2kg
              escale = kg2ug
!           case ("DUEM")
!             offset = -1
!             iscale = ug2kg
!             escale = kg2ug
            case ("SS")
              offset = p_seas_1 - 1
              iscale = ug2kg
              escale = kg2ug
!           case ("SSEM")
!             offset = p_dust_1 - p_seas_1
            case ("SSSD")
              offset = -p_seas_d
              idiagn = 1
            case ("SSDP")
              offset = -p_seas_d
              idiagn = 2
            case ("SSWT")
              offset = -p_seas_d
              idiagn = 3
            case ("SSSV")
              offset = -p_seas_d
              idiagn = 4
            case ("CAphobicCA.bc")
              offset = p_bc_1 - 1
              iscale = ug2kg
              escale = kg2ug
            case ("CAphilicCA.bc")
              offset = p_bc_2 - 1
              iscale = ug2kg
              escale = kg2ug
            case ("CAphobicCA.oc")
              offset = p_oc_1 - 1
              iscale = ug2kg
              escale = kg2ug
            case ("CAphilicCA.oc")
              offset = p_oc_2 - 1
              iscale = ug2kg
              escale = kg2ug
            case ("NH3")
              offset = p_nh3 - 1
              iscale = ug2kg
              escale = kg2ug
            case ("NH4a")
              offset = p_nh4a - 1
              iscale = ug2kg
              escale = kg2ug
            case ("NO3an1")
              offset = p_no3an1 - 1
              iscale = ug2kg
              escale = kg2ug
            case ("NO3an2")
              offset = p_no3an2 - 1
              iscale = ug2kg
              escale = kg2ug
            case ("NO3an3")
              offset = p_no3an3 - 1
              iscale = ug2kg
              escale = kg2ug
            case default
              nullify(fp3d, fp4d)
          end select

          if (stateintent == ESMF_STATEINTENT_IMPORT) then
            if (offset < 0) nullify(fp3d,fp4d)
            if (associated(fp3d)) then
              v = offset + 1
              do k = 1, nk
                kk = nk - k + 1
                do j = 1, nj
                  do i = 1, ni
                    fp3d(i,j,kk) = iscale * max(0._ESMF_KIND_R8, q(i,j,k,v))
                  end do
                end do
              end do
            else if (associated(fp4d)) then
              do n = 1, size(fp4d, 4)
                v = offset + n
                do k = 1, nk
                  kk = nk - k + 1
                  do j = 1, nj
                    do i = 1, ni
                      fp4d(i,j,kk,n) = iscale * max(0._ESMF_KIND_R8, q(i,j,k,v))
                    end do
                  end do
                end do
              end do
            end if
          else if (stateintent == ESMF_STATEINTENT_EXPORT) then
            if (associated(fp3d)) then
              if (offset < 0) then
                ! -- diagnostics
                kk = size(fp3d, dim=3)
                if (idiagn > 0) then
                  if (associated(qd)) then
                    do k = 1, kk
                      do j = 1, nj
                        do i = 1, ni
                          qd(i,j,k-offset,idiagn) = escale * fp3d(i,j,k)
                        end do
                      end do
                    end do
                  end if
                else
                  if (associated(qu)) then
                    do k = 1, kk
                      do j = 1, nj
                        do i = 1, ni
                          qu(i,j,k-offset) = escale * fp3d(i,j,k)
                        end do
                      end do
                    end do
                  end if
                end if
              else
                v = offset + 1
                do k = 1, nk
                  kk = nk - k + 1
                  do j = 1, nj
                    do i = 1, ni
                      q(i,j,k,v) = escale * max(0._ESMF_KIND_R4, fp3d(i,j,kk))
                    end do
                  end do
                end do
              end if
            else if (associated(fp4d)) then
              do n = 1, size(fp4d, 4)
                v = offset + n
                do k = 1, nk
                  kk = nk - k + 1
                  do j = 1, nj
                    do i = 1, ni
                      q(i,j,k,v) = escale * max(0._ESMF_KIND_R4, fp4d(i,j,kk,n))
                    end do
                  end do
                end do
              end do
            end if
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

  subroutine AerosolDiagUpdate(model, cap, rc)
    type(ESMF_GridComp)            :: model
    type(MAPL_Cap)                 :: cap
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, j, k, kk, ni, nj, nk
    real(ESMF_KIND_R4) :: pm25, pm10
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer :: rho
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: q

    type(ESMF_Field) :: field
    type(ESMF_State) :: state

    ! -- local parameters
    real(ESMF_KIND_R8), parameter :: w_du2  = log(1.250_ESMF_KIND_R8) / log(1.8_ESMF_KIND_R8)
    real(ESMF_KIND_R8), parameter :: w_du4  = log(1.667_ESMF_KIND_R8) / log(2.0_ESMF_KIND_R8)
    real(ESMF_KIND_R8), parameter :: w_ss3  = log(2.50_ESMF_KIND_R8) / log(3._ESMF_KIND_R8)
    real(ESMF_KIND_R8), parameter :: w_so4  = 132.14_ESMF_KIND_R8 / 96.06_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: w_no3  = 80.043_ESMF_KIND_R8 / 62.0_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: w25_no3an2 = 0.138_ESMF_KIND_R8 * w_no3
    real(ESMF_KIND_R8), parameter :: w10_no3an2 = 0.808_ESMF_KIND_R8 * w_no3
    real(ESMF_KIND_R8), parameter :: w10_no3an3 = 0.164_ESMF_KIND_R8 * w_no3

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- retrieve export state
    call ESMF_GridCompGet(model, exportState=state, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- retrieve tracers
    nullify(q)
    call AerosolGetPtr(state, "inst_tracer_mass_frac", q, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- retrieve air density
    call ESMF_StateGet(cap % cap_gc % import_state, "AIRDENS", field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(rho)
    call ESMF_FieldGet(field, farrayPtr=rho, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- compute pm
    ni = size(q,1)
    nj = size(q,2)
    nk = size(q,3)

    do k = 1, nk
      kk = nk - k + 1
      do j = 1, nj
        do i = 1, ni
          ! -- compute partial PM2.5
          pm25 = q(i,j,k,p_bc_1) + q(i,j,k,p_bc_2)     &
               + q(i,j,k,p_oc_1) + q(i,j,k,p_oc_2)     &
               + q(i,j,k,p_dust_1)                     &
               + q(i,j,k,p_seas_1) + q(i,j,k,p_seas_2) &
               + w_so4 * q(i,j,k,p_sulf)               &
               + w_no3 * q(i,j,k,p_no3an1)

          ! -- compute PM10
          pm10 = pm25 &
               + q(i,j,k,p_dust_2) + q(i,j,k,p_dust_3) + w_du4 * q(i,j,k,p_dust_4) &
               + w10_no3an2 * q(i,j,k,p_no3an2) + w10_no3an3 * q(i,j,k,p_no3an3)   &
               + q(i,j,k,p_seas_3) + q(i,j,k,p_seas_4)

          ! -- complete PM2.5
          pm25 = pm25 &
               + w_du2      * q(i,j,k,p_dust_2) &
               + w25_no3an2 * q(i,j,k,p_no3an2) &
               + w_ss3      * q(i,j,k,p_seas_3)

          ! -- convert from ug/kg to ug/m3
          q(i,j,k,p_pm25) = pm25 * rho(i,j,kk)
          q(i,j,k,p_pm10) = pm10 * rho(i,j,kk)
        end do
      end do
    end do

  end subroutine AerosolDiagUpdate

end module Aerosol_Comp_Mod
