module Aerosol_Comp_Mod

  use ESMF
  use NUOPC
  use MAPL

  implicit none

  integer, parameter :: fieldMapSize = 18
  character(len=*), dimension(fieldMapSize, 2), parameter :: &
    fieldMap = reshape((/ &
      "FROCEAN                         ", "ocean_fraction                  ", &
      "FRACI                           ", "ice_fraction                    ", &
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
      "PFI_LSAN                        ", "inst_liq_nonconv_tendency_levels", &
      "PFL_LSAN                        ", "inst_liq_nonconv_tendency_levels", &
      "FCLD                            ", "inst_cloud_frac_levels          "  &
      /), (/fieldMapSize, 2/), order=(/2,1/))

  ! gravity (m/s2)
  real(ESMF_KIND_R8), parameter :: con_g = 9.80665e+0_ESMF_KIND_R8
  ! inverse gravity (m/s2)
  real(ESMF_KIND_R8), parameter :: onebg = 1._ESMF_KIND_R8 / con_g
  ! unit conversion
  real(ESMF_KIND_R8), parameter :: kg2ug   = 1.e+09_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: ug2kg   = 1.e-09_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: one2ppm = 1.e+06_ESMF_KIND_R8
  real(ESMF_KIND_R8), parameter :: ppm2one = 1.e-06_ESMF_KIND_R8


  interface AerosolGetPtr
    module procedure AerosolGetPtr2D
    module procedure AerosolGetPtr3D
    module procedure AerosolGetPtr4D
  end interface

  private

  public :: AerosolStateUpdate
  public :: AerosolTracerReroute

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
    integer :: i,j, k, kk, n, ni, nj, nk, nlev, offset, v
    integer, dimension(1) :: plb, pub, rlb, rub
    real(ESMF_KIND_R4) :: blkevap, blkesat
    real(ESMF_KIND_R8) :: dt, tmp
    real(ESMF_KIND_R4), dimension(:,:),     pointer :: fp2dr
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer :: fp3dr
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: fp4dr
    real(ESMF_KIND_R8), dimension(:,:),     pointer :: fp2dp
    real(ESMF_KIND_R8), dimension(:,:),     pointer :: rain, rainc, zorl
    real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: fp3dp
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
                        fp3dr(i,j,kk) = prsl(i,j,k) / (rv * temp(i,j,k) * ( 1._ESMF_KIND_R8 + fv * q(i,j,k,1) ) )
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

    ! -- local parameters
    integer, parameter :: p_o3     = 7
    integer, parameter :: p_so2    = p_o3 + 1
    integer, parameter :: p_sulf   = p_o3 + 2
    integer, parameter :: p_dms    = p_o3 + 3
    integer, parameter :: p_msa    = p_o3 + 4
    integer, parameter :: p_bc_1   = p_o3 + 6
    integer, parameter :: p_bc_2   = p_o3 + 7
    integer, parameter :: p_oc_1   = p_o3 + 8
    integer, parameter :: p_oc_2   = p_o3 + 9
    integer, parameter :: p_dust_1 = p_o3 + 10
    integer, parameter :: p_seas_1 = p_o3 + 15
    integer, parameter :: p_nh3    = p_o3 + 20
    integer, parameter :: p_nh4a   = p_o3 + 21
    integer, parameter :: p_no3an1 = p_o3 + 22
    integer, parameter :: p_no3an2 = p_o3 + 23
    integer, parameter :: p_no3an3 = p_o3 + 24
    integer, parameter :: p_seas_d = p_seas_1 - p_so2 - 3

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

  subroutine AerosolGetPtr2D(state, fieldName, farrayPtr, localDe, rc)
    type(ESMF_State)                :: state
    character(len=*),   intent(in)  :: fieldName
    real(ESMF_KIND_R8), pointer     :: farrayPtr(:,:)
    integer, optional,  intent(in)  :: localDe
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer          :: localrc
    type(ESMF_Field) :: field

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    nullify(farrayPtr)

    call ESMF_StateGet(state, fieldName, field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=farrayPtr, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolGetPtr2D

  subroutine AerosolGetPtr3D(state, fieldName, farrayPtr, localDe, rc)
    type(ESMF_State)                :: state
    character(len=*),   intent(in)  :: fieldName
    real(ESMF_KIND_R8), pointer     :: farrayPtr(:,:,:)
    integer, optional,  intent(in)  :: localDe
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer          :: localrc
    type(ESMF_Field) :: field

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    nullify(farrayPtr)

    call ESMF_StateGet(state, fieldName, field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=farrayPtr, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolGetPtr3D

  subroutine AerosolGetPtr4D(state, fieldName, farrayPtr, localDe, rc)
    type(ESMF_State)                :: state
    character(len=*),   intent(in)  :: fieldName
    real(ESMF_KIND_R8), pointer     :: farrayPtr(:,:,:,:)
    integer, optional,  intent(in)  :: localDe
    integer, optional,  intent(out) :: rc

    ! -- local variables
    integer          :: localrc
    type(ESMF_Field) :: field

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    nullify(farrayPtr)

    call ESMF_StateGet(state, fieldName, field, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_FieldGet(field, localDe=localDe, farrayPtr=farrayPtr, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

  end subroutine AerosolGetPtr4D

  ! -- Additional tools --

  subroutine AerosolTracerReroute(model, rc)
    type(ESMF_GridComp)            :: model
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: localDe, localDeCount
    type(ESMF_Field) :: iField, eField
    type(ESMF_State) :: importState, exportState
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: fpImp, fpExp

    logical, save :: first = .true.

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- query the Component for its VM and importState
    call ESMF_GridCompGet(model, importState=importState, &
      exportState=exportState, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- retrieve
    call ESMF_StateGet(importState, "inst_tracer_mass_frac", iField, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(exportState, "inst_tracer_mass_frac", eField, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    if (first) then
      call ESMF_FieldFill(eField, dataFillScheme="one", const1=1.e-05_ESMF_KIND_R8, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out
      first = .false.
      return
    end if

    call ESMF_FieldGet(eField, localDeCount=localDeCount, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    do localDe = 0, localDeCount - 1
      nullify(fpExp, fpImp)
      call ESMF_FieldGet(eField, localDe=localDe, farrayPtr=fpImp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return  ! bail out
      call ESMF_FieldGet(iField, localDe=localDe, farrayPtr=fpExp, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return  ! bail out
      fpExp = fpImp
    end do

  end subroutine AerosolTracerReroute

  subroutine AerosolFieldDiagnostics(model, rc)
    type(ESMF_GridComp)            :: model
    integer, optional, intent(out) :: rc

    ! -- local variables
    logical :: isConnected
    integer :: localrc
    integer :: item
    integer :: localDe, localDeCount, rank
    integer :: fieldCount, maxLength
    character(len=ESMF_MAXSTR)          :: msgString
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR), pointer :: connectedList(:)
    character(len=ESMF_MAXSTR), pointer :: fieldNames(:)
    real(ESMF_KIND_R8)                              :: minValue, maxValue
    real(ESMF_KIND_R8), dimension(:),       pointer :: localMin, localMax, globalMin, globalMax
    real(ESMF_KIND_R8), dimension(:),       pointer :: fp1d
    real(ESMF_KIND_R8), dimension(:,:),     pointer :: fp2d
    real(ESMF_KIND_R8), dimension(:,:,:),   pointer :: fp3d
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: fp4d
    type(ESMF_Field), pointer :: fieldList(:)
    type(ESMF_State)  :: importState
    type(ESMF_VM)     :: vm

    ! -- local parameters
    character(len=*), parameter :: rName = "diagnostics"

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- get component information
    call NUOPC_CompGet(model, name=name, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- query the Component for its VM and importState
    call ESMF_GridCompGet(model, vm=vm, importState=importState, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(fieldNames, connectedList, fieldList)
    call NUOPC_GetStateMemberLists(importState, StandardNameList=fieldNames, &
      ConnectedList=connectedList, fieldList=fieldList, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    isConnected = .false.
    if (associated(fieldNames)) isConnected = any(connectedList == "true")

    nullify(localMin, localMax, globalMin, globalMax)

    ! -- check values of imported fields, if requested
    if (isConnected) then

      fieldCount = size(fieldList)

      allocate(localMin(fieldCount), localMax(fieldCount), &
        globalMin(fieldCount), globalMax(fieldCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

      localMin = huge(0._ESMF_KIND_R8)
      localMax = -localMin
      globalMin = 0._ESMF_KIND_R8
      globalMax = 0._ESMF_KIND_R8

      ! -- find longest field name for formatting purpose
      maxLength = 0
      do item = 1, fieldCount
        maxLength = max(maxLength, len_trim(fieldNames(item)))
      end do

      do item = 1, fieldCount
        if (connectedList(item) == "true") then
          ! --- get field data
          call ESMF_FieldGet(fieldList(item), rank=rank, &
            localDeCount=localDeCount, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out

          do localDe = 0, localDeCount - 1
            minValue = localMin(item)
            maxValue = localMax(item)
            select case(rank)
              case(1)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp1d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp1d)
                maxValue = maxval(fp1d)
              case(2)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp2d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp2d)
                maxValue = maxval(fp2d)
              case(3)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp3d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp3d)
                maxValue = maxval(fp3d)
              case(4)
                call ESMF_FieldGet(fieldList(item), localDe=localDe, farrayPtr=fp4d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp4d)
                maxValue = maxval(fp4d)
              case default
                call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="Field rank not implemented.", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return ! bail out
            end select
            localMin(item) = min(minValue, localMin(item))
            localMax(item) = max(maxValue, localMax(item))
          end do
        end if
      end do

      ! -- compute global min and max values
      call ESMF_VMAllReduce(vm, localMin, globalMin, fieldCount, ESMF_REDUCE_MIN, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call ESMF_VMAllReduce(vm, localMax, globalMax, fieldCount, ESMF_REDUCE_MAX, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      ! -- log results
      do item = 1, fieldCount
        if (connectedList(item) == "true") then
          write(msgString,'(a,": ",a,": ",a,"[",i0,"]: local/global min/max =",4g20.8)') &
            trim(name), rName, fieldNames(item)(1:maxLength), item, localMin(item), &
            localMax(item), globalMin(item), globalMax(item)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
        end if
      end do

      deallocate(localMin, localMax, globalMin, globalMax, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

      if (associated(fieldNames)) then
        deallocate(fieldNames, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if
      if (associated(connectedList)) then
        deallocate(connectedList, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if
      if (associated(fieldList)) then
        deallocate(fieldList, stat=localrc)
        if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
          line=__LINE__, &
          file=__FILE__, &
          rcToReturn=rc)) return
      end if
    end if

  end subroutine AerosolFieldDiagnostics

  subroutine MAPLFieldDiagnostics(model, state, label, rc)
    type(ESMF_GridComp)            :: model
    type(ESMF_State)               :: state
    character(len=*),  intent(in)  :: label
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, item
    integer :: localDe, localDeCount, rank
    integer :: fieldCount, maxLength
    character(len=ESMF_MAXSTR)          :: msgString
    character(len=ESMF_MAXSTR)          :: name
    character(len=ESMF_MAXSTR), pointer :: shapeList(:)
    character(len=ESMF_MAXSTR), pointer :: fieldNames(:)
    real(ESMF_KIND_R4)                              :: minValue, maxValue
    real(ESMF_KIND_R4), dimension(:),       pointer :: localMin, localMax, globalMin, globalMax
    real(ESMF_KIND_R4), dimension(:),       pointer :: fp1d
    real(ESMF_KIND_R4), dimension(:,:),     pointer :: fp2d
    real(ESMF_KIND_R4), dimension(:,:,:),   pointer :: fp3d
    real(ESMF_KIND_R4), dimension(:,:,:,:), pointer :: fp4d
    type(ESMF_Field) :: field
    type(ESMF_VM)    :: vm
    type(ESMF_FieldStatus_Flag) :: fieldStatus
    type(ESMF_StateItem_Flag), pointer :: itemTypeList(:)

    ! -- local parameters
    character(len=*), parameter :: rName = "MAPL: diagnostic"

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    ! -- get component information
    call NUOPC_CompGet(model, name=name, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    ! -- query the Component for its VM and importState
    call ESMF_GridCompGet(model, vm=vm, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    call ESMF_StateGet(state, itemCount=fieldCount, nestedFlag=.true., rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    write(msgString,'(a,": ",a,": ",a,": found ",i0," items")') &
      trim(name), rName, trim(label), fieldCount
    call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    allocate(fieldNames(fieldCount), itemTypeList(fieldCount), shapeList(fieldCount))
    fieldNames = ""
    shapeList = ""

    call ESMF_StateGet(state, itemNameList=fieldNames, &
      itemTypeList=itemTypeList, nestedFlag=.true., rc=localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    nullify(localMin, localMax, globalMin, globalMax)

    ! -- check values of imported fields, if requested
    if (fieldCount > 0) then

      allocate(localMin(fieldCount), localMax(fieldCount), &
        globalMin(fieldCount), globalMax(fieldCount), stat=localrc)
      if (ESMF_LogFoundAllocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return

      localMin = huge(0._ESMF_KIND_R4)
      localMax = -localMin
      globalMin = 0._ESMF_KIND_R4
      globalMax = 0._ESMF_KIND_R4

      ! -- find longest field name for formatting purpose
      maxLength = 0
      do item = 1, fieldCount
        maxLength = max(maxLength, len_trim(fieldNames(item)))
      end do

      do item = 1, fieldCount
        if (itemTypeList(item) == ESMF_STATEITEM_FIELD) then
          ! --- get field data
          call ESMF_StateGet(state, trim(fieldNames(item)), field, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          localDeCount = 0
          rank = 0
          call ESMF_FieldGet(field, status=fieldStatus, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
          if (fieldStatus == ESMF_FIELDSTATUS_COMPLETE) then
            call ESMF_FieldGet(field, rank=rank, &
              localDeCount=localDeCount, rc=localrc)
            if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
              line=__LINE__,  &
              file=__FILE__,  &
              rcToReturn=rc)) return  ! bail out
          else
            shapeList(item) = "incomplete"
          end if

          do localDe = 0, localDeCount - 1
            minValue = localMin(item)
            maxValue = localMax(item)
            select case(rank)
              case(1)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp1d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp1d)
                maxValue = maxval(fp1d)
                write(shapeList(item), '(i0,":",i0)') lbound(fp1d,1), ubound(fp1d,1)
              case(2)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp2d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp2d)
                maxValue = maxval(fp2d)
                write(shapeList(item), '(i0,":",i0,",",i0,":",i0)') (lbound(fp2d,i), ubound(fp2d,i), i=1,2)
              case(3)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp3d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp3d)
                maxValue = maxval(fp3d)
                write(shapeList(item), '(i0,":",i0,2(",",i0,":",i0))') (lbound(fp3d,i), ubound(fp3d,i), i=1,3)
              case(4)
                call ESMF_FieldGet(field, localDe=localDe, farrayPtr=fp4d, rc=localrc)
                if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)) return  ! bail out
                minValue = minval(fp4d)
                maxValue = maxval(fp4d)
                write(shapeList(item), '(i0,":",i0,3(",",i0,":",i0))') (lbound(fp4d,i), ubound(fp4d,i), i=1,4)
              case default
                call ESMF_LogSetError(ESMF_RC_NOT_IMPL, &
                  msg="Field rank not implemented.", &
                  line=__LINE__, &
                  file=__FILE__, &
                  rcToReturn=rc)
                return ! bail out
            end select
            localMin(item) = min(minValue, localMin(item))
            localMax(item) = max(maxValue, localMax(item))
          end do
        end if
      end do

      ! -- compute global min and max values
      call ESMF_VMAllReduce(vm, localMin, globalMin, fieldCount, ESMF_REDUCE_MIN, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      call ESMF_VMAllReduce(vm, localMax, globalMax, fieldCount, ESMF_REDUCE_MAX, rc=localrc)
      if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__,  &
        file=__FILE__,  &
        rcToReturn=rc)) return  ! bail out

      ! -- log results
      do item = 1, fieldCount
        if (itemTypeList(item) == ESMF_STATEITEM_FIELD) then
          write(msgString,'(a,": ",a,": ",a,"[",i0,"]: [",a,"] local/global min/max =",4g20.8)') &
            trim(name), rName, fieldNames(item)(1:maxLength), item, trim(shapeList(item)), &
            localMin(item), localMax(item), globalMin(item), globalMax(item)
          call ESMF_LogWrite(trim(msgString), ESMF_LOGMSG_INFO, rc=localrc)
          if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__,  &
            file=__FILE__,  &
            rcToReturn=rc)) return  ! bail out
        end if
      end do

      deallocate(localMin, localMax, globalMin, globalMax, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

    if (associated(fieldNames)) then
      deallocate(fieldNames, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if
    if (associated(itemTypeList)) then
      deallocate(itemTypeList, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if
    if (associated(shapeList)) then
      deallocate(shapeList, stat=localrc)
      if (ESMF_LogFoundDeallocError(statusToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
        line=__LINE__, &
        file=__FILE__, &
        rcToReturn=rc)) return
    end if

  end subroutine MAPLFieldDiagnostics

end module Aerosol_Comp_Mod
