module Aerosol_Diag_Mod

  use ESMF
  use gFTL_StringIntegerMap
  use MAPL, only: MAPL_Cap

  use Aerosol_Internal_Mod, only: Aerosol_InternalState_T
  use Aerosol_Shared_Mod,   only: AerosolGetPtr
  use Aerosol_Tracer_Mod,   only: Aerosol_Tracer_T

  implicit none

  interface AerosolDiagUpdate
    module procedure ComputePM
  end interface

  private

  ! -- methods
  public :: AerosolDiagUpdate

contains

  subroutine ComputePM(model, rc)
    type(ESMF_GridComp)            :: model
    integer, optional, intent(out) :: rc

    ! -- local variables
    integer :: localrc
    integer :: i, j, k, kk, n, ni, nj, nk
    integer :: bin, nbins, s
    integer, pointer                :: idx, pm
    real(ESMF_KIND_R8)              :: w(4)
    real(ESMF_KIND_R4), pointer     :: rho(:,:,:)
    real(ESMF_KIND_R8), pointer     :: q(:,:,:,:)
    type(ESMF_Field)                :: field
    type(ESMF_State)                :: state
    type(MAPL_Cap), pointer         :: cap
    type(Aerosol_InternalState_T)   :: is
    type(Aerosol_Tracer_T), pointer :: trp
    type(StringIntegerMapIterator)  :: iter

    ! -- local parameters
    character(len=*), parameter :: pm_size(2) = [ 'PM25', 'PM10' ]

    ! -- begin
    if (present(rc)) rc = ESMF_SUCCESS

    nullify(cap, trp)

    ! -- retrieve model internal state
    call ESMF_GridCompGetInternalState(model, is, localrc)
    if (ESMF_LogFoundError(rcToCheck=localrc, msg=ESMF_LOGERR_PASSTHRU, &
      line=__LINE__,  &
      file=__FILE__,  &
      rcToReturn=rc)) return  ! bail out

    cap => is % wrap % maplCap
    if (.not.associated(cap)) return
    trp => is % wrap % tracers
    if (.not.associated(trp)) return

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

    do n = 1, size(pm_size)
      pm => trp % indexMap % at(pm_size(n))
      if (associated(pm)) then
        q(:,:,:,pm) = 0._ESMF_KIND_R8
        iter = trp % indexMap % begin()
        do while (iter /= trp % indexMap % end())
          idx => iter % value()
          nbins = PMGetTracerWeight(iter % key(), pm_size(n), w)
          s = idx
          do bin = 1, nbins
            do k = 1, nk
              kk = nk - k + 1
              do j = 1, nj
                do i = 1, ni
                  ! -- compute PM (ug/m3)
                  q(i,j,k,pm) = q(i,j,k,pm) + w(bin) * q(i,j,k,s) * rho(i,j,kk)
                end do
              end do
            end do
            s = s + 1
          end do
          call iter % next()
        end do
      end if
    end do

  end subroutine ComputePM

  function PMGetTracerWeight(name, pm, w) result(nbins)
    character(len=*),   intent(in)  :: name
    character(len=*),   intent(in)  :: pm
    real(ESMF_KIND_R8), intent(out) :: w(4)

    ! -- local variables
    integer :: nbins

    ! -- local parameters
    real(ESMF_KIND_R8), parameter :: zero = 0.0_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: one  = 1.0_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: w25_du2    = log(1.250_ESMF_KIND_R8) / log(1.8_ESMF_KIND_R8)
    real(ESMF_KIND_R8), parameter :: w_du4      = log(1.667_ESMF_KIND_R8) / log(2.0_ESMF_KIND_R8)
    real(ESMF_KIND_R8), parameter :: w_so4      = 132.14_ESMF_KIND_R8 / 96.06_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: w_no3      = 80.043_ESMF_KIND_R8 / 62.0_ESMF_KIND_R8
    real(ESMF_KIND_R8), parameter :: w25_ss3    = log(2.50_ESMF_KIND_R8) / log(3._ESMF_KIND_R8)
    real(ESMF_KIND_R8), parameter :: w10_no3an2 = 0.808_ESMF_KIND_R8 * w_no3
    real(ESMF_KIND_R8), parameter :: w25_no3an2 = 0.138_ESMF_KIND_R8 * w_no3
    real(ESMF_KIND_R8), parameter :: w10_no3an3 = 0.164_ESMF_KIND_R8 * w_no3

    ! -- begin
    w = one

    nbins = 1
    select case (trim(name))
      case ('CAphilicCA.bc', &
            'CAphilicCA.oc', &
            'CAphobicCA.bc', &
            'CAphobicCA.oc')
       ! -- set to 1 by default
      case ('DU')
        select case (trim(pm))
          case ('PM10')
            nbins = 4
            w(4) = w_du4
          case ('PM25')
            nbins = 2
            w(2) = w25_du2
        end select
      case ('NO3an1')
        w(1) = w_no3
      case ('NO3an2')
        select case (trim(pm))
          case ('PM10')
            w(1) = w10_no3an2
          case ('PM25')
            w(1) = w25_no3an2
        end select
      case ('NO3an3')
        w(1) = w10_no3an3
      case ('SO4')
        w(1) = w_so4
      case ('SS')
        select case (trim(pm))
          case ('PM10')
            nbins = 4
          case ('PM25')
            nbins = 3
            w(3) = w25_ss3
        end select
      case default
        nbins = 0
    end select

  end function PMGetTracerWeight

end module Aerosol_Diag_Mod
