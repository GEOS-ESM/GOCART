module Mediator_InterpolatorMod
  use ESMF
  use MAPL_Mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  implicit none

  private
  public Mediator_Interpolator

  type Mediator_Interpolator
     type(ESMF_Time) :: t0, t1 ! interpolation window
     type(ESMF_State) :: state0, state1 !states to interpolate between
     type(ESMF_TimeInterval) :: interval
   contains
     procedure :: update
     procedure :: interpolate
  end type Mediator_Interpolator

  interface Mediator_Interpolator
     module procedure new_Mediator_Interpolator
  end interface Mediator_Interpolator

  interface lerp
     module procedure lerp_r4, lerp_r8
  end interface lerp

contains

  function new_Mediator_Interpolator(start_time, interp_interval) result(interpolator)
    type(mediator_interpolator) :: interpolator
    type(ESMF_Time), intent(in) :: start_time
    type(ESMF_TimeInterval), intent(in) :: interp_interval
    interpolator%interval = interp_interval
    interpolator%t0 = start_time
    interpolator%t1 = start_time + interpolator%interval
  end function new_mediator_interpolator

  subroutine update(this)
    class(mediator_interpolator) :: this
    ! update the time window for interpolation
    this%t0 = this%t0 + this%interval
    this%t1 = this%t0 + this%interval
  end subroutine update


  subroutine interpolate(this, current_time, in_field0, in_field1, out_field)
    class(mediator_interpolator) :: this
    type(ESMF_Time) :: current_time

    type(ESMF_Field) :: in_field0, in_field1, out_field
    integer :: field_rank
    type(ESMF_TypeKind_Flag) :: field_type_kind
    integer :: rc
    real(real32) :: t

    real(real32), pointer :: in0_ptr2d_r4(:,:), in1_ptr2d_r4(:,:), out_ptr2d_r4(:,:)
    real(real32), pointer :: in0_ptr3d_r4(:,:,:), in1_ptr3d_r4(:,:,:), out_ptr3d_r4(:,:,:)

    real(real64), pointer :: in0_ptr2d_r8(:,:), in1_ptr2d_r8(:,:), out_ptr2d_r8(:,:)
    real(real64), pointer :: in0_ptr3d_r8(:,:,:), in1_ptr3d_r8(:,:,:), out_ptr3d_r8(:,:,:)

    character(len=ESMF_MAXSTR) :: field_name

    if (current_time < this%t0 .or. current_time > this%t1) then
       print *, "interpolation time outside interpolation range"
       stop
    end if

    call ESMF_FieldGet(in_field0, name = field_name, rc = rc)

    ! add checks that all fields match
    call ESMF_FieldGet(in_field0, rank = field_rank, typekind = field_type_kind, rc = rc)
    if (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, &
         file=__FILE__)) &
         return  ! bail out

    t = (current_time - this%t0) / (this%t1 - this%t0)

    if (field_rank == 2 .and. field_type_kind == ESMF_TYPEKIND_R4) then
       call ESMF_FieldGet(in_field0, farrayPtr = in0_ptr2d_r4, rc = rc)
       call ESMF_FieldGet(in_field1, farrayPtr = in1_ptr2d_r4, rc = rc)
       call ESMF_FieldGet(out_field, farrayPtr = out_ptr2d_r4, rc = rc)
       call lerp(in0_ptr2d_r4, in1_ptr2d_r4, t, out_ptr2d_r4)
    else if(field_rank == 3 .and. field_type_kind == ESMF_TYPEKIND_R4) then
       call ESMF_FieldGet(in_field0, farrayPtr = in0_ptr3d_r4, rc = rc)
       call ESMF_FieldGet(in_field1, farrayPtr = in1_ptr3d_r4, rc = rc)
       call ESMF_FieldGet(out_field, farrayPtr = out_ptr3d_r4, rc = rc)
       call lerp(in0_ptr3d_r4, in1_ptr3d_r4, t, out_ptr3d_r4)
    else if (field_rank == 2 .and. field_type_kind == ESMF_TYPEKIND_R8) then
       call ESMF_FieldGet(in_field0, farrayPtr = in0_ptr2d_r8, rc = rc)
       call ESMF_FieldGet(in_field1, farrayPtr = in1_ptr2d_r8, rc = rc)
       call ESMF_FieldGet(out_field, farrayPtr = out_ptr2d_r8, rc = rc)
       call lerp(in0_ptr2d_r8, in1_ptr2d_r8, t, out_ptr2d_r8)
    else if(field_rank == 3 .and. field_type_kind == ESMF_TYPEKIND_R8) then
       call ESMF_FieldGet(in_field0, farrayPtr = in0_ptr3d_r8, rc = rc)
       call ESMF_FieldGet(in_field1, farrayPtr = in1_ptr3d_r8, rc = rc)
       call ESMF_FieldGet(out_field, farrayPtr = out_ptr3d_r8, rc = rc)
       call lerp(in0_ptr3d_r8, in1_ptr3d_r8, t, out_ptr3d_r8)
    end if

  end subroutine interpolate

  elemental subroutine lerp_r4(v0, v1, t, v_out)
    real(real32), intent(in) :: v0, v1, t      
    real(real32), intent(out) :: v_out
    v_out = (1.0 - t) * v0 + t * v1
  end subroutine lerp_r4

  elemental subroutine lerp_r8(v0, v1, t, v_out)
    real(real64), intent(in) :: v0, v1
    real(real32), intent(in) :: t
    real(real64), intent(out) :: v_out
    v_out = (1.0 - t) * v0 + t * v1
  end subroutine lerp_r8

end module Mediator_InterpolatorMod
