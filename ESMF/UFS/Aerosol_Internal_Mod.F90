module Aerosol_Internal_Mod

  use MAPL
  use gFTL_StringIntegerMap
  use gFTL_StringStringMap

  use Aerosol_Tracer_Mod, only: Aerosol_Tracer_T

  implicit none

  type Aerosol_InternalData_T
    type(MAPL_Cap),         pointer :: maplCap => null()
    type(Aerosol_Tracer_T), pointer :: tracers => null()
  end type Aerosol_InternalData_T

  type Aerosol_InternalState_T
    type(Aerosol_InternalData_T), pointer :: wrap => null()
  end type Aerosol_InternalState_T

  private

  public :: Aerosol_InternalState_T
  public :: Aerosol_InternalData_T

end module Aerosol_Internal_Mod
