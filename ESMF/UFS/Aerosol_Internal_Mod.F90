module Aerosol_Internal_Mod

  use MAPL
  use gFTL_StringIntegerMap

  implicit none

  type Aerosol_InternalData_T
    type(MAPL_Cap), pointer :: maplCap
    type(StringIntegerMap)  :: tracerMap
  end type Aerosol_InternalData_T

  type Aerosol_InternalState_T
    type(Aerosol_InternalData_T), pointer :: wrap => null()
  end type Aerosol_InternalState_T

  private

  public :: Aerosol_InternalState_T
  public :: Aerosol_InternalData_T

end module Aerosol_Internal_Mod
