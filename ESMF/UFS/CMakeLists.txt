esma_set_this (OVERRIDE UFS_Aerosols)

set (srcs
  Aerosol_Cap.F90
  Aerosol_Comp_Mod.F90
  Aerosol_Diag_Mod.F90
  Aerosol_Internal_Mod.F90
  Aerosol_Logger_Mod.F90
  Aerosol_Shared_Mod.F90
  Aerosol_Tracer_Mod.F90
  )

if(UFS_TRACING)
set (dependencies Aerosol_GridComp MAPL ESMF::ESMF ufs_tracing)
else()
set (dependencies Aerosol_GridComp MAPL ESMF::ESMF)
endif()
esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES ${dependencies})
