esma_set_this (OVERRIDE UFS_Aerosols)

set (srcs
  Aerosol_Cap.F90
  Aerosol_Comp_Mod.F90
  )

set (dependencies Aerosol_GridComp MAPL esmf)
esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES ${dependencies})
