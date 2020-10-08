esma_set_this ()

set (srcs
  Aerosol_Cap.F90
  )

set (dependencies GOCART2G_GridComp esmf)
esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES ${dependencies})
