esma_set_this ()

set (srcs
  GOCART2G_MieMod.F90
  GOCART2G_Process.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES esmf NetCDF::NetCDF_Fortran
  )
