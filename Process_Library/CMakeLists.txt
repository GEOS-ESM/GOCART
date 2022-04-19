esma_set_this ()

set (srcs
  Chem_MieTableMod2G.F90
  GOCART2G_Process.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES esmf NetCDF::NetCDF_Fortran
  )
if (OMP)
  target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif()
