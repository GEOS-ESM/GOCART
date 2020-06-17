esma_set_this ()

set (srcs
  Chem_MieTableMod2G.F90
  Chem_AeroGeneric.F90
  GOCART2G_Process.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES MAPL GMAO_mpeu)

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

#include_directories (${esma_include}/${this})

