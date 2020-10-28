esma_set_this (OVERRIDE Chem_Shared2G)

set (srcs
  Chem_AeroGeneric.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES MAPL)

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/@GSW )
    set (gsw_ECBUILD_SYSTEM_INCLUDED TRUE)
endif()

esma_add_subdirectories (
        MAPL
        GMAO_Shared
        Configuration
)
