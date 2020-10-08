esma_set_this (OVERRIDE Chem_Shared2G)

set (srcs
  Chem_AeroGeneric.F90
  )

esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES MAPL)

if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/@GSW )
    set (gsw_ECBUILD_SYSTEM_INCLUDED TRUE)
endif()

if ( IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/MAPL@ )
  esma_add_subdirectories (
    MAPL
    )
endif ()

if ( IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/GMAO_Shared@ )
  esma_add_subdirectories (
    GMAO_Shared
    )
endif ()
