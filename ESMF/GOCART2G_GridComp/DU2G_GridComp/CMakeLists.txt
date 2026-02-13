esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES GA_Environment MAPL MAPL.generic3g MAPL.utilities Chem_Shared2G Process_Library ESMF::ESMF NetCDF::NetCDF_Fortran
  TYPE SHARED)

mapl_acg (
  ${this} DU2G_StateSpecs.rc
  IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS GET_POINTERS DECLARE_POINTERS
  3g
)

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()

# CMake has an OpenMP issue with NAG Fortran: https://gitlab.kitware.com/cmake/cmake/-/issues/21280
if (NOT CMAKE_Fortran_COMPILER_ID MATCHES "NAG")
   target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()
