esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES MAPL GA_Environment Chem_Shared2G Process_Library esmf NetCDF::NetCDF_Fortran)

mapl_acg (${this}   CA2G_StateSpecs.rc 
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS 
          GET_POINTERS DECLARE_POINTERS)

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc *.yaml)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()


