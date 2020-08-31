esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES MAPL Chem_Shared2G Chem_Shared)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

#esma_generate_gocart_code (${this} "-B\;-E\;-N\;DU2G")
mapl_acg (${this}   DU2G_StateSpecs.rc 
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS 
          GET_POINTERS DECLARE_POINTERS)

file (GLOB_RECURSE rc_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
