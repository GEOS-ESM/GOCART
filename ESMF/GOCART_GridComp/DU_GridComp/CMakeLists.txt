esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES Chem_Shared Chem_Base GMAO_mpeu)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

esma_generate_gocart_code (${this} "-B\;-E\;-C\;-N\;GOCART")

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
