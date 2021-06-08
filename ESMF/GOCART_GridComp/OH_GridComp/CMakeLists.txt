esma_set_this()

add_definitions (-DGEOS5)

set (srcs 
 calc_oh.F90
 cblock_cf_mod.F90
 cblock_CO_budget_mod.F90
 cblock_CO_mod.F90
 cblock_OH_mod.F90
 cblock_size_mod.F90
 computeParameters_mod.F90
 computeReflectivity_mod.F90
 correctOH.F90
 donefunctions01.F90
 donefunctions02.F90
 donefunctions03.F90
 donefunctions04.F90
 donefunctions05.F90
 donefunctions06.F90
 donefunctions07.F90
 donefunctions08.F90
 donefunctions09.F90
 donefunctions10.F90
 donefunctions11.F90
 donefunctions12.F90
 donefunctions13.F90
 donefunctions14.F90
 donefunctions15.F90
 donefunctions16.F90
 donefunctions17.F90
 donefunctions18.F90
 donefunctions19.F90
 donefunctions20.F90
 donefunctions21.F90
 donefunctions22.F90
 donefunctions23.F90
 donefunctions24.F90
 donefunctions25.F90
 donefunctions26.F90
 donefunctions27.F90
 donefunctions28.F90
 donefunctions29.F90
 donefunctions30.F90
 donefunctions31.F90
 donefunctions32.F90
 donefunctions33.F90
 donefunctions34.F90
 donefunctions35.F90
 donefunctions36.F90
 func_ifs.F90
 getinfo_mod.F90
 getoh_mod.F90
 GNUmakefile
 identification_mod.F90
 inputFileReading_mod.F90
 interpOH.F90
 OH_GridCompMod.F90
 ohParameterizationMethod_mod.F90
 polint.F90
 solarZenithAngle_mod.F90
 )

esma_add_library (${this}
  SRCS ${srcs}
  DEPENDENCIES MAPL Chem_Shared Chem_Base GMAO_mpeu MAPL_cfio_r4)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

#esma_generate_gocart_code (${this} "-B\;-E\;-C\;-N\;GOCART")

file (GLOB_RECURSE rc_files CONFIGURE_DEPENDS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.rc)
foreach ( file ${rc_files} )
   get_filename_component( dir ${file} DIRECTORY )
   install( FILES ${file} DESTINATION etc/${dir} )
endforeach()
