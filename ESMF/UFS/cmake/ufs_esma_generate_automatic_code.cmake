# UFS/GOCART interface
#
# UFS porting of ESMA_cmake package.
# See original version at: https://github.com/GEOS-ESM/ESMA_cmake (tag: v3.4.2)

set (acg_flags -v)

macro (new_ufs_esma_generate_automatic_code
    target registry headers rcs rcs_destination flags)

  find_file (generator
    NAME mapl_acg.pl
    )

  add_custom_command (
    OUTPUT ${rcs}
    BYPRODUCTS ${headers}
    COMMAND ${generator} ${acg_flags} ${flags} ${CMAKE_CURRENT_SOURCE_DIR}/${registry}
    COMMAND ${CMAKE_COMMAND} -E copy ${rcs} ${rcs_destination}
    COMMAND touch foo
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/${registry}
    DEPENDS ${generator}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating automated code for ${registry}"
    )
  add_custom_target (phony_${target} DEPENDS ${rcs})
  add_dependencies (${target} phony_${target})
  install(FILES ${rcs_destination}/${rcs} DESTINATION etc)

endmacro ()

macro (esma_generate_gocart_code target flags)

  string (REPLACE "_GridComp" "" name ${target})

  set (automatic_headers
    ${name}_ExportSpec___.h
    ${name}_GetPointer___.h
    )
  set (automatic_rc
    ${name}_History___.rc
    )

  set (registry ${name}_Registry.rc)

  new_ufs_esma_generate_automatic_code (
    ${target} ${registry}
    "${automatic_headers}" "${automatic_rc}"
    ${_ufs_esma_etc}
    ${flags}
  )

endmacro ()
