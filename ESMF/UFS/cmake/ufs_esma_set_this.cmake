# UFS/GOCART interface
#
# UFS porting of ESMA_cmake package.
# See original version at: https://github.com/GEOS-ESM/ESMA_cmake (tag: v3.4.2)

set (_ufs_esma_include  ${CMAKE_BINARY_DIR}/include CACHE PATH "include directory")
set (_ufs_esma_etc      ${CMAKE_BINARY_DIR}/etc     CACHE PATH "etc directory")

file (MAKE_DIRECTORY ${_ufs_esma_include})
file (MAKE_DIRECTORY ${_ufs_esma_etc})

# set for compatibility
set (esma_include ${_ufs_esma_include})

macro (esma_set_this)

  set (options OPTIONAL EXCLUDE_FROM_ALL)
  set (oneValueArgs OVERRIDE)
  set (multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (ARGS_UNPARSED_ARGUMENTS)
    message (FATAL_ERROR "esma_set_this unparsed arguments: ${ARGS_UNPARSED_ARGUMENTS}")
  endif ()

  # (re) define variable "this"
  if (ARGS_OVERRIDE)
    set (this ${ARGS_OVERRIDE})
  else ()
    get_filename_component (dir "${CMAKE_CURRENT_BINARY_DIR}" NAME)
    string(REPLACE "@" "" this ${dir})
  endif ()

  set (include_${this} ${_ufs_esma_include}/${this})

  file (MAKE_DIRECTORY ${_ufs_esma_include}/${this})
  file (MAKE_DIRECTORY ${_ufs_esma_etc}/${this})

endmacro ()

