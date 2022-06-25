# UFS/GOCART interface
#
# UFS porting of ESMA_cmake package.
# See original version at: https://github.com/GEOS-ESM/ESMA_cmake (tag: v3.4.2)

macro (esma_add_library this)

  set (options EXCLUDE_FROM_ALL NOINSTALL)
  set (oneValueArgs
    # shared with ecbuild
    TYPE)
  set (multiValueArgs
    # esma unique
    SUBCOMPONENTS SUBDIRS NEVER_STUB PRIVATE_DEFINITIONS PUBLIC_DEFINITIONS
    # shared with ecbuild (and not deprecated)
    SOURCES DEPENDS PUBLIC_LIBS
    # deprecated in esma (produces explicit warnings)
    SRCS INCLUDES DEPENDENCIES
    # deprecated in ecbuild
    PUBLIC_INCLUDES
    )
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (ARGS_UNPARSED_ARGUMENTS)
      message (WARNING "Unrecognized keyword arguments passed to esma_add_library: ${ARGS_UNPARSED_ARGUMENTS}")
  endif ()

  # Subdirs must exist and should be configured prior to subcomponents.
  foreach (subdir ${ARGS_SUBDIRS})
    add_subdirectory(${subdir})
  endforeach()

  # Handle deprecated
  if (ARGS_SRCS)
    set (ARGS_SOURCES ${ARGS_SRCS})
  endif ()
  if (ARGS_INCLUDES)
    set (ARGS_PUBLIC_INCLUDES ${ARGS_INCLUDES})
  endif ()
  if (ARGS_DEPENDENCIES)
    set (ARGS_PUBLIC_LIBS ${ARGS_DEPENDENCIES})
  endif ()

  # Configure subcomponents.  These can be stubbed and may have a
  # different name than the directory they reside in.  (Most
  # unfortunate.)
  set (non_stubbed)
  foreach (subdir ${ARGS_SUBCOMPONENTS})

    string(REPLACE "@" "" mod_name ${subdir})

    if (NOT rename_${subdir}) # usual case
      set (module_name ${mod_name})
    else ()
      set(module_name ${rename_${mod_name}})
    endif ()

    esma_add_subdirectory(${mod_name} FOUND found)
    if (found)
      list (APPEND non_stubbed ${mod_name})
    else ()
      message(WARNING "sub-component ${module_name} NOT FOUND")
    endif ()

  endforeach ()

  # This library depends on all DEPENDENCIES and _non-stubbed_ subcomponents.
  set (all_dependencies ${ARGS_PUBLIC_LIBS} ${non_stubbed})
  if (ARGS_TYPE)
    set(ARGS_TYPE TYPE ${ARGS_TYPE})
  endif()
  if (ARGS_NOINSTALL)
    set(NOINSTALL_ "NOINSTALL")
  endif()

  add_library(${this} ${ARGS_SOURCES})
  target_link_libraries (${this} PUBLIC ${all_dependencies} )
  target_include_directories (${this} PUBLIC ${ARGS_PUBLIC_INCLUDES} )

  set (CPP_DEBUG_${this} "" CACHE STRING "List of files to pass -DDEBUG")
  foreach (file ${CPP_DEBUG_${this}})
    message (STATUS "setting debug option for ${file}")
    set_source_files_properties (${file} PROPERTIES COMPILE_DEFINITIONS -DDEBUG)
  endforeach()

  set_target_properties (${this} PROPERTIES EXCLUDE_FROM_ALL ${ARGS_EXCLUDE_FROM_ALL})
  set_target_properties (${this} PROPERTIES Fortran_MODULE_DIRECTORY ${_ufs_esma_include}/${this})

  set (install_dir include/${this})
  # Export target  include directories for other targets
  target_include_directories(${this} PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}> # stubs
# modules and copied *.h, *.inc
    $<BUILD_INTERFACE:${_ufs_esma_include}/${this}>
    $<INSTALL_INTERFACE:${install_dir}>
    )

  if (ARGS_PUBLIC_INCLUDES)
    target_include_directories(${this} PUBLIC $<BUILD_INTERFACE:${ARGS_PUBLIC_INCLUDES}>)
  endif ()

  if (ARGS_PRIVATE_DEFINITIONS)
    target_compile_definitions(${this} PRIVATE ${ARGS_PRIVATE_DEFINITIONS})
  endif ()
  if (ARGS_PUBLIC_DEFINITIONS)
    target_compile_definitions(${this} PUBLIC ${ARGS_PUBLIC_DEFINITIONS})
  endif ()

  export (TARGETS ${this} APPEND FILE "${PROJECT_TARGETS_FILE}" )

  install ( TARGETS ${this}
            EXPORT  ${PROJECT_NAME}-targets
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} )

endmacro ()
