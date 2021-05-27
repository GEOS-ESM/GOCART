# UFS/GOCART interface
#
# UFS porting of ESMA_cmake package.
# See original version at: https://github.com/GEOS-ESM/ESMA_cmake (tag: v3.4.2)

function (esma_add_subdirectory dir)

  set (options)
  set (oneValueArgs FOUND)
  set (multiValueArgs)
  cmake_parse_arguments(ARGS "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

  if (ARGS_FOUND)
    set (${ARGS_FOUND} FALSE PARENT_SCOPE)
  endif ()

  if (IS_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${dir})
    add_subdirectory (${dir})
    if (ARGS_FOUND)
      set (${ARGS_FOUND} TRUE PARENT_SCOPE)
    endif ()
  else ()
    message(STATUS "Directory not found ${dir} (possibly sparse checkout)")
  endif ()

endfunction ()

function (esma_add_subdirectories dirs)

  set (dirs_ ${dirs} ${ARGN})
  message (DEBUG "esma_add_subdirectories:  ${dirs}")

  foreach (subdir ${dirs_})
    esma_add_subdirectory (${subdir})
  endforeach()

endfunction ()
