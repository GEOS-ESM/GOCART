cmake_minimum_required (VERSION 3.17)
cmake_policy (SET CMP0053 NEW)
cmake_policy (SET CMP0054 NEW)

# Let's look for the ESMA cmake directory in a few places
# to see if we are building standalone
set (ESMA_CMAKE_DIRS
  cmake
  @cmake
  cmake@
  )

foreach (dir IN LISTS ESMA_CMAKE_DIRS)
  if (EXISTS ${CMAKE_CURRENT_LIST_DIR}/${dir})
    list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/${dir}")
    set (ESMA_CMAKE_PATH "${CMAKE_CURRENT_LIST_DIR}/${dir}" CACHE PATH "Path to ESMA_cmake code")
    set(GOCART_STANDALONE TRUE)
  endif ()
endforeach ()

if (GOCART_STANDALONE)
  project (
          GOCART
          VERSION 2.3.0
          LANGUAGES Fortran CXX C)  # Note - CXX is required for ESMF

  if ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
      message(SEND_ERROR "In-source builds are disabled. Please
            issue cmake command in separate build directory.")
  endif ("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

  # Set the default build type to release
  if (NOT CMAKE_BUILD_TYPE)
      message (STATUS "Setting build type to 'Release' as none was specified.")
      set (CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
      # Set the possible values of build type for cmake-gui
      set_property (CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
              "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
  endif ()
endif()

# Set build options
option (UFS_GOCART "Build GOCART component for UFS" OFF)

set (DOING_GEOS5 YES)

# Should find a better place for this - used in Chem component
set (ACG_FLAGS -v)

if (UFS_GOCART)
# Ensure we build as 32-bit
  message ("Force 32-bit build for GOCART")
  if (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    string (REPLACE "-real-size 64" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    string (REPLACE "-fdefault-real-8" "" CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}")
  endif()
  list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/ESMF/UFS/cmake")
endif()

if (NOT COMMAND esma)
    include (esma)
endif()

if(GOCART_STANDALONE)
  ecbuild_declare_project()
endif()

if (NOT Baselibs_FOUND)
  # Find dependencies
  if (NOT NetCDF_Fortran_FOUND)
    find_package (NetCDF REQUIRED Fortran)
  endif ()
  add_definitions(-DHAS_NETCDF3)

  if (NOT ESMF_FOUND)
    find_package (ESMF REQUIRED)
  endif ()
endif ()

if (UFS_GOCART)
  find_package (GFTL_SHARED REQUIRED)
  # Dom Heinzeller 2023/08/30 - workaround until https://github.com/GEOS-ESM/MAPL/pull/2320
  # is merged and finds its way into the ufs-weather-model dependency tree
  find_package (YAFYAML QUIET)
  find_package (FARGPARSE QUIET)
  find_package (PFLOGGER QUIET)
  #
  find_package (MAPL REQUIRED)
  include(mapl_acg)
elseif (IS_DIRECTORY "${PROJECT_SOURCE_DIR}/ESMF/Shared/MAPL@")
  message (STATUS "Using local MAPL source tree: ${PROJECT_SOURCE_DIR}/ESMF/Shared/MAPL@")
endif ()


# This is a KLUDGE to emulate paths used in GEOSgcm so that GOCART Legacy code generation works
file (MAKE_DIRECTORY ${esma_include}/GEOSchem_GridComp)

# Generic DFLAGS
# These should be relocated and/or eliminated.
add_definitions(-Dsys${CMAKE_SYSTEM_NAME} -DESMA64)
add_definitions(${MPI_Fortran_COMPILE_FLAGS})
include_directories(${MPI_Fortran_INCLUDE_PATH})

# Recursively build source tree
add_subdirectory (ESMF)
add_subdirectory (Process_Library)

if(GOCART_STANDALONE)
  ecbuild_install_project (NAME GOCART)
endif()

# https://www.scivision.dev/cmake-auto-gitignore-build-dir/
# --- auto-ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()

# Piggyback that file into install
install(
   FILES ${PROJECT_BINARY_DIR}/.gitignore
   DESTINATION ${CMAKE_INSTALL_PREFIX}
   )
