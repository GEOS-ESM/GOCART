cmake_minimum_required (VERSION 3.17)
cmake_policy (SET CMP0053 NEW)
cmake_policy (SET CMP0054 NEW)

project (
        GOCART
        VERSION 2.0.7
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
    list (APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake@")
    include (esma)
endif()

ecbuild_declare_project()

if (NOT Baselibs_FOUND)
  # Find dependencies
  if (NOT NetCDF_Fortran_FOUND)
    find_package (NetCDF REQUIRED Fortran)
  endif ()
  add_definitions(-DHAS_NETCDF3)

  if (NOT ESMF_FOUND)
    find_package (ESMF REQUIRED)
  endif ()

  find_package (GFTL REQUIRED)
  find_package (GFTL_SHARED REQUIRED)
  find_package (PFLOGGER QUIET)
  find_package (YAFYAML REQUIRED)
endif ()

if (UFS_GOCART)
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

ecbuild_install_project (NAME GOCART)

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
