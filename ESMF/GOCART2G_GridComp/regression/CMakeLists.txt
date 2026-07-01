file(STRINGS "cases.txt" TEST_CASES)

if(APPLE)
  set(LD_PATH "DYLD_LIBRARY_PATH")
else()
  set(LD_PATH "LD_LIBRARY_PATH")
  # OpenMPI flags
  string(REPLACE " " ";" MPI_Fortran_LIBRARY_VERSION_LIST ${MPI_Fortran_LIBRARY_VERSION_STRING})
  list(GET MPI_Fortran_LIBRARY_VERSION_LIST 0 MPI_Fortran_LIBRARY_VERSION_FIRSTWORD)
  if(MPI_Fortran_LIBRARY_VERSION_FIRSTWORD MATCHES "Open")
    set(MPIEXEC_PREFLAGS "--oversubscribe;--bind-to;core;--map-by;core" CACHE STRING "Flags for mpiexec" FORCE)
  endif()
endif()

set(GOCART2G_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/GOCART2G_GridComp)
set(REGRESSION_DATA_DIR ${CMAKE_BINARY_DIR}/regression-data/${GOCART2G_PATH})

set(SYNC_WORK_DIR "${CMAKE_BINARY_DIR}/regression-data/sync_work/gocart2g")

# Locate the cmake script used internally by esma_sync_data()
find_file(ESMA_SYNC_DATA_SCRIPT
  NAMES esma_sync_data_script.cmake
  PATHS ${CMAKE_MODULE_PATH}
  NO_DEFAULT_PATH
  REQUIRED
)

# Fixture test: sync regression data from S3 at ctest (NOT configure) time
add_test(
  NAME gocart2g/sync_data
  COMMAND ${CMAKE_COMMAND}
    "-DLOCAL_DIR=${REGRESSION_DATA_DIR}"
    "-DWORK_DIR=${SYNC_WORK_DIR}"
    "-DESMA_SYNC_DATA_SCRIPT=${ESMA_SYNC_DATA_SCRIPT}"
    -P ${CMAKE_CURRENT_SOURCE_DIR}/sync_data.cmake
)
set_tests_properties(gocart2g/sync_data PROPERTIES
  FIXTURES_SETUP gocart2g_regression_data
)

# Regression tests
foreach(TEST_CASE ${TEST_CASES})
  add_test(
    NAME "${TEST_CASE}"
    COMMAND ${CMAKE_COMMAND}
      -DTEST_CASE=${TEST_CASE}
      -DMPIEXEC_EXECUTABLE=${MPIEXEC_EXECUTABLE}
      -DMPIEXEC_NUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}
      -DMY_BINARY_DIR=${CMAKE_BINARY_DIR}/bin
      "-DMPIEXEC_PREFLAGS=${MPIEXEC_PREFLAGS}"
      -DREGRESSION_DATA_DIR=${REGRESSION_DATA_DIR}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/run_case.cmake
  )
  set_tests_properties("${TEST_CASE}" PROPERTIES
    LABELS "REGRESSION"
    ENVIRONMENT "${LD_PATH}=${CMAKE_BINARY_DIR}/lib:$ENV{${LD_PATH}}"
    FIXTURES_REQUIRED gocart2g_regression_data
  )
endforeach()
