# Standalone shared library (loaded via dso: in YAML, like NI2G_GridComp.so),
# used only by the nitrate-sa standalone test harness below. Not part of the
# production GOCART2G_GridComp library, so it is built here rather than via
# the parent CMakeLists.txt's 'alldirs'/SUBCOMPONENTS list.
add_subdirectory(DuSsStub_GridComp)

file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME gocart2g
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/GOCART2G_GridComp
  TEST_CASES ${TEST_CASES}
)
