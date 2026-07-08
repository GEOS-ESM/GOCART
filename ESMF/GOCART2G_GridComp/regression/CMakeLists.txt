file(STRINGS "cases.txt" TEST_CASES)

include(esma_regression_tests)
esma_add_regression_tests(
  NAME      gocart2g
  DATA_PATH GEOSagcm_GridComp/GEOSphysics_GridComp/GEOSchem_GridComp/GOCART2G_GridComp
  TEST_CASES ${TEST_CASES}
)
