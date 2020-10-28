add_definitions(-DUSE_MPI)

set (TEST_SRCS
        test_LinearFields.pf
        test_ConvertUnits_private.pf
        test_MissingFields.pf
)

add_library (cinderella_extras
  MAPL_Initialize.F90
  )

add_pfunit_ctest(Cinderella_private_tests.x
        TEST_SOURCES ${TEST_SRCS}
        LINK_LIBRARIES MAPL MAPL.pfio cinderella_extras MAPL_pFUnit Cinderella
        EXTRA_INITIALIZE Initialize
        EXTRA_USE MAPL_pFUnit_Initialize
        MAX_PES 6
        )

add_dependencies(build-tests cinderella_extras Cinderella_private_tests.x)
target_link_libraries (cinderella_extras PUBLIC ${ESMF_LIBRARIES} Cinderella MAPL_pFUnit)


set_tests_properties(Cinderella_private_tests.x PROPERTIES LABELS "Cinderella")

configure_file(test_ConvertUnits.rc   test_ConvertUnits.rc   COPYONLY)
