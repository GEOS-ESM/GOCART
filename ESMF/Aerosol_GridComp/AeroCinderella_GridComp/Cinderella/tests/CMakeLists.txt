add_definitions(-DUSE_MPI)

set (TEST_SRCS
        test_LinearFields.pf
        test_ConvertUnits_private.pf
        test_NOAA_TracerMap.pf
        test_NOAA_GOCARTtracers.pf
        test_MissingFields.pf
)

add_pfunit_ctest(Cinderella_private_tests.x
        TEST_SOURCES ${TEST_SRCS}
        LINK_LIBRARIES MAPL MAPL.pfio base_extras MAPL_pFUnit Cinderella
        EXTRA_INITIALIZE Initialize
        EXTRA_USE MAPL_pFUnit_Initialize
        MAX_PES 6
        )

add_dependencies(build-tests Cinderella_private_tests.x)

configure_file(test_GOCART_Tracers.rc test_GOCART_Tracers.rc COPYONLY)
configure_file(field_table            field_table            COPYONLY)
configure_file(GOCART_Tracers.yaml    GOCART_Tracers.yaml    COPYONLY)
configure_file(test_ConvertUnits.rc   test_ConvertUnits.rc   COPYONLY)
