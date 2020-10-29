add_definitions(-DUSE_MPI)

set (TEST_SRCS
        test_NOAA_TracerMap.pf
        test_NOAA_GOCARTtracers.pf
        test_NOAA_MAPLconfig.pf
)

add_library (configuration_extras
  MAPL_Initialize.F90
  )

add_pfunit_ctest(Configuration_tests.x
        TEST_SOURCES ${TEST_SRCS}
        LINK_LIBRARIES MAPL MAPL.pfio configuration_extras MAPL_pFUnit Cinderella
        EXTRA_INITIALIZE Initialize
        EXTRA_USE MAPL_pFUnit_Initialize
        MAX_PES 6
        )


add_dependencies(build-tests configuration_extras Configuration_tests.x)
target_link_libraries (configuration_extras PUBLIC ${ESMF_LIBRARIES} Configuration MAPL_pFUnit)

set_tests_properties(Configuration_tests.x PROPERTIES LABELS "Configuration")

configure_file(test_GOCART_Tracers.rc test_GOCART_Tracers.rc COPYONLY)
configure_file(field_table            field_table            COPYONLY)
configure_file(GOCART_Tracers.yaml    GOCART_Tracers.yaml    COPYONLY)
configure_file(NOAA_Config.yaml       NOAA_Config.yaml       COPYONLY)
