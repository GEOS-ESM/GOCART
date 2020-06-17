add_definitions(-DUSE_MPI)

set (TEST_SRCS
        test_LinearFields.pf
        test_ConvertUnits_private.pf
        test_ChemTracers_private.pf
        test_FieldDictEntry.pf
        test_FieldDict.pf
)

add_pfunit_ctest(Cinderella_private_tests
        TEST_SOURCES ${TEST_SRCS}
        LINK_LIBRARIES MAPL MAPL.pfio base_extras MAPL_pFUnit Cinderella
        EXTRA_INITIALIZE Initialize
        EXTRA_USE MAPL_pFUnit_Initialize
        MAX_PES 6
        )

add_dependencies(build-tests Cinderella_private_tests)

configure_file(test_ChemTracers.rc test_ChemTracers.rc COPYONLY)
configure_file(test_ConvertUnits.rc test_ConvertUnits.rc COPYONLY)
configure_file(test_FieldDict.yaml test_FieldDict.yaml COPYONLY)
