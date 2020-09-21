add_definitions(-DUSE_MPI)

set (TEST_SRCS
    test_GOCART2G_Process.pf
    )

set(SRCS
    )


add_pfunit_ctest(Process_Library_tests.x
    TEST_SOURCES ${TEST_SRCS}
    OTHER_SOURCES ${SRCS}
    LINK_LIBRARIES Process_Library
    MAX_PES 4
    )

add_dependencies(build-tests Process_Library_tests.x)

set_tests_properties(Process_Library_tests.x PROPERTIES LABELS "Process_Library")
