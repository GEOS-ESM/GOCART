if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/@GSW )
    set (gsw_ECBUILD_SYSTEM_INCLUDED TRUE)
endif()

esma_add_subdirectories (
        @MAPL
)
