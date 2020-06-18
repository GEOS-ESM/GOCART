if( EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/@GSW )
    set (gsw_ECBUILD_SYSTEM_INCLUDED TRUE)
endif()

esma_add_subdirectories (
        MAPL
        NCEP_Shared
        GMAO_Shared
        Chem_Base
        Chem_Shared
)
