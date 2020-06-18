esma_set_this()

set (alldirs
        Cinderella
        )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/AeroCinderella_GridComp.F90)

    esma_add_library (${this}
            SRCS AeroCinderella_GridComp.F90
            SUBCOMPONENTS ${alldirs}
            DEPENDENCIES MAPL GEOS_Shared Chem_Shared
            INCLUDES ${INC_ESMF})

    mapl_acg (${this} AeroCinderella_StateSpecs.rc
            IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
            GET_POINTERS DECLARE_POINTERS)

else ()

    esma_add_subdirectories (${alldirs})

endif()
