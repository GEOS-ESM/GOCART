esma_set_this()

set (alldirs
        AeroCinderella_GridComp
        )

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Aerosol_GridComp.F90)

    esma_add_library (${this}
            SRCS Aerosol_GridComp.F90
            SUBCOMPONENTS ${alldirs}
            DEPENDENCIES MAPL GEOS_Shared GOCART2G_GridComp
            INCLUDES ${INC_ESMF})

    mapl_acg (${this} Aerosol_StateSpecs.rc
            IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
            GET_POINTERS DECLARE_POINTERS)

else ()

    esma_add_subdirectories (${alldirs})

endif()
