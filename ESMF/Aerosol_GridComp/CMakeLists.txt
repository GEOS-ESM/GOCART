esma_set_this()

if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Aerosol_GridComp.F90)

    esma_add_library (${this}
            SRCS Aerosol_GridComp.F90
            DEPENDENCIES GOCART2G_GridComp MAPL esmf)

    mapl_acg (${this} Aerosol_StateSpecs.rc
            IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
            GET_POINTERS DECLARE_POINTERS)

endif()
