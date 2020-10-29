esma_set_this()

set (srcs
        NOAA_Tracers.F90
        NOAA_MAPLconfig.F90
        NOAA_GOCARTtracers.F90
        )

esma_add_library(${this}
        SRCS ${srcs}
        DEPENDENCIES MAPL)

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (PFUNIT_FOUND)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif ()
