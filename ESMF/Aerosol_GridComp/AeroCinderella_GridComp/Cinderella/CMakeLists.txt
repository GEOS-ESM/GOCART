esma_set_this()

set (srcs
        LinearFields.F90
        ConvertUnits_private.F90
        ConvertUnits.F90
        NOAA_TracerMap.F90
        NOAA_GOCARTtracers.F90
#       MissingFields.F90
        )

esma_add_library(${this}
        SRCS ${srcs}
        DEPENDENCIES MAPL esmf NetCDF::NetCDF_Fortran gftl-shared yafyaml pfunit pflogger
        )

if (PFUNIT_FOUND)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif ()
