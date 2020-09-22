esma_set_this ()

set (srcs
  Chem_MieTableMod2G.F90
  GOCART2G_Process.F90
  EmissionReaderMod.F90
  DustEmissionsMod.F90
  )

esma_add_library(${this}
  SRCS ${srcs}
  )

target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

if (PFUNIT_FOUND)
    add_subdirectory(tests EXCLUDE_FROM_ALL)
endif ()
