set (HEMCO_EXTERNAL_CONFIG TRUE)
set (MAPL_ESMF TRUE)
set (BUILD_GEOS_INTERFACE TRUE)

if (UFS_GOCART)
  esma_add_subdirectories(
    Aerosol_GridComp
    GOCART2G_GridComp
    Shared
    UFS
  )
else ()
  esma_add_subdirectories(
    Shared
    Aerosol_GridComp
    GOCART_GridComp
    GOCART2G_GridComp
    HEMCO_GridComp
  )
endif ()
