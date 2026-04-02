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
    # Apps              # not yet ported to MAPL3 - restore when ported
    Aerosol_GridComp
    # GOCART_GridComp   # not yet ported to MAPL3 - restore when ported
    GOCART2G_GridComp
  )
endif ()
