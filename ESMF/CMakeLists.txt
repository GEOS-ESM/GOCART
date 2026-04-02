set (HEMCO_EXTERNAL_CONFIG TRUE)
set (MAPL_ESMF TRUE)
set (BUILD_GEOS_INTERFACE TRUE)

esma_add_subdirectories(
  Shared
  # Apps              # not yet ported to MAPL3 - restore when ported
  # Aerosol_GridComp  # not yet ported to MAPL3 - restore when ported
  # GOCART_GridComp   # not yet ported to MAPL3 - restore when ported
  GOCART2G_GridComp
)

if (UFS_GOCART)
  esma_add_subdirectories(
    # UFS             # not yet ported to MAPL3 - restore when ported
  )
endif ()
