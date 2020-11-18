set (HEMCO_EXTERNAL_CONFIG TRUE)
set (MAPL_ESMF TRUE)
set (BUILD_GEOS_INTERFACE TRUE)

esma_add_subdirectories(
        Apps
        Shared
        Aerosol_GridComp
        GOCART2G_GridComp
        HEMCO_GridComp
        UFS
)
