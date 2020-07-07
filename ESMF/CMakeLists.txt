set (HEMCO_EXTERNAL_CONFIG TRUE)
set (MAPL_ESMF TRUE)

esma_add_subdirectories(
        Apps
        include
        Shared
        Aerosol_GridComp
        GOCART_GridComp
        HEMCO_GridComp
)
