ecbuild_add_executable(
        TARGET AerosolApp.x
        SOURCES AerosolMain.F90 AerosolDriver.F90 UFS_TestingCap.F90
        LIBS gftl-shared MAPL Aerosol_GridComp
        INCLUDES ${esma_include}/NUOPC ${esma_include}/MAPL
)

ecbuild_add_executable(
        TARGET MAPL-NUOPCApp.x
        SOURCES MaplNuopcMain.F90 MaplNuopcDriver.F90 MaplProvide_GridComp.F90 MaplRecieve_GridComp.F90
        LIBS gftl-shared MAPL
        INCLUDES ${esma_include}/NUOPC ${esma_include}/MAPL
)
