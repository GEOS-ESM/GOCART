ecbuild_add_executable(
        TARGET AerosolApp.x
        SOURCES AerosolMain.F90 AerosolDriver.F90
        LIBS gftl-shared MAPL Aerosol_GridComp
        INCLUDES ${esma_include}/NUOPC
)