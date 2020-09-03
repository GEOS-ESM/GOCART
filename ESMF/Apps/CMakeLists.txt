ecbuild_add_executable(
        TARGET SyntheticApp.x
        SOURCES SyntheticApp.F90 Synthetic_Driver.F90 Provider_GridComp.F90 Reciever_GridComp.F90 UFS_TestingCap.F90
        LIBS gftl-shared MAPL MAPL.NUOPC
        INCLUDES ${esma_include}/MAPL
)

ecbuild_add_executable(
        TARGET AerosolApp.x
        SOURCES AerosolMain.F90 AerosolDriver.F90
        LIBS gftl-shared MAPL MAPL.NUOPC Aerosol_GridComp
        INCLUDES ${esma_include}/MAPL
)
