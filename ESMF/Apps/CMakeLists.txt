ecbuild_add_executable(
        TARGET AerosolApp.x
        SOURCES AerosolMain.F90
        LIBS gftl-shared MAPL
        INCLUDES ${esma_include}/NUOPC
)