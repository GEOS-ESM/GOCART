add_executable(SyntheticApp.x
  SyntheticApp.F90 Synthetic_Driver.F90 Provider_GridComp.F90 Reciever_GridComp.F90 UFS_TestingCap.F90
  )
target_link_libraries(SyntheticApp.x PRIVATE
  gftl gftl-shared yafyaml pflogger MAPL MAPL.NUOPC
  )

add_executable(AerosolApp.x
  AerosolMain.F90 AerosolDriver.F90
  )
target_link_libraries(AerosolApp.x PRIVATE
  Aerosol_GridComp
  )
