esma_set_this ()

set (alldirs
  O3_GridComp
  CO_GridComp
  CO2_GridComp
  CFC_GridComp
  Rn_GridComp
  CH4_GridComp
  )

set (srcs
  GOCART_GridCompMod.F90
  Aero_GridCompMod.F90
  )

set (resource_files
   GOCARTdata_AerRegistry.rc
   )

install( FILES ${resource_files}
   DESTINATION etc
   )

set (dependencies Chem_Base Chem_Shared MAPL GMAO_mpeu esmf)
esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES ${dependencies})

