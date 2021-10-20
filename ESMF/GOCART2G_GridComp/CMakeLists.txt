esma_set_this ()

set (alldirs
  GA_Environment
  DU2G_GridComp
  SS2G_GridComp
  CA2G_GridComp
  NI2G_GridComp
  SU2G_GridComp
  )

set (srcs
  GOCART2G_GridCompMod.F90
  )

set (resource_files
   GOCART2G_GridComp.rc
   ) 

install( FILES ${resource_files} 
   DESTINATION etc
   )

set (dependencies MAPL Chem_Shared2G Process_Library esmf)
esma_add_library (${this}
  SRCS ${srcs}
  SUBCOMPONENTS ${alldirs}
  DEPENDENCIES ${dependencies}
  )

mapl_acg (${this}   GOCART2G_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)
