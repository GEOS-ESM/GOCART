esma_set_this ()

add_subdirectory (GA_Environment)
add_subdirectory (DU2G_GridComp)
add_subdirectory (SS2G_GridComp)
add_subdirectory (CA2G_GridComp
add_subdirectory (NI2G_GridComp
add_subdirectory (SU2G_GridComp

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
  DEPENDENCIES ${dependencies}
  )

mapl_acg (${this}   GOCART2G_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)
