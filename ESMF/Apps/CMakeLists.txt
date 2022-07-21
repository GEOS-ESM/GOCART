esma_set_this ()

set (srcs
  GOCART2G_SimpleBundleMod.F90 
  aop_calculator.F90
  GOCART2G_AopMod.F90
  )

foreach (exe aop_calculator.x)
  get_filename_component (basename ${exe} NAME_WE) # base file name
  ecbuild_add_executable (
    TARGET ${exe}
    SOURCES ${basename}.F90
    LIBS ${this} MAPL)
endforeach ()


set (resource_files
   aop_calculator.rc
   )
install( FILES ${resource_files} 
   DESTINATION etc
   )
esma_add_library(${this} 
  SRCS ${srcs} 
  DEPENDENCIES MAPL Process_Library esmf NetCDF::NetCDF_Fortran
  )
