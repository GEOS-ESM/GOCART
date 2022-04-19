esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES MAPL Chem_Shared2G Process_Library
  )

if (OMP)
   target_link_libraries(${this} PRIVATE OpenMP::OpenMP_Fortran)
endif ()
