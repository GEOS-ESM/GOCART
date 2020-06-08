esma_set_this ()

esma_add_library (${this}
  SRCS ${this}Mod.F90
  DEPENDENCIES Chem_Shared Chem_Base GMAO_mpeu)
target_include_directories (${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})

