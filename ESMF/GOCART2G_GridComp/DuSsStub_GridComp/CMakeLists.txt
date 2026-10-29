esma_set_this (OVERRIDE DuSsStub_GridComp)

esma_add_library(${this}
  SRCS DuSsStub_GridCompMod.F90
  DEPENDENCIES MAPL.generic MAPL ESMF::ESMF
  TYPE SHARED)
