esma_set_this()

set (srcs
  Chem_ConstMod.F90
  Chem_SettlingMod.F90
  DryDepositionMod.F90
  DustEmissionMod.F90
  WetRemovalMod.F90
  SeasaltEmissionMod.F90
  NitrateChemDriverMod.F90
  SulfateChemDriverMod.F90
  ConvectionMod.F90
  SmokeEmissionMod.F90
  VegLaiMod.F90
  GmiPrintError_mod.F90
  GmiTimeControl_mod.F90
  GmiSolar_mod.F90
  GmiResistance.F90
  GmiEmissUtils.F90
  Henrys_Law_cts.F90
  Chem_GroupMod.F90
  Bundle_IncrementMod.F90
  )

if (DOING_GEOS5)
  list (APPEND srcs mod_diag.F90 diaglist.F90)
  list (APPEND srcs Chem_UtilMod.F90 Chem_StateMod.F90)
else ()
  list (APPEND srcs Chem_UtilMod.F90)
  list (APPEND srcs Chem_UtilMod_geos4.F90)
endif ()

esma_add_library(${this} 
   SRCS ${srcs} 
   DEPENDENCIES GMAO_mpeu MAPL Chem_Base)
target_include_directories(${this} PUBLIC ${INC_ESMF} ${INC_NETCDF})
#add_definitions(-DDEBUG)
add_definitions(-DSPMD -DFVCHEM)

file (COPY gmi_drydep_data.h DESTINATION ${include_${this}})
file (COPY gmi_emiss_constants.h DESTINATION ${include_${this}})
file (COPY gmi_phys_constants.h DESTINATION ${include_${this}})
file (COPY gmi_time_constants.h DESTINATION ${include_${this}})

install(
  PROGRAMS bundleParser.py
  DESTINATION bin
)
  
  
